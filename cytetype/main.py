import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Literal, get_args

import anndata
import numpy as np
from natsort import natsorted

from .api import submit_annotation_job, wait_for_completion
from .api.client import (
    upload_obs_duckdb as upload_obs_duckdb_file,
)
from .api.client import (
    upload_vars_h5 as upload_vars_h5_file,
)
from .api.exceptions import AuthenticationError, AuthorizationError
from .api.schemas import JobVisibility
from .config import get_default_api_url, load_credentials, logger, validate_api_url
from .core.artifacts import (
    _is_integer_valued,
    save_features_matrix,
)
from .core.artifacts import (
    save_obs_duckdb as save_obs_duckdb_file,
)
from .core.payload import build_annotation_payload, save_query_to_file
from .core.results import (
    fetch_remote_results,
    load_local_results,
    store_annotations,
    store_job_details,
)
from .preprocessing import (
    aggregate_cluster_metadata,
    aggregate_expression_percentages,
    extract_marker_genes,
    extract_visualization_coordinates,
    resolve_gene_symbols_column,
    validate_adata,
)
from .preprocessing.validation import (
    _generate_unique_na_label,
    materialize_canonical_gene_symbols_column,
)

__all__ = ["CyteType"]

# The server gives API keys without an organization license this prefix and refuses
# organization visibility for them. Organization keys are "cyt_" and random characters,
# so about 1 in 4096 of them also starts with it: the prefix alone is not proof.
_PERSONAL_API_KEY_PREFIX = "cyt_p_"


def _get_cytetype_version() -> str | None:
    try:
        return version("cytetype")
    except PackageNotFoundError:
        return None


def _unreadable_credentials_message(error: Exception) -> str:
    return (
        f"Saved CyteType credentials could not be read ({error}). Pass `auth_token`, "
        "or make sure XDG_CONFIG_HOME or HOME points to a readable configuration "
        "directory."
    )


def _remove_partial_artifact(*paths: str) -> None:
    """Delete files left behind by a failed artifact build."""
    for path in paths:
        try:
            Path(path).unlink(missing_ok=True)
        except OSError as exc:
            logger.warning(f"Failed to delete partial artifact {path}: {exc}")


class CyteType:
    """CyteType class for characterizing clusters from single-cell RNA-seq data.

    This class provides an object-oriented interface for cluster characterization using the CyteType API.
    The expensive data preparation steps (validation, expression percentage calculation, and marker
    gene extraction) are performed during initialization, allowing for efficient reuse when making
    multiple requests with different parameters.

    Initialization also writes the `vars.h5` and `obs.duckdb` artifact files and adds a temporary
    `__cytetype_gene_symbols` column (with a numeric suffix if that name is taken) to `adata.var`
    holding the cleaned gene symbols. Call `cleanup()` after `run()` to delete the files and remove
    the column; a file whose build failed is deleted right away. After `cleanup()`, create a new
    `CyteType` object to run another job with artifacts.
    """

    # Type annotations for instance attributes
    adata: anndata.AnnData
    group_key: str
    rank_key: str
    gene_symbols_column: str | None
    n_top_genes: int
    pcent_batch_size: int
    coordinates_key: str | None
    cluster_map: dict[str, str]
    clusters: list[str]
    expression_percentages: dict[str, dict[str, float]]
    marker_genes: dict[str, list[str]]
    group_metadata: dict[str, dict[str, dict[str, int]]]
    visualization_data: dict[str, Any]
    _auth_token_api_url: str | None
    __version__: str | None = _get_cytetype_version()

    def __init__(
        self,
        adata: anndata.AnnData,
        group_key: str,
        rank_key: str = "rank_genes_groups",
        gene_symbols_column: str | None = None,
        n_top_genes: int = 50,
        aggregate_metadata: bool = True,
        min_percentage: int = 10,
        pcent_batch_size: int = 5000,
        coordinates_key: str = "X_umap",
        max_cells_per_group: int = 1000,
        vars_h5_path: str = "vars.h5",
        obs_duckdb_path: str = "obs.duckdb",
        max_metadata_categories: int = 500,
        api_url: str | None = None,
        auth_token: str | None = None,
        label_na: bool = False,
    ) -> None:
        """Initialize CyteType with AnnData object and perform data preparation.

        Args:
            adata (anndata.AnnData): The AnnData object to annotate. Must contain log1p-normalized
                gene expression data in `adata.X` and gene names in `adata.var_names`.
            group_key (str): The key in `adata.obs` containing the cluster labels.
                These clusters will receive cell type annotations.
            rank_key (str, optional): The key in `adata.uns` containing differential expression
                results from `sc.tl.rank_genes_groups`. Must use the same `groupby` as `group_key`.
                Defaults to "rank_genes_groups".
            gene_symbols_column (str | None, optional): Name of the column in `adata.var` that
                contains gene symbols. When None (default), auto-detects by checking well-known
                column names (feature_name, gene_symbols, etc.), then adata.var_names, then
                heuristic scan of all var columns. Defaults to None.
            n_top_genes (int, optional): Number of top marker genes per cluster to extract during
                initialization. Higher values may improve annotation quality but increase memory usage.
                Defaults to 50.
            aggregate_metadata (bool, optional): Whether to aggregate metadata from the AnnData object.
                Defaults to True.
            min_percentage (int, optional): Minimum percentage of cells in a group to include in the
                cluster context. Defaults to 10.
            pcent_batch_size (int, optional): Number of cells to process per chunk when
                calculating expression percentages. Defaults to 5000.
            coordinates_key (str, optional): Key in `adata.obsm` with 2D coordinates for the
                report's visualization. This key is tried first, then "X_umap", "X_tsne" and
                "X_pca", all matched case-insensitively; the first array with `adata.n_obs`
                rows and at least 2 columns is used (only its first 2 columns). If none
                qualifies, visualization is disabled. Defaults to "X_umap".
            max_cells_per_group (int, optional): Maximum number of cells to sample per group
                for visualization. If a group has more cells than this limit, a random sample
                will be taken. Defaults to 1000.
            vars_h5_path (str, optional): Path where the expression artifact (normalized and,
                when available, raw counts) is written during initialization and uploaded by
                `run()`. A relative path is written to the current working directory.
                `cleanup()` deletes the file. Defaults to "vars.h5".
            obs_duckdb_path (str, optional): Path where the cell metadata artifact (`adata.obs`
                and the visualization coordinates) is written during initialization and
                uploaded by `run()`. A relative path is written to the current working
                directory. `cleanup()` deletes the file. Defaults to "obs.duckdb".
            max_metadata_categories (int, optional): Maximum number of unique values a categorical
                obs column may have to be included in cluster metadata aggregation. Columns with
                more unique values (e.g. cell barcodes, per-cell IDs) are skipped to avoid
                excessive memory usage. Defaults to 500.
            api_url (str | None, optional): URL for the CyteType API endpoint. Defaults to
                `CYTETYPE_API_URL` when set, otherwise "https://cytetype.nygen.io".
            auth_token (str | None, optional): API key sent as a Bearer token. When None, the key
                saved by `cytetype setup` (or `cytetype login`) for `api_url` is read each time it
                is needed, and initialization logs a note if there is none. Defaults to None.
            label_na (bool, optional): If True, cells with NaN values in the
                ``group_key`` column are assigned an ``'Unknown'`` cluster label
                (or ``'Unknown 2'``, etc. if that label already exists) in a new
                AnnData object, which ``run()`` returns; the labels in the original
                object's ``obs`` are not changed. If False (default), a ``ValueError``
                is raised instead.

        Raises:
            KeyError: If the required keys are missing in `adata.obs` or `adata.uns`
            ValueError: If the data format is incorrect, there are validation errors, or
                `api_url` is not a valid server origin
        """
        self.adata = adata
        self.group_key = group_key
        self.rank_key = rank_key
        self.n_top_genes = n_top_genes
        self.pcent_batch_size = pcent_batch_size
        self.coordinates_key = coordinates_key
        self.max_cells_per_group = max_cells_per_group
        self.api_url = validate_api_url(
            api_url if api_url is not None else get_default_api_url()
        )
        self.auth_token = auth_token
        self._auth_token_api_url = self.api_url if auth_token else None
        self._artifact_build_errors: list[tuple[str, Exception]] = []
        self._vars_h5_path: str | None = None
        self._obs_duckdb_path: str | None = None
        self._original_gene_symbols_column: str | None = None
        self._temporary_gene_symbols_column: str | None = None
        self._artifacts_cleaned_up = False

        # Report a missing sign-in before the slow preparation, not only at run().
        if not auth_token:
            self._check_saved_credentials()

        try:
            self.gene_symbols_column = resolve_gene_symbols_column(
                adata, gene_symbols_column
            )
            self._original_gene_symbols_column = self.gene_symbols_column

            self.coordinates_key = validate_adata(
                adata,
                group_key,
                rank_key,
                self.gene_symbols_column,
                coordinates_key,
                label_na=label_na,
            )

            if label_na:
                nan_mask = adata.obs[group_key].isna()
                if nan_mask.any():
                    n_nan = int(nan_mask.sum())
                    pct = round(100 * n_nan / adata.n_obs, 1)
                    existing_labels = set(
                        str(v) for v in adata.obs[group_key].dropna().unique()
                    )
                    na_label = _generate_unique_na_label(existing_labels)
                    logger.warning(
                        f"Relabeling {n_nan} cells ({pct}%) with NaN values "
                        f"in '{group_key}' as '{na_label}'."
                    )
                    adata = anndata.AnnData(
                        X=adata.X,
                        obs=adata.obs.copy(),
                        var=adata.var,
                        uns=adata.uns,
                        obsm=adata.obsm,
                        varm=adata.varm,
                        layers=adata.layers,
                        obsp=adata.obsp,
                        varp=adata.varp,
                    )
                    col = adata.obs[group_key]
                    if hasattr(col, "cat"):
                        col = col.cat.add_categories(na_label)
                    adata.obs[group_key] = col.fillna(na_label)
                    self.adata = adata

            (
                self.gene_symbols_column,
                self._original_gene_symbols_column,
            ) = materialize_canonical_gene_symbols_column(
                adata, self.gene_symbols_column
            )
            self._temporary_gene_symbols_column = self.gene_symbols_column

            # Use original labels as IDs if all are short (<=3 chars), otherwise enumerate
            _unique_group_categories: list[str | int] = natsorted(
                adata.obs[group_key].unique().tolist()
            )
            _short_ids = all(len(str(x)) <= 3 for x in _unique_group_categories)
            self.cluster_map = {
                str(x): str(x) if _short_ids else str(n)
                for n, x in enumerate(_unique_group_categories)
            }
            self.clusters = [
                self.cluster_map[str(x)] for x in adata.obs[group_key].values.tolist()
            ]

            gene_names = adata.var[self.gene_symbols_column].tolist()
            self.expression_percentages = aggregate_expression_percentages(
                adata=adata,
                clusters=self.clusters,
                gene_names=gene_names,
                cell_batch_size=pcent_batch_size,
            )

            logger.info("Extracting marker genes...")
            self.marker_genes = extract_marker_genes(
                adata=self.adata,
                cell_group_key=self.group_key,
                rank_genes_key=self.rank_key,
                cluster_map=self.cluster_map,
                n_top_genes=n_top_genes,
                gene_symbols_col=self.gene_symbols_column,
            )

            if aggregate_metadata:
                logger.info("Aggregating cluster metadata...")
                self.group_metadata = aggregate_cluster_metadata(
                    adata=self.adata,
                    group_key=self.group_key,
                    min_percentage=min_percentage,
                    max_categories=max_metadata_categories,
                )
                # Replace keys in group_metadata using cluster_map
                self.group_metadata = {
                    self.cluster_map.get(str(key), str(key)): value
                    for key, value in self.group_metadata.items()
                }
                self.group_metadata = {
                    k: self.group_metadata[k]
                    for k in sorted(self.group_metadata.keys())
                }
            else:
                self.group_metadata = {}

            # Prepare visualization data with sampling
            sampled_coordinates, sampled_cluster_labels = (
                extract_visualization_coordinates(
                    adata=adata,
                    coordinates_key=self.coordinates_key,
                    group_key=self.group_key,
                    cluster_map=self.cluster_map,
                    max_cells_per_group=self.max_cells_per_group,
                )
            )

            self.visualization_data = {
                "coordinates": sampled_coordinates,
                "clusters": sampled_cluster_labels,
            }

            # Resolve raw counts once and cache
            self._raw_counts_result = self._resolve_raw_counts()
            if self._raw_counts_result is None:
                logger.warning(
                    "No integer raw counts found in adata.layers['counts'], "
                    "adata.raw.X, or adata.X. Skipping raw counts in vars.h5."
                )

            # Build vars.h5
            try:
                raw_mat, raw_col_indices = (
                    self._raw_counts_result
                    if self._raw_counts_result is not None
                    else (None, None)
                )
                save_features_matrix(
                    out_file=vars_h5_path,
                    mat=self.adata.X,
                    var_df=self.adata.var,
                    var_names=self.adata.var_names,
                    raw_mat=raw_mat,
                    raw_col_indices=raw_col_indices,
                    gene_symbols_column=self.gene_symbols_column,
                )
                sys.stderr.flush()
                self._vars_h5_path = vars_h5_path
            except Exception as exc:
                logger.warning(f"vars.h5 artifact failed during build: {exc}")
                self._artifact_build_errors.append(("vars_h5", exc))
                _remove_partial_artifact(vars_h5_path)

            # Build obs.duckdb
            try:
                logger.info("Writing obs data to duckdb artifact...")
                obsm_coordinates = (
                    self.adata.obsm[self.coordinates_key]
                    if self.coordinates_key and self.coordinates_key in self.adata.obsm
                    else None
                )
                save_obs_duckdb_file(
                    out_file=obs_duckdb_path,
                    obs_df=self.adata.obs,
                    obsm_coordinates=obsm_coordinates,
                    coordinates_key=self.coordinates_key,
                )
                sys.stderr.flush()
                self._obs_duckdb_path = obs_duckdb_path
            except Exception as exc:
                logger.warning(f"obs.duckdb artifact failed during build: {exc}")
                self._artifact_build_errors.append(("obs_duckdb", exc))
                _remove_partial_artifact(obs_duckdb_path, f"{obs_duckdb_path}.wal")

            logger.info("Data preparation completed. Ready for submitting jobs.")
        except Exception:
            self._cleanup_temporary_gene_symbols_column()
            raise

    def _cleanup_temporary_gene_symbols_column(self) -> None:
        temp_column = self._temporary_gene_symbols_column
        if temp_column is None:
            return

        if temp_column in self.adata.var.columns:
            del self.adata.var[temp_column]
            logger.info(
                f"Deleted temporary canonical gene-symbol column '{temp_column}'."
            )

        self.gene_symbols_column = self._original_gene_symbols_column
        self._temporary_gene_symbols_column = None

    def _resolve_raw_counts(
        self,
    ) -> "tuple[Any, np.ndarray | None] | None":
        if "counts" in self.adata.layers:
            mat = self.adata.layers["counts"]
            if _is_integer_valued(mat):
                return mat, None

        if self.adata.raw is not None:
            raw_mat = self.adata.raw.X
            col_indices = self.adata.raw.var_names.get_indexer(self.adata.var_names)
            if (col_indices == -1).any():
                logger.warning(
                    "Some var_names not found in adata.raw — skipping adata.raw.X as raw counts source."
                )
            elif _is_integer_valued(raw_mat):
                return raw_mat, col_indices.astype(np.intp)

        if _is_integer_valued(self.adata.X):
            return self.adata.X, None

        return None

    def _upload_artifacts(
        self,
        auth_token: str,
        upload_timeout_seconds: int,
        upload_max_workers: int = 4,
    ) -> tuple[dict[str, str], list[tuple[str, Exception]]]:
        """Upload pre-built artifact files to the server.

        Returns (uploaded_ids, errors) so the caller can decide whether
        partial success is acceptable.
        """
        uploaded: dict[str, str] = {}
        errors: list[tuple[str, Exception]] = list(self._artifact_build_errors)
        timeout = (60.0, float(upload_timeout_seconds))

        # --- vars.h5 upload ---
        if self._vars_h5_path is not None:
            try:
                logger.info("Uploading vars.h5 artifact...")
                vars_upload = upload_vars_h5_file(
                    self.api_url,
                    auth_token,
                    self._vars_h5_path,
                    timeout=timeout,
                    max_workers=upload_max_workers,
                )
                if vars_upload.file_kind != "vars_h5":
                    raise ValueError(
                        f"Unexpected upload file_kind for vars artifact: {vars_upload.file_kind}"
                    )
                uploaded["vars_h5"] = vars_upload.upload_id
            except Exception as exc:
                logger.warning(f"vars.h5 upload failed: {exc}")
                errors.append(("vars_h5", exc))

        print()

        # --- obs.duckdb upload ---
        if self._obs_duckdb_path is not None:
            try:
                logger.info("Uploading obs.duckdb artifact...")
                obs_upload = upload_obs_duckdb_file(
                    self.api_url,
                    auth_token,
                    self._obs_duckdb_path,
                    timeout=timeout,
                    max_workers=upload_max_workers,
                )
                if obs_upload.file_kind != "obs_duckdb":
                    raise ValueError(
                        f"Unexpected upload file_kind for obs artifact: {obs_upload.file_kind}"
                    )
                uploaded["obs_duckdb"] = obs_upload.upload_id
            except Exception as exc:
                logger.warning(f"obs.duckdb upload failed: {exc}")
                errors.append(("obs_duckdb", exc))

        return uploaded, errors

    def cleanup(self) -> None:
        """Delete the artifact files built during initialization.

        Call this after run() completes to remove the vars.h5 and obs.duckdb
        files from disk and drop the temporary canonical gene-symbol column.
        Paths are cleared so repeated calls are safe. A later `run()` raises
        unless `require_artifacts=False`; create a new `CyteType` object to
        rebuild the artifacts.
        """
        for attr, path in [
            ("_vars_h5_path", self._vars_h5_path),
            ("_obs_duckdb_path", self._obs_duckdb_path),
        ]:
            if path is not None:
                try:
                    Path(path).unlink(missing_ok=True)
                    logger.info(f"Deleted artifact file: {path}")
                except OSError as exc:
                    logger.warning(f"Failed to delete artifact {path}: {exc}")
                setattr(self, attr, None)
                self._artifacts_cleaned_up = True

        self._cleanup_temporary_gene_symbols_column()

    def _check_saved_credentials(self) -> None:
        """Report early when run() will not find a saved API key for ``self.api_url``.

        Only checks: it never raises and does not keep the token, which
        ``_resolve_auth_token`` reads when it is needed.
        """
        try:
            credentials = load_credentials(self.api_url)
        except ValueError as error:
            logger.warning(
                f"{error}. Unless you pass auth_token to run(), run `cytetype logout`, "
                "then `cytetype setup`, before calling run()."
            )
            return
        except (OSError, RuntimeError) as error:
            # For example an unreadable file, or no home directory in a container.
            logger.warning(_unreadable_credentials_message(error))
            return

        if credentials is None:
            # Not a problem when auth_token is passed to run(), as in CI.
            logger.info(
                f"No saved CyteType sign-in found for {self.api_url}. Unless you pass "
                "auth_token to run(), run `cytetype setup` (or `cytetype login`) "
                "before calling run()."
            )

    def _resolve_auth_token(
        self,
        api_url: str,
        auth_token: str | None = None,
    ) -> str:
        """Return the API key to send to ``api_url``.

        A token given to ``__init__`` or ``run()`` is kept for its API URL. Saved
        credentials are read on every call and never stored on the instance, so a
        key replaced with ``cytetype setup --force`` is used right away.
        """
        normalized_api_url = validate_api_url(api_url)
        if auth_token:
            self.auth_token = auth_token
            self._auth_token_api_url = normalized_api_url
            return auth_token

        if self.auth_token and self._auth_token_api_url == normalized_api_url:
            return self.auth_token

        try:
            credentials = load_credentials(normalized_api_url)
        except ValueError as error:
            raise AuthenticationError(
                f"{error}. Run `cytetype logout`, then `cytetype setup`.",
                error_code="AUTHENTICATION_REQUIRED",
            ) from error
        except (OSError, RuntimeError) as error:
            raise AuthenticationError(
                _unreadable_credentials_message(error),
                error_code="AUTHENTICATION_REQUIRED",
            ) from error

        if credentials is None:
            raise AuthenticationError(
                "CyteType sign-in is required. Run `cytetype setup` first.",
                error_code="AUTHENTICATION_REQUIRED",
            )
        return credentials.apiToken

    def run(
        self,
        study_context: str,
        llm_configs: list[dict[str, Any]] | None = None,
        metadata: dict[str, Any] | None = None,
        n_parallel_clusters: int = 2,
        results_prefix: str = "cytetype",
        poll_interval_seconds: int = 10,
        timeout_seconds: int = 7200,
        api_url: str | None = None,
        auth_token: str | None = None,
        save_query: bool = True,
        query_filename: str = "query.json",
        upload_timeout_seconds: int = 3600,
        upload_max_workers: int = 4,
        require_artifacts: bool = True,
        show_progress: bool = True,
        override_existing_results: bool = False,
        visibility: Literal["private", "org", "public"] | None = None,
    ) -> anndata.AnnData:
        """Perform cluster characterization using the CyteType API.

        Args:
            study_context (str): Biological context for the experimental setup. Required.
                For example, include information about 'organisms', 'tissues', 'diseases', 'developmental_stages',
                'single_cell_methods', 'experimental_conditions'.
            llm_configs (list[dict[str, Any]] | None, optional): Configuration for the large language
                models to be used. Each dict must match the LLMModelConfig schema with required 'provider' and 'name' fields.
                Either 'apiKey' or all AWS credentials ('awsAccessKeyId', 'awsSecretAccessKey', 'awsDefaultRegion') must be provided.
                Defaults to None, using the API's default model.
            metadata (dict[str, Any] | None, optional): Custom metadata tags to include in the report header.
                Values must be strings; values that look like URLs will be made clickable in the report.
                A "title" key (case-insensitive) also names the job in the dashboard. Defaults to None.
            n_parallel_clusters (int, optional): Number of parallel requests to make to the model (1-50).
                Note that high values can lead to rate limit errors. Ignored when `llm_configs` is not
                provided; the server then uses its own setting. Defaults to 2.
            results_prefix (str, optional): Prefix for the keys added to `adata.obs` and `adata.uns`
                to store results (see Returns). Defaults to "cytetype".
            poll_interval_seconds (int, optional): How often (in seconds) to check for results from
                the API. Defaults to 10.
            timeout_seconds (int, optional): Maximum time (in seconds) to wait for API results before
                raising a timeout error. Defaults to 7200.
            api_url (str, optional): URL for the CyteType API endpoint. If not provided, the API URL
                specified in the initialization will be used. If provided, it will override the API URL
                specified in the initialization.
            auth_token (str | None, optional): API key sent as a Bearer token. If provided, it overrides
                the token specified in the initialization and is kept for later calls to the same API
                URL. If neither is provided, the key saved by `cytetype setup` (or `cytetype login`)
                for the API URL is read at each call. Defaults to None.
            save_query (bool, optional): Whether to save the query to a file. Defaults to True.
            query_filename (str, optional): Filename for saving the query when save_query is True.
                Defaults to "query.json".
            upload_timeout_seconds (int, optional): Socket read timeout used for each artifact upload.
                Defaults to 3600.
            upload_max_workers (int, optional): Number of parallel threads used to upload file
                chunks. Each worker holds one chunk in memory; the chunk size is set by the server.
                Defaults to 4.
            require_artifacts (bool, optional): Whether to raise an error if artifact building or
                uploading fails, or if `cleanup()` already deleted the artifacts. When True (default),
                any artifact failure stops the run. Set to False to skip artifacts and continue with
                annotation only. Authentication and authorization errors always stop the run.
                Defaults to True.
            show_progress (bool, optional): Whether to display progress updates with spinner and
                cluster status. Set to False to disable all visual progress output. Defaults to True.
            override_existing_results (bool, optional): Whether to override existing results with the
                same results_prefix. If False (default) and results already exist, will raise an error
                to prevent accidental overwriting. Set to True to explicitly allow overwriting.
                Defaults to False.
            visibility (Literal["private", "org", "public"] | None, optional): Who can view the
                report. Defaults to None, which uses the server default ("private"). "org" shares
                the report with your organization and needs an organization API key. With a key
                that looks like a personal key (starting with `cyt_p_`), a warning is logged
                before the artifacts are uploaded. The server rejects "org" for any key without
                an organization license, such as a personal key, when the job is submitted,
                after the upload. "public" makes the report viewable by anyone with the link, and
                any signed-in CyteType user can find it in Chat's report search.

        Returns:
            anndata.AnnData: The AnnData object the results were added to. This is the object passed
            to `CyteType()`, except when `label_na=True` relabeled NaN cluster labels: then it is
            a new object and the obs columns exist only there, so use the returned object. It
            gets the following additions:

            - `adata.obs[f"{results_prefix}_annotation_{group_key}"]`: Cell type annotations
            - `adata.obs[f"{results_prefix}_cellOntologyTerm_{group_key}"]`: Cell Ontology term names
            - `adata.obs[f"{results_prefix}_cellOntologyTermID_{group_key}"]`: Cell Ontology term IDs
            - `adata.obs[f"{results_prefix}_cellState_{group_key}"]`: Cell states
            - `adata.obs[f"{results_prefix}_reviewerConfidence_{group_key}"]`: Reviewer confidence
              ("High", "Moderate" or "Low"; "Unknown" when missing)
            - `adata.uns[f"{results_prefix}_results"]`: The job ID and the normalized results
              (see `get_results()`) as a JSON string under "result"
            - `adata.uns[f"{results_prefix}_jobDetails"]`: The job ID, report URL and API URL,
              plus the requested visibility when given

            The obs columns are categorical. Cells in clusters without a result get "Unknown"
            ("" for the cell state).

        Raises:
            ValueError: If results with the same prefix already exist and override_existing_results
                is False, if `api_url` is not a valid server origin, if `visibility` is not a valid
                value, if `upload_timeout_seconds` is not positive, or if `cleanup()` deleted the
                artifacts and `require_artifacts` is True.
            AuthenticationError: If no API key is available (including saved credentials that cannot
                be read) or the server rejects it.
            pydantic.ValidationError: If `llm_configs`, `metadata`, `n_parallel_clusters` or the
                prepared input data are invalid.
            RateLimitError: If the server's rate limit for job submissions is exceeded.
            QuotaExceededError: If the account's quota is exceeded.
            AuthorizationError: If the account may not submit jobs or upload files (for example, an
                expired license).
            LLMValidationError: If a model in `llm_configs` fails the server's validation.
            APIError: For other errors returned by the server.
            NetworkError: If the server cannot be reached or returns an unexpected response.
            cytetype.api.exceptions.TimeoutError: If a request times out or the job does not complete
                within `timeout_seconds` (also importable as `cytetype.CyteTypeTimeoutError`).
            JobFailedError: If the job fails on the server.
            Exception: With require_artifacts=True, the original exception raised while building or
                uploading an artifact.

        """
        # Check for existing results
        job_details_key = f"{results_prefix}_jobDetails"
        if job_details_key in self.adata.uns and not override_existing_results:
            existing_job_id = self.adata.uns[job_details_key].get("job_id", "unknown")
            raise ValueError(
                f"Results with prefix '{results_prefix}' already exist in this AnnData object "
                f"(job_id: {existing_job_id}). To prevent accidental overwriting, please either:\n"
                f"  1. Use a different results_prefix (e.g., results_prefix='cytetype_v2')\n"
                f"  2. Set override_existing_results=True to explicitly overwrite existing results\n"
                f"  3. Use annotator.get_results(results_prefix='{results_prefix}') to retrieve existing results"
            )

        if visibility is not None and visibility not in get_args(JobVisibility):
            raise ValueError(
                "visibility must be 'private', 'org', 'public' or None, "
                f"got {visibility!r}"
            )

        if self._artifacts_cleaned_up:
            if require_artifacts:
                raise ValueError(
                    "The artifact files were deleted by cleanup(). Create a new "
                    "CyteType object to rebuild them, or pass require_artifacts=False "
                    "to run without them."
                )
            logger.warning(
                "The artifact files were deleted by cleanup(). Continuing without them."
            )

        if api_url is not None:
            self.api_url = validate_api_url(api_url)
        token = self._resolve_auth_token(self.api_url, auth_token)
        # The server refuses a personal key only at submission, after the upload, so warn
        # now. The prefix is not proof, so the run goes on and the server decides.
        if visibility == "org" and token.startswith(_PERSONAL_API_KEY_PREFIX):
            logger.warning(
                "visibility='org' needs an organization API key, and this key looks like "
                "a personal key (cyt_p_...). If it is, the server rejects the job after "
                "the upload. Use an organization key, or choose 'private' or 'public'."
            )
        if upload_timeout_seconds <= 0:
            raise ValueError("upload_timeout_seconds must be greater than 0")

        # Validate the base payload before uploading.
        payload = build_annotation_payload(
            study_context,
            metadata,
            self.cluster_map,
            self.group_metadata,
            self.marker_genes,
            self.visualization_data,
            self.expression_percentages,
            n_parallel_clusters,
            llm_configs,
            visibility=visibility,
        )

        uploaded_file_refs, artifact_errors = self._upload_artifacts(
            token,
            upload_timeout_seconds=upload_timeout_seconds,
            upload_max_workers=upload_max_workers,
        )
        if uploaded_file_refs:
            payload["uploaded_files"] = uploaded_file_refs

        if artifact_errors:
            # Credential problems are not artifact bugs, and the job submission would
            # fail the same way, so raise them with the fix instead of the issue link.
            for _, error in artifact_errors:
                if isinstance(error, AuthenticationError):
                    raise AuthenticationError(
                        f"Artifact upload was rejected: {error.message.rstrip('.')}. "
                        "Run `cytetype setup --force` to sign in again (or `cytetype "
                        "login` to save a new key), or pass a valid `auth_token`.",
                        error_code=error.error_code,
                    ) from error
                if isinstance(error, AuthorizationError):
                    raise AuthorizationError(
                        f"Artifact upload was not authorized: {error.message.rstrip('.')}. "
                        f"Check your account and license at {self.api_url}/dashboard.",
                        error_code=error.error_code,
                    ) from error

            failed_names = ", ".join(name for name, _ in artifact_errors)
            if require_artifacts:
                logger.error(
                    f"Artifact build/upload failed for: {failed_names}. "
                    "Rerun with `require_artifacts=False` to skip this error.\n"
                    "Please report the error below in a new issue at "
                    "https://github.com/NygenAnalytics/CyteType\n"
                    f"({type(artifact_errors[0][1]).__name__}: {str(artifact_errors[0][1]).strip()})"
                )
                raise artifact_errors[0][1]
            logger.warning(
                f"Artifact build/upload failed for: {failed_names}. "
                "Continuing without those artifacts. "
                "Set `require_artifacts=True` to see the full traceback."
            )

        # Save query if requested
        if save_query:
            save_query_to_file(payload["input_data"], query_filename)

        # Submit job and store details
        print()
        job_id = submit_annotation_job(self.api_url, token, payload)
        store_job_details(
            self.adata, job_id, self.api_url, results_prefix, visibility=visibility
        )

        # Wait for completion
        result = wait_for_completion(
            self.api_url,
            token,
            job_id,
            poll_interval_seconds,
            timeout_seconds,
            show_progress,
        )

        # Store results
        store_annotations(
            self.adata,
            result,
            job_id,
            results_prefix,
            self.group_key,
            self.clusters,
            check_unannotated=True,
        )

        report_url = self.adata.uns[f"{results_prefix}_jobDetails"]["report_url"]
        logger.info(f"Report:\n{report_url}")

        return self.adata

    def get_results(
        self,
        results_prefix: str = "cytetype",
    ) -> dict[str, Any] | None:
        """Retrieve the CyteType results for this AnnData object.

        Results stored in `adata.uns[f"{results_prefix}_results"]` are returned first. They are
        re-normalized, so results saved by older versions also get the filled fields and the
        `reviewerConfidence`, `isHeterogeneous` and `heterogeneityType` keys (obs columns are not
        added in this case). Otherwise, if `adata.uns[f"{results_prefix}_jobDetails"]` holds a job
        ID, the API is asked for the job status and, if the job has completed, for its results
        (up to two requests, no polling). Fetched results are written the same way as by `run()`,
        to the AnnData object that `run()` returns (`self.adata`; a new object if `label_na=True`
        relabeled NaN cluster labels).

        Args:
            results_prefix: The prefix used when storing results. Defaults to "cytetype".

        Returns:
            dict[str, Any] | None: The normalized results with "annotations", "summary",
            "clusterCategories", "studyContext" and "raw_annotations" (the API's per-cluster
            results, unchanged). Each annotation has "clusterId", "annotation", "ontologyTerm",
            "ontologyTermID", "granularAnnotation", "cellState", "justification",
            "supportingMarkers", "conflictingMarkers", "missingExpression",
            "unexpectedExpression", "reviewerConfidence" ("High", "Moderate", "Low" or None),
            "isHeterogeneous" (bool or None) and "heterogeneityType" (str or None). None if no
            results are stored and they could not be fetched (no job details, job not
            completed, or the request failed).

        Raises:
            ValueError: If the API URL stored in the job details is not a valid server origin.
            AuthenticationError: If results must be fetched and no API key is available for the
                job's API URL (including saved credentials that cannot be read).
        """
        logger.info(f"Retrieving results with prefix: '{results_prefix}'")

        # Try loading from local storage first
        local_result = load_local_results(self.adata, results_prefix)
        if local_result:
            return local_result

        # Try fetching from API if job details exist
        job_details_key = f"{results_prefix}_jobDetails"
        if job_details_key not in self.adata.uns:
            logger.info(
                "No results found locally and no job details available for retrieval."
            )
            return None

        job_details = self.adata.uns[job_details_key]
        job_id = job_details.get("job_id")

        if not job_id:
            logger.error("Job details found but missing job_id.")
            return None

        stored_api_url = job_details.get("api_url")
        job_api_url = validate_api_url(
            self.api_url if stored_api_url is None else stored_api_url
        )
        auth_token = self._resolve_auth_token(job_api_url)
        return fetch_remote_results(
            self.adata,
            job_id,
            job_api_url,
            auth_token,
            results_prefix,
            self.group_key,
            self.clusters,
        )
