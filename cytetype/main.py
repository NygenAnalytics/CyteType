from pathlib import Path
from typing import Any
from importlib.metadata import PackageNotFoundError, version

import anndata
from natsort import natsorted

from .config import logger
from .api import submit_annotation_job, wait_for_completion
from .api.client import (
    upload_obs_duckdb as upload_obs_duckdb_file,
    upload_vars_h5 as upload_vars_h5_file,
)
from .preprocessing import (
    validate_adata,
    aggregate_expression_percentages,
    extract_marker_genes,
    aggregate_cluster_metadata,
    extract_visualization_coordinates,
)
from .core.payload import build_annotation_payload, save_query_to_file
from .core.artifacts import (
    save_features_matrix,
    save_obs_duckdb as save_obs_duckdb_file,
)
from .core.results import (
    store_job_details,
    store_annotations,
    load_local_results,
    fetch_remote_results,
)

__all__ = ["CyteType"]


def _get_cytetype_version() -> str | None:
    try:
        return version("cytetype")
    except PackageNotFoundError:
        return None


class CyteType:
    """CyteType class for characterizing clusters from single-cell RNA-seq data.

    This class provides an object-oriented interface for cluster characterization using the CyteType API.
    The expensive data preparation steps (validation, expression percentage calculation, and marker
    gene extraction) are performed during initialization, allowing for efficient reuse when making
    multiple requests with different parameters.
    """

    # Type annotations for instance attributes
    adata: anndata.AnnData
    group_key: str
    rank_key: str
    gene_symbols_column: str
    n_top_genes: int
    pcent_batch_size: int
    coordinates_key: str | None
    cluster_map: dict[str, str]
    clusters: list[str]
    expression_percentages: dict[str, dict[str, float]]
    marker_genes: dict[str, list[str]]
    group_metadata: dict[str, dict[str, dict[str, int]]]
    visualization_data: dict[str, Any]
    __version__: str | None = _get_cytetype_version()

    def __init__(
        self,
        adata: anndata.AnnData,
        group_key: str,
        rank_key: str = "rank_genes_groups",
        gene_symbols_column: str = "gene_symbols",
        n_top_genes: int = 50,
        aggregate_metadata: bool = True,
        min_percentage: int = 10,
        pcent_batch_size: int = 2000,
        coordinates_key: str = "X_umap",
        max_cells_per_group: int = 1000,
        api_url: str = "https://prod.cytetype.nygen.io",
        auth_token: str | None = None,
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
            gene_symbols_column (str, optional): Name of the column in `adata.var` that contains
                the gene symbols. Defaults to "gene_symbols".
            n_top_genes (int, optional): Number of top marker genes per cluster to extract during
                initialization. Higher values may improve annotation quality but increase memory usage.
                Defaults to 50.
            aggregate_metadata (bool, optional): Whether to aggregate metadata from the AnnData object.
                Defaults to True.
            min_percentage (int, optional): Minimum percentage of cells in a group to include in the
                cluster context. Defaults to 10.
            pcent_batch_size (int, optional): Batch size for calculating expression percentages to
                optimize memory usage. Defaults to 2000.
            coordinates_key (str, optional): Key in adata.obsm containing 2D coordinates for
                visualization. Must be a 2D array with same number of elements as clusters.
                Defaults to "X_umap".
            max_cells_per_group (int, optional): Maximum number of cells to sample per group
                for visualization. If a group has more cells than this limit, a random sample
                will be taken. Defaults to 1000.
            api_url (str, optional): URL for the CyteType API endpoint. Only change if using a custom
                deployment. Defaults to "https://prod.cytetype.nygen.io".
            auth_token (str | None, optional): Bearer token for API authentication. If provided,
                will be included in the Authorization header as "Bearer {auth_token}". Defaults to None.

        Raises:
            KeyError: If the required keys are missing in `adata.obs` or `adata.uns`
            ValueError: If the data format is incorrect or there are validation errors
        """
        self.adata = adata
        self.group_key = group_key
        self.rank_key = rank_key
        self.gene_symbols_column = gene_symbols_column
        self.n_top_genes = n_top_genes
        self.pcent_batch_size = pcent_batch_size
        self.coordinates_key = coordinates_key
        self.max_cells_per_group = max_cells_per_group
        self.api_url = api_url
        self.auth_token = auth_token

        # Validate data and get the best available coordinates key
        self.coordinates_key = validate_adata(
            adata, group_key, rank_key, gene_symbols_column, coordinates_key
        )

        self.cluster_map = {
            str(x): str(n + 1)
            for n, x in enumerate(natsorted(adata.obs[group_key].unique().tolist()))
        }
        self.clusters = [
            self.cluster_map[str(x)] for x in adata.obs[group_key].values.tolist()
        ]

        logger.info("Calculating expression percentages...")
        self.expression_percentages = aggregate_expression_percentages(
            adata=adata,
            clusters=self.clusters,
            batch_size=pcent_batch_size,
            gene_names=adata.var[self.gene_symbols_column].tolist(),
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
            )
            # Replace keys in group_metadata using cluster_map
            self.group_metadata = {
                self.cluster_map.get(str(key), str(key)): value
                for key, value in self.group_metadata.items()
            }
            self.group_metadata = {
                k: self.group_metadata[k] for k in sorted(self.group_metadata.keys())
            }
        else:
            self.group_metadata = {}

        # Prepare visualization data with sampling
        logger.info("Extracting sampled visualization coordinates...")
        sampled_coordinates, sampled_cluster_labels = extract_visualization_coordinates(
            adata=adata,
            coordinates_key=self.coordinates_key,
            group_key=self.group_key,
            cluster_map=self.cluster_map,
            max_cells_per_group=self.max_cells_per_group,
        )

        self.visualization_data = {
            "coordinates": sampled_coordinates,
            "clusters": sampled_cluster_labels,
        }

        logger.info("Data preparation completed. Ready for submitting jobs.")

    def _build_and_upload_artifacts(
        self,
        vars_h5_path: str,
        obs_duckdb_path: str,
        upload_timeout_seconds: int,
    ) -> dict[str, str]:
        """Build local artifacts and upload them before annotate."""
        logger.info("Saving vars.h5 artifact from normalized counts...")
        save_features_matrix(
            out_file=vars_h5_path,
            mat=self.adata.X,
            var_df=self.adata.var,
            var_names=self.adata.var_names,
        )

        logger.info("Saving obs.duckdb artifact from observation metadata...")
        save_obs_duckdb_file(
            out_file=obs_duckdb_path,
            obs_df=self.adata.obs,
        )

        logger.info("Uploading obs.duckdb artifact...")
        obs_upload = upload_obs_duckdb_file(
            self.api_url,
            self.auth_token,
            obs_duckdb_path,
            timeout=(30.0, float(upload_timeout_seconds)),
        )
        if obs_upload.file_kind != "obs_duckdb":
            raise ValueError(
                f"Unexpected upload file_kind for obs artifact: {obs_upload.file_kind}"
            )

        logger.info("Uploading vars.h5 artifact...")
        vars_upload = upload_vars_h5_file(
            self.api_url,
            self.auth_token,
            vars_h5_path,
            timeout=(30.0, float(upload_timeout_seconds)),
        )
        if vars_upload.file_kind != "vars_h5":
            raise ValueError(
                f"Unexpected upload file_kind for vars artifact: {vars_upload.file_kind}"
            )

        return {
            "obs_duckdb": obs_upload.upload_id,
            "vars_h5": vars_upload.upload_id,
        }

    @staticmethod
    def _cleanup_artifact_files(paths: list[str]) -> None:
        for artifact_path in paths:
            try:
                Path(artifact_path).unlink(missing_ok=True)
            except OSError as exc:
                logger.warning(f"Failed to cleanup artifact {artifact_path}: {exc}")

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
        vars_h5_path: str = "vars.h5",
        obs_duckdb_path: str = "obs.duckdb",
        upload_timeout_seconds: int = 3600,
        cleanup_artifacts: bool = False,
        require_artifacts: bool = True,
        show_progress: bool = True,
        override_existing_results: bool = False,
    ) -> anndata.AnnData:
        """Perform cluster characterization using the CyteType API.

        Args:
            study_context (str, optional): Biological context for the experimental setup.
                For example, include information about 'organisms', 'tissues', 'diseases', 'developmental_stages',
                'single_cell_methods', 'experimental_conditions'. Defaults to None.
            llm_configs (list[dict[str, Any]] | None, optional): Configuration for the large language
                models to be used. Each dict must match the LLMModelConfig schema with required 'provider' and 'name' fields.
                Either 'apiKey' or all AWS credentials ('awsAccessKeyId', 'awsSecretAccessKey', 'awsDefaultRegion') must be provided.
                Defaults to None, using the API's default model.
            metadata (dict[str, Any] | None, optional): Custom metadata tags to include in the report header.
                Values that look like URLs will be made clickable in the report.
                Defaults to None.
            n_parallel_clusters (int, optional): Number of parallel requests to make to the model. Maximum is 50. Note than high values can lead to rate limit errors.
                Defaults to 2.
            results_prefix (str, optional): Prefix for keys added to `adata.obs` and `adata.uns` to
                store results. The final annotation column will be
                `adata.obs[f"{results_key}_{group_key}"]`. Defaults to "cytetype".
            poll_interval_seconds (int, optional): How often (in seconds) to check for results from
                the API. Defaults to 10.
            timeout_seconds (int, optional): Maximum time (in seconds) to wait for API results before
                raising a timeout error. Defaults to 7200.
            api_url (str, optional): URL for the CyteType API endpoint. If not provided, the API URL
                specified in the initialization will be used. If provided, it will override the API URL
                specified in the initialization.
            auth_token (str | None, optional): Bearer token for API authentication. If not provided, the auth token
                specified in the initialization will be used. If provided, it will override the auth token
                specified in the initialization.
            save_query (bool, optional): Whether to save the query to a file. Defaults to True.
            query_filename (str, optional): Filename for saving the query when save_query is True.
                Defaults to "query.json".
            vars_h5_path (str, optional): Local output path for generated vars.h5 artifact.
                Defaults to "vars.h5".
            obs_duckdb_path (str, optional): Local output path for generated obs.duckdb artifact.
                Defaults to "obs.duckdb".
            upload_timeout_seconds (int, optional): Socket read timeout used for each artifact upload.
                Defaults to 3600.
            cleanup_artifacts (bool, optional): Whether to delete generated artifact files after run
                completes or fails. Defaults to False.
            require_artifacts (bool, optional): Whether to raise an error if artifact building or
                uploading fails. When True (default), any artifact failure stops the run. Set to
                False to skip artifacts and continue with annotation only. Defaults to True.
            show_progress (bool, optional): Whether to display progress updates with spinner and
                cluster status. Set to False to disable all visual progress output. Defaults to True.
            override_existing_results (bool, optional): Whether to override existing results with the
                same results_prefix. If False (default) and results already exist, will raise an error
                to prevent accidental overwriting. Set to True to explicitly allow overwriting.
                Defaults to False.

        Returns:
            anndata.AnnData: The input AnnData object, modified in place with the following additions:
                - `adata.obs[f"{results_prefix}_{group_key}"]`: Cell type annotations as categorical values
                - `adata.uns[f"{results_prefix}_results"]`: Complete API response data and job tracking info

        Raises:
            ValueError: If results with the same prefix already exist and override_existing_results is False
            CyteTypeAPIError: If the API request fails or returns invalid data
            CyteTypeTimeoutError: If the API does not return results within the specified timeout period

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

        if api_url:
            self.api_url = api_url.strip("/")
        if auth_token:
            self.auth_token = auth_token
        if upload_timeout_seconds <= 0:
            raise ValueError("upload_timeout_seconds must be greater than 0")

        # Validate the base payload before doing potentially expensive uploads.
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
        )

        artifact_paths = [vars_h5_path, obs_duckdb_path]
        try:
            try:
                uploaded_file_refs = self._build_and_upload_artifacts(
                    vars_h5_path=vars_h5_path,
                    obs_duckdb_path=obs_duckdb_path,
                    upload_timeout_seconds=upload_timeout_seconds,
                )
                payload["uploaded_files"] = uploaded_file_refs
            except Exception as exc:
                if require_artifacts:
                    logger.error(
                        "Artifact build/upload failed. "
                        "Rerun with `require_artifacts=False` to skip this error.\n"
                        "Please report the error below in a new issue at "
                        "https://github.com/NygenAnalytics/CyteType\n"
                        f"({type(exc).__name__}: {exc})"
                    )
                    raise
                logger.warning(
                    "Artifact build/upload failed. Continuing without artifacts. "
                    "Set `require_artifacts=True` to see the full traceback."
                )

            # Save query if requested
            if save_query:
                save_query_to_file(payload["input_data"], query_filename)

            # Submit job and store details
            job_id = submit_annotation_job(self.api_url, self.auth_token, payload)
            store_job_details(self.adata, job_id, self.api_url, results_prefix)

            # Wait for completion
            result = wait_for_completion(
                self.api_url,
                self.auth_token,
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

            return self.adata
        finally:
            if cleanup_artifacts:
                self._cleanup_artifact_files(artifact_paths)

    def get_results(
        self,
        results_prefix: str = "cytetype",
    ) -> dict[str, Any] | None:
        """Retrieve the CyteType results from the AnnData object.

        If results are not available locally but job details exist, attempts to retrieve
        results from the API with a single request (no polling).

        Args:
            results_prefix: The prefix used when storing results. Defaults to "cytetype".

        Returns:
            dict[str, Any] | None: The original result dictionary from the API, or None if not found.
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

        return fetch_remote_results(
            self.adata,
            job_id,
            self.api_url,
            self.auth_token,
            results_prefix,
            self.group_key,
            self.clusters,
        )
