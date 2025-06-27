from typing import Any
import json

import anndata
import pandas as pd
from natsort import natsorted
from pydantic import BaseModel, Field, ConfigDict

from .config import logger, DEFAULT_API_URL, DEFAULT_POLL_INTERVAL, DEFAULT_TIMEOUT
from .client import submit_job, poll_for_results, check_job_status
from .anndata_helpers import (
    _validate_adata,
    _calculate_pcent,
    _get_markers,
    _aggregate_metadata,
    _extract_sampled_coordinates,
)

__all__ = ["CyteType", "BioContext", "ModelConfig", "RunConfig"]


class BioContext(BaseModel):
    """Biological context information for the data."""

    model_config = ConfigDict(populate_by_name=True)

    studyContext: str = Field(default="")
    clusterContext: dict[str, dict[str, dict[str, int]]] = Field(default_factory=dict)


class ModelConfig(BaseModel):
    """Configuration for the large language model to be used."""

    model_config = ConfigDict(populate_by_name=True)

    provider: str
    name: str | None = None
    apiKey: str | None = Field(default=None)
    baseUrl: str | None = Field(default=None)
    modelSettings: dict[str, Any] | None = Field(default=None)


class RunConfig(BaseModel):
    """Configuration for the annotation run."""

    model_config = ConfigDict(populate_by_name=True)

    concurrentClusters: int | None = Field(default=None)
    maxAnnotationRevisions: int | None = Field(default=None)
    maxLLMRequests: int | None = Field(default=None)


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

        # Validate data and get the best available coordinates key
        validated_coordinates_key = _validate_adata(
            adata, group_key, rank_key, gene_symbols_column, coordinates_key
        )
        self.coordinates_key = validated_coordinates_key

        self.cluster_map = {
            str(x): str(n + 1)
            for n, x in enumerate(natsorted(adata.obs[group_key].unique().tolist()))
        }
        self.clusters = [
            self.cluster_map[str(x)] for x in adata.obs[group_key].values.tolist()
        ]

        logger.info("Calculating expression percentages.")
        self.expression_percentages = _calculate_pcent(
            adata=adata,
            clusters=self.clusters,
            batch_size=pcent_batch_size,
            gene_names=adata.var[self.gene_symbols_column].tolist(),
        )

        logger.info("Extracting marker genes.")
        self.marker_genes = _get_markers(
            adata=self.adata,
            cell_group_key=self.group_key,
            rank_genes_key=self.rank_key,
            ct_map=self.cluster_map,
            n_top_genes=n_top_genes,
            gene_symbols_col=self.gene_symbols_column,
        )

        if aggregate_metadata:
            self.group_metadata = _aggregate_metadata(
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
        logger.info("Preparing visualization data with sampling.")
        sampled_coordinates, sampled_cluster_labels = _extract_sampled_coordinates(
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

    def _store_results_and_annotations(
        self,
        result_data: dict[str, Any],
        job_id: str,
        results_prefix: str,
        check_unannotated: bool = True,
    ) -> None:
        """Store API results and update annotations in the AnnData object.

        Args:
            result_data: The result dictionary from the API
            job_id: The job ID
            results_prefix: Prefix for storing results
            check_unannotated: Whether to check and warn about unannotated clusters
        """
        # Store results in AnnData object (excluding marker_genes and visualization_data)
        filtered_result_data = {
            k: v
            for k, v in result_data.items()
            if k not in ["marker_genes", "visualization_data"]
        }
        self.adata.uns[f"{results_prefix}_results"] = {
            "job_id": job_id,
            "result": json.dumps(
                filtered_result_data
            ),  # Convert to JSON string for HDF5 compatibility
        }

        # Update annotations
        annotation_map = {
            item["clusterId"]: item["annotation"]
            for item in result_data.get("annotations", [])
        }
        self.adata.obs[f"{results_prefix}_annotation_{self.group_key}"] = pd.Series(
            [annotation_map.get(cluster_id, "Unknown") for cluster_id in self.clusters],
            index=self.adata.obs.index,
        ).astype("category")

        # Update ontology terms
        ontology_map = {
            item["clusterId"]: item["ontologyTerm"]
            for item in result_data.get("annotations", [])
        }
        self.adata.obs[f"{results_prefix}_cellOntologyTerm_{self.group_key}"] = (
            pd.Series(
                [
                    ontology_map.get(cluster_id, "Unknown")
                    for cluster_id in self.clusters
                ],
                index=self.adata.obs.index,
            ).astype("category")
        )

        # Check for unannotated clusters if requested
        if check_unannotated:
            unannotated_clusters = set(
                [
                    cluster_id
                    for cluster_id in self.clusters
                    if cluster_id not in annotation_map
                ]
            )

            if unannotated_clusters:
                logger.warning(
                    f"No annotations received from API for cluster IDs: {unannotated_clusters}. "
                    f"Corresponding cells marked as 'Unknown Annotation'."
                )

        # Log success message
        logger.success(
            f"Annotations successfully added to `adata.obs['{results_prefix}_annotation_{self.group_key}']` "
            f", ontology term added to `adata.obs['{results_prefix}_cellOntologyTerm_{self.group_key}']` "
            f"and, full results added to `adata.uns['{results_prefix}_results']`."
        )

    def run(
        self,
        study_context: str,
        model_config: list[dict[str, Any]] | None = None,
        run_config: dict[str, Any] | None = None,
        metadata: dict[str, Any] | None = None,
        results_prefix: str = "cytetype",
        poll_interval_seconds: int = DEFAULT_POLL_INTERVAL,
        timeout_seconds: int = DEFAULT_TIMEOUT,
        api_url: str = DEFAULT_API_URL,
        save_query: bool = True,
        query_filename: str = "query.json",
        auth_token: str | None = None,
    ) -> anndata.AnnData:
        """Perform cluster characterization using the CyteType API.

        Args:
            study_context (str, optional): Biological context for the experimental setup.
                For example, include information about 'organisms', 'tissues', 'diseases', 'developmental_stages',
                'single_cell_methods', 'experimental_conditions'. Defaults to None.
            model_config (list[dict[str, Any]] | None, optional): Configuration for the large language
                models to be used. Each dict must include 'provider', 'name', 'apiKey', 'baseUrl' (optional), 'modelSettings' (optional).
                Defaults to None, using the API's default model.
            run_config (dict[str, Any] | None, optional): Configuration for the annotation run.
                Can include 'maxAnnotationRevisions'.
                Defaults to None, using the API's default settings.
            metadata (dict[str, Any] | None, optional): Custom metadata to send with the API request.
                Can include information like run labels, experiment names, or any other user-defined data.
                This metadata will be sent to the API for tracking purposes but not stored locally.
                Defaults to None.
            results_prefix (str, optional): Prefix for keys added to `adata.obs` and `adata.uns` to
                store results. The final annotation column will be
                `adata.obs[f"{results_key}_{group_key}"]`. Defaults to "cytetype".
            poll_interval_seconds (int, optional): How often (in seconds) to check for results from
                the API. Defaults to DEFAULT_POLL_INTERVAL.
            timeout_seconds (int, optional): Maximum time (in seconds) to wait for API results before
                raising a timeout error. Defaults to DEFAULT_TIMEOUT.
            api_url (str, optional): URL for the CyteType API endpoint. Only change if using a custom
                deployment. Defaults to DEFAULT_API_URL.
            save_query (bool, optional): Whether to save the query to a file. Defaults to True.
            query_filename (str, optional): Filename for saving the query when save_query is True.
                Defaults to "query.json".
            auth_token (str | None, optional): Bearer token for API authentication. If provided,
                will be included in the Authorization header as "Bearer {auth_token}". Defaults to None.

        Returns:
            anndata.AnnData: The input AnnData object, modified in place with the following additions:
                - `adata.obs[f"{results_prefix}_{group_key}"]`: Cell type annotations as categorical values
                - `adata.uns[f"{results_prefix}_results"]`: Complete API response data and job tracking info

        Raises:
            CyteTypeAPIError: If the API request fails or returns invalid data
            CyteTypeTimeoutError: If the API does not return results within the specified timeout period

        """
        api_url = api_url.rstrip("/")

        bio_context = BioContext(
            studyContext=study_context, clusterContext=self.group_metadata
        ).model_dump()

        # Process model config using Pydantic model
        if model_config is not None:
            model_config_list = [
                ModelConfig(**config).model_dump() for config in model_config
            ]
        else:
            model_config_list = []

        # Process run config using Pydantic model
        if run_config is not None:
            run_config_dict = RunConfig(**run_config).model_dump()
            # Remove None values
            run_config_dict = {
                k: v for k, v in run_config_dict.items() if v is not None
            }
        else:
            run_config_dict = {}

        # Prepare API query
        query: dict[str, Any] = {
            "bioContext": bio_context,
            "markerGenes": self.marker_genes,
            "expressionData": self.expression_percentages,
            "modelConfig": model_config_list,
            "runConfig": run_config_dict,
            "visualizationData": self.visualization_data,
            "clusterLabels": {v: k for k, v in self.cluster_map.items()},
        }

        # Add metadata if provided
        if metadata is not None:
            query["metadata"] = metadata

        if save_query:
            with open(query_filename, "w") as f:
                json.dump(query, f)

        # Submit job and poll for results
        job_id = submit_job(
            query,
            api_url,
            auth_token=auth_token,
        )

        # Store job details (job_id, report link) for potential later retrieval
        report_url = f"{api_url}/report/{job_id}"

        self.adata.uns[f"{results_prefix}_jobDetails"] = {
            "job_id": job_id,
            "report_url": report_url,
            "api_url": api_url,
            "auth_token": auth_token if auth_token else None,
        }

        logger.info(f"Waiting for results for job ID: {job_id}")

        # Log the report URL that updates automatically
        logger.info(
            f"View the automatically updating visualization report at: {report_url}"
        )

        result = poll_for_results(
            job_id,
            api_url,
            poll_interval_seconds,
            timeout_seconds,
            auth_token=auth_token,
        )

        self._store_results_and_annotations(
            result,
            job_id,
            results_prefix,
            check_unannotated=True,
        )

        return self.adata

    def get_results(self, results_prefix: str = "cytetype") -> dict[str, Any] | None:
        """Retrieve the CyteType results from the AnnData object.

        If results are not available locally but job details exist, attempts to retrieve
        results from the API with a single request (no polling).

        Args:
            results_prefix (str): The prefix used when storing results. Defaults to "cytetype".

        Returns:
            dict[str, Any] | None: The original result dictionary from the API, or None if not found.
        """
        results_key = f"{results_prefix}_results"

        # First check if results already exist locally
        if results_key in self.adata.uns:
            stored_results = self.adata.uns[results_key]
            if "result" in stored_results:
                # The result is stored as a JSON string for HDF5 compatibility
                try:
                    result = json.loads(stored_results["result"])
                    if isinstance(result, dict):
                        return result
                    else:
                        logger.warning(
                            f"Expected dict from stored result, got {type(result)}"
                        )
                except (json.JSONDecodeError, TypeError):
                    # Fallback for cases where result might still be a dict (backwards compatibility)
                    result = stored_results["result"]
                    if isinstance(result, dict):
                        return result
                    else:
                        logger.warning(
                            f"Expected dict from stored result fallback, got {type(result)}"
                        )

        # If no results found locally, try to retrieve using stored job details
        job_details_key = f"{results_prefix}_jobDetails"
        if job_details_key not in self.adata.uns:
            logger.info(
                "No results found locally and no job details available for retrieval."
            )
            return None

        job_details = self.adata.uns[job_details_key]
        job_id = job_details.get("job_id")
        api_url = job_details.get("api_url", DEFAULT_API_URL)
        auth_token = job_details.get("auth_token")

        if not job_id:
            logger.error("Job details found but missing job_id.")
            return None

        logger.info(
            f"No results found locally. Attempting to retrieve results for job_id: {job_id}"
        )

        try:
            # Use the client function to check job status
            api_url = api_url.rstrip("/")
            status_response = check_job_status(job_id, api_url, auth_token)

            status = status_response["status"]

            if status == "completed":
                logger.info(
                    f"Job {job_id} completed successfully. Storing results locally."
                )
                result_data = status_response["result"]

                # Ensure we have a proper dict result instead of Any
                if not isinstance(result_data, dict):
                    logger.error(
                        f"Expected dict result from API, got {type(result_data)}"
                    )
                    return None

                # Store the retrieved results locally
                self.adata.uns[results_key] = {
                    "job_id": job_id,
                    "result": json.dumps(result_data),
                }

                # Update annotations
                self._store_results_and_annotations(
                    result_data,
                    job_id,
                    results_prefix,
                    check_unannotated=False,
                )

                return result_data

            elif status == "error":
                logger.error(f"Job {job_id} failed: {status_response['message']}")
                return None

            elif status in ["processing", "pending"]:
                logger.info(
                    f"Job {job_id} is still {status}. Results not yet available."
                )
                logger.info(
                    f"Check progress at: {job_details.get('report_url', 'N/A')}"
                )
                return None

            elif status == "not_found":
                logger.info(f"Job {job_id} results not yet available (404).")
                return None

            else:
                logger.warning(f"Job {job_id} has unknown status: '{status}'.")
                return None

        except Exception as e:
            logger.error(f"Failed to retrieve results for job {job_id}: {e}")
            return None
