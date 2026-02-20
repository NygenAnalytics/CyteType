import json
import anndata
import pandas as pd
from typing import Any, Callable

from ..config import logger
from ..api import get_job_status, fetch_job_results


def store_job_details(
    adata: anndata.AnnData,
    job_id: str,
    api_url: str,
    results_prefix: str,
) -> None:
    """Store job details in adata.uns for later retrieval.

    Args:
        adata: AnnData object to store details in
        job_id: Job ID from API
        api_url: API base URL
        results_prefix: Prefix for storage keys
    """
    report_url = f"{api_url}/report/{job_id}"
    adata.uns[f"{results_prefix}_jobDetails"] = {
        "job_id": job_id,
        "report_url": report_url,
        "api_url": api_url,
    }


def store_annotations(
    adata: anndata.AnnData,
    result_data: dict[str, Any],
    job_id: str,
    results_prefix: str,
    group_key: str,
    clusters: list[str],
    check_unannotated: bool = True,
) -> None:
    """Store API results and update annotations in the AnnData object.

    Args:
        adata: AnnData object to update
        result_data: Result dictionary from API
        job_id: Job ID
        results_prefix: Prefix for storage keys
        group_key: Column name in adata.obs for grouping
        clusters: List of cluster assignments
        check_unannotated: Whether to check and warn about unannotated clusters
    """
    # Store raw results
    _store_raw_results(adata, result_data, job_id, results_prefix)

    # Update annotation columns using DRY helper
    _add_annotation_column(
        adata,
        result_data,
        "annotation",
        f"{results_prefix}_annotation_{group_key}",
        clusters,
        "Unknown",
    )
    _add_annotation_column(
        adata,
        result_data,
        "ontologyTerm",
        f"{results_prefix}_cellOntologyTerm_{group_key}",
        clusters,
        "Unknown",
    )
    _add_annotation_column(
        adata,
        result_data,
        "ontologyTermID",
        f"{results_prefix}_cellOntologyTermID_{group_key}",
        clusters,
        "Unknown",
    )
    _add_annotation_column(
        adata,
        result_data,
        "cellState",
        f"{results_prefix}_cellState_{group_key}",
        clusters,
        "",
        field_getter=lambda item: item.get("cellState", ""),
    )

    # Check for unannotated clusters
    if check_unannotated:
        _check_unannotated_clusters(result_data, clusters)

    # Log success
    logger.info(
        f"     Annotation labels → adata.obs['{results_prefix}_annotation_{group_key}']\n"
        f"   Cell Ontology terms → adata.obs['{results_prefix}_cellOntologyTerm_{group_key}']\n"
        f"Cell Ontology term IDs → adata.obs['{results_prefix}_ontologyTermID_{group_key}']\n"
        f"           Cell states → adata.obs['{results_prefix}_cellState_{group_key}']\n"
        f"          Full results → adata.uns['{results_prefix}_results']"
    )


def load_local_results(
    adata: anndata.AnnData,
    results_prefix: str,
) -> dict[str, Any] | None:
    """Load results from adata.uns if they exist.

    Args:
        adata: AnnData object to load results from
        results_prefix: Prefix used when storing results

    Returns:
        Result dictionary if found, None otherwise
    """
    results_key = f"{results_prefix}_results"

    if results_key not in adata.uns:
        return None

    stored_results = adata.uns[results_key]
    if "result" not in stored_results:
        return None

    # Try JSON string format (new format for HDF5 compatibility)
    try:
        result = json.loads(stored_results["result"])
        if isinstance(result, dict):
            return result
        else:
            logger.warning(f"Expected dict from stored result, got {type(result)}")
    except (json.JSONDecodeError, TypeError):
        # Fallback to dict format (backwards compatibility)
        result = stored_results["result"]
        if isinstance(result, dict):
            return result
        else:
            logger.warning(
                f"Expected dict from stored result fallback, got {type(result)}"
            )

    return None


def fetch_remote_results(
    adata: anndata.AnnData,
    job_id: str,
    api_url: str,
    auth_token: str | None,
    results_prefix: str,
    group_key: str,
    clusters: list[str],
) -> dict[str, Any] | None:
    """Fetch results from API for a specific job.

    Args:
        adata: AnnData object to store results in
        job_id: Job ID to fetch results for
        api_url: API base URL
        auth_token: Authentication token for API
        results_prefix: Prefix for storage keys
        group_key: Column name in adata.obs for grouping
        clusters: List of cluster assignments

    Returns:
        Result dictionary if fetch successful, None otherwise
    """
    logger.info(
        f"No results found locally. Attempting to retrieve results for job_id: {job_id}"
    )

    report_url = f"{api_url}/report/{job_id}"

    try:
        api_url = api_url.rstrip("/")
        status_response = get_job_status(api_url, auth_token, job_id)
        status = status_response["jobStatus"]

        if status == "completed":
            result_data = fetch_job_results(api_url, auth_token, job_id)

            if not isinstance(result_data, dict):
                logger.error(f"Expected dict result from API, got {type(result_data)}")
                return None

            # Store retrieved results locally
            _store_raw_results(adata, result_data, job_id, results_prefix)

            # Update annotations
            store_annotations(
                adata,
                result_data,
                job_id,
                results_prefix,
                group_key,
                clusters,
                check_unannotated=False,
            )

            logger.info(f"Interactive report available at: {report_url}")

            return result_data

        elif status == "failed":
            logger.error(f"Job {job_id} failed")
            return None

        elif status in ["processing", "pending"]:
            logger.info(f"Job {job_id} is still {status}. Results not yet available.")
            logger.info(f"Interactive report available at: {report_url}")
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


def _store_raw_results(
    adata: anndata.AnnData,
    result_data: dict[str, Any],
    job_id: str,
    results_prefix: str,
) -> None:
    """Store raw results in adata.uns (excluding large fields)."""
    filtered_result_data = {
        k: v
        for k, v in result_data.items()
        if k not in ["marker_genes", "visualization_data"]
    }
    adata.uns[f"{results_prefix}_results"] = {
        "job_id": job_id,
        "result": json.dumps(
            filtered_result_data
        ),  # JSON string for HDF5 compatibility
    }


def _add_annotation_column(
    adata: anndata.AnnData,
    result_data: dict[str, Any],
    field_key: str,
    column_name: str,
    clusters: list[str],
    default_value: str,
    field_getter: Callable[[dict[str, Any]], Any] | None = None,
) -> None:
    """Add a single annotation column to adata.obs (DRY helper).

    Args:
        adata: AnnData object to update
        result_data: Result dictionary containing annotations
        field_key: Key in annotation dict to extract (e.g., "annotation", "ontologyTerm")
        column_name: Name for the new column in adata.obs
        clusters: List of cluster assignments
        default_value: Default value for missing annotations
        field_getter: Optional custom function to extract field value
    """
    if field_getter is None:

        def field_getter(item: dict[str, Any]) -> Any:
            return item.get(field_key, default_value)

    field_map = {
        item["clusterId"]: field_getter(item)
        for item in result_data.get("annotations", [])
    }

    adata.obs[column_name] = pd.Series(
        [field_map.get(cluster_id, default_value) for cluster_id in clusters],
        index=adata.obs.index,
    ).astype("category")


def _check_unannotated_clusters(
    result_data: dict[str, Any],
    clusters: list[str],
) -> None:
    """Check for and warn about unannotated clusters."""
    annotation_map = {
        item["clusterId"]: item["annotation"]
        for item in result_data.get("annotations", [])
    }

    unannotated_clusters = set(
        cluster_id for cluster_id in clusters if cluster_id not in annotation_map
    )

    if unannotated_clusters:
        logger.warning(
            f"No annotations received from API for cluster IDs: {unannotated_clusters}. "
            f"Corresponding cells marked as 'Unknown Annotation'."
        )
