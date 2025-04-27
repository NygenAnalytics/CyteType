import anndata
import numpy as np
import pandas as pd

from .config import logger, DEFAULT_API_URL, DEFAULT_POLL_INTERVAL, DEFAULT_TIMEOUT
from .client import submit_annotation_job, poll_for_results
from .anndata_helpers import _validate_adata, _calculate_pcent, _get_markers


def annotate_anndata(
    adata: anndata.AnnData,
    cell_group_key: str,
    rank_genes_key: str = "rank_genes_groups",
    results_key_added: str = "CyteType",
    organism: str = "Homo sapiens",
    tissues: list[str] | None = None,
    diseases: list[str] | None = None,
    developmental_stages: list[str] | None = None,
    single_cell_methods: list[str] | None = None,
    experimental_conditions: list[str] | None = None,
    n_top_genes: int = 50,
    pcent_batch_size: int = 2000,
    poll_interval_seconds: int = DEFAULT_POLL_INTERVAL,
    timeout_seconds: int = DEFAULT_TIMEOUT,
    api_url: str = DEFAULT_API_URL,
) -> anndata.AnnData:
    """
    Annotates cell types in an AnnData object using the CyteType API.

    Args:
        adata: The AnnData object to annotate. Requires raw counts in `adata.raw` and
               log1p normalized data in `adata.X`. Gene names should be in `adata.var_names`.
        cell_group_key: The key in `adata.obs` containing the preliminary cell clusters/groups.
        results_key_added: Prefix for keys added to `adata.obs` and `adata.uns`.
                           Defaults to 'CyteType'.
        organism: The organism name (e.g., "Homo sapiens").
        tissues: List of tissues relevant to the dataset.
        diseases: List of diseases relevant to the dataset.
        developmental_stages: List of developmental stages.
        single_cell_methods: List of single-cell methods used.
        experimental_conditions: List of experimental conditions.
        n_top_genes: Number of top marker genes per cluster to send to the API.
        pcent_batch_size: Batch size for calculating expression percentages.
        poll_interval_seconds: How often (in seconds) to check for results from the API.
                               Defaults to `cytetype.config.DEFAULT_POLL_INTERVAL`.
        timeout_seconds: Maximum time (in seconds) to wait for API results.
                         Defaults to `cytetype.config.DEFAULT_TIMEOUT`.
        api_url: Optional URL for the CyteType API endpoint. If None, uses the default
                 from `cytetype.config.DEFAULT_API_URL`.

    Returns:
        The input AnnData object, modified in place with annotation results.
    """
    job_id = None

    _validate_adata(adata, cell_group_key, rank_genes_key)

    ct_map = {
        str(x): n + 1 for n, x in enumerate(sorted(adata.obs[cell_group_key].unique()))
    }
    clusters_int = np.array([ct_map[str(x)] for x in adata.obs[cell_group_key]])

    logger.info("Calculating expression percentages.")
    pcent = _calculate_pcent(adata, clusters_int, pcent_batch_size)
    logger.info("Extracting marker genes.")
    markers = _get_markers(
        adata=adata,
        cell_group_key=cell_group_key,
        rank_genes_key=rank_genes_key,
        ct_map=ct_map,
        n_top_genes=n_top_genes,
    )

    bio_context = {
        "organisms": [organism] if organism else ["Unknown"],
        "tissues": tissues if tissues is not None else ["Unknown"],
        "diseases": diseases if diseases is not None else ["Unknown"],
        "developmentalStages": developmental_stages
        if developmental_stages is not None
        else ["Unknown"],
        "singleCellMethods": single_cell_methods
        if single_cell_methods is not None
        else ["Unknown"],
        "experimentalConditions": experimental_conditions
        if experimental_conditions is not None
        else ["Unknown"],
    }

    query = {
        "bioContext": bio_context,
        "markerGenes": markers,
        "expressionData": pcent,
    }

    logger.info("Submitting job to CyteType API.")
    job_id = submit_annotation_job(query, api_url)
    logger.info(f"Polling for results for job ID: {job_id}")
    annotation_results = poll_for_results(
        job_id, api_url, poll_interval_seconds, timeout_seconds
    )

    adata.uns[f"{results_key_added}_{cell_group_key}"] = {
        "job_id": job_id,
        "result": annotation_results,
    }

    anno_map = {
        int(i["clusterId"]): i["annotation"]
        for i in annotation_results.get("annotations", [])
    }
    adata.obs[f"{results_key_added}_{cell_group_key}"] = pd.Series(
        [anno_map.get(x, "Unknown Annotation") for x in clusters_int],
        index=adata.obs.index,
    ).astype("category")

    unannotated_clusters = {
        clust_int for clust_int in clusters_int if clust_int not in anno_map
    }
    if unannotated_clusters:
        logger.warning(
            f"No annotations received from API for cluster IDs: {unannotated_clusters}. Corresponding cells marked as 'Unknown Annotation'."
        )

    logger.info(
        f"Annotations successfully added to `adata.obs['{results_key_added}_{cell_group_key}']`."
    )

    return adata
