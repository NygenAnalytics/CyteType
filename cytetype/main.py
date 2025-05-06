import anndata
import numpy as np
import pandas as pd
from typing import Any

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
    gene_symbols_column_name: str = "gene_symbols",
    pcent_batch_size: int = 2000,
    poll_interval_seconds: int = DEFAULT_POLL_INTERVAL,
    timeout_seconds: int = DEFAULT_TIMEOUT,
    api_url: str = DEFAULT_API_URL,
    model_config: list[dict[str, Any]] | None = None,
) -> anndata.AnnData:
    """Annotate cell types in an AnnData object using the CyteType API service.

    This function takes a pre-processed AnnData object with clustered cells and differential expression
    results, submits the data to the CyteType API, and integrates the resulting cell type annotations
    back into the AnnData object.

    Requirements:
        - AnnData object with log1p-normalized data in `adata.X`
        - Cell clusters/groups in `adata.obs[cell_group_key]`
        - Results from `sc.tl.rank_genes_groups` in `adata.uns[rank_genes_key]`
        - Internet connection to reach the CyteType API service

    Args:
        adata (anndata.AnnData): The AnnData object to annotate. Must contain log1p-normalized
            gene expression data in `adata.X` and gene names in `adata.var_names`.
        cell_group_key (str): The key in `adata.obs` containing the preliminary cell cluster labels.
            These clusters will receive cell type annotations.
        rank_genes_key (str, optional): The key in `adata.uns` containing differential expression
            results from `sc.tl.rank_genes_groups`. Must use the same `groupby` as `cell_group_key`.
            Defaults to "rank_genes_groups".
        results_key_added (str, optional): Prefix for keys added to `adata.obs` and `adata.uns` to
            store results. The final annotation column will be
            `adata.obs[f"{results_key_added}_{cell_group_key}"]`. Defaults to "CyteType".
        organism (str, optional): The scientific name of the organism being analyzed.
            Defaults to "Homo sapiens".
        tissues (list[str] | None, optional): List of tissues/organs relevant to the dataset.
            Helps improve annotation specificity. Defaults to None.
        diseases (list[str] | None, optional): List of diseases/conditions relevant to the dataset.
            Helps improve annotation specificity. Defaults to None.
        developmental_stages (list[str] | None, optional): List of developmental stages relevant
            to the dataset. Helps improve annotation specificity. Defaults to None.
        single_cell_methods (list[str] | None, optional): List of single-cell methods used to
            generate the data. Defaults to None.
        experimental_conditions (list[str] | None, optional): List of experimental conditions or
            treatments applied to the cells. Defaults to None.
        n_top_genes (int, optional): Number of top marker genes per cluster to send to the API.
            Higher values may improve annotation quality but increase API request size.
            Defaults to 50.
        gene_symbols_column_name (str, optional): Name of the column in `adata.var` that contains
            the gene symbols. Defaults to "gene_symbols".
        pcent_batch_size (int, optional): Batch size for calculating expression percentages to
            optimize memory usage. Defaults to 2000.
        poll_interval_seconds (int, optional): How often (in seconds) to check for results from
            the API. Defaults to DEFAULT_POLL_INTERVAL.
        timeout_seconds (int, optional): Maximum time (in seconds) to wait for API results before
            raising a timeout error. Defaults to DEFAULT_TIMEOUT.
        api_url (str, optional): URL for the CyteType API endpoint. Only change if using a custom
            deployment. Defaults to DEFAULT_API_URL.
        model_config (list[dict[str, Any]] | None, optional): Configuration for the large language
            models used for annotation. Allows specifying provider, model name, API key, and base URL.
            Defaults to None, using the API's default model.

    Returns:
        anndata.AnnData: The input AnnData object, modified in place with the following additions:
            - `adata.obs[f"{results_key_added}_{cell_group_key}"]`: Cell type annotations as categorical values
            - `adata.uns[f"{results_key_added}_results"]`: Complete API response data and job tracking info

    Raises:
        KeyError: If the required keys are missing in `adata.obs` or `adata.uns`
        ValueError: If the data format is incorrect or there are validation errors
        CyteTypeAPIError: If the API request fails or returns invalid data
        CyteTypeTimeoutError: If the API does not return results within the specified timeout period

    Examples:
        >>> import scanpy as sc
        >>> import cytetype
        >>>
        >>> # Load and preprocess data (example with standard scanpy workflow)
        >>> adata = sc.read_h5ad("my_dataset.h5ad")
        >>> sc.pp.normalize_total(adata, target_sum=1e4)
        >>> sc.pp.log1p(adata)
        >>> sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        >>> sc.pp.pca(adata, use_highly_variable=True)
        >>> sc.pp.neighbors(adata)
        >>> sc.tl.leiden(adata, resolution=0.8)
        >>> sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
        >>>
        >>> # Annotate cell types with CyteType
        >>> adata = cytetype.annotate_anndata(
        >>>     adata,
        >>>     cell_group_key='leiden',
        >>>     tissues=["Brain", "Nervous system"],
        >>>     n_top_genes=100
        >>> )
        >>>
        >>> # Access annotations
        >>> sc.pl.umap(adata, color='CyteType_leiden', legend_loc='on data')
    """
    job_id = None

    _validate_adata(adata, cell_group_key, rank_genes_key, gene_symbols_column_name)

    ct_map = {
        str(x): str(n + 1)
        for n, x in enumerate(sorted(adata.obs[cell_group_key].unique().tolist()))
    }
    clusters = [ct_map[str(x)] for x in adata.obs[cell_group_key].values.tolist()]

    logger.info("Calculating expression percentages.")
    pcent = _calculate_pcent(
        adata=adata,
        clusters=clusters,
        batch_size=pcent_batch_size,
        gene_names=adata.var_names,
    )
    logger.info("Extracting marker genes.")
    markers = _get_markers(
        adata=adata,
        cell_group_key=cell_group_key,
        rank_genes_key=rank_genes_key,
        ct_map=ct_map,
        n_top_genes=n_top_genes,
        gene_symbols_col=gene_symbols_column_name,
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

    job_id = submit_annotation_job(query, api_url, model_config=model_config)
    logger.info(f"Waiting for results for job ID: {job_id}")
    annotation_results = poll_for_results(
        job_id, api_url, poll_interval_seconds, timeout_seconds
    )

    adata.uns[f"{results_key_added}_results"] = {
        "job_id": job_id,
        "result": annotation_results,
    }

    anno_map = {
        i["clusterId"]: i["annotation"]
        for i in annotation_results.get("annotations", [])
    }
    adata.obs[f"{results_key_added}_{cell_group_key}"] = pd.Series(
        [anno_map.get(x, "Unknown Annotation") for x in clusters],
        index=adata.obs.index,
    ).astype("category")

    unannotated_clusters = {
        x for x in clusters if x not in anno_map
    }
    if unannotated_clusters:
        logger.warning(
            f"No annotations received from API for cluster IDs: {unannotated_clusters}. Corresponding cells marked as 'Unknown Annotation'."
        )

    logger.info(
        f"Annotations successfully added to `adata.obs['{results_key_added}_{cell_group_key}']` and `adata.uns['{results_key_added}_results']`."
    )

    return adata
