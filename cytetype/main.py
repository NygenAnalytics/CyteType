from typing import Any
import json

import anndata
import pandas as pd
from pydantic import BaseModel, Field, ConfigDict

from .config import logger, DEFAULT_API_URL, DEFAULT_POLL_INTERVAL, DEFAULT_TIMEOUT
from .client import submit_annotation_job, poll_for_results
from .anndata_helpers import _validate_adata, _calculate_pcent, _get_markers


class BioContext(BaseModel):
    """Biological context information for cell type annotation."""
    model_config = ConfigDict(populate_by_name=True)
    
    organisms: list[str] = Field(default=["Unknown"])
    tissues: list[str] = Field(default=["Unknown"])
    diseases: list[str] = Field(default=["Unknown"])
    developmental_stages: list[str] = Field(default=["Unknown"], alias="developmentalStages")
    single_cell_methods: list[str] = Field(default=["Unknown"], alias="singleCellMethods")
    experimental_conditions: list[str] = Field(default=["Unknown"], alias="experimentalConditions")


class ModelConfig(BaseModel):
    """Configuration for the large language model used for annotation."""
    model_config = ConfigDict(populate_by_name=True)
    
    provider: str
    name: str | None = None
    api_key: str | None = Field(default=None, alias="apiKey")
    base_url: str | None = Field(default=None, alias="baseUrl")


class CyteType:
    """CyteType cell type annotation class.
    
    This class provides an object-oriented interface for cell type annotation using the CyteType API.
    The expensive data preparation steps (validation, expression percentage calculation, and marker
    gene extraction) are performed during initialization, allowing for efficient reuse when making
    multiple annotation requests with different parameters.
    """

    def __init__(
        self,
        adata: anndata.AnnData,
        cell_group_key: str,
        rank_genes_key: str = "rank_genes_groups",
        results_key_added: str = "CyteType",
        gene_symbols_column_name: str = "gene_symbols",
        n_top_genes: int = 50,
        pcent_batch_size: int = 2000,
    ):
        """Initialize CyteType with AnnData object and perform data preparation.

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
            gene_symbols_column_name (str, optional): Name of the column in `adata.var` that contains
                the gene symbols. Defaults to "gene_symbols".
            n_top_genes (int, optional): Number of top marker genes per cluster to extract during
                initialization. Higher values may improve annotation quality but increase memory usage.
                Defaults to 50.
            pcent_batch_size (int, optional): Batch size for calculating expression percentages to
                optimize memory usage. Defaults to 2000.
            
        Raises:
            KeyError: If the required keys are missing in `adata.obs` or `adata.uns`
            ValueError: If the data format is incorrect or there are validation errors
        """
        self.adata = adata
        self.cell_group_key = cell_group_key
        self.rank_genes_key = rank_genes_key
        self.results_key_added = results_key_added
        self.gene_symbols_column_name = gene_symbols_column_name
        self.pcent_batch_size = pcent_batch_size
        self.n_top_genes = n_top_genes

        # Validate the AnnData object
        _validate_adata(adata, cell_group_key, rank_genes_key, gene_symbols_column_name)

        # Create cluster mapping and cluster list
        self.cluster_map = {
            str(x): str(n + 1)
            for n, x in enumerate(sorted(adata.obs[cell_group_key].unique().tolist()))
        }
        self.clusters = [self.cluster_map[str(x)] for x in adata.obs[cell_group_key].values.tolist()]

        # Perform expensive computations during initialization
        logger.info("Calculating expression percentages.")
        self.expression_percentages = _calculate_pcent(
            adata=adata,
            clusters=self.clusters,
            batch_size=pcent_batch_size,
            gene_names=adata.var[self.gene_symbols_column_name].tolist(),
        )

        # Extract marker genes during initialization
        logger.info("Extracting marker genes.")
        self.marker_genes = _get_markers(
            adata=self.adata,
            cell_group_key=self.cell_group_key,
            rank_genes_key=self.rank_genes_key,
            ct_map=self.cluster_map,
            n_top_genes=n_top_genes,
            gene_symbols_col=self.gene_symbols_column_name,
        )

        logger.info("Data preparation completed. Ready for annotation.")

    def run(
        self,
        bio_context: dict[str, Any] | None = None,
        model_config: list[dict[str, Any]] | None = None,
        poll_interval_seconds: int = DEFAULT_POLL_INTERVAL,
        timeout_seconds: int = DEFAULT_TIMEOUT,
        api_url: str = DEFAULT_API_URL,
    ) -> anndata.AnnData:
        """Perform cell type annotation using the CyteType API.

        Args:
            bio_context (dict[str, Any] | None, optional): Biological context information.
                Can include keys: 'organisms', 'tissues', 'diseases', 'developmental_stages',
                'single_cell_methods', 'experimental_conditions'. Defaults to None.
            model_config (list[dict[str, Any]] | None, optional): Configuration for the large language
                models used for annotation. Each dict can include 'provider', 'name', 'apiKey', 'baseUrl'.
                Defaults to None, using the API's default model.
            poll_interval_seconds (int, optional): How often (in seconds) to check for results from
                the API. Defaults to DEFAULT_POLL_INTERVAL.
            timeout_seconds (int, optional): Maximum time (in seconds) to wait for API results before
                raising a timeout error. Defaults to DEFAULT_TIMEOUT.
            api_url (str, optional): URL for the CyteType API endpoint. Only change if using a custom
                deployment. Defaults to DEFAULT_API_URL.

        Returns:
            anndata.AnnData: The input AnnData object, modified in place with the following additions:
                - `adata.obs[f"{results_key_added}_{cell_group_key}"]`: Cell type annotations as categorical values
                - `adata.uns[f"{results_key_added}_results"]`: Complete API response data and job tracking info

        Raises:
            CyteTypeAPIError: If the API request fails or returns invalid data
            CyteTypeTimeoutError: If the API does not return results within the specified timeout period

        """
        api_url = api_url.rstrip("/")

        # Process bio context using Pydantic model
        if bio_context is None:
            bio_context = {}
        
        bio_context_model = BioContext(**bio_context)
        bio_context_dict = bio_context_model.model_dump(by_alias=True)

        # Process model config using Pydantic model
        model_config_list = None
        if model_config is not None:
            model_config_list = [ModelConfig(**config).model_dump(by_alias=True) for config in model_config]

        # Prepare API query
        query: dict[str, Any] = {
            "bioContext": bio_context_dict,
            "markerGenes": self.marker_genes,
            "expressionData": self.expression_percentages,
        }

        # Save query for debugging (optional)
        with open("query.json", "w") as f:
            json.dump(query, f)

        # Submit job and poll for results
        job_id = submit_annotation_job(query, api_url, model_config=model_config_list)
        logger.info(f"Waiting for results for job ID: {job_id}")
        
        # Log the report URL that updates automatically
        report_url = f"{api_url}/report/{job_id}"
        logger.info(f"View the automatically updating visualization report at: {report_url}")
        
        annotation_results = poll_for_results(
            job_id, api_url, poll_interval_seconds, timeout_seconds
        )

        # Store results in AnnData object
        self.adata.uns[f"{self.results_key_added}_results"] = {
            "job_id": job_id,
            "result": annotation_results,
        }

        # Create annotation mapping and add to observations
        annotation_map = {
            item["clusterId"]: item["annotation"]
            for item in annotation_results.get("annotations", [])
        }
        self.adata.obs[f"{self.results_key_added}_{self.cell_group_key}"] = pd.Series(
            [annotation_map.get(cluster_id, "Unknown Annotation") for cluster_id in self.clusters],
            index=self.adata.obs.index,
        ).astype("category")

        # Check for unannotated clusters
        unannotated_clusters = {cluster_id for cluster_id in self.clusters if cluster_id not in annotation_map}
        if unannotated_clusters:
            logger.warning(
                f"No annotations received from API for cluster IDs: {unannotated_clusters}. "
                f"Corresponding cells marked as 'Unknown Annotation'."
            )

        logger.info(
            f"Annotations successfully added to `adata.obs['{self.results_key_added}_{self.cell_group_key}']` "
            f"and `adata.uns['{self.results_key_added}_results']`."
        )

        return self.adata
