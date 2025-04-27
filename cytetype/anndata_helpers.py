import anndata
import numpy as np
import pandas as pd

from .config import logger


def _validate_adata(
    adata: anndata.AnnData, cell_group_key: str, rank_genes_key: str
) -> None:
    """Validate the AnnData object structure."""

    if cell_group_key not in adata.obs:
        raise KeyError(f"Cell group key '{cell_group_key}' not found in `adata.obs`.")
    if adata.X is None:
        raise ValueError(
            "`adata.X` is required for ranking genes. Please ensure it contains log1p normalized data."
        )
    if len(adata.var_names) != adata.shape[1]:
        raise ValueError("`adata.var_names` is not same size as `adata.X`")
    if rank_genes_key not in adata.uns:
        raise KeyError(
            f"'{rank_genes_key}' not found in `adata.uns`. Run `sc.tl.rank_genes_groups` first."
        )
    if adata.uns[rank_genes_key]["params"]["groupby"] != cell_group_key:
        raise ValueError(
            f"`rank_genes_groups` run with groupby='{adata.uns[rank_genes_key]['params']['groupby']}', expected '{cell_group_key}'."
        )
    if "names" not in adata.uns[rank_genes_key] or not hasattr(
        adata.uns[rank_genes_key]["names"], "dtype"
    ):
        raise ValueError(
            f"'names' field in `adata.uns['{rank_genes_key}']` is missing or invalid."
        )


def _calculate_pcent(
    adata: anndata.AnnData, clusters: np.ndarray, batch_size: int
) -> dict[str, dict[int, float]]:
    """Calculate percentage of cells expressing each gene within clusters."""

    pcent = {}
    n_genes = adata.shape[1]
    gene_names = adata.var_names

    for s in range(0, n_genes, batch_size):
        e = min(s + batch_size, n_genes)
        batch_data = adata.X[:, s:e]
        if hasattr(batch_data, "toarray"):
            batch_data = batch_data.toarray()
        elif isinstance(batch_data, np.ndarray):
            pass
        else:
            raise TypeError(
                f"Unexpected data type in `adata.raw.X` slice: {type(batch_data)}"
            )

        df = pd.DataFrame(batch_data > 0, columns=gene_names[s:e]) * 100
        df["clusters"] = clusters
        pcent.update(df.groupby("clusters").mean().round(2).to_dict())
        del df, batch_data
    return pcent


def _get_markers(
    adata: anndata.AnnData,
    cell_group_key: str,
    rank_genes_key: str,
    ct_map: dict[str, int],
    n_top_genes: int,
) -> dict[str, list[str]]:
    """Extract top marker genes from rank_genes_groups results."""
    try:
        mdf = pd.DataFrame(adata.uns[rank_genes_key]["names"])
    except ValueError:
        logger.warning(
            "Could not directly convert `rank_genes_groups['names']` to DataFrame. Attempting alternative."
        )
        try:
            names_rec = adata.uns[rank_genes_key]["names"]
            mdf = pd.DataFrame(
                {field: names_rec[field] for field in names_rec.dtype.names}
            )
        except Exception as e:
            raise ValueError(
                f"Failed to extract marker gene names from `rank_genes_groups`. Error: {e}"
            )

    markers = {}
    # Use var_names directly if features attribute doesn't exist
    if hasattr(adata.var, "features"):
        gene_ids_to_name = adata.var.features.to_dict()
    else:
        # Create a simple mapping from index (gene_id) to itself (gene_name)
        gene_ids_to_name = {gene_id: gene_id for gene_id in adata.var_names}

    for group_name in list(mdf.columns):
        api_cluster_id = None
        for original_label, mapped_id in ct_map.items():
            if str(original_label) == str(group_name):
                api_cluster_id = str(mapped_id)
                break

        if api_cluster_id is None:
            raise ValueError(
                f"Internal inconsistency: Group name '{group_name}' from rank_genes_groups results "
                f"was not found in the mapping generated from adata.obs['{cell_group_key}']. "
                f"Ensure rank_genes_groups was run on the same cell grouping."
            )

        top_genes = mdf[group_name].values[: min(n_top_genes, len(mdf))]
        # Filter out genes that aren't in the gene_ids_to_name dictionary
        valid_genes = [gene for gene in top_genes if gene in gene_ids_to_name]
        markers[api_cluster_id] = [gene_ids_to_name[gene] for gene in valid_genes]

    return markers
