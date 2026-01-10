import anndata
import numpy as np
import pandas as pd


def aggregate_expression_percentages(
    adata: anndata.AnnData, clusters: list[str], batch_size: int, gene_names: list[str]
) -> dict[str, dict[str, float]]:
    """Aggregate gene expression percentages per cluster.

    Args:
        adata: AnnData object containing expression data
        clusters: List of cluster assignments for each cell
        batch_size: Number of genes to process per batch (for memory efficiency)
        gene_names: List of gene names corresponding to columns in adata.X

    Returns:
        Dictionary mapping gene names to cluster-level expression percentages
    """
    pcent = {}
    n_genes = adata.shape[1]

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


def aggregate_cluster_metadata(
    adata: anndata.AnnData,
    group_key: str,
    min_percentage: int = 10,
) -> dict[str, dict[str, dict[str, int]]]:
    """Aggregate categorical metadata per cluster.

    For each categorical column (excluding the group_key), calculates the percentage
    distribution of values within each group and returns only values that represent
    more than min_percentage of cells in that group.

    Args:
        adata: AnnData object containing single-cell data
        group_key: Column name in adata.obs to group cells by
        min_percentage: Minimum percentage of cells in a group to include

    Returns:
        Nested dictionary structure:
        {group_name: {column_name: {value: percentage}}}
        where percentage is the percentage of cells in that group having that value
        (only values >min_percentage are included)
    """
    grouped_data = adata.obs.groupby(group_key, observed=False)
    column_distributions: dict[str, dict[str, dict[str, int]]] = {}

    # Process each column in adata.obs
    for column_name in adata.obs.columns:
        if column_name == group_key:
            continue

        column_dtype = adata.obs[column_name].dtype
        if column_dtype in ["object", "category", "string"]:
            # Calculate value counts for each group
            value_counts_df = grouped_data[column_name].value_counts().unstack().T

            # Convert to percentages and filter for values >min_percentage
            percentage_df = (
                (100 * value_counts_df / value_counts_df.sum())
                .fillna(0)
                .astype(int)
                .T.stack()
            )
            significant_values = percentage_df[percentage_df > min_percentage].to_dict()

            # Reorganize into nested dictionary structure
            group_value_percentages: dict[str, dict[str, int]] = {}
            for (group_name, value), percentage in significant_values.items():
                group_name = str(group_name)
                value = str(value)
                if group_name not in group_value_percentages:
                    group_value_percentages[group_name] = {}
                group_value_percentages[group_name][value] = percentage

            column_distributions[column_name] = group_value_percentages

    # Reorganize final structure: {group_name: {column_name: {value: percentage}}}
    result: dict[str, dict[str, dict[str, int]]] = {
        str(group_name): {} for group_name in grouped_data.groups.keys()
    }

    for column_name in column_distributions:
        for group_name in column_distributions[column_name]:
            result[group_name][column_name] = column_distributions[column_name][
                group_name
            ]

    return result
