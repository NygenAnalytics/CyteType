import anndata
import numpy as np
import pandas as pd

from .marker_detection import _accumulate_group_stats


def aggregate_expression_percentages(
    adata: anndata.AnnData,
    clusters: list[str],
    gene_names: list[str],
    cell_batch_size: int = 5000,
) -> dict[str, dict[str, float]]:
    """Aggregate gene expression percentages per cluster.

    Uses a single-pass row-batched accumulation (fast for CSR / backed data).

    Args:
        adata: AnnData object containing expression data
        clusters: List of cluster assignments for each cell
        gene_names: List of gene names corresponding to columns in adata.X
        cell_batch_size: Number of cells to process per chunk

    Returns:
        Dictionary mapping gene names to cluster-level expression percentages
    """
    unique_clusters = sorted(set(clusters))
    n_groups = len(unique_clusters)
    cluster_to_idx = {c: i for i, c in enumerate(unique_clusters)}
    cell_group_indices = np.array([cluster_to_idx[c] for c in clusters])

    stats = _accumulate_group_stats(
        adata.X,
        cell_group_indices,
        n_groups,
        adata.shape[1],
        cell_batch_size=cell_batch_size,
        compute_nnz=True,
        progress_desc="expression_percentages",
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        pct_matrix = np.round(stats.nnz / stats.n[:, None] * 100, 2)
    pct_matrix = np.nan_to_num(pct_matrix, nan=0.0)

    pcent: dict[str, dict[str, float]] = {}
    for gene_idx, name in enumerate(gene_names):
        pcent[name] = {
            unique_clusters[g_idx]: float(pct_matrix[g_idx, gene_idx])
            for g_idx in range(n_groups)
        }
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
