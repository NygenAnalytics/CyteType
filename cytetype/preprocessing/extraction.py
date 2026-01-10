import anndata
import pandas as pd

from ..config import logger


def extract_marker_genes(
    adata: anndata.AnnData,
    cell_group_key: str,
    rank_genes_key: str,
    cluster_map: dict[str, str],
    n_top_genes: int,
    gene_symbols_col: str,
) -> dict[str, list[str]]:
    """Extract top marker genes from rank_genes_groups results.

    Args:
        adata: AnnData object containing rank_genes_groups results
        cell_group_key: Key in adata.obs containing cluster labels
        rank_genes_key: Key in adata.uns containing rank_genes_groups results
        cluster_map: Dictionary mapping original labels to cluster IDs
        n_top_genes: Number of top genes to extract per cluster
        gene_symbols_col: Column in adata.var containing gene symbols

    Returns:
        Dictionary mapping cluster IDs to lists of marker gene symbols

    Raises:
        ValueError: If marker genes cannot be extracted or validated
    """
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

    gene_ids_to_name = adata.var[gene_symbols_col].to_dict()
    markers = {}
    any_genes_found = False

    for group_name in mdf.columns.tolist():
        cluster_id = cluster_map.get(str(group_name), "")
        if not cluster_id:
            raise ValueError(
                f"Internal inconsistency: Group name '{group_name}' from rank_genes_groups results "
                f"was not found in the mapping generated from adata.obs['{cell_group_key}']. "
                f"Ensure rank_genes_groups was run on the same cell grouping."
            )
        top_genes = mdf[group_name].values[: min(n_top_genes, len(mdf))]
        if len(top_genes) == 0:
            logger.warning(
                f"No top genes found for group '{group_name}' (cluster '{cluster_id}')"
            )
        else:
            any_genes_found = True

        markers[cluster_id] = [
            gene_ids_to_name[gene] for gene in top_genes if gene in gene_ids_to_name
        ]

    if not any_genes_found:
        raise ValueError(
            "No marker genes found for any group. This could indicate issues with the "
            "rank_genes_groups analysis or that all groups have insufficient marker genes."
        )

    return markers


def extract_visualization_coordinates(
    adata: anndata.AnnData,
    coordinates_key: str | None,
    group_key: str,
    cluster_map: dict[str, str],
    max_cells_per_group: int = 1000,
    random_state: int = 42,
) -> tuple[list[list[float]] | None, list[str]]:
    """Extract and sample coordinates for visualization.

    Args:
        adata: AnnData object containing single-cell data
        coordinates_key: Key in adata.obsm containing coordinates
        group_key: Column name in adata.obs to group cells by
        cluster_map: Dictionary mapping original cluster labels to new cluster IDs
        max_cells_per_group: Maximum number of cells to sample per group
        random_state: Random seed for reproducible sampling

    Returns:
        tuple: (sampled_coordinates, sampled_cluster_labels)
            - sampled_coordinates: List of [x, y] coordinate pairs, or None if no coordinates
            - sampled_cluster_labels: List of cluster labels corresponding to sampled coordinates
    """
    if coordinates_key is None:
        logger.warning("No coordinates key provided, returning None coordinates.")
        return None, []

    coordinates = adata.obsm[coordinates_key]

    # Take only the first 2 dimensions for visualization
    if coordinates.shape[1] > 2:
        coordinates = coordinates[:, :2]
        logger.info(
            f"Using first 2 dimensions of '{coordinates_key}' for visualization."
        )

    # Create DataFrame with coordinates and group labels
    coord_df = pd.DataFrame(
        {
            "x": coordinates[:, 0],
            "y": coordinates[:, 1],
            "group": adata.obs[group_key].values,
        }
    )

    # Sample cells from each group using pandas
    sampled_coords = []
    for group_label in coord_df["group"].unique():
        group_mask = coord_df["group"] == group_label
        group_size = group_mask.sum()
        sample_size = min(max_cells_per_group, group_size)

        sampled_group = coord_df[group_mask].sample(
            n=sample_size, random_state=random_state
        )
        sampled_coords.append(sampled_group)

        if group_size > max_cells_per_group:
            logger.info(
                f"Sampled {sample_size} cells from group '{group_label}' "
                f"(originally {group_size} cells)"
            )

    # Concatenate all sampled groups
    sampled_coord_df: pd.DataFrame = pd.concat(sampled_coords, ignore_index=True)

    # Extract coordinates and labels
    sampled_coordinates = sampled_coord_df[["x", "y"]].values.tolist()

    # Map original cluster labels to new cluster IDs
    sampled_cluster_labels = [
        cluster_map.get(str(label), str(label))
        for label in sampled_coord_df["group"].values
    ]

    logger.info(
        f"Extracted {len(sampled_coordinates)} coordinate points "
        f"(sampled from {len(coordinates)} total cells)"
    )

    return sampled_coordinates, sampled_cluster_labels
