import anndata
import numpy as np
import pandas as pd
import re

from .config import logger


def _is_gene_id_like(value: str) -> bool:
    """Check if a value looks like a gene ID rather than a gene symbol.

    Common gene ID patterns:
    - Ensembl: ENSG00000000003, ENSMUSG00000000001, etc.
    - RefSeq: NM_000001, XM_000001, etc.
    - Numeric IDs: just numbers
    - Other database IDs with similar patterns

    Args:
        value: String value to check

    Returns:
        bool: True if the value looks like a gene ID, False if it looks like a gene symbol
    """
    if not isinstance(value, str) or not value.strip():
        return False

    value = value.strip()

    # Ensembl IDs (human, mouse, etc.)
    if re.match(r"^ENS[A-Z]*G\d{11}$", value, re.IGNORECASE):
        return True

    # RefSeq IDs
    if re.match(r"^[NX][MR]_\d+$", value):
        return True

    # Purely numeric IDs
    if re.match(r"^\d+$", value):
        return True

    # Other common ID patterns (long alphanumeric with underscores/dots)
    if re.match(r"^[A-Z0-9]+[._][A-Z0-9._]+$", value) and len(value) > 10:
        return True

    return False


def _validate_gene_symbols_column(
    adata: anndata.AnnData, gene_symbols_col: str
) -> None:
    """Validate that the gene_symbols_col contains gene symbols rather than gene IDs.

    Args:
        adata: AnnData object
        gene_symbols_col: Column name in adata.var that should contain gene symbols

    Raises:
        ValueError: If the column appears to contain gene IDs instead of gene symbols
    """
    gene_values = adata.var[gene_symbols_col].dropna().astype(str)

    if len(gene_values) == 0:
        logger.warning(
            f"Column '{gene_symbols_col}' is empty or contains only NaN values."
        )
        return

    # Sample a subset for efficiency (check up to 1000 non-null values)
    sample_size = min(1000, len(gene_values))
    sample_values = gene_values.sample(n=sample_size)

    # Count how many look like gene IDs vs gene symbols
    id_like_count = sum(1 for value in sample_values if _is_gene_id_like(value))
    id_like_percentage = (id_like_count / len(sample_values)) * 100

    if id_like_percentage > 50:
        example_ids = [
            value for value in sample_values.iloc[:5] if _is_gene_id_like(value)
        ]
        logger.warning(
            f"Column '{gene_symbols_col}' appears to contain gene IDs rather than gene symbols. "
            f"{id_like_percentage:.1f}% of values look like gene IDs (e.g., {example_ids[:3]}). "
            f"The annotation might not be accurate. Consider using a column that contains "
            f"human-readable gene symbols (e.g., 'TSPAN6', 'DPM1', 'SCYL3') instead of database identifiers."
        )
    elif id_like_percentage > 20:
        logger.warning(
            f"Column '{gene_symbols_col}' contains {id_like_percentage:.1f}% values that look like gene IDs. "
            f"Please verify this column contains gene symbols rather than gene identifiers."
        )


def _validate_adata(
    adata: anndata.AnnData,
    cell_group_key: str,
    rank_genes_key: str,
    gene_symbols_col: str,
    coordinates_key: str,
) -> str | None:
    """Validate the AnnData object structure and return the best available coordinates key.

    Returns:
        str | None: The coordinates key that was found and validated, or None if no suitable coordinates found.
    """

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
    if hasattr(adata.var, gene_symbols_col) is False:
        raise KeyError(f"Column '{gene_symbols_col}' not found in `adata.var`.")
    _validate_gene_symbols_column(adata, gene_symbols_col)

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

    # Validate coordinates with fallback options (case-insensitive matching)
    common_coordinate_keys = [coordinates_key, "X_umap", "X_tsne", "X_pca"]
    found_coordinates_key: str | None = None

    # Create a case-insensitive lookup for available keys
    available_keys = list(adata.obsm.keys())
    key_lookup = {key.lower(): key for key in available_keys}

    for key in common_coordinate_keys:
        # Try case-insensitive match
        actual_key = key_lookup.get(key.lower())
        if actual_key is not None:
            coordinates = adata.obsm[actual_key]
            if coordinates.shape[0] == adata.shape[0]:
                if coordinates.shape[1] >= 2:
                    found_coordinates_key = actual_key
                    if actual_key != key:
                        logger.info(
                            f"Using coordinates from '{actual_key}' (matched '{key}' case-insensitively) for visualization."
                        )
                    else:
                        logger.info(
                            f"Using coordinates from '{actual_key}' for visualization."
                        )
                    break
                else:
                    logger.warning(
                        f"Coordinates in '{actual_key}' have shape {coordinates.shape}, need at least 2 dimensions."
                    )
            else:
                logger.warning(
                    f"Coordinates in '{actual_key}' have {coordinates.shape[0]} rows, expected {adata.shape[0]}."
                )

    if found_coordinates_key is None:
        logger.warning(
            f"No suitable 2D coordinates found in adata.obsm. "
            f"Looked for: {common_coordinate_keys} (case-insensitive). "
            f"Available keys: {available_keys}. "
            f"Visualization will be disabled."
        )

    return found_coordinates_key


def _extract_sampled_coordinates(
    adata: anndata.AnnData,
    coordinates_key: str | None,
    group_key: str,
    cluster_map: dict[str, str],
    max_cells_per_group: int = 1000,
    random_state: int = 42,
) -> tuple[list[list[float]] | None, list[str]]:
    """Extract coordinates with sampling to limit the number of points per group.

    Args:
        adata: AnnData object containing single-cell data
        coordinates_key: Key in adata.obsm containing coordinates
        group_key: Column name in adata.obs to group cells by
        cluster_map: Dictionary mapping original cluster labels to new cluster IDs
        max_cells_per_group: Maximum number of cells to sample per group (default: 1000)
        random_state: Random seed for reproducible sampling (default: 42)

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


def _calculate_pcent(
    adata: anndata.AnnData, clusters: list[str], batch_size: int, gene_names: list[str]
) -> dict[str, dict[str, float]]:
    """Calculate percentage of cells expressing each gene within clusters."""

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


def _get_markers(
    adata: anndata.AnnData,
    cell_group_key: str,
    rank_genes_key: str,
    ct_map: dict[str, str],
    n_top_genes: int,
    gene_symbols_col: str,
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

    gene_ids_to_name = adata.var[gene_symbols_col].to_dict()
    markers = {}
    any_genes_found = False

    for group_name in mdf.columns.tolist():
        cluster_id = ct_map.get(str(group_name), "")
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


def _aggregate_metadata(
    adata: anndata.AnnData,
    group_key: str,
    min_percentage: int = 10,
) -> dict[str, dict[str, dict[str, int]]]:
    """
    Build group metadata by analyzing categorical/string columns in adata.obs.

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
