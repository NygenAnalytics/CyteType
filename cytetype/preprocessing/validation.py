import re
import anndata

from ..config import logger


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


def validate_adata(
    adata: anndata.AnnData,
    cell_group_key: str,
    rank_genes_key: str,
    gene_symbols_col: str,
    coordinates_key: str,
) -> str | None:
    """Validate the AnnData object structure and return the best available coordinates key.

    Args:
        adata: AnnData object to validate
        cell_group_key: Key in adata.obs containing cluster labels
        rank_genes_key: Key in adata.uns containing rank_genes_groups results
        gene_symbols_col: Column in adata.var containing gene symbols
        coordinates_key: Preferred key in adata.obsm for coordinates

    Returns:
        str | None: The coordinates key that was found and validated, or None if no suitable coordinates found.

    Raises:
        KeyError: If required keys are missing
        ValueError: If data format is incorrect
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
