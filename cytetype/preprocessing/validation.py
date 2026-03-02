import re
import anndata

from ..config import logger

_KNOWN_GENE_SYMBOL_COLUMNS = [
    "feature_name",
    "gene_symbols",
    "gene_symbol",
    "gene_short_name",
    "gene_name",
    "symbol",
]


def _is_gene_id_like(value: str) -> bool:
    if not isinstance(value, str) or not value.strip():
        return False

    value = value.strip()

    if re.match(r"^ENS[A-Z]*G\d{11}$", value, re.IGNORECASE):
        return True

    if re.match(r"^[NX][MR]_\d+$", value):
        return True

    if re.match(r"^\d+$", value):
        return True

    if re.match(r"^[A-Z0-9]+[._][A-Z0-9._]+$", value) and len(value) > 10:
        return True

    return False


def _has_composite_gene_values(values: list[str]) -> bool:
    """Detect values like 'TSPAN6_ENSG00000000003' (symbol_id or id_symbol)."""
    composite_count = 0
    for v in values[:200]:
        parts = re.split(r"[_|]", v, maxsplit=1)
        if len(parts) == 2:
            id_flags = [_is_gene_id_like(p) for p in parts]
            if id_flags[0] != id_flags[1]:
                composite_count += 1
    return len(values) > 0 and (composite_count / min(200, len(values))) > 0.5


def _id_like_percentage(values: list[str]) -> float:
    if not values:
        return 100.0
    n = min(500, len(values))
    sample = values[:n]
    return sum(1 for v in sample if _is_gene_id_like(v)) / n * 100


def _validate_gene_symbols_column(
    adata: anndata.AnnData, gene_symbols_col: str
) -> None:
    gene_values = adata.var[gene_symbols_col].dropna().astype(str)

    if len(gene_values) == 0:
        logger.warning(
            f"Column '{gene_symbols_col}' is empty or contains only NaN values."
        )
        return

    values_list = gene_values.tolist()
    pct = _id_like_percentage(values_list)

    if pct > 50:
        example_ids = [v for v in values_list[:20] if _is_gene_id_like(v)][:3]
        logger.warning(
            f"Column '{gene_symbols_col}' appears to contain gene IDs rather than gene symbols. "
            f"{pct:.1f}% of values look like gene IDs (e.g., {example_ids}). "
            f"The annotation might not be accurate. Consider using a column that contains "
            f"human-readable gene symbols (e.g., 'TSPAN6', 'DPM1', 'SCYL3') instead of database identifiers."
        )
    elif pct > 20:
        logger.warning(
            f"Column '{gene_symbols_col}' contains {pct:.1f}% values that look like gene IDs. "
            f"Please verify this column contains gene symbols rather than gene identifiers."
        )


def resolve_gene_symbols_column(
    adata: anndata.AnnData, gene_symbols_column: str | None
) -> str | None:
    """Resolve which source contains gene symbols.

    Returns the column name in adata.var, or None if var_names should be used directly.
    """
    if gene_symbols_column is not None:
        if gene_symbols_column not in adata.var.columns:
            raise KeyError(
                f"Column '{gene_symbols_column}' not found in `adata.var`. "
                f"Available columns: {list(adata.var.columns)}. "
                f"Set gene_symbols_column=None for auto-detection."
            )
        _validate_gene_symbols_column(adata, gene_symbols_column)
        logger.info(f"Using gene symbols from column '{gene_symbols_column}'.")
        return gene_symbols_column

    # --- Auto-detection: score all candidates, then pick the best ---
    # Each candidate: (column_name | None, id_like_pct, unique_ratio, priority)
    #   column_name=None → use var_names.
    #   priority: 0 = known column, 1 = var_names, 2 = other column.
    # Sorted by (id_like_pct ↑, priority ↑, unique_ratio ↓ with 1.0 penalized)
    # so the lowest ID-like % wins; ties broken by known columns first, then
    # by higher unique ratio (gene names have high cardinality, unlike
    # categorical metadata — but exactly 1.0 is slightly penalized because
    # IDs are always unique while gene symbols occasionally have duplicates).
    _KNOWN_SET = set(_KNOWN_GENE_SYMBOL_COLUMNS)
    candidates: list[tuple[str | None, float, float, int]] = []

    for col in _KNOWN_GENE_SYMBOL_COLUMNS:
        if col not in adata.var.columns:
            continue
        values = adata.var[col].dropna().astype(str).tolist()
        if not values:
            continue
        if _has_composite_gene_values(values):
            logger.warning(
                f"Column '{col}' appears to contain composite gene name/ID values "
                f"(e.g., '{values[0]}'). Skipping."
            )
            continue
        pct = _id_like_percentage(values)
        unique_ratio = len(set(values)) / len(values)
        candidates.append((col, pct, unique_ratio, 0))

    var_names_list = adata.var_names.astype(str).tolist()
    if var_names_list:
        var_id_pct = _id_like_percentage(var_names_list)
        var_unique_ratio = len(set(var_names_list)) / len(var_names_list)
        candidates.append((None, var_id_pct, var_unique_ratio, 1))

    for col in adata.var.columns:
        if col in _KNOWN_SET:
            continue
        try:
            values = adata.var[col].dropna().astype(str).tolist()
        except (TypeError, ValueError):
            continue
        if not values:
            continue
        if _has_composite_gene_values(values):
            continue
        n_unique = len(set(values))
        if n_unique < max(10, len(values) * 0.05):
            continue
        pct = _id_like_percentage(values)
        unique_ratio = n_unique / len(values)
        candidates.append((col, pct, unique_ratio, 2))

    viable = [c for c in candidates if c[1] < 50]

    def _ur_sort_key(ur: float) -> float:
        return ur if ur < 1.0 else ur - 0.02

    if viable:
        viable.sort(key=lambda c: (c[1], c[3], -_ur_sort_key(c[2])))
        best_col, best_pct, best_ur, _ = viable[0]

        if best_col is not None:
            source = f"column '{best_col}'"
            _validate_gene_symbols_column(adata, best_col)
        else:
            source = "adata.var_names (index)"

        logger.info(
            f"Auto-detected gene symbols in {source} "
            f"({best_pct:.0f}% ID-like, {best_ur:.0%} unique)."
        )
        return best_col

    # No viable candidate: fall back to var_names with warning
    fallback_pct = var_id_pct if var_names_list else 100.0
    logger.warning(
        "Could not find a column containing gene symbols in adata.var. "
        "Falling back to adata.var_names, but they appear to contain gene IDs "
        f"({fallback_pct:.0f}% ID-like). Annotation quality may be affected. "
        "Consider providing gene_symbols_column explicitly."
    )
    return None


def validate_adata(
    adata: anndata.AnnData,
    cell_group_key: str,
    rank_genes_key: str,
    gene_symbols_col: str | None,
    coordinates_key: str,
) -> str | None:
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
