import random
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
_CANONICAL_GENE_SYMBOLS_COLUMN = "__cytetype_gene_symbols"


def _is_gene_id_like(value: str) -> bool:
    if not isinstance(value, str) or not value.strip():
        return False

    value = value.strip()

    if re.match(r"^ENS[A-Z]*G\d{11}(\.\d+)?$", value, re.IGNORECASE):
        return True

    if re.match(r"^[NX][MR]_\d+$", value):
        return True

    if re.match(r"^\d+(?:\.0)?$", value):
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


def _extract_symbol_from_composite(value: str) -> str:
    parts = re.split(r"[_|]", value, maxsplit=1)
    if len(parts) != 2:
        return value
    id_flags = [_is_gene_id_like(p) for p in parts]
    if id_flags[0] and not id_flags[1]:
        return parts[1]
    if not id_flags[0] and id_flags[1]:
        return parts[0]
    return value


def clean_gene_names(names: list[str]) -> list[str]:
    """Extract gene symbols from composite gene name/ID values.

    If >50% of values are composite (e.g. ``TSPAN6_ENSG00000000003``),
    splits each value and returns the gene-symbol part.  Non-composite
    lists are returned unchanged.
    """
    if not _has_composite_gene_values(names):
        return names
    cleaned = [_extract_symbol_from_composite(n) for n in names]
    logger.info(
        f"Cleaned {len(cleaned)} composite gene values "
        f"(e.g., '{names[0]}' -> '{cleaned[0]}')."
    )
    return cleaned


def _temporary_gene_symbols_column_name(adata: anndata.AnnData) -> str:
    candidate = _CANONICAL_GENE_SYMBOLS_COLUMN
    suffix = 1
    while candidate in adata.var.columns:
        candidate = f"{_CANONICAL_GENE_SYMBOLS_COLUMN}_{suffix}"
        suffix += 1
    return candidate


def materialize_canonical_gene_symbols_column(
    adata: anndata.AnnData, gene_symbols_column: str | None
) -> tuple[str, str | None]:
    """Create a temporary canonical gene-symbol column in ``adata.var``."""
    if gene_symbols_column is None:
        source_values = adata.var_names.astype(str).tolist()
        source_name = "adata.var_names"
    else:
        source_values = [
            str(value)
            for value in adata.var[gene_symbols_column].astype("string").fillna("")
        ]
        source_name = f"column '{gene_symbols_column}'"

    canonical_column = _temporary_gene_symbols_column_name(adata)
    cleaned = clean_gene_names(source_values)

    id_pct = _id_like_percentage(cleaned)
    if id_pct > 49:
        raise ValueError(
            f"\n\nGene Symbol Detection Error\n"
            f"{'─' * 50}\n"
            f"CyteType requires human-readable gene symbols (e.g., TSPAN6, DPM1, SCYL3)\n"
            f"To fix this, either:\n"
            f"  1. Set gene_symbols_column to a column in adata.var that contains gene symbols\n"
            f"  2. Convert your gene identifiers to symbols before running CyteType\n"
        )

    adata.var[canonical_column] = cleaned
    logger.info(
        f"Materialized canonical gene symbols in temporary column '{canonical_column}' "
        f"from {source_name}."
    )
    return canonical_column, gene_symbols_column


def _id_like_percentage(values: list[str], seed: int = 42) -> float:
    if not values:
        return 100.0
    n = min(2000, len(values))
    if n < len(values):
        rng = random.Random(seed)
        sample = rng.sample(values, n)
    else:
        sample = values
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
        values = adata.var[gene_symbols_column].dropna().astype(str).tolist()
        if _has_composite_gene_values(values):
            logger.info(
                f"Column '{gene_symbols_column}' contains composite gene name/ID values "
                f"(e.g., '{values[0]}'). Gene symbols will be extracted automatically."
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
        score_values = (
            [_extract_symbol_from_composite(v) for v in values]
            if _has_composite_gene_values(values)
            else values
        )
        pct = _id_like_percentage(score_values)
        unique_ratio = len(set(score_values)) / len(score_values)
        candidates.append((col, pct, unique_ratio, 0))

    var_names_list = adata.var_names.astype(str).tolist()
    if var_names_list:
        var_score_values = (
            [_extract_symbol_from_composite(v) for v in var_names_list]
            if _has_composite_gene_values(var_names_list)
            else var_names_list
        )
        var_id_pct = _id_like_percentage(var_score_values)
        var_unique_ratio = len(set(var_score_values)) / len(var_score_values)
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
        score_values = (
            [_extract_symbol_from_composite(v) for v in values]
            if _has_composite_gene_values(values)
            else values
        )
        n_unique = len(set(score_values))
        if n_unique < max(10, len(score_values) * 0.05):
            continue
        pct = _id_like_percentage(score_values)
        unique_ratio = n_unique / len(score_values)
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


def _generate_unique_na_label(existing_labels: set[str]) -> str:
    label = "Unknown"
    if label not in existing_labels:
        return label
    n = 2
    while f"{label} {n}" in existing_labels:
        n += 1
    return f"{label} {n}"


def validate_adata(
    adata: anndata.AnnData,
    cell_group_key: str,
    rank_genes_key: str,
    gene_symbols_col: str | None,
    coordinates_key: str,
    label_na: bool = False,
) -> str | None:
    if cell_group_key not in adata.obs:
        raise KeyError(f"Cell group key '{cell_group_key}' not found in `adata.obs`.")

    nan_mask = adata.obs[cell_group_key].isna()
    n_nan = int(nan_mask.sum())
    if n_nan > 0:
        pct = round(100 * n_nan / adata.n_obs, 1)
        if n_nan == adata.n_obs:
            raise ValueError(
                f"All {n_nan} cells have NaN values in '{cell_group_key}'. "
                f"Cannot proceed with annotation."
            )
        if not label_na:
            raise ValueError(
                f"{n_nan} cells ({pct}%) have NaN values in '{cell_group_key}'. "
                f"Either fix the data or set label_na=True to assign these cells "
                f"an 'Unknown' cluster label."
            )

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

    rank_names = adata.uns[rank_genes_key]["names"]
    try:
        sample_genes = [
            str(g)
            for field in rank_names.dtype.names[:3]
            for g in rank_names[field][:20]
        ]
    except Exception:
        sample_genes = []

    if sample_genes:
        id_pct = _id_like_percentage(sample_genes)
        if id_pct > 50:
            examples = [g for g in sample_genes if _is_gene_id_like(g)][:3]
            logger.warning(
                f"rank_genes_groups results contain gene IDs rather than gene symbols "
                f"(e.g. {examples}). This typically happens when var_names were Ensembl "
                f"IDs at the time rank_genes_groups was run but have since been replaced "
                f"with gene symbols. Marker gene extraction may fail or produce empty "
                f"results. Consider re-running sc.tl.rank_genes_groups on the current adata."
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
