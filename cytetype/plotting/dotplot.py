from typing import Any

import anndata

from ..config import logger
from ..core.results import load_local_results


def _gene_symbols(value: object) -> list[str]:
    if not isinstance(value, list):
        return []
    return [gene for gene in value if isinstance(gene, str) and gene.strip()]


def _plot_gene_names(
    adata: anndata.AnnData,
    gene_symbols: str | None,
    use_raw: bool | None,
    layer: str | None,
) -> tuple[list[str], str]:
    """Return the gene names ``sc.pl.dotplot`` looks up, and where they come from."""
    # Same default as scanpy: use .raw when it exists and no layer is given.
    if use_raw is None:
        use_raw = layer is None and adata.raw is not None
    raw = adata.raw if use_raw else None
    var = raw.var if raw is not None else adata.var
    source = "adata.raw.var" if raw is not None else "adata.var"
    if gene_symbols is None:
        return [str(name) for name in var.index], f"{source}_names"
    if gene_symbols not in var.columns:
        raise KeyError(f"Column '{gene_symbols}' not found in {source}.")
    names = [str(value) for value in var[gene_symbols].dropna()]
    return names, f"{source}['{gene_symbols}']"


def marker_dotplot(
    adata: anndata.AnnData,
    group_key: str,
    results_prefix: str = "cytetype",
    n_top_markers: int = 3,
    gene_symbols: str | None = None,
    **kwargs: Any,
) -> Any:
    """Dotplot of top marker genes grouped by CyteType cluster categories.

    Reads stored CyteType results from ``adata.uns`` and builds a
    category-grouped marker dict suitable for ``sc.pl.dotplot``.

    Parameters
    ----------
    adata
        AnnData object that has been annotated by CyteType
        (results stored in ``adata.uns``).
    group_key
        The original cluster key used during annotation (e.g. ``"leiden"``).
    results_prefix
        Prefix used when the results were stored. Must match the
        ``results_prefix`` passed to :meth:`CyteType.run`.
    n_top_markers
        Number of top supporting genes to display per cluster.
    gene_symbols
        Column in ``adata.var`` containing gene symbols.  Forwarded to
        ``sc.pl.dotplot`` via the ``gene_symbols`` parameter; marker genes
        are matched to its values instead of ``var_names``.
    **kwargs
        Additional keyword arguments forwarded to ``sc.pl.dotplot``
        (e.g. ``cmap``, ``use_raw``).

    Returns
    -------
    The return value of ``sc.pl.dotplot``.

    Notes
    -----
    Rows follow the order of the cluster categories. Labels without a
    category (e.g. ``"Unknown"``) are placed last. Genes come from the
    structured output's ``keySupportingGenes``, falling back to
    ``supportingMarkers`` for jobs without it. They are matched
    case-insensitively to the names ``sc.pl.dotplot`` looks up (the server
    reports upper-case symbols, e.g. ``CD8A`` for mouse ``Cd8a``). Genes that
    are not found are skipped with a warning, and categories left without
    genes are dropped.
    """
    results = load_local_results(adata, results_prefix)
    if results is None:
        raise KeyError(
            f"No CyteType results found in adata.uns with prefix '{results_prefix}'. "
            "Run CyteType annotation first or check the results_prefix."
        )

    cluster_categories = results.get("clusterCategories", [])
    if not cluster_categories:
        raise ValueError(
            "No cluster categories found in CyteType results. "
            "The API response may not include category groupings for this run."
        )

    raw_annotations = results.get("raw_annotations", {})
    if not raw_annotations:
        raise ValueError("No raw_annotations found in CyteType results.")

    groupby = f"{results_prefix}_annotation_{group_key}"
    if groupby not in adata.obs:
        raise KeyError(
            f"Annotation column '{groupby}' not found in adata.obs. "
            "Check group_key and results_prefix."
        )
    # scanpy logs an error and does not plot unless categories_order holds
    # exactly the groupby categories, so derive it from them.
    obs_categories = list(adata.obs[groupby].astype("category").cat.categories)
    label_lookup = set(obs_categories)

    gene_names, gene_source = _plot_gene_names(
        adata, gene_symbols, kwargs.get("use_raw"), kwargs.get("layer")
    )
    exact_names = set(gene_names)
    names_by_folded_case: dict[str, str] = {}
    for name in gene_names:
        names_by_folded_case.setdefault(name.casefold(), name)

    markers: dict[str, list[str]] = {}
    categories_order: list[Any] = []
    unmatched_genes: list[str] = []

    for category in cluster_categories:
        category_name = category["categoryName"]
        markers[category_name] = []

        for cluster_id in category["clusterIds"]:
            cluster_data = raw_annotations.get(cluster_id)
            latest = (
                cluster_data.get("latest") if isinstance(cluster_data, dict) else None
            )
            annotation = latest.get("annotation") if isinstance(latest, dict) else None
            if not isinstance(annotation, dict):
                logger.warning(
                    "Cluster '{}' listed in clusterCategories but has no annotation in raw_annotations, skipping.",
                    cluster_id,
                )
                continue

            full_output = annotation.get("fullOutput")
            cell_type = (
                full_output.get("cellType") if isinstance(full_output, dict) else None
            )
            if not isinstance(cell_type, dict):
                cell_type = {}

            # The groupby column holds annotation.annotation, which the server may
            # have refined; cellType.label is the fallback.
            for label in (annotation.get("annotation"), cell_type.get("label")):
                if isinstance(label, str) and label in label_lookup:
                    if label not in categories_order:
                        categories_order.append(label)
                    break

            genes = _gene_symbols(cell_type.get("keySupportingGenes")) or _gene_symbols(
                annotation.get("supportingMarkers")
            )
            for gene in genes[:n_top_markers]:
                matched = (
                    gene
                    if gene in exact_names
                    else names_by_folded_case.get(gene.casefold())
                )
                if matched is None:
                    unmatched_genes.append(gene)
                else:
                    markers[category_name].append(matched)

        markers[category_name] = sorted(set(markers[category_name]))

    if unmatched_genes:
        missing = list(dict.fromkeys(unmatched_genes))
        logger.warning(
            f"{len(missing)} marker genes are not shown because they were not found "
            f"in {gene_source}: {', '.join(missing[:5])}"
            f"{', ...' if len(missing) > 5 else ''}."
        )
    # scanpy would still draw a bracket and label for an empty gene group.
    markers = {name: genes for name, genes in markers.items() if genes}
    if not markers:
        raise ValueError(
            f"None of the marker genes were found in {gene_source}. Check "
            "gene_symbols (and use_raw), or whether the data contains these genes."
        )

    categories_order += [c for c in obs_categories if c not in categories_order]

    try:
        import scanpy as sc
    except ImportError:
        raise ImportError(
            "scanpy is required for plotting. Install it with: pip install scanpy"
        ) from None

    return sc.pl.dotplot(
        adata,
        markers,
        groupby=groupby,
        gene_symbols=gene_symbols,
        categories_order=categories_order,
        **kwargs,
    )
