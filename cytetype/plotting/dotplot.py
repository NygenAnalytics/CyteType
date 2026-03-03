from typing import Any

import anndata

from ..config import logger
from ..core.results import load_local_results


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
        ``sc.pl.dotplot`` via the ``gene_symbols`` parameter.
    **kwargs
        Additional keyword arguments forwarded to ``sc.pl.dotplot``
        (e.g. ``cmap``, ``use_raw``).

    Returns
    -------
    The return value of ``sc.pl.dotplot``.
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

    markers: dict[str, list[str]] = {}
    categories_order: list[str] = []

    for category in cluster_categories:
        category_name = category["categoryName"]
        markers[category_name] = []

        for cluster_id in category["clusterIds"]:
            cluster_data = raw_annotations.get(cluster_id)
            if cluster_data is None:
                logger.warning(
                    "Cluster '{}' listed in clusterCategories but missing from raw_annotations, skipping.",
                    cluster_id,
                )
                continue

            full_output = cluster_data["latest"]["annotation"]["fullOutput"]
            cell_type = full_output["cellType"]

            categories_order.append(cell_type["label"])
            markers[category_name].extend(
                cell_type["keySupportingGenes"][:n_top_markers]
            )

        markers[category_name] = sorted(set(markers[category_name]))

    try:
        import scanpy as sc
    except ImportError:
        raise ImportError(
            "scanpy is required for plotting. Install it with: pip install scanpy"
        ) from None

    groupby = f"{results_prefix}_annotation_{group_key}"

    return sc.pl.dotplot(
        adata,
        markers,
        groupby=groupby,
        gene_symbols=gene_symbols,
        categories_order=categories_order,
        **kwargs,
    )
