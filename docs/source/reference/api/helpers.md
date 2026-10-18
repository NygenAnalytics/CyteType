---
myst:
  html_meta:
    description: Reference for the CyteType helper functions marker_dotplot, rank_genes_groups_backed and subsample_by_group.
---

# Helpers

The `cytetype` package also exports three helper functions.

## Plot the supporting marker genes

`marker_dotplot` reads the results that `CyteType.run()` stored in the AnnData object and
draws a Scanpy dot plot of the genes that support each label, with the clusters grouped by
the categories in the results. It needs Scanpy, which is not a dependency of `cytetype`.

```python
from cytetype import marker_dotplot

marker_dotplot(adata, group_key="leiden")
```

```{eval-rst}
.. autofunction:: cytetype.marker_dotplot
```

## Find marker genes in backed AnnData

`rank_genes_groups_backed` computes marker genes for an AnnData object opened in backed
mode, such as with `anndata.read_h5ad(path, backed="r")`. It reads the cells in chunks of
`cell_batch_size` instead of loading the whole matrix, compares each group with all other
cells using Welch's t-test, and writes Scanpy-compatible results to `adata.uns[key_added]`.
Pass the same key to `CyteType` as `rank_key`. [Large datasets](../guides/large-datasets)
shows where it fits in a workflow.

```{eval-rst}
.. autofunction:: cytetype.rank_genes_groups_backed
```

## Subsample cells per cluster

`subsample_by_group` keeps at most `max_cells_per_group` randomly chosen cells from each
group and returns them as a new in-memory AnnData object. Groups with fewer cells are kept
whole, and `random_state` makes the selection reproducible.

```{eval-rst}
.. autofunction:: cytetype.subsample_by_group
```
