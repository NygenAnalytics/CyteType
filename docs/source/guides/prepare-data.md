---
description: Get an AnnData or Seurat object ready for CyteType, from normalized expression and cluster labels to marker genes, gene symbols, coordinates, raw counts and cell metadata.
---

(prepare-data)=
# Prepare your data

CyteType annotates the clusters you already have, starting from the marker genes you computed for them. Before a job is submitted, the client checks your object and summarizes every cluster on your computer. This page explains what it needs, what it does with each input, and which files it writes.

## Checklist

| Input | AnnData (Python) | Seurat (R) |
|---|---|---|
| Normalized expression | Log-normalized values in `adata.X` | `NormalizeData()` on the `RNA` assay |
| Clusters | A column of `adata.obs`, named by `group_key`, with a label for every cell | A column of `obj@meta.data`, named by `group_key` (default `seurat_clusters`), with a label for every cell |
| Marker genes | `sc.tl.rank_genes_groups()` results for the same column, stored under `rank_key` | A `FindAllMarkers()` table with the columns `cluster`, `gene` and `avg_log2FC` |
| Gene symbols | In `adata.var_names` or a column of `adata.var` | As the row names of the object |
| Coordinates (optional) | A 2D embedding in `adata.obsm`, by default `X_umap` | A reduction, by default `umap` |
| Raw counts (optional) | Integer counts in `adata.layers["counts"]`, `adata.raw.X` or `adata.X` | Integer counts in the `counts` layer of the `RNA` assay |
| Cell metadata (optional) | Categorical columns of `adata.obs` | Factor or character columns of `obj@meta.data` |

A minimal preparation looks like this:

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
import scanpy as sc
from cytetype import CyteType

# adata.X holds log-normalized expression; adata.obs["leiden"] holds the clusters.
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", key_added="rank_genes_leiden")

# Check the inputs, summarize the clusters and write the upload files.
annotator = CyteType(adata, group_key="leiden", rank_key="rank_genes_leiden")
```
:::
:::{tab-item} R
:sync: r
```r
library(Seurat)
library(CyteTypeR)

# The RNA assay holds log-normalized data; meta.data$seurat_clusters holds the clusters.
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(obj, only.pos = TRUE)

# Check the inputs, summarize the clusters and write the upload files.
prepped_data <- PrepareCyteTypeR(obj, markers, group_key = "seurat_clusters")
```
:::
::::

The checks of your data and the cluster summaries run in this step, not when you submit the job, so problems with your data show up here. On a large dataset the step takes a while; see [Large datasets](large-datasets.md). Every argument is described in the [Python API reference](../reference/api/cytetype.html#cytetype.CyteType), and in R with `?PrepareCyteTypeR`.

## Normalized expression

CyteType uses the normalized values in two ways: it computes, for every gene, the percentage of cells in each cluster with a non-zero value, and it uploads the matrix so that the server can compute further evidence from your cells.

::::{tab-set}
:::{tab-item} Python
:sync: python
`adata.X` must hold log-normalized expression, for example the result of `sc.pp.normalize_total()` followed by `sc.pp.log1p()`.

- **Not scaled values.** After `sc.pp.scale()`, nearly every value is non-zero, so every gene would look expressed in every cell.
- **All genes, not only highly variable ones.** The expression summaries and the uploaded matrix cover only the genes in the object. Genes you leave out cannot be checked during annotation or searched in the report.

If you scaled `adata.X` or reduced it to highly variable genes, and the log-normalized values for all genes are in `adata.raw`, prepare the job from those:

```python
# Keeps obs, obsm and uns, so the clusters, the UMAP and the marker results come along.
adata_all = adata.raw.to_adata()
```

`to_adata()` does not copy layers, so a counts layer is not carried over.
:::
:::{tab-item} R
:sync: r
CyteTypeR reads the `data` layer of the `RNA` assay, which `NormalizeData()` fills with log-normalized values. The gene names and the raw counts come from the same assay. Without a `data` layer, `PrepareCyteTypeR()` stops with `The RNA assay has no normalized 'data' layer. Run NormalizeData() first.`

- **Any default assay.** After `SCTransform()` or integration, another assay, such as `SCT`, is usually the default. CyteTypeR still reads the `RNA` assay and warns once that it does, for example `DefaultAssay(seurat_obj) is 'SCT'. Using assay='RNA' instead.` An object without an `RNA` assay stops with `This workflow requires an assay named 'RNA'`.
- **Split layers.** If the `RNA` assay of a Seurat v5 object is split into layers, for example by sample before integration, CyteTypeR joins them on a copy. Your object is not changed, but the copy takes memory (see [Large datasets](large-datasets.md#large-datasets-memory)).
:::
::::

(prepare-data-clusters)=
## Clusters

`group_key` names the column that holds your cluster labels. These can be Leiden, Louvain or Seurat cluster numbers, or named groups such as an earlier annotation. CyteType annotates each cluster as a whole, and every cell in a cluster receives the cluster's result.

Keep the cluster column categorical or text, as Scanpy's Leiden and Seurat's `seurat_clusters` columns are. The server finds the column in the uploaded cell metadata by its values and does not recognize a numeric column; the evidence that CyteType computes from your cells is then missing, without an error.

### Missing labels

Every cell needs a label.

::::{tab-set}
:::{tab-item} Python
:sync: python
If any cell has a missing (NaN) label, `CyteType()` stops with an error such as `10 cells (1.7%) have NaN values in 'leiden'`. You have two options:

1. **Label those cells yourself** before you compute marker genes. They then form a cluster that is annotated like the others:

   ```python
   adata.obs["leiden"] = adata.obs["leiden"].cat.add_categories("Unassigned").fillna("Unassigned")
   ```

2. **Pass `label_na=True`.** CyteType puts the cells into a cluster named `Unknown` (or `Unknown 2`, and so on, if that label exists) in a new AnnData object, which `run()` returns. Your original object is not changed, so always use the return value of `run()`. Because the marker results do not include this cluster, it is not annotated, and its cells get `Unknown` in the result columns.

In general, only clusters that appear in the marker results are annotated.
:::
:::{tab-item} R
:sync: r
If any cell has no label (`NA`), `PrepareCyteTypeR()` stops with an error such as `3 cells (3.8%) have NA values in 'cell_group'`. CyteTypeR has no option like `label_na`, so give those cells a label, or remove them, before you compute the marker genes:

```r
obj$cell_group <- as.character(obj$cell_group)
obj$cell_group[is.na(obj$cell_group)] <- "Unassigned"
```
:::
::::

### Cluster IDs and what you see in the report

CyteType gives each cluster an ID, which can differ from your label. The Python client keeps short labels, such as Leiden numbers, as IDs, while CyteTypeR always numbers the clusters from `1`, so Seurat cluster `0` becomes cluster `1`. [Results in AnnData and Seurat](../results/results-in-your-data.md#results-cluster-ids) gives the exact rules and shows how to map IDs to labels.

The report lists clusters by ID and shows your label for each cluster as its **User Label**; you can also switch the cluster list to your labels. When your labels are names rather than numbers or generic terms such as `Unknown`, the reviewer compares CyteType's label with yours, and the report shows the match as a percentage. The result columns in your object are mapped back to your cells, so you only need the IDs for per-cluster tables.

(prepare-data-markers)=
## Marker genes

Your marker genes are the starting evidence for each cluster. CyteType does not compute them, so choose the method and thresholds that suit your study.

::::{tab-set}
:::{tab-item} Python
:sync: python
Run `sc.tl.rank_genes_groups()` on the cluster column and pass its key to `CyteType`:

```python
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", key_added="rank_genes_leiden")
annotator = CyteType(adata, group_key="leiden", rank_key="rank_genes_leiden")
```

- **`rank_key` must be the `key_added` you used.** Without `key_added`, Scanpy stores the results under `"rank_genes_groups"`, which is also the default `rank_key`. If the key is missing, `CyteType()` stops and asks you to run `sc.tl.rank_genes_groups` first.
- **The `groupby` must equal `group_key`.** Otherwise `CyteType()` stops with an error such as `` `rank_genes_groups` run with groupby='louvain', expected 'leiden'. ``
- **Marker names must be in `adata.var_names`.** Markers that are not found are dropped with a warning. This happens when the ranking used `adata.raw` and `adata.raw` holds more genes than `adata.X`.
- **`n_top_genes`** (default 50) sets how many genes per cluster CyteType takes, in the order of the ranking.

Any ranking method works. For an AnnData opened in backed mode, use [`rank_genes_groups_backed()`](large-datasets.md#large-datasets-backed).
:::
:::{tab-item} R
:sync: r
Pass the table from `FindAllMarkers()`. It must follow these rules:

- **Columns.** It needs `cluster`, `gene` and `avg_log2FC`, as `FindAllMarkers()` returns them.
- **Exactly the same clusters.** The `cluster` values must be the same set as the labels in `group_key`. Compute the markers with `Idents(obj)` set to that column, and do not filter out a whole cluster. Otherwise `PrepareCyteTypeR()` stops with `Please check if cluster labels are consistent between marker table and seurat obj!`.
- **Top genes by fold change.** For each cluster, CyteTypeR keeps the `n_top_genes` genes (default 50) with the highest `avg_log2FC`. It does not filter by p-value, so filter the table first if you want a threshold. `only.pos = TRUE` keeps the table to genes that are higher in the cluster.
- **At least five markers per cluster** after that step. Otherwise `PrepareCyteTypeR()` stops with `Invalid marker genes, some clusters have fewer than 5 markers`. Keep `n_top_genes` at 5 or more, and check that strict filters leave five genes in every cluster.

```r
markers <- FindAllMarkers(obj, only.pos = TRUE)
markers <- markers[markers$p_val_adj < 0.05, ]
table(markers$cluster)  # every cluster needs at least five rows
```
:::
::::

The server uses at most 250 marker genes per cluster; any beyond that are ignored.

(prepare-data-gene-symbols)=
## Gene symbols

CyteType reasons about genes by name, so it needs gene symbols such as `CD3E` or `Cd3e`, not database identifiers such as Ensembl IDs.

::::{tab-set}
:::{tab-item} Python
:sync: python
With `gene_symbols_column=None` (the default), CyteType finds the symbols itself. It scores these candidates:

1. the columns `feature_name`, `gene_symbols`, `gene_symbol`, `gene_short_name`, `gene_name` and `symbol` of `adata.var`, when present;
2. `adata.var_names`;
3. any other column of `adata.var` with enough distinct values.

It picks the candidate with the lowest share of values that look like Ensembl, RefSeq, Entrez or similar accession IDs; on a tie, the named columns win, then `adata.var_names`. The log tells you which one it used, for example `Auto-detected gene symbols in column 'gene_name'`. To choose yourself, pass the column name, for example `gene_symbols_column="gene_name"`.

- **Combined names are split.** When most values combine a symbol and an ID, such as `TSPAN6_ENSG00000000003` (joined by `_` or `|`), CyteType keeps the symbol.
- **IDs stop the preparation.** If about half or more of the chosen names still look like IDs, `CyteType()` stops with `Gene Symbol Detection Error`. Point `gene_symbols_column` at a column with symbols, or convert the IDs to symbols first.
- **A temporary column.** CyteType stores the cleaned symbols in a column named `__cytetype_gene_symbols` in `adata.var`. `cleanup()` removes it.
:::
:::{tab-item} R
:sync: r
CyteTypeR uses the row names of the `RNA` assay as gene names: the expression summaries and the uploaded matrix use them, and so does the marker table that `FindAllMarkers()` returns. If your row names are IDs, convert them to symbols before you build the object. `Read10X()` reads symbols by default.

CyteTypeR checks the genes but only warns: when more than half of the genes in the marker table, or of the gene names of the `RNA` assay, look like database IDs (Ensembl IDs with or without a version suffix, RefSeq accessions or numeric IDs), it logs a warning such as `100.0% of the genes in the marker table look like gene IDs` and continues.

The `gene_symbols` argument (default `"gene_symbols"`) names a column of the `RNA` assay's feature metadata to check instead of the row names, and changes nothing else. When the column does not exist, you see the message `Gene symbols not found in meta.features, using rownames`, which is expected.
:::
::::

## Coordinates

Coordinates are optional. They draw the embedding in the report.

::::{tab-set}
:::{tab-item} Python
:sync: python
`coordinates_key` (default `"X_umap"`) names an embedding in `adata.obsm`. If it is missing or unusable, CyteType tries `X_umap`, `X_tsne` and `X_pca` in turn, ignoring case. It uses the first array that has one row per cell and at least two columns, and takes its first two columns. If none qualifies, the job still runs, but the report has no embedding.
:::
:::{tab-item} R
:sync: r
`coordinates_key` (default `"umap"`) names a reduction of the Seurat object, such as `"umap"` or `"tsne"`. If it is missing or unusable, CyteTypeR tries `umap`, `tsne` and `pca` in turn, ignoring case. It uses the first reduction that has one row per cell and at least two dimensions, and takes its first two. If none qualifies, the job still runs, but the report has no embedding.
:::
::::

`max_cells_per_group` (default 1000) sets how many randomly chosen cells per cluster are drawn in the report's embedding. The coordinates of all cells are also saved in `obs.duckdb`.

(prepare-data-raw-counts)=
## Raw counts

Raw counts are optional. When they are found, they are uploaded in `vars.h5` next to the normalized values, and the server can use them for analyses that work on counts. [Target analysis](../platform/target-analysis.md) requires them.

::::{tab-set}
:::{tab-item} Python
:sync: python
CyteType looks in this order and uses the first source that holds whole numbers:

1. `adata.layers["counts"]`;
2. `adata.raw.X`, but only if `adata.raw` has every gene in `adata.var_names`;
3. `adata.X`.

If none qualifies, CyteType logs a warning and uploads the normalized values only. To keep the counts, copy them into a layer before you normalize:

```python
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```
:::
:::{tab-item} R
:sync: r
CyteTypeR uploads the `counts` layer of the `RNA` assay, which Seurat keeps when you normalize, if every value in it is a whole number. Otherwise it logs `The raw counts are not integer-valued. Skipping raw counts in vars.h5.` and uploads the normalized values only.
:::
::::

(prepare-data-metadata)=
## Cell metadata

CyteType summarizes your categorical cell metadata per cluster, for example which conditions, donors or tissue regions a cluster's cells come from. [Study context and metadata](study-context.md#study-context-cell-metadata) describes which columns help.

- **Which columns.** In Python, columns of type category, object or string; in R, factor and character columns. Numeric and logical columns are not summarized, and neither is the cluster column itself.
- **`min_percentage`** (default 10). For each cluster, only the values held by more than this percentage of its cells (in whole percentages) are kept.
- **`max_metadata_categories`** (default 500). Columns with more distinct values than this, such as cell barcodes, are skipped.
- **`aggregate_metadata=False`** turns the summaries off.

The summaries are separate from `obs.duckdb`, which always contains your complete cell metadata table, numeric columns included. Turning the summaries off does not remove anything from that file. If some columns must not leave your computer, for example donor identifiers, remove them from a copy of the object before you prepare the job. [Security and privacy](../support/security-and-privacy.md) describes what is uploaded and who can access it.

## Species

CyteType supports human and mouse data. Name the organism in your [study context](study-context.md): CyteType identifies it from that text and the marker genes. Other organisms are not supported: their results have not been validated, and a job can fail when CyteType's reference does not cover the organism.

(prepare-data-files)=
## Files CyteType writes

Preparing a job writes two files, which `run()` and `CyteTypeR()` upload with the job, and submitting it writes a copy of the request:

| File | Contents | Written by | Location option |
|---|---|---|---|
| `vars.h5` | The normalized expression matrix, raw counts when found, and the gene metadata (`adata.var`, or the assay's feature metadata in R). | `CyteType()`, `PrepareCyteTypeR()` | `vars_h5_path` |
| `obs.duckdb` | Every column of your cell metadata, plus the 2D coordinates of all cells. | `CyteType()`, `PrepareCyteTypeR()` | `obs_duckdb_path` |
| `query.json` | The request: study context, metadata tags, cluster labels, marker genes, metadata summaries, the sampled coordinates and the expression percentages. | `run()`, `CyteTypeR()` | `query_filename`; turn off with `save_query=False` |

- **Location.** All three default to the working directory; relative paths are resolved against it. Set the options above to write somewhere else, for example a scratch disk (see [Automation, HPC and CI](automation.md#automation-files)).
- **Keys in `query.json`.** The Python client never writes your LLM keys to `query.json`. In CyteTypeR the file also contains `llm_configs`, including any API keys, so use `save_query = FALSE` when you pass keys (see [Use your own LLM](llm-configuration.md#llm-configuration-keys)).
- **Removing the files.** The upload files stay on disk until you delete them, because every job you submit from the same prepared object uploads them again. Delete them with `annotator.cleanup()` in Python or `CleanUpArtifacts(prepped_data)` in R once you no longer need them. After that, a new job needs a new preparation (see [Run, monitor and rerun jobs](running-jobs.md#running-jobs-rerun)).

## Next steps

- [Study context and metadata](study-context.md): describe your study and choose the metadata to pass with the job.
- [Run, monitor and rerun jobs](running-jobs.md): submit the job and follow it.
- [Large datasets](large-datasets.md): prepare data that does not fit in memory or takes long to upload.
