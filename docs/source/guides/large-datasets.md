---
description: Annotate large single-cell datasets with CyteType by working from backed AnnData, ranking marker genes in chunks, capping cells per cluster, and planning memory, disk space and uploads.
---

(large-datasets)=
# Large datasets

The annotation runs on CyteType's servers, but the preparation runs in your session. It reads the whole expression matrix to summarize the clusters, writes `vars.h5` and `obs.duckdb` to disk, and uploads them with every job. For a large dataset, plan for memory, local disk space and upload size.

| Concern | Python | R |
|---|---|---|
| The matrix does not fit in memory | [Open the AnnData in backed mode](#large-datasets-backed) and rank markers with `rank_genes_groups_backed()` | [Cap the cells per cluster](#large-datasets-subsample) |
| The upload files are too large | Cap the cells per cluster with `subsample_by_group()`; drop metadata columns you do not need | Cap the cells per cluster with Seurat; drop metadata columns |
| Memory peaks during preparation | Lower `pcent_batch_size` | Lower `pcent_batch_size` |
| Slow or unreliable uploads | `upload_max_workers`, `upload_timeout_seconds` | `upload_max_workers`, `upload_timeout_seconds` |

(large-datasets-backed)=
## Keep the matrix on disk

::::{tab-set}
:::{tab-item} Python
:sync: python
When you open an `.h5ad` file with `backed="r"`, `adata.X` stays on disk. CyteType then reads it in chunks of cells, both for the expression summaries and when it writes `vars.h5`, so the matrix is never loaded whole.

```python
import scanpy as sc
import cytetype

# adata.X in the file must already hold log-normalized values.
adata = sc.read_h5ad("atlas.h5ad", backed="r")

# Rank marker genes by streaming the matrix in chunks of cells.
cytetype.rank_genes_groups_backed(adata, groupby="leiden", key_added="rank_genes_leiden", n_genes=100)

annotator = cytetype.CyteType(adata, group_key="leiden", rank_key="rank_genes_leiden")
adata = annotator.run(study_context="...")
```

Keep in mind:

- **Only `adata.X` and `adata.raw.X` stay on disk.** AnnData reads everything else into memory, including layers. A large `adata.layers["counts"]` is therefore loaded in full; raw counts kept in `adata.raw` stay on disk.
- **Every cell needs a cluster label,** and every cluster at least two cells, for `rank_genes_groups_backed()`.
- **The results are added in memory.** `run()` adds the result columns to `adata.obs` and the results to `adata.uns`, but the file you opened is not changed. Write a new file, or save only the new columns:

  ```python
  annotator.cleanup()
  adata.write_h5ad("atlas_cytetype.h5ad")

  # Or keep only the per-cell result columns.
  columns = [c for c in adata.obs.columns if c.startswith("cytetype_")]
  adata.obs[columns].to_csv("atlas_cytetype_labels.csv")
  ```
:::
:::{tab-item} R
:sync: r
CyteTypeR works on Seurat objects in memory and cannot read the expression matrix from disk in chunks. To reduce the memory it needs, [cap the number of cells per cluster](#large-datasets-subsample) before you prepare the job, and see [Memory](#large-datasets-memory).
:::
::::

### `rank_genes_groups_backed()` (Python)

[`cytetype.rank_genes_groups_backed()`](../reference/api/helpers.html#cytetype.rank_genes_groups_backed) computes the same ranking as Scanpy's `sc.tl.rank_genes_groups(method="t-test")` (Welch's t-test of each cluster against all other cells) in a single pass over the matrix. It reads a chunk of cells at a time, so it works on backed data without loading the matrix. It writes the results to `adata.uns[key_added]` in Scanpy's format, ready for `CyteType`.

| Argument | Default | Meaning |
|---|---|---|
| `adata` | | The AnnData object, backed or in memory. |
| `groupby` | | The cluster column. Use the same column as `group_key`. |
| `use_raw` | `False` | Rank `adata.raw.X` instead of `adata.X`. |
| `layer` | `None` | Rank a layer instead of `adata.X`. |
| `n_genes` | `None` | Genes to keep per cluster; `None` keeps all genes. |
| `key_added` | `"rank_genes_groups"` | Key in `adata.uns`. Pass the same value to `CyteType` as `rank_key`. |
| `cell_batch_size` | `5000` | Cells read per chunk. Lower it to reduce memory. |
| `corr_method` | `"benjamini-hochberg"` | Multiple-testing correction: `"benjamini-hochberg"` or `"bonferroni"`. |
| `pts` | `False` | Also store the fraction of cells in and outside each cluster that express each gene. |
| `copy` | `False` | Return a modified copy instead of writing into `adata`. |

With `n_genes=None`, the results hold a ranking of every gene for every cluster, which can take a lot of memory with many clusters. CyteType uses only the top `n_top_genes` genes per cluster (50 by default, at most 250 on the server), so `n_genes` between 100 and 250 is enough.

(large-datasets-subsample)=
## Cap the number of cells per cluster

Because CyteType annotates clusters, not cells, you can run it on a sample that keeps every cluster but caps the number of cells in each. This reduces memory, preparation time and the size of the upload files. In exchange, the expression summaries, the uploaded files and the report describe only the sampled cells.

::::{tab-set}
:::{tab-item} Python
:sync: python
`cytetype.subsample_by_group()` returns a new AnnData in memory with at most `max_cells_per_group` cells per cluster. Smaller clusters are kept whole, and cells without a cluster label are left out. It also works on a backed AnnData, and then reads only the sampled cells.

| Argument | Default | Meaning |
|---|---|---|
| `adata` | | The AnnData object, backed or in memory. |
| `group_key` | | The cluster column. |
| `max_cells_per_group` | `1000` | The maximum number of cells kept per cluster. |
| `random_state` | `0` | Seed for the random sample. |

The sample does not include `adata.uns`, so copy the marker results from the full data, or rank the markers again on the sample:

```python
small = cytetype.subsample_by_group(adata, group_key="leiden", max_cells_per_group=2000)

# Reuse the markers ranked on all cells.
small.uns["rank_genes_leiden"] = adata.uns["rank_genes_leiden"]

annotator = cytetype.CyteType(small, group_key="leiden", rank_key="rank_genes_leiden")
small = annotator.run(study_context="...")
```

The labels apply to whole clusters, so you can copy them to the full dataset through the cluster column:

```python
labels = dict(zip(small.obs["leiden"].astype(str), small.obs["cytetype_annotation_leiden"].astype(str)))
adata.obs["cytetype_annotation_leiden"] = adata.obs["leiden"].astype(str).map(labels).astype("category")
```
:::
:::{tab-item} R
:sync: r
CyteTypeR works on Seurat objects in memory and has no sampling function. Seurat's `subset()` with `downsample` keeps at most that many cells per identity class, so set the identities to your cluster column first. Marker genes computed on the full object can be reused, because the sample keeps every cluster:

```r
Idents(obj) <- "seurat_clusters"
small <- subset(obj, downsample = 2000)

prepped_data <- PrepareCyteTypeR(small, markers, group_key = "seurat_clusters")
small <- CyteTypeR(small, prepped_data, study_context = "...")
```

The labels apply to whole clusters, so you can copy them to the full object through the cluster column:

```r
labels <- unique(small@meta.data[, c("seurat_clusters", "cytetype_annotation_seurat_clusters")])
obj$cytetype_annotation_seurat_clusters <- labels$cytetype_annotation_seurat_clusters[
  match(as.character(obj$seurat_clusters), as.character(labels$seurat_clusters))
]
```
:::
::::

(large-datasets-memory)=
## Memory

::::{tab-set}
:::{tab-item} Python
:sync: python
- **`pcent_batch_size`** (default 5000) is the number of cells read per chunk while CyteType computes the expression percentages. With a sparse `adata.X` the chunks are small. With a dense `adata.X`, each chunk is expanded into dense 8-byte arrays, so 5,000 cells with 30,000 genes need a few gigabytes; lower `pcent_batch_size`, or store `adata.X` as a sparse matrix.
- **Writing `vars.h5`** from a backed sparse matrix stored by rows (CSR, the usual layout in `.h5ad` files) uses up to about 4 GB of working memory. From a matrix in memory, CyteType copies batches of genes while it writes, on top of the matrix itself.
- **Writing `obs.duckdb`** uses up to 4 GB of memory and spills to temporary files beyond that (see [Disk space](#large-datasets-disk)).
:::
:::{tab-item} R
:sync: r
- **The Seurat object is in memory.** While it writes `vars.h5`, CyteTypeR also makes a transposed copy of the normalized matrix, and when the `RNA` assay is split into layers, it joins them into a copy first. Plan for more memory than the object itself. If that is too much, [cap the cells per cluster](#large-datasets-subsample).
- **`pcent_batch_size`** (default 5000) is, in CyteTypeR, the number of genes processed per batch while it computes the expression percentages. Lower it to reduce the memory of each batch.
- **Writing `obs.duckdb`** uses up to 4 GB of memory and spills to temporary files beyond that (see [Disk space](#large-datasets-disk)).
:::
::::

(large-datasets-disk)=
## Disk space

- **Both files exist at the same time,** from preparation until you call `cleanup()` (Python) or `CleanUpArtifacts()` (R). Write them to a disk with room by setting `vars_h5_path` and `obs_duckdb_path` when you prepare the job.
- **`vars.h5`** stores each non-zero value of the normalized matrix in 8 bytes (the value and its position), and the same again for the raw counts when they are included, before LZ4 compression. For example, 1 million cells with 2,000 detected genes each hold about 2 billion non-zero values, or about 16 GB per matrix before compression.
- **`obs.duckdb`** holds one row per cell with every metadata column, plus two coordinate columns. Many long text columns make it larger.
- **Temporary files.** When the cell metadata table does not fit in the memory set aside for it, DuckDB spills to a temporary directory while it writes `obs.duckdb`. Both clients create a private directory under the system temporary directory (in R, under the session's `tempdir()`) and delete it afterwards. To move it, for example to a scratch disk, set `TMPDIR` before you start Python or R.

(large-datasets-upload)=
## Upload limits

| File | Largest upload |
|---|---|
| `obs.duckdb` | 2 GiB |
| `vars.h5` | 50 GiB |

The client checks the size of each file before uploading it. A file over the limit stops the run. With `require_artifacts=False`, the job runs without that file instead; CyteTypeR then runs without either file (see [Run, monitor and rerun jobs](running-jobs.md#running-jobs-options)). To bring the files under the limits:

- [cap the cells per cluster](#large-datasets-subsample);
- remove cell metadata columns you do not need, especially long text columns, which shrinks `obs.duckdb`;
- in Python, leave out the raw counts by removing the source CyteType would find (see [Raw counts](prepare-data.md#prepare-data-raw-counts)), which roughly halves `vars.h5`. Raw counts are optional for annotation, but [target analysis](../platform/target-analysis.md) and some cell-level evidence need them.

## Upload settings

Each file is uploaded in chunks, several at a time. The files are uploaded again for every job, including a rerun from the same prepared object.

| Option | Python default | R default | What it does |
|---|---|---|---|
| `upload_max_workers` | `4` | `6` | The number of chunks uploaded in parallel. Each worker holds one chunk in memory; the server sets the chunk size, currently up to 100 MB. |
| `upload_timeout_seconds` | `3600` | `3600` | The time limit, in seconds, for each chunk request. |

- **Retries.** A chunk that fails with a network error, a timeout or a server error is retried up to three times, after 1, 5 and 20 seconds.
- **Slow or shared connections.** Use fewer workers. The default timeout is generous; raise it only if chunk uploads time out.
- **Parallel uploads in R** need the furrr package, and also progressr when `show_progress = TRUE` (see [Installation](../getting-started/installation.md#installation-r-notes)). Each worker is a background R session. Without them, CyteTypeR uploads one chunk at a time.

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
adata = annotator.run(
    study_context="...",
    upload_max_workers=8,
    upload_timeout_seconds=7200,
)
```
:::
:::{tab-item} R
:sync: r
```r
obj <- CyteTypeR(
  obj, prepped_data,
  study_context = "...",
  upload_max_workers = 8L,
  upload_timeout_seconds = 7200L
)
```
:::
::::

Large jobs also take longer on the server. [Run, monitor and rerun jobs](running-jobs.md) explains how long the client waits and how to fetch the results later, and [Automation, HPC and CI](automation.md) covers batch jobs with a time limit.
