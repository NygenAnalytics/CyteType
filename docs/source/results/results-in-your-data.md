---
description: What CyteType adds to your AnnData or Seurat object, what get_results() and GetResults() return, and how to save, reload, plot and export the results.
---

(results-in-your-data)=
# Results in AnnData and Seurat

When a job finishes, the client writes the results into your object: one label per cell in the cell metadata, and the per-cluster results next to it. The [interactive report](report.md) holds the full evidence for every cluster; this page covers what lands in your data and how to use it.

The examples use the cluster column `leiden` in Python and `seurat_clusters` in R, with the default `results_prefix = "cytetype"`. Column and key names change with both.

(results-return-value)=
## Keep the returned object

::::{tab-set}
:::{tab-item} Python
:sync: python
`run()` returns the annotated AnnData. It is the object you passed to `CyteType()`, except when `label_na=True` relabeled missing cluster labels: then it is a new object, and only that object carries the results. Always assign the return value.

```python
adata = annotator.run(study_context="Healthy human PBMCs, 10x Genomics 3' scRNA-seq")
```
:::
:::{tab-item} R
:sync: r
`CyteTypeR()` returns a modified copy of your Seurat object; R does not change the object you pass in. Assign the result. If you don't, the labels are not stored in your session, although the report still has them and `GetResults(job_id = ...)` can fetch them.

```r
obj <- CyteTypeR(obj, prepped_data, study_context = "Healthy human PBMCs, 10x Genomics 3' scRNA-seq")
```

If `CyteTypeR()` stops before it returns, because `timeout_seconds` ran out or you interrupted it, the job keeps running on the server but your object gets no results. Fetch them with `GetResults(job_id = "<job_id>")`, using the job ID from the message (see [Fetch results later](../guides/running-jobs.md#running-jobs-fetch-later)).
:::
::::

(results-cell-columns)=
## Cell-level columns

Annotation happens per cluster, so every cell in a cluster gets the same values.

::::{tab-set}
:::{tab-item} Python
:sync: python
`run()` adds five categorical columns to `adata.obs`, named `{results_prefix}_<field>_{group_key}`:

| Column | Content |
|---|---|
| `cytetype_annotation_leiden` | The cluster's cell type label, for example `Classical Monocyte`. |
| `cytetype_cellOntologyTerm_leiden` | The Cell Ontology term name for the label, for example `classical monocyte`. Can be empty when no term matched. |
| `cytetype_cellOntologyTermID_leiden` | The Cell Ontology ID, for example `CL:0000860`. Can be empty when no term matched. |
| `cytetype_cellState_leiden` | The cell state, for example `proliferating`. Empty when no state was called. |
| `cytetype_reviewerConfidence_leiden` | The reviewer's confidence: `High`, `Moderate` or `Low`, or `Unknown` when the result has none. |

The categories of the confidence column are listed as `High`, `Moderate`, `Low`, `Unknown` (only those that occur), so legends follow that order; the categorical is not ordered, so comparisons such as `>` do not apply. Cells in clusters without a result get `Unknown` (an empty value for the cell state).
:::
:::{tab-item} R
:sync: r
`CyteTypeR()` adds factor columns to `obj@meta.data`. The first five use the same names as the Python client; the last three are legacy copies kept for scripts written against earlier versions:

| Column | Content |
|---|---|
| `cytetype_annotation_seurat_clusters` | The cluster's cell type label. |
| `cytetype_cellOntologyTerm_seurat_clusters` | The Cell Ontology term name. |
| `cytetype_cellOntologyTermID_seurat_clusters` | The Cell Ontology ID. |
| `cytetype_cellState_seurat_clusters` | The cell state. Empty when no state was called. |
| `cytetype_reviewerConfidence_seurat_clusters` | The reviewer's confidence: `High`, `Moderate` or `Low`, or `Unknown` when the result has none. |
| `cytetype_seurat_clusters` | Legacy copy of the label column. |
| `cytetype_ontologyTerm_seurat_clusters` | Legacy copy of the term name column. |
| `cytetype_ontologyID_seurat_clusters` | Legacy copy of the term ID column. |

Prefer the new names. The levels of the confidence column are `High`, `Moderate`, `Low`, `Unknown` (only those that occur), so legends follow that order; the factor is not ordered, so comparisons such as `>` do not apply. Cells in clusters without a result get `Unknown` (an empty value for the cell state) in the first five columns and `NA` in the legacy columns.
:::
::::

Cell Ontology IDs use the `CL:0000860` form (occasionally `PCL:`, the Provisional Cell Ontology). Results from older jobs can use the `CL_0000860` form.

(results-cluster-ids)=
## Cluster IDs and your labels

The cell-level columns are already mapped to cells. Per-cluster results, the report and its CSV export identify clusters by CyteType's cluster ID, which can differ from the labels in your cluster column. In the report and the CSV, `cluster_id` is CyteType's ID and `author_label` is your label.

::::{tab-set}
:::{tab-item} Python
:sync: python
If every label in the cluster column has at most three characters (for example Leiden numbers), CyteType uses the labels themselves as IDs. Otherwise it numbers the clusters `0`, `1`, `2`, ... in natural sort order of the labels. With `label_na=True`, the added label `Unknown` is longer than three characters, so all clusters are then numbered. `annotator.cluster_map` maps each of your labels to its ID:

```python
id_to_label = {cluster_id: label for label, cluster_id in annotator.cluster_map.items()}
```
:::
:::{tab-item} R
:sync: r
CyteTypeR always numbers clusters `1`, `2`, `3`, ... in the sorted order of your labels (level order for a factor such as `seurat_clusters`), so Seurat cluster `0` is CyteType cluster `1`. `prepped_data$clusterLabels` and `obj@misc$cytetype_jobDetails$cluster_labels` map each ID to your label. The table in `obj@misc$cytetype_results` already uses your labels.

```r
id_to_label <- unlist(prepped_data$clusterLabels)
id_to_label[["1"]]  # your label for CyteType cluster 1
```
:::
::::

(results-per-cluster)=
## Per-cluster results

::::{tab-set}
:::{tab-item} Python
:sync: python
`get_results()` returns a dictionary:

```python
results = annotator.get_results()
```

| Key | Content |
|---|---|
| `annotations` | A list with one dictionary per annotated cluster (fields below). |
| `clusterCategories` | Display groups: a list of dictionaries with `categoryName`, `clusterIds` and `color`. They define the report's Biological Order and the grouping in `marker_dotplot()`. |
| `studyContext` | CyteType's description of the study, shown as Study Overview in the report. |
| `summary` | A study-level assessment from earlier pipeline versions; empty for current jobs. |
| `raw_annotations` | The server's per-cluster results, unchanged (see [Raw results](#results-raw)). |

Each entry in `annotations` has these keys:

| Key | Type | Content |
|---|---|---|
| `clusterId` | str | CyteType's cluster ID. |
| `annotation` | str | The cell type label (as in the obs annotation column). |
| `ontologyTerm` | str | The Cell Ontology term name. |
| `ontologyTermID` | str | The Cell Ontology ID. |
| `granularAnnotation` | str | A one-sentence description that combines the cell type, its state and key genes. |
| `cellState` | str | The cell state; empty when none was called. |
| `justification` | str | A summary of the evidence for the label. |
| `supportingMarkers` | list of str | Genes whose expression supports the label. |
| `conflictingMarkers` | list of str | Genes expressed in the cluster that are inconsistent with the label. |
| `missingExpression` | str | Genes expected for the label that are weak or absent, with an explanation. |
| `unexpectedExpression` | str | An explanation of the inconsistent genes. |
| `reviewerConfidence` | str or None | `High`, `Moderate` or `Low`. |
| `isHeterogeneous` | bool or None | Whether the reviewer judged that the cluster contains distinct populations. |
| `heterogeneityType` | str or None | `distinct_populations`, `mixed_lineage`, `doublet_artifact`, `transitional_state` or `unresolved`. None for homogeneous clusters and when the type was not classified. |

[Interpret and review results](interpreting-results.md) explains how to read the confidence and heterogeneity values.

`adata.uns` holds two entries:

| Key | Content |
|---|---|
| `cytetype_results` | `{"job_id": ..., "result": ...}`, where `result` is the `get_results()` dictionary stored as a JSON string so that the AnnData can be written to h5ad. |
| `cytetype_jobDetails` | `job_id`, `report_url` and `api_url`, plus `visibility` when you passed it to `run()`. |
:::
:::{tab-item} R
:sync: r
`CyteTypeR()` stores two entries in `obj@misc`:

| Entry | Content |
|---|---|
| `obj@misc$cytetype_results` | A data.frame with one row per annotated cluster (columns below). The row names and `clusterId` are your cluster labels. Its attributes `"job_id"` and `"result"` hold the job ID and the full results. |
| `obj@misc$cytetype_jobDetails` | `job_id`, `report_url`, `api_url`, `group_key` and `cluster_labels` (CyteType ID to your label), plus `visibility` when you passed it to `CyteTypeR()`. |

```r
res <- obj@misc$cytetype_results
full <- attr(res, "result")
```

The data.frame has these columns:

| Column | Type | Content |
|---|---|---|
| `clusterId` | character | Your cluster label. |
| `annotation` | character | The cell type label (as in the meta.data annotation column). |
| `ontologyTerm` | character | The Cell Ontology term name. |
| `ontologyTermID` | character | The Cell Ontology ID. |
| `granularAnnotation` | character | A one-sentence description that combines the cell type, its state and key genes. |
| `cellState` | character | The cell state; empty when none was called. |
| `justification` | character | A summary of the evidence for the label. |
| `supportingMarkers` | character | Genes whose expression supports the label, separated by `"; "`. |
| `conflictingMarkers` | character | Genes expressed in the cluster that are inconsistent with the label, separated by `"; "`. |
| `missingExpression` | character | Genes expected for the label that are weak or absent, with an explanation. |
| `unexpectedExpression` | character | An explanation of the inconsistent genes. |
| `reviewerConfidence` | character | `High`, `Moderate` or `Low`; `NA` when missing. |
| `isHeterogeneous` | logical | Whether the reviewer judged that the cluster contains distinct populations; `NA` when missing. |
| `heterogeneityType` | character | `distinct_populations`, `mixed_lineage`, `doublet_artifact`, `transitional_state` or `unresolved`. `NA` for homogeneous clusters and when the type was not classified. |

Clusters without a result have no row. [Interpret and review results](interpreting-results.md) explains how to read the confidence and heterogeneity values.

The full results are a list with `annotations` (one list per annotated cluster with the fields of the table, where `clusterId` is CyteType's ID and `supportingMarkers` and `conflictingMarkers` are lists of genes), `raw_annotations` (the server's per-cluster results, unchanged; see [Raw results](#results-raw)) and the server's `studyContext`, `summary`, `semanticOrder`, `clusterCategories` and `displayLevel`.

A run with another prefix, for example `results_prefix = "run2"`, stores `obj@misc$run2_results` (a list with `job_id` and `result`, the full results) and `obj@misc$run2_jobDetails`, and leaves `obj@misc$cytetype_results` as it is.

`GetResults()` returns different things depending on how you call it:

| Call | Returns |
|---|---|
| `GetResults(obj)` after a run with the default prefix | The data.frame in `obj@misc$cytetype_results`. |
| `GetResults(obj, results_prefix = "run2")` after a run with that prefix | The full results. |
| `GetResults(job_id = "<job_id>")` | For a completed job, a data.frame with the same columns and attributes as `obj@misc$cytetype_results`, where `clusterId` is CyteType's ID. It also saves the server's response as `cytetypeR_results_<job_id>.json` in the working directory, unless you pass `save_json = FALSE`. While the job is still running, and for a job ID the server does not know, it returns `NULL` with a message; a failed job is an error. |
:::
::::

(results-evidence-fields)=
:::{note}
**Evidence fields for current jobs.** The server stores the evidence for current jobs in a structured form inside the raw results. The Python client 0.21.0 and CyteTypeR 0.6.0 fill the five evidence fields from it: `justification` from the conclusion for the label (the **Summary** in the report's Justification section), `supportingMarkers` from the key supporting genes, `conflictingMarkers` from the unexpected genes, and `missingExpression` and `unexpectedExpression` from the explanations (or from the gene lists when there is no explanation). Values the server stored are never overwritten. The server applies the same rules when it returns results, so earlier Python versions also receive these fields filled.
:::

`get_results()` returns the results stored in `adata.uns` without contacting the server. If there are none but `adata.uns` holds the job details, for example because `run()` stopped waiting, it checks the job once and, if the job has finished, downloads the results and adds the columns as `run()` would. It returns `None` when nothing is available yet, and raises `AuthenticationError` when it needs the server but finds no API key. [Run, monitor and rerun jobs](../guides/running-jobs.md) covers recovering results this way.

(results-raw)=
## Raw results

:::{warning}
The raw results follow the server's internal format, which can change between releases. Prefer the documented fields above in pipelines, and pin the client version when your code depends on a raw field.
:::

Some details exist only in the raw results, such as the annotation levels L1 to L5 (see [Identity, levels and state](interpreting-results.md#interpret-levels)), the classification and the reviewer's reasoning. They are in `raw_annotations` of the full results: `results["raw_annotations"]` in Python, a dictionary keyed by CyteType cluster ID, and `attr(obj@misc$cytetype_results, "result")$raw_annotations` in R, a named list. Each value has `latest` (the cluster's newest result) and `history` (earlier results for the same cluster, usually empty). The fields people use most:

| Path under `raw_annotations[cluster_id]` | Content |
|---|---|
| `latest.annotation.fullOutput.cellType.annotationLevels` | The hierarchy from `L1` (broadest) to `L5` (most specific). Each level has a `label`; levels can be missing. |
| `latest.annotation.fullOutput.cellType.classification` | How the label relates to the reference hierarchy: `catalogue`, `outOfCatalog` or `unresolved` for current jobs (see [Classification](interpreting-results.md#interpret-classification)). |
| `latest.annotation.fullOutput.cellType.keySupportingGenes`, `missingGenes`, `unexpectedGenes` | The gene lists behind the evidence fields. |
| `latest.annotation.fullOutput.cellStateGeneGroups` | Gene sets that support the cell state. |
| `latest.review.fullOutput.confidenceAssessment` | The reviewer's confidence `label` and its `reasoning`. |
| `latest.review.fullOutput.heterogeneityAssessment` | `label`, `reasoning`, `heterogeneityType` and `heterogeneityTypeReasoning`. |
| `latest.review.labelComparisonScore` | The similarity between CyteType's label and your cluster label, 0 to 100. `labelComparisonLlmScore` holds the fallback estimate when this is empty. |

::::{tab-set}
:::{tab-item} Python
:sync: python
Build a table of the annotation levels:

```python
import pandas as pd

rows = []
for cluster_id, cluster in results["raw_annotations"].items():
    annotation = (cluster.get("latest") or {}).get("annotation") or {}
    cell_type = (annotation.get("fullOutput") or {}).get("cellType") or {}
    levels = cell_type.get("annotationLevels") or {}
    rows.append({
        "clusterId": cluster_id,
        **{name: (levels.get(name) or {}).get("label", "") for name in ["L1", "L2", "L3", "L4", "L5"]},
    })
levels_table = pd.DataFrame(rows)
```
:::
:::{tab-item} R
:sync: r
Build a table of the annotation levels:

```r
raw <- attr(obj@misc$cytetype_results, "result")$raw_annotations

# The labels of the levels L1 to L5 of one cluster ("" when a level is missing).
level_labels <- function(cluster) {
  levels <- cluster$latest$annotation$fullOutput$cellType$annotationLevels
  vapply(paste0("L", 1:5), function(name) {
    label <- levels[[name]]$label
    if (is.null(label)) "" else label
  }, character(1))
}
levels_table <- data.frame(clusterId = names(raw), t(vapply(raw, level_labels, character(5))))
```

The names of `raw` are CyteType's cluster IDs. For a job fetched with `GetResults(job_id = ...)`, use `attr(res, "result")$raw_annotations`.
:::
::::

(results-save)=
## Save and reload

::::{tab-set}
:::{tab-item} Python
:sync: python
The results live in `adata.obs` and `adata.uns`, so they are saved with the AnnData. Call `annotator.cleanup()` first if you want to remove the upload files and the temporary gene symbol column that CyteType added to `adata.var`.

```python
import json
import anndata as ad

adata.write_h5ad("annotated.h5ad")

adata = ad.read_h5ad("annotated.h5ad")
results = json.loads(adata.uns["cytetype_results"]["result"])
report_url = adata.uns["cytetype_jobDetails"]["report_url"]
```

Parsing the JSON returns the results exactly as they were stored. Files written by 0.21.0 already contain `reviewerConfidence`, `isHeterogeneous` and `heterogeneityType` and the filled evidence fields. For files written by older versions, `get_results()` and `marker_dotplot()` in 0.21.0 normalize the stored results when they read them: they fill the evidence fields and add the three review keys (as `None` for files from version 0.7.0 and older, which have no raw results). The `reviewerConfidence` column is not added to `adata.obs` of older files. `get_results()` is a method of the `CyteType` object, so in a new session parsing the JSON is the simplest way to read the results.
:::
:::{tab-item} R
:sync: r
The columns and the `misc` entries are saved with the Seurat object:

```r
saveRDS(obj, "annotated.rds")

obj <- readRDS("annotated.rds")
res <- GetResults(obj)
obj@misc$cytetype_jobDetails$report_url
```
:::
::::

(results-plot)=
## Plot the results

::::{tab-set}
:::{tab-item} Python
:sync: python
Color the embedding by the label or by the reviewer's confidence:

```python
import scanpy as sc

sc.pl.umap(adata, color="cytetype_annotation_leiden")

confidence_colors = {"High": "#1a9850", "Moderate": "#fdae61", "Low": "#d73027", "Unknown": "#bdbdbd"}
sc.pl.umap(adata, color="cytetype_reviewerConfidence_leiden", palette=confidence_colors)
```

`marker_dotplot()` draws a Scanpy dot plot of the genes that support each label, with the labels ordered and the genes grouped by display group (the `clusterCategories`):

```python
from cytetype import marker_dotplot

marker_dotplot(adata, group_key="leiden")
marker_dotplot(adata, group_key="leiden", n_top_markers=5, standard_scale="var")
```

| Argument | Default | Meaning |
|---|---|---|
| `adata` | | The annotated AnnData; the results are read from `adata.uns`. |
| `group_key` | | The cluster column used for the run, for example `"leiden"`. |
| `results_prefix` | `"cytetype"` | Must match the `results_prefix` of the run. |
| `n_top_markers` | `3` | Number of key supporting genes shown per cluster. |
| `gene_symbols` | `None` | A column of `adata.var` with gene symbols, if `var_names` are IDs. It is passed to `sc.pl.dotplot`. |
| `**kwargs` | | Passed to `sc.pl.dotplot`, for example `standard_scale`, `cmap`, `use_raw`, `layer`, `show`, `save` or `return_fig`. |

Gene names are matched without regard to case, because CyteType reports upper-case symbols (for example `CD8A` for mouse `Cd8a`). Genes that are not found are skipped with a warning. Like Scanpy, the function uses `adata.raw` when it exists and you pass no `layer`. It needs Scanpy, and it raises an error when the results have no display groups (some older jobs) or no raw results (files from version 0.7.0 and older).
:::
:::{tab-item} R
:sync: r
Color the embedding by the label or by the reviewer's confidence:

```r
library(Seurat)

DimPlot(obj, group.by = "cytetype_annotation_seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()

confidence_colors <- c(High = "#1a9850", Moderate = "#fdae61", Low = "#d73027", Unknown = "#bdbdbd")
DimPlot(obj, group.by = "cytetype_reviewerConfidence_seurat_clusters", cols = confidence_colors)
```

CyteTypeR has no dot plot helper, but the `supportingMarkers` column of the results table lists the genes that support each label, so you can plot them with Seurat's `DotPlot()`:

```r
res <- obj@misc$cytetype_results
genes <- unique(unlist(lapply(strsplit(res$supportingMarkers, "; ", fixed = TRUE), head, 3)))

# CyteType reports upper-case symbols; match them to the case used in your object
genes <- rownames(obj)[match(toupper(genes), toupper(rownames(obj)))]
genes <- unique(genes[!is.na(genes)])

DotPlot(obj, features = genes, group.by = "cytetype_annotation_seurat_clusters") + RotatedAxis()
```

:::
::::

(results-export)=
## Export a per-cluster table

The report's **Download as CSV** button exports one row per cluster with the label, confidence, heterogeneity, supporting genes and annotation levels (see [Download as CSV](report.md#report-csv)). To build a table from your own object:

::::{tab-set}
:::{tab-item} Python
:sync: python
From `get_results()`, with your cluster labels added:

```python
import pandas as pd

results = annotator.get_results()
table = pd.DataFrame(results["annotations"])
for column in ["supportingMarkers", "conflictingMarkers"]:
    table[column] = table[column].apply(lambda genes: ";".join(genes) if isinstance(genes, list) else "")

id_to_label = {cluster_id: label for label, cluster_id in annotator.cluster_map.items()}
table.insert(1, "leiden", table["clusterId"].map(id_to_label))
table.to_csv("cytetype_clusters.csv", index=False)
```

From the obs columns, which also works on a reloaded file without the `CyteType` object:

```python
group_key = "leiden"
fields = ["annotation", "cellOntologyTerm", "cellOntologyTermID", "cellState", "reviewerConfidence"]
columns = [f"cytetype_{field}_{group_key}" for field in fields]

per_cluster = adata.obs.groupby(group_key, observed=True)[columns].first()
per_cluster["n_cells"] = adata.obs[group_key].value_counts()
per_cluster.to_csv("cytetype_clusters.csv")
```
:::
:::{tab-item} R
:sync: r
From the results table, which already uses your cluster labels:

```r
res <- obj@misc$cytetype_results
write.csv(res, "cytetype_clusters.csv", row.names = FALSE)
```
:::
::::

## See also

- [The CyteType report](report.md): the evidence behind each label.
- [Interpret and review results](interpreting-results.md): how to read confidence, heterogeneity and classifications.
- [Run, monitor and rerun jobs](../guides/running-jobs.md): fetch results after a disconnect and rerun with new settings.
- [Python API reference](../reference/api/index.html) and [R function reference](../reference/r-api.md).
