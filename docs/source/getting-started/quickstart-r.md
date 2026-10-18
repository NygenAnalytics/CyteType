---
description: Annotate the clusters of the 10x PBMC 3k dataset with CyteTypeR. Preprocess with Seurat, run CyteType, and explore the labels, the evidence and the report.
---

(quickstart-r)=
# Quickstart (R)

This page annotates the clusters of a small public dataset, 2,700 peripheral blood mononuclear cells (PBMCs) from 10x Genomics, with Seurat and CyteTypeR. The preprocessing runs on your computer; the annotation itself runs on CyteType's servers. The code on this page is not run when the site is built, so no outputs are shown.

Before you start, you need:

- R 4.2 or newer with CyteTypeR installed (see [Installation](installation.md)).
- A CyteType sign-in on this computer. Run `SetupCyteTypeR()` once (see [Sign in and API keys](sign-in.md)).

Working in Python? The same workflow with Scanpy is in the [Python quickstart](quickstart.ipynb).

## Load and cluster the data

CyteType annotates the clusters you define, using the marker genes you compute for them. These are the standard Seurat steps from the PBMC 3k tutorial: quality control, normalization, clustering and marker detection.

```r
# Load Seurat for preprocessing, CyteTypeR for annotation and dplyr for tables.
library(Seurat)
library(CyteTypeR)
library(dplyr)
```

```r
# Download and unpack the 10x Genomics PBMC 3k filtered matrices (about 7.6 MB).
url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
download.file(url, destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz", mode = "wb")
untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")

# Read the counts; the row names are gene symbols.
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
```

```r
# Drop cells with unusual numbers of genes or a high fraction of mitochondrial reads.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Log-normalize; CyteTypeR uploads the normalized values together with the raw counts.
pbmc <- NormalizeData(pbmc)
```

```r
# Find variable genes, scale, reduce dimensions and build the neighbor graph.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Cluster the cells and compute a UMAP for plotting and for the report.
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
```

```r
# Find the marker genes of every cluster; CyteTypeR reads them from this table.
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
```

## Prepare the annotation

`PrepareCyteTypeR()` checks the inputs and summarizes each cluster: its top marker genes, how widely each gene is expressed, and the cell metadata. It also writes two files to the working directory, `vars.h5` (expression) and `obs.duckdb` (cell metadata and UMAP coordinates), which are uploaded with the job. [Prepare your data](../guides/prepare-data.md) explains each requirement.

```r
# Check the inputs and build the files that are uploaded with the job.
prepped_data <- PrepareCyteTypeR(
  obj = pbmc,
  marker_table = pbmc.markers,
  group_key = "seurat_clusters",
  coordinates_key = "umap"
)
```

A few things to know about this step:

- For each cluster, CyteTypeR keeps the 50 marker genes with the highest `avg_log2FC` (set `n_top_genes` to change this). It does not filter by p-value, so filter `pbmc.markers` first if you want a threshold.
- Every cluster needs at least five marker genes, and the marker table must contain exactly the clusters in `seurat_clusters`. Otherwise `PrepareCyteTypeR()` stops with an error.
- The message `Gene symbols not found in meta.features, using rownames` is expected here, because the gene symbols are the row names of the object.

## Run the annotation

`CyteTypeR()` uploads the files, submits the job and waits for the results. The report link is printed as soon as the job is submitted, so you can follow its progress in the browser. If your R session ends before the job finishes, the job keeps running on the server (see [Run, monitor and rerun jobs](../guides/running-jobs.md)).

The study context tells CyteType what the sample is. Include the organism, tissue and condition; [Study context and metadata](../guides/study-context.md) has more examples. The `title` entry in `metadata` names the job in your dashboard.

```r
# Describe the study, submit the job and wait for it to finish.
pbmc <- CyteTypeR(
  obj = pbmc,
  prepped_data = prepped_data,
  study_context = paste(
    "Peripheral blood mononuclear cells (PBMCs) from a healthy human donor,",
    "profiled with 10x Genomics 3' single-cell RNA-seq."
  ),
  metadata = list(title = "PBMC 3k quickstart (R)")
)
```

Always assign the result, as in `pbmc <- CyteTypeR(pbmc, ...)`. `CyteTypeR()` returns a copy of the object with the results added and leaves the object you passed in unchanged.

## Explore the results

Every cell now carries its cluster's annotation in `pbmc@meta.data`:

| Column | Content |
|---|---|
| `cytetype_annotation_seurat_clusters` | Cell type label |
| `cytetype_cellOntologyTerm_seurat_clusters` | Cell Ontology term name |
| `cytetype_cellOntologyTermID_seurat_clusters` | Cell Ontology term ID, such as `CL:0000236` |
| `cytetype_cellState_seurat_clusters` | Cell state (empty when none was called) |
| `cytetype_reviewerConfidence_seurat_clusters` | The reviewer's confidence: `High`, `Moderate` or `Low` |

CyteTypeR also adds three columns under older names, with the same label, term and ID: `cytetype_seurat_clusters`, `cytetype_ontologyTerm_seurat_clusters` and `cytetype_ontologyID_seurat_clusters`.

```r
# Show one row per cluster with its label, Cell Ontology term ID, cell state and reviewer confidence.
pbmc@meta.data |>
  distinct(
    seurat_clusters,
    cytetype_annotation_seurat_clusters,
    cytetype_cellOntologyTermID_seurat_clusters,
    cytetype_cellState_seurat_clusters,
    cytetype_reviewerConfidence_seurat_clusters
  ) |>
  arrange(seurat_clusters)
```

```r
# Plot the CyteType labels on the UMAP.
DimPlot(pbmc, group.by = "cytetype_annotation_seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
```

The Python client's `marker_dotplot()` shows the genes that support each label. In R, the `supportingMarkers` column of the results table (see below) lists them, so you can plot them with Seurat's `DotPlot()`:

```r
# Plot the first three genes that support each label.
supporting <- strsplit(pbmc@misc$cytetype_results$supportingMarkers, "; ", fixed = TRUE)
genes <- intersect(unlist(lapply(supporting, head, 3)), rownames(pbmc))
DotPlot(pbmc, features = genes, group.by = "cytetype_annotation_seurat_clusters") + RotatedAxis()
```

## Read the results table

`pbmc@misc$cytetype_results` has one row per cluster. `clusterId` holds your Seurat cluster labels, and `GetResults(pbmc)` returns the same table.

```r
# Show the label, Cell Ontology term, granular label and cell state of every cluster.
results <- pbmc@misc$cytetype_results
results[, c("clusterId", "annotation", "ontologyTerm", "ontologyTermID", "granularAnnotation", "cellState")]
```

The table also holds the evidence behind each label (`justification`, `supportingMarkers`, `conflictingMarkers`, `missingExpression` and `unexpectedExpression`) and the reviewer's assessment (`reviewerConfidence`, `isHeterogeneous` and `heterogeneityType`). [Results in AnnData and Seurat](../results/results-in-your-data.md#results-per-cluster) describes every column.

```r
# Look at the evidence recorded for the first cluster.
as.list(results[1, c("clusterId", "annotation", "cellState", "reviewerConfidence", "supportingMarkers")])
```

## Open the report

The interactive report holds the full reasoning for every cluster: the candidate cell types that were considered, the marker evidence, the reviewer's confidence and quality-control checks, and pathway context. Your reports are private unless you share them (see [Visibility and sharing](../platform/sharing.md)).

```r
# Get the link to the interactive report.
pbmc@misc$cytetype_jobDetails$report_url

# Or open the report in your browser; you are asked to sign in if needed.
ViewCyteTypeJob(pbmc@misc$cytetype_jobDetails$job_id)
```

[The CyteType report](../results/report.md) describes every section.

## Save and clean up

The results are stored in `pbmc@meta.data` and `pbmc@misc`, so they are saved with the Seurat object. `CleanUpArtifacts()` deletes the two upload files.

```r
# Remove vars.h5 and obs.duckdb from the working directory.
CleanUpArtifacts(prepped_data)

# Save the annotated object; the CyteType results are saved with it.
saveRDS(pbmc, "pbmc3k_cytetype.rds")
```

The working directory also holds `query.json`, a copy of the request that was sent. You do not need it afterwards.

## Next steps

- [Prepare your data](../guides/prepare-data.md): requirements for clusters, marker genes and gene symbols.
- [The CyteType report](../results/report.md): every section of the report.
- [Results in AnnData and Seurat](../results/results-in-your-data.md): all columns and keys.
- [Interpret and review results](../results/interpreting-results.md): confidence, heterogeneity and limits.
- [Use your own LLM](../guides/llm-configuration.md): run CyteType with your own model provider.
