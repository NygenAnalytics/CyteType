---
description: CyteType annotates the clusters of single-cell RNA-seq data from Python or R and returns a cell type, a Cell Ontology term, a cell state and a confidence rating for each cluster, with a report of the evidence.
---

(home)=
# CyteType

CyteType annotates the clusters of single-cell RNA-seq data. You provide the clusters and
marker genes you have computed, from Python with Scanpy and AnnData or from R with Seurat,
and CyteType returns for each cluster a cell type label, a Cell Ontology term, a cell state
and a confidence rating from an independent review step. The annotation runs on CyteType's
servers, and an interactive report shows the evidence behind every label.

## A first run

Install the client, sign in once, then annotate a clustered dataset. Using CyteType
requires an API key: `cytetype setup` or `SetupCyteTypeR()` signs you in with your browser
and saves your key on this computer.

::::{tab-set}
:::{tab-item} Python
:sync: python
```bash
pip install cytetype
cytetype setup
```

```python
import scanpy as sc
from cytetype import CyteType

# adata holds log-normalized expression and Leiden clusters in adata.obs["leiden"].
sc.tl.rank_genes_groups(adata, groupby="leiden", key_added="rank_genes_leiden")

annotator = CyteType(adata, group_key="leiden", rank_key="rank_genes_leiden")
adata = annotator.run(
    study_context="Peripheral blood mononuclear cells from a healthy human donor."
)
```
:::
:::{tab-item} R
:sync: r
```r
install.packages("devtools")
devtools::install_github("NygenAnalytics/CyteTypeR")
CyteTypeR::SetupCyteTypeR()
```

```r
library(Seurat)
library(CyteTypeR)

# obj holds normalized data in its RNA assay, a UMAP and clusters from
# FindClusters().
markers <- FindAllMarkers(obj, only.pos = TRUE)
prepped <- PrepareCyteTypeR(obj, markers, group_key = "seurat_clusters")
obj <- CyteTypeR(
  obj,
  prepped,
  study_context = "Peripheral blood mononuclear cells from a healthy human donor."
)
```
:::
::::

The link to the report is printed as soon as the job starts. The
[Python quickstart](getting-started/quickstart.ipynb) and the
[R quickstart](getting-started/quickstart-r.md) walk through a complete example.

## Where to go next

- **Getting started:** [install the client](getting-started/installation.md),
  [sign in](getting-started/sign-in.md), and read
  [how CyteType works](getting-started/how-cytetype-works.md).
- **Prepare and run:** [prepare your data](guides/prepare-data.md), work with
  [large datasets](guides/large-datasets.md), write a
  [study context](guides/study-context.md),
  [run, monitor and rerun jobs](guides/running-jobs.md),
  [use your own LLM](guides/llm-configuration.md), and run CyteType from
  [scripts, HPC and CI](guides/automation.md).
- **Results:** find the [results in AnnData and Seurat](results/results-in-your-data.md),
  read [the CyteType report](results/report.md), and
  [interpret and review the results](results/interpreting-results.md).
- **Platform:** manage your jobs in the [dashboard](platform/dashboard.md), control
  [visibility and sharing](platform/sharing.md), ask questions about your results in
  [Chat](platform/chat.md), run a [target analysis](platform/target-analysis.md),
  [connect AI assistants](platform/mcp.md) through MCP, and check the
  [usage limits](platform/limits.md).
- **Reference:** the [Python API](reference/api/index.html), the
  [R functions](reference/r-api.md), the [command line](reference/cli.md) and the
  [error codes](reference/errors.md).
- **Support:** [troubleshooting](support/troubleshooting.md), the [FAQ](support/faq.md),
  [security and privacy](support/security-and-privacy.md) and
  [how to cite CyteType](support/citing.md).

## Source code and issues

The Python client is developed in
[NygenAnalytics/CyteType](https://github.com/NygenAnalytics/CyteType) and the R package in
[NygenAnalytics/CyteTypeR](https://github.com/NygenAnalytics/CyteTypeR). Report problems in
the [CyteType issues](https://github.com/NygenAnalytics/CyteType/issues) or the
[CyteTypeR issues](https://github.com/NygenAnalytics/CyteTypeR/issues). The
[documentation index for LLMs](https://docs.nygen.io/CyteType/llms.txt) lists every page in
plain text for assistants that read this site.
