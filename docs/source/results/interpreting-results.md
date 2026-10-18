---
description: How to read CyteType labels, reviewer confidence, heterogeneity flags and classifications, what the results cannot tell you, what to record for reproducibility, and a checklist for reviewing a run.
---

(interpreting-results)=
# Interpret and review results

CyteType proposes a label for each cluster, together with the evidence for it and an independent review. Treat the labels as reviewed proposals: they are most useful when you check the evidence before you rely on them. This page explains each signal, what CyteType cannot tell you, and how to review a run. [How CyteType works](../getting-started/how-cytetype-works.md) describes the pipeline that produces them.

(interpret-label)=
## What a label is

- **One label per cluster.** CyteType annotates clusters, not cells. Every cell in a cluster gets the cluster's label, so a small population inside a larger cluster is not labeled on its own.
- **Identity, state and description are separate.** The label names the cell's identity. The cell state, such as proliferating or activated, is a separate field. The granular annotation is a one-sentence description that combines the two with key genes.
- **The label is chosen, then reviewed.** CyteType weighs competing candidate identities over several rounds and places the winner in its reference hierarchy of cell types, which is linked to the Cell Ontology. Once all clusters are annotated, it compares similar clusters and harmonizes their labels. An independent reviewer then assesses each final label.
- **The Cell Ontology term can be broader than the label.** The term is matched to the identity before label refinement, which can add a qualifier to tell apart clusters that share an identity. The qualifier is not part of the term.

(interpret-levels)=
## Identity, levels and state

The reference hierarchy runs from L1, the broadest level (a lineage), to L5, the most specific. The report shows the path as the **Annotation hierarchy** in Justification, the CSV export has L1 to L4, and the raw results have all five levels (see [Raw results](results-in-your-data.md#results-raw)). Not every cluster has every level: placement stops where the evidence stops, and some branches of the hierarchy are shallower than others.

### Why not deeper?

When CyteType stops at a level, the report explains why under **Why this level was retained**, in the Selection Rounds of the Justification section. The usual reasons:

| Reason given in the report (shortened) | What it means for you |
|---|---|
| More specific descendants were investigated, but their evidence was not strong enough. | The data does not separate the subtypes. More informative markers, a finer clustering or a more specific study context can help. |
| The identity is already the most specific one in the reference hierarchy. | There is nothing narrower to test. |
| The evidence consistently supported this identity, and no narrower distinction was needed. | A confident stop. |
| The evidence did not support a narrower label reliably within the permitted rounds. | Treat the depth of the label as uncertain. |
| The identity did not map reliably to a reference identity. | It was kept outside the hierarchy (see [Classification](#interpret-classification)). |
| The evidence supported a combined identity, but its components share no reliable deeper reference category. | The cluster appears to contain two identities; check its heterogeneity and consider subclustering. |

### Cell state

The cell state describes what the cells are doing rather than what they are, for example proliferation, exhaustion or senescence. It is supported by the gene sets under Cell State Evidence in the report. An unexpected state, such as proliferation in a cell type that is usually quiescent in your tissue, deserves a look: it can be real biology or a sign that the cluster mixes populations.

### Display groups

The report's Biological Order and the CSV's `coarse_cell_type` group clusters at one level of the hierarchy, chosen for each study so that there is a manageable number of groups. Unresolved clusters form a group of their own, and each identity outside the hierarchy is listed under its own label. Use the groups for orientation; a cluster's own label and levels are what count.

(interpret-classification)=
## Classification

Each label has a classification that relates it to the reference hierarchy. It is in the raw results as `classification` (see [Raw results](results-in-your-data.md#results-raw)), and the report shows a chip in the Justification header for some values. Current jobs use three values:

| Value | In the report | Meaning | What to do |
|---|---|---|---|
| `catalogue` | No chip | The identity was placed on an entry of the reference hierarchy. | Review the evidence as usual. |
| `outOfCatalog` | No chip; the hierarchy shows only the identity | A coherent identity with no matching entry in the hierarchy. The Justification section names the broader biological group it was placed in. | Read the evidence closely and confirm the identity with known markers or the literature. Make sure it is not an artifact, such as a population of doublets. |
| `unresolved` | **Unresolved** | No identity could be assigned, and the label is `Unknown`. The reviewer does not run, so the confidence is Low. | Check the cluster's quality (counts, mitochondrial fraction, doublets) and its marker genes, and consider merging, splitting or removing it. |

Reports of older jobs can also show these values:

| Value | Chip | Meaning | What to do |
|---|---|---|---|
| `broaderParent` | **Broader parent** | Only a broader entry of the hierarchy is supported. | Treat the label as a parent-level call, and subcluster if you need more detail. |
| `mixedDoublet` | **Mixed / doublet** | Signals of more than one cell type, such as doublets or a mixed population. | Check for doublets and subcluster (see [Heterogeneity](#interpret-heterogeneity)). |
| `malignantAberrant` | **Malignant / aberrant** | Malignant or otherwise aberrant cells. | Confirm with independent evidence, such as copy-number inference or clinical annotations. |

In current jobs, a mixed population shows up as a heterogeneity flag, or as a combined identity in the Justification section.

(interpret-confidence)=
## Reviewer confidence

After annotation, an independent reviewer assesses each final label. A panel of expert perspectives generated by CyteType's models examines the marker evidence. The reviewer weighs their agreement together with the cell classifier and the genes that set the cluster apart from its neighbors, and assigns one of three levels with a short rationale:

| Level | The reviewer's criterion |
|---|---|
| High | Strong agreement across the panel, with compelling evidence. |
| Moderate | Mixed assessments or notable concerns. |
| Low | Fundamental problems identified by most of the panel. |

Read the confidence as a categorical triage signal:

- **It is not a calibrated probability.** High does not mean that the label is correct with a known probability, and the three levels are not evenly spaced. Use them to decide where to look first, not as an acceptance threshold.
- **It is not any of the other numbers in the report.** The Review confidence and Reference match confidence percentages in the Justification section rank candidates inside the annotation step. The Match score compares CyteType's label with yours. The cell classifier's accuracy measures how separable the cluster is. None of them is the reviewer's confidence.
- **Low has more than one cause.** The reviewer may have found real problems. The cluster may also be unresolved, in which case no review runs, or the review may have failed. In those two cases the Quality Control section shows "Confidence details not available." instead of a rationale. Read the rationale before you act on a Low.

`reviewerConfidence` is missing (`None` in Python, `NA` in the R results table) when a stored result has no confidence at all, for example results saved by cytetype 0.7.0 and older. The cell-level confidence column shows `Unknown` for cells whose cluster has no confidence, including clusters without a result.

What to do with each level:

- **High:** spot-check the supporting genes and the heterogeneity call.
- **Moderate:** read the rationale and the panel's weaknesses and alternatives, and look at the missing and unexpected genes.
- **Low:** treat the label as unconfirmed. Check whether the review ran, then the marker genes, the heterogeneity and the cluster's quality. A finer or cleaner clustering, or a more specific study context, followed by a new run often helps.

To list the clusters that need a closer look:

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
import pandas as pd

table = pd.DataFrame(annotator.get_results()["annotations"])
needs_review = table[(table["reviewerConfidence"] != "High") | table["isHeterogeneous"].eq(True)]
needs_review[["clusterId", "annotation", "reviewerConfidence", "isHeterogeneous", "heterogeneityType"]]
```
:::
:::{tab-item} R
:sync: r
```r
res <- obj@misc$cytetype_results
needs_review <- subset(res, !(reviewerConfidence %in% "High") | isHeterogeneous %in% TRUE)
needs_review[, c("clusterId", "annotation", "reviewerConfidence", "isHeterogeneous", "heterogeneityType")]
```
:::
::::

(interpret-heterogeneity)=
## Heterogeneity

The reviewer also judges whether the cluster is one population (**Homogeneous**) or contains distinct populations (**Heterogeneous**). It weighs the panel's assessments, conflicts between lineage markers, the cell classifier's breakdown of the cluster's cells and the differences from neighboring clusters. Variation in state within one cell type, such as activation or maturation, does not count as heterogeneity.

When a cluster is heterogeneous and cell-level data was uploaded with the job, a second step classifies the type of heterogeneity from how pairs of conflicting genes are co-expressed in single cells (the Co-expression Evidence in the report):

| Type | Meaning | What to do |
|---|---|---|
| `distinct_populations` | Two or more identities in separate subsets of cells. | Subcluster and annotate the subclusters. |
| `mixed_lineage` | Markers of mutually exclusive lineages occur together in the same cells, in a way that doublets do not explain. | Check doublet scores, then consider an aberrant or malignant program and validate it independently. |
| `doublet_artifact` | Conflicting markers occur together in the same cells in a pattern that suggests two cells captured in one droplet. | Run doublet detection, remove the doublets and recluster. |
| `transitional_state` | Cells form a continuum along a developmental or activation axis. | Consider a trajectory or pseudotime analysis rather than splitting the cluster. |
| `unresolved` | The co-expression evidence is too weak to choose a type. | Inspect the Co-expression Evidence, and subcluster to test whether a split holds. |

The type is missing (`None` in Python, `NA` in R) for homogeneous clusters, and also when it could not be classified, for example because the job ran without uploaded cell-level data (`require_artifacts=False`). The report shows the heterogeneity call and its reasoning but not the type: both clients return it as `heterogeneityType`, and the raw results also hold the reasoning for it.

To subcluster a heterogeneous cluster, for example cluster `5`:

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
import scanpy as sc

sc.tl.leiden(
    adata,
    restrict_to=("leiden", ["5"]),
    resolution=0.5,
    flavor="igraph",
    n_iterations=2,
    key_added="leiden_sub",
)
```
:::
:::{tab-item} R
:sync: r
```r
library(Seurat)

# Splits cluster "5" of the active identities using the shared nearest-neighbor graph
obj <- FindSubCluster(obj, cluster = "5", graph.name = "RNA_snn",
                      subcluster.name = "sub_cluster", resolution = 0.5)
```
:::
::::

Then compute marker genes for the new column and run CyteType on it with a new `results_prefix` (see [Prepare your data](../guides/prepare-data.md) and [Run, monitor and rerun jobs](../guides/running-jobs.md)). For doublets, use a dedicated detector such as Scrublet (`sc.pp.scrublet` in Scanpy) or scDblFinder in R.

(interpret-compare)=
## Compare with your existing labels

When the labels in your cluster column are names rather than numbers, the reviewer scores how similar CyteType's label is to yours, from 0 to 100. The score appears as **Match** in the report, as `author_label_similarity_score` in the CSV export and as `labelComparisonScore` in the raw results.

- **How it is calculated.** The score comes from CyteOnto, an ontology-based measure of the semantic similarity between cell type labels. When CyteOnto returns no score, the reviewer estimates one instead (`labelComparisonLlmScore` in the raw results) from the relationship between the two labels: 100 for the same type, 95 for synonyms, 70 when one is a broader or narrower version of the other, 40 for sibling types, 20 for distant relatives in the same lineage and 0 for different lineages, minus a penalty when the labels differ in state, such as activation or differentiation.
- **Bands.** The report's filter groups scores into High (80 or more), Moderate (40 to 79) and Low (below 40).
- **When there is no score.** Labels that are cluster numbers or generic names, such as `0`, `C1`, `cluster 3`, `unknown` or `NA`, are not compared.

How to read it:

- A high score means that the two labels agree, not that either is correct.
- A low score marks a disagreement worth reading, not an error. Many disagreements are about granularity (your `T cell` against CyteType's `CD8+ T Cell`), naming conventions or state. Read both labels and the evidence before you decide.
- In the report, switch **Labels** to **User Label** and use the **Match** filter to review the low matches first (see [Compare with your own labels](report.md#report-compare)).

(interpret-limits)=
## Known limits

- **Cluster-level only.** Every cell in a cluster gets the same label, so the results depend on your clustering. A cluster that merges two cell types gets one label (often flagged as heterogeneous), and a cell type split across clusters gets several labels, which may differ.
- **Not deterministic.** Rerunning the same data can change the wording, the reasoning, the confidence and sometimes the label. Compare runs by their biology rather than by exact strings, and keep the job ID of the run you report.
- **Gene expression only.** CyteType uses the expression matrix (in R, the RNA assay), the marker genes, the cell metadata and your study context. Other modalities in your object, such as surface protein or chromatin accessibility, are not used.
- **Human and mouse.** CyteType is validated on human and mouse data. Results for other species have not been validated.
- **Gene symbols.** CyteType reasons about gene symbols, so data indexed by Ensembl or other IDs needs a symbol column (see [Prepare your data](../guides/prepare-data.md)).
- **Input quality matters.** Weak or contaminated marker genes, poorly normalized data, or a vague or wrong study context lead to weaker labels.
- **Model-generated review.** The reviewer and its panel are generated by CyteType's models. They make the evidence easier to check; they do not replace review by a scientist who knows the tissue.
- **No detection of personal data.** CyteType does not detect or remove identifiers or protected health information. Your cell metadata is uploaded with the job, so remove identifying information before you run it (see [Security and privacy](../support/security-and-privacy.md)).
- **Display sampling.** The report's embedding shows up to `max_cells_per_group` cells per cluster (1,000 by default). This affects only the display.

(interpret-reproducibility)=
## Reproducibility: what to record

The annotation runs on CyteType's servers, and the server workflow and hosted models can change between runs. Neither the results nor the report record a server release, so record enough to identify each run:

- **The job ID.** It identifies the run and its report.
- **The client and its version**, for example `cytetype` 0.21.0 or `CyteTypeR` 0.6.0. The report also shows it as **Client** under Study Information.
- **The date of the run.**
- **`query.json`.** `run()` and `CyteTypeR()` write it to the working directory by default (`save_query` and `query_filename`). It holds what was submitted: the study context, metadata tags, cluster labels, marker genes, cluster metadata summaries, expression summaries and the client version. The Python client never writes LLM API keys to it. CyteTypeR writes the complete request, including your `llm_configs` and their API keys, so pass `save_query = FALSE` when you use your own LLM keys in R, or remove the keys before you archive the file.
- **Your LLM configuration**, if you used your own models: the provider and model names.

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
import cytetype

details = adata.uns["cytetype_jobDetails"]
print(details["job_id"], details["report_url"], cytetype.__version__)
```
:::
:::{tab-item} R
:sync: r
```r
details <- obj@misc$cytetype_jobDetails
cat(details$job_id, details$report_url, as.character(packageVersion("CyteTypeR")), "\n")
```
:::
::::

To cite CyteType and report the version you used, see [Cite CyteType and report versions](../support/citing.md).

(interpret-checklist)=
## Review checklist

Work through the report with this list before you use the labels:

1. **Completeness.** The job finished and every cluster has a label. Note any Unresolved clusters.
2. **Overview.** In Biological Order, the display groups and their sizes make sense for your tissue. An unexpected lineage, such as hepatocytes in blood, points to contamination, ambient RNA or a wrong study context.
3. **Low and Moderate confidence.** Filter the cluster list to them. For each, read the Quality Control rationale; if there is none, the review did not run or failed.
4. **Heterogeneity.** For each heterogeneous cluster, read the reasoning and the type, and decide whether to subcluster or remove doublets.
5. **Marker evidence.** The supporting genes are ones you would expect for the label, and the missing and unexpected genes have a plausible explanation. Color the embedding by one or two key genes to check.
6. **Depth.** Where a label is broader than you need, read **Why this level was retained** before you assign a finer type yourself.
7. **Outside the hierarchy.** Identities outside the reference hierarchy get extra scrutiny.
8. **Cell states.** The states are plausible for your tissue and condition.
9. **Your labels.** If you compared CyteType with existing labels, read every Low match.
10. **Record and export.** Save the job ID, client version, date and `query.json`, then download the CSV or save your object.

If you change the clustering, the marker genes or the study context, run CyteType again and compare the runs (see [Run, monitor and rerun jobs](../guides/running-jobs.md)).

## See also

- [The CyteType report](report.md): where each signal appears.
- [Results in AnnData and Seurat](results-in-your-data.md): the same signals in your object.
- [How CyteType works](../getting-started/how-cytetype-works.md): the annotation pipeline.
