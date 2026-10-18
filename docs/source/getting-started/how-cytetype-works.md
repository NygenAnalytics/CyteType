---
description: What CyteType does with your clusters, the stages of an annotation job, what you get back, what is uploaded, and how to treat the results.
---

(how-cytetype-works)=
# How CyteType works

CyteType annotates the clusters of a single-cell RNA-seq dataset. For each cluster it proposes a cell type, maps it to the Cell Ontology, calls a cell state, and records the evidence behind each decision. The work is done on CyteType's servers by a team of agents: each one handles a part of the problem, such as understanding the study, annotating a cluster or reviewing an annotation, and uses analysis tools on your data along the way. An independent reviewer checks every call before the results come back to you.

## What you provide

| Input | Python | R | What it is for |
|---|---|---|---|
| Clusters | a column of `adata.obs`, named by `group_key` | a column of `obj@meta.data`, named by `group_key` | The unit of annotation. Every cell in a cluster receives the cluster's label. |
| Marker genes | `sc.tl.rank_genes_groups()` results | a `FindAllMarkers()` table | The starting evidence for each cluster (the top 50 genes per cluster by default). |
| Study context | `study_context` text | `study_context` text | What the sample is: organism, tissue, condition, assay and anything else that frames the study. |
| Report metadata (optional) | `metadata` dictionary | `metadata` list | Tags shown in the report header. A `title` entry names the job in the dashboard. |
| Coordinates (optional) | a 2D embedding in `adata.obsm`, such as `X_umap` | a reduction, such as `umap` | The embedding shown in the report. |

From these inputs, the client prepares the rest on your computer:

- **Cluster summaries:** the percentage of cells in each cluster that express each gene, and the make-up of each cluster for your categorical cell metadata, such as sample, donor or condition.
- **A plotting sample:** up to 1,000 cells per cluster with their coordinates.
- **Two data files**, uploaded with the job: `vars.h5` holds normalized expression for all genes, raw counts when the client finds them, and gene metadata; `obs.duckdb` holds your cell metadata table and the 2D coordinates.

The uploaded files let CyteType compute evidence from your individual cells on the server, such as differential expression between related clusters. Chat and target analysis use them later as well. Without them, CyteType still annotates from the marker genes and cluster summaries, with less evidence to draw on.

[Prepare your data](../guides/prepare-data.md) covers each requirement, and [Study context and metadata](../guides/study-context.md) shows how to describe your study.

## What happens during a job

While a job runs, the report shows its stage: **Preparing**, **Understanding context**, **Annotating clusters**, **Refining annotations**, **Reviewing** and **Finalizing**. After the short preparation step, the stages do the following.

1. **Understanding the study.** CyteType reads your study context together with the data. It determines the organism, infers the likely tissue from gene expression, writes a short study summary with keywords, and chooses the pathway libraries that fit the study. When files were uploaded, it also works out their layout, for example which column holds the clusters and which columns describe samples, donors or conditions.
2. **Annotating each cluster.** Clusters are processed in parallel. For each one, CyteType first builds its context: pathway enrichment, which metadata values the cluster is enriched in or depleted of, and how much of it comes from its largest donor. It then weighs competing candidate identities over rounds of evidence review, looking at supporting markers, expected markers that are weak or absent, and unexpected markers. A supported identity is placed in a curated reference hierarchy of cell types, from broad lineage (L1) to specific identity (L5), linked to the Cell Ontology. CyteType moves to a more specific level only when the evidence supports it, and otherwise keeps the broader label. The cell state, such as activation or proliferation, is called separately from the identity.
3. **Refining across clusters.** CyteType compares clusters with their closest relatives, using differential expression between them (computed from your cells when files were uploaded). It can add a short qualifier to tell similar clusters apart, harmonizes the wording of labels, and sorts identities that are not in the reference hierarchy into biological groups.
4. **Independent review.** A separate reviewer checks each annotation against independent evidence: assessments from a panel of expert perspectives, marker genes from external databases (SingleCellBase and CellGuide), and differential genes against neighboring clusters. When files were uploaded, it also uses a classifier trained on your cells, which tests how well the key genes separate the cluster from the others, and co-expression checks for clusters that look mixed. The reviewer then:
   - assigns a confidence of **High**, **Moderate** or **Low**;
   - decides whether the cluster is heterogeneous and, if so, what kind (distinct populations, mixed lineage, doublet artifact, transitional state, or unresolved);
   - writes a one-sentence granular label with the cell type, state and key genes;
   - when your cluster names are labels rather than numbers or generic names, scores how closely CyteType's label matches your label for the cluster, from 0 to 100.
5. **Finalizing.** CyteType builds the study summary: display groups and colors, a marker gene dot plot and, when coordinates were provided, cluster outlines on the embedding. Then it saves the results and the report is complete.

By default, the agents use models that CyteType hosts. You can run them on your own model provider instead (see [Use your own LLM](../guides/llm-configuration.md)).

## What you get back

| Result | In your AnnData or Seurat object | In the report |
|---|---|---|
| Cell type label | `{prefix}_annotation_{group_key}` column | Cluster card and the embedding |
| Cell Ontology term and ID | `{prefix}_cellOntologyTerm_{group_key}` and `{prefix}_cellOntologyTermID_{group_key}` columns | Cluster card |
| Cell state | `{prefix}_cellState_{group_key}` column | Cluster card and the Cell State Evidence section |
| Reviewer confidence | `{prefix}_reviewerConfidence_{group_key}` column | Cluster card and the Quality Control section |
| Heterogeneity | `isHeterogeneous` and `heterogeneityType` in the per-cluster results | Cluster card and the Quality Control section |
| Granular label | `granularAnnotation` in the per-cluster results | Cluster card |
| L1 to L5 hierarchy | `raw_annotations` in the full results | Justification section; L1 to L4 also in the CSV download |
| Evidence: supporting, missing and unexpected genes, justification | The per-cluster results | Cell Type Evidence and Justification sections |
| Match with your own labels | `raw_annotations` in the full results | Filter in the cluster list and the CSV download |

`{prefix}` is `cytetype` unless you set `results_prefix`. In Python, `get_results()` returns the full results, with the per-cluster results in `annotations`. In R, `obj@misc$cytetype_results` is a table of the per-cluster results, with the full results in its `"result"` attribute.

Every job also has a report at `https://cytetype.nygen.io/report/<job_id>`. Its link is printed when the job starts and stored in your object. See [Results in AnnData and Seurat](../results/results-in-your-data.md) for every column and key, and [The CyteType report](../results/report.md) for every section of the report.

## What is uploaded

Each job sends:

- the request: your study context, report metadata, cluster names, marker genes, cluster summaries and the plotting sample, plus the model settings and API keys when you [use your own LLM](../guides/llm-configuration.md);
- the two files `vars.h5` and `obs.duckdb`, described above.

Your original AnnData or Seurat file is not uploaded. Note that `obs.duckdb` contains every column of your cell metadata, so remove columns you do not want to share before you prepare the job (before `CyteType()` in Python, `PrepareCyteTypeR()` in R). Turning off the metadata summaries with `aggregate_metadata` does not remove columns from the file. [Security and privacy](../support/security-and-privacy.md) describes how uploaded data is stored and who can access it.

## Treat the results as reviewed proposals

- **Annotation is per cluster, not per cell.** All cells in a cluster share its label. A cluster that mixes cell types is flagged as heterogeneous instead of being split. To resolve it, subcluster those cells and annotate them in a new job.
- **Runs are not deterministic.** Running the same data again can give different wording, a different level of detail and, for ambiguous clusters, occasionally a different call. Keep the job ID and report link of the run you report on (see [Cite CyteType and report versions](../support/citing.md)).
- **Confidence is a categorical judgment.** High, Moderate and Low summarize the reviewer's assessment; they are not calibrated probabilities.
- **Review the calls that matter to your study.** Start with clusters that have low confidence, are flagged as heterogeneous, or disagree with your own labels, and check their evidence in the report against the markers you know.

[Interpret and review results](../results/interpreting-results.md) explains each signal and the limits to keep in mind.

## After annotation

- **The report** shows every cluster with its evidence, quality-control checks and context, and can be shared ([The CyteType report](../results/report.md), [Visibility and sharing](../platform/sharing.md)).
- **Chat** answers questions about one or more completed runs and can compute new results from your uploaded data ([Chat](../platform/chat.md)).
- **The MCP connector** lets AI assistants read your jobs ([Connect AI assistants (MCP)](../platform/mcp.md)).
- **Target analysis** turns a completed job into cell-type-resolved target hypotheses, each with its evidence. It is available on licensed plans ([Target analysis](../platform/target-analysis.md)).
