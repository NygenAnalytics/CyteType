---
description: A guide to the interactive CyteType report, from opening it and following a running job to reading every evidence section, comparing with your own labels, sharing and exporting a CSV.
---

(report)=
# The CyteType report

Every job has an interactive report at `https://cytetype.nygen.io/report/<job_id>`. It shows your embedding, the label for each cluster and the evidence behind it: the candidates CyteType considered, the marker genes, the reviewer's assessment and how the label fits the study. The labels are also copied into your object (see [Results in AnnData and Seurat](results-in-your-data.md)); the report is where you check them.

(report-open)=
## Open the report

The client prints the report link as soon as the job starts, so you can follow the job before it finishes. The link is also stored with the results.

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
adata.uns["cytetype_jobDetails"]["report_url"]
```

To open a report by its job ID from a terminal:

```bash
cytetype view <JOB_ID>
```
:::
:::{tab-item} R
:sync: r
```r
obj@misc$cytetype_jobDetails$report_url

ViewCyteTypeJob("<JOB_ID>")
```

If you installed the command-line launcher with `InstallCyteTypeRCli()`, `cytetyper view <JOB_ID>` does the same from a terminal.
:::
::::

`cytetype view` and `ViewCyteTypeJob()` open the report through the sign-in page of the server you signed in to; if your browser is already signed in, you go straight to the report (see the [command-line reference](../reference/cli.md)). You can also open any of your reports from the Reports list in the [dashboard](../platform/dashboard.md) with **View**.

(report-access)=
## Who can open a report

The report uses your browser sign-in (the emailed one-time code at `https://cytetype.nygen.io/login`), not your API key. Who can open it depends on the job's visibility:

| Visibility | Who can open the report |
|---|---|
| Private (the default) | You, and the people you shared it with by email. |
| Organization | Members of your organization's license, and the people you shared it with. |
| Public | Anyone with the link, without signing in. Any signed-in user can also find it in Chat's report search. |

If you are not signed in, the report shows **Authentication Required** with a **Sign In to View Report** button, and after signing in you return to the report. An account without access sees the same screen even when it is signed in: ask the owner to share the report with that email address. A deleted report shows **Report Deleted**, and a job ID that does not exist shows **Error Loading Report**. [Visibility and sharing](../platform/sharing.md) explains how to change the visibility and manage shares, and what a public report exposes.

(report-live)=
## Follow a running job

The report opens while the job is still running and updates itself, so you don't need to reload the page.

- With no cluster selected, **Study Information** shows a progress tracker with the job's steps: Preparing, Understanding context, Annotating clusters, Refining annotations, Reviewing and Finalizing. During Annotating clusters and Reviewing it also counts how many clusters have passed that step.
- In the cluster list, a spinner marks each cluster while the job runs (a clock before it starts). Selecting one shows **Annotation Pending** or **Annotation in Progress** with its current step, for example Annotating or Refining.
- The embedding is available from the start. The labels and the evidence appear when the job reaches its final step.
- Gene search becomes available once the report's gene index is ready; until then the search box reads "Loading gene index…".

(report-layout)=
## Page layout

The report has three columns:

- **Left sidebar:** the cluster list, with sorting, filters, a label switch and search.
- **Center:** gene search, the embedding, and a card for the selected cluster (or the study information).
- **Right:** a toolbar (Share, Open Chat, Download as CSV) above the evidence for the selected cluster, or the study summary when no cluster is selected.

Select a cluster in the list or click it on the embedding. The address bar then includes `?cluster=<id>`, so you can bookmark or send a link to a specific cluster, and the browser's back and forward buttons move between the clusters you viewed. The up and down arrow keys move to the previous or next cluster in the list.

:::{note}
Reports of jobs run with earlier CyteType versions can show a tabbed panel (Markers, Justification, Review, Clinical Info, Cluster Context) instead of the sections described on this page.
:::

(report-sidebar)=
## Left sidebar

The toolbar at the top of the sidebar has four controls:

| Control | Options |
|---|---|
| Sort | **Biological Order** groups clusters under display-group headers. **Alphabetical** sorts by label. **Cluster ID** sorts by CyteType's cluster ID. |
| Filter | **Confidence**: High, Moderate, Low. **Match**: High (≥80%), Moderate (40-79%), Low (<40%), No Original Label. |
| Labels | **CyteType Label** or **User Label** (the labels in your cluster column). |
| Search clusters... | Filters the list by cluster ID or label as you type. Esc clears it. |

Each row shows CyteType's cluster ID, the label and three dots. The dots show the reviewer's confidence, or the Match score when you display user labels; hover over them to see the value. Clusters that are still running are always listed. Clusters without a confidence value are hidden as soon as you untick any of the confidence boxes.

In Biological Order, each header names a display group and counts its clusters. CyteType picks the level of its reference hierarchy for these groups so that the study has a manageable number of them, puts unresolved clusters in a group of their own, and lists each identity outside the hierarchy under its own label. The headers appear once the job has finished. While a cluster is selected, hover over a header to highlight its group on the embedding.

**Dataset Overview**, below the toolbar, deselects the cluster and shows the [study summary](#report-summary).

The footer links to your reports in the **Dashboard**, the **Citation** (the CyteType preprint), a light and dark theme switch, and the CyteType repository on GitHub.

(report-center)=
## Center column

(report-gene-search)=
### Gene search

Type a gene symbol in the search box ("Search gene expression"), then pick a suggestion or press Enter. The embedding colors each cluster by the percentage of its cells that express the gene, on a scale from 0% to 100%; hover over a cluster to see its value in the legend. Only the genes in your data can be found. Click the gene name in the legend for a short description of the gene from MyGene.info (the top match, not filtered by species), and click × to clear the coloring.

Gene names in the right panel work the same way: click a gene to color the embedding by it, and click it again to clear.

### Embedding

The embedding shows the 2D coordinates sent with the job (UMAP by default) for a sample of up to 1,000 cells per cluster (`max_cells_per_group`), and at most 50,000 points in all. This sampling affects only the display. The embedding cannot be zoomed or panned. Hover over a cluster to see its label, and click it to select it. When cluster shapes are available, a button below the embedding cycles the view between Scatter (cells only), Hybrid (cells over the cluster shapes) and Hull (shapes only), and another one shows abbreviated cluster labels, which you can drag into place; **Reset label positions** puts them back. The shapes are created at the end of the job, so reload a report that you opened while the job was running. If the job had no coordinates, the embedding area reads "UMAP Visualization Not Available".

(report-cluster-card)=
### Cluster card

When a cluster is selected, a card below the embedding summarizes it. Click the title, the subtitle, the cell state or the confidence to jump to the matching section of the right panel.

| Field | Content |
|---|---|
| Title | The CyteType label. Opens Justification. |
| Subtitle | The granular annotation, a one-sentence description of the cluster. Opens Cell Type Evidence. |
| Cell State | The cell state, if one was called. Opens Cell State Evidence. |
| Confidence | The reviewer's confidence, with a **Heterogeneous** badge when the reviewer found distinct populations in the cluster. Opens Quality Control. |
| Ontology | The Cell Ontology ID and term name, linked to the term in the EBI Ontology Lookup Service. |
| User Label | Your label for the cluster, with its **Match** score when there is one. |

### Study information

With no cluster selected, the card shows **Study Information**: the job **Status** (with the progress tracker while the job runs), a **Completed** row with the time the job was submitted, the **User-Provided Metadata** you passed as `metadata` (values that are URLs become links), and the **Client** that submitted the job with its version, for example `anndata v0.21.0` or `seurat v0.6.0`. [Study context and metadata](../guides/study-context.md) explains which metadata to pass.

(report-toolbar)=
## Toolbar

The toolbar sits above the right panel.

| Button | What it does |
|---|---|
| Share | Shown only to the report's owner. Invites one person by email; they get view access. **Manage all shares in dashboard** opens your reports in the dashboard. |
| Open Chat | Opens [Chat](../platform/chat.md) in a new tab with this run attached, and with the selected cluster if there is one. Chat needs you to be signed in; the button reads **Sign in required** when you are not. |
| Download as CSV | Downloads the cluster table described in [Download as CSV](#report-csv). |

People you share a report with by email can open it after they sign in with that address, including people who don't have an account yet. See [Visibility and sharing](../platform/sharing.md).

(report-evidence)=
## Right panel: the evidence for a cluster

Select a cluster to see its evidence in up to five sections: Justification, Cell Type Evidence, Cell State Evidence, Quality Control and Context Fit. A navigation bar at the top jumps to each section and adds a small summary: the number of genes for Cell Type Evidence, the number of gene sets for Cell State Evidence, and a colored dot for the confidence in Quality Control. Sections, and most items inside them, collapse and expand when you click their header.

Gene chips are colored by how often the gene is detected in this cluster compared with the average of the other clusters: red when it is detected more often, blue when less often. Hover over a chip to see the log2 ratio and the percentage of the cluster's cells that express the gene.

(report-justification)=
### Justification

How CyteType arrived at the label.

- **Classification chip.** The header shows **Unresolved** when no identity could be assigned. Reports of older jobs can also show **Broader parent**, **Mixed / doublet** or **Malignant / aberrant**. [Interpret and review results](interpreting-results.md#interpret-classification) explains each.
- **Annotation hierarchy.** The label's place in CyteType's reference hierarchy of cell types, as a path that starts at L1 (broadest). For an identity outside the hierarchy it shows only the identity. When the label combines two identities, each has its own row.
- **Summary.** The conclusion for the label. Python returns the same text as `justification`.
- **Candidates Investigated.** The identities that were considered and rejected, marked with ×. Each shows **Evidence against** and **Why considered**, then **Evidence details** and **Miscellaneous**.
- **Selection Rounds.** The identity selected in each round, marked with a **Round** chip and ✓. CyteType evaluates candidates in rounds and moves to more specific identities while the evidence supports them. Each round shows **Evidence supporting**, then either **How the annotation was refined** (a more specific identity went on to the next round) or **Why this level was retained** (why the annotation stopped at this level), followed by the details.
- **Label Refinement.** When several clusters share an identity, CyteType compares them and can add a qualifier that tells them apart. This item shows the **Distinguishing feature** and the **Reasoning**; the CSV export has the qualifier in its own column.
- **Biological order group.** For an identity outside the hierarchy, the broader biological group it was placed in and why.

The **Evidence details** of a candidate can include:

| Item | Content |
|---|---|
| Evidence supporting, Evidence against, Why considered | The arguments for and against the identity, and why it was proposed. |
| Evidence synthesis | A balanced assessment of the supporting, missing and unexpected genes. |
| Pathways | Pathways that support the identity. |
| Gene evidence | Genes listed as **Evidence supporting**, **Expected but weak or absent** and **Unexpected expression**. |
| Reference | The **Reference hierarchy** path and the **Reference match** reasoning for placing the identity in the hierarchy. |
| Population consistency | Whether the evidence points to one population or to two distinct identities. |

**Miscellaneous** holds **Review confidence** and **Reference match confidence** percentages, the **Reference node identifier** and the **Reference ranking rationale**. A candidate with a **combined identity** chip merges two identities, for example when the cluster appears to contain both; its **Merged identity** item lists them.

:::{important}
The percentages under Miscellaneous rank candidates inside the annotation step. They are not the reviewer's confidence shown in Quality Control, and they are not probabilities that the label is correct.
:::

Reports of older jobs show **Other candidates investigated** instead of the rounds, with ×××, ×× or × for a strong, moderate or weak rejection, a **Consideration reason** and a **Rejection reason**, and the **Initial candidates** under **Supplementary information**.

(report-cell-type-evidence)=
### Cell Type Evidence

The genes behind the label, and how the cluster differs from the clusters closest to it.

- **Overview graphic:** the number of **Supporting**, **Missing** and **Unexpected** genes, drawn as a bar.
- **Supporting:** genes whose expression supports the label, with an explanation.
- **Missing:** genes expected for the label that are weak or absent in the cluster.
- **Unexpected:** genes expressed in the cluster that are inconsistent with the label.
- **Differentially Expressed vs Nearest Clusters:** one item per nearby cluster, with a volcano plot (log2FC against -log10 pAdj) and a table of the genes with pAdj below 0.05 and an absolute log2FC above 1 (columns **Gene**, **log2FC** and **pAdj**). Hover over a gene in the table to find it in the plot. A chip names the method:
  - **Wilcoxon:** a cell-level test that treats every cell as independent, so p-values can overstate significance.
  - **DESeq2:** a pseudobulk test in which samples, not cells, are the replicates.
  - **Log2FC:** the fold change of the percentage of cells expressing each gene. No test is run, so there is no volcano plot and the table has no p-values.
- **Co-expression Evidence:** for pairs of genes, the percentage of cells that express one gene only, both genes or neither, with the phi coefficient (φ, from -1 for mutually exclusive to +1 for always co-expressed) and a chi-squared p-value. It appears, for up to three pairs of a supporting and an unexpected gene, when the reviewer found the cluster heterogeneous and cell-level data was uploaded with the job.

(report-cell-state-evidence)=
### Cell State Evidence

The cell state, kept separate from the cell's identity. Chips give the state (**State:**) and the number of gene sets. **Summary** explains the state, and **Supporting gene sets** lists named gene programs with their genes and reasoning. The section is absent when no state was called.

(report-quality-control)=
### Quality Control

The independent reviewer's assessment of the final label. Chips show **Confidence** (High, Moderate or Low) and **Heterogeneity** (Heterogeneous or Homogeneous).

| Item | Content |
|---|---|
| Confidence | The reviewer's reasoning for the confidence level. "Confidence details not available." means the review did not run, or failed, for this cluster (see [Reviewer confidence](interpreting-results.md#interpret-confidence)). |
| Heterogeneity | The reasoning for the heterogeneity call. |
| Cell classifier | Results of a classifier trained on the marker genes of all clusters and tested on cells it did not see during training. The **accuracy** chip shows how well the genes separate this cluster's cells from the others (green from 85%, amber from 65%, red below). The **Entropy** chip is the mean uncertainty of the per-cell predictions; lower means the cells look more alike. The bars show the labels the classifier gave to this cluster's cells. It needs cell-level data uploaded with the job. |
| Reviewer Panel | Several reviewer perspectives, each with a specialty (for example a T cell immunologist), a one-line conclusion, and a **Strength**, **Weakness**, **Biological interpretation** and **Alternatives**. These assessments are generated by CyteType's models, not written by human experts. |
| Markers from External Databases | How the cluster's genes compare with reference marker lists, and the **Database markers** from SingleCellBase, CellGuide Canonical and CellGuide Computational. |

(report-context-fit)=
### Context Fit

How the label fits the study and the cluster.

- **Tissue:** how the identity fits the tissue and the study context.
- **Cluster:** how it fits this cluster's own context. **Cluster metadata** has a **Summary**, **Keywords**, and the **Full Metadata Breakdown (%)**: the percentage of the cluster's cells in each category of your categorical cell metadata columns.
- **Pathways Enriched/Deleted:** the **Supporting pathways** for the identity with an explanation and a **Pathway summary**, followed by **Detailed pathway data** for each pathway library. Each table lists the **Pathway** and its normalized enrichment score (**NES**). The **Single sample** tab ranks genes by the share of this cluster's cells that express them; the **Differential** tab ranks them against the average of all clusters. A positive NES means the pathway's genes rank high, a negative NES that they rank low. Only pathways with a false discovery rate below 0.01 are listed, up to 50 per library and tab.

(report-summary)=
## Study summary

With no cluster selected (click **Dataset Overview**), the right panel shows:

- **Study Overview:** CyteType's description of the study, built from your study context and the data. Python returns the same text as `studyContext` from `get_results()`.
- **Marker Gene Expression:** a dot plot with one column per cluster, labeled with its cell type and ordered by display group. The rows are up to three key supporting genes per cluster, grouped by display group. Dot size is the percentage of cells that express the gene; color is the mean expression when cell-level data was uploaded with the job, and otherwise the percentage. Hover over a dot to see the values. The plot appears once the job has finished.

The embedding colors every cluster in this view, and the center column shows the study information.

(report-compare)=
## Compare with your own labels

The report treats the labels in your cluster column (`group_key`) as your labels. Switch **Labels** to **User Label** to show them in the cluster list and on the embedding; the dots then show each cluster's **Match** score, and the cluster card shows the score next to your label.

The Match score is the similarity between your label and CyteType's label, from 0% to 100%. The filter bands are High (80% or more), Moderate (40% to 79%) and Low (below 40%). Use the **Match** filter to focus on the clusters where you and CyteType disagree.

There is no score when your labels are cluster numbers or generic names (such as `0`, `C1`, `cluster 3` or `unknown`); such clusters count as No Original Label. To compare CyteType with an existing annotation, run it with that annotation column as `group_key`, with marker genes computed for the same column (see [Prepare your data](../guides/prepare-data.md)). [Interpret and review results](interpreting-results.md#interpret-compare) explains how the score is calculated and how to read it.

(report-csv)=
## Download as CSV

**Download as CSV** saves `cytetype_annotations_<job_id>.csv`, with one row per finished cluster. Clusters that are still running are left out.

| Column | Content |
|---|---|
| `cluster_id` | CyteType's cluster ID. |
| `cell_type` | The CyteType label. |
| `coarse_cell_type` | The display group (the Biological Order header). |
| `cell_state` | The cell state. |
| `cl_term` | The Cell Ontology term name. |
| `cl_term_id` | The Cell Ontology ID. |
| `confidence` | The reviewer's confidence: High, Moderate or Low (Unknown when there is none). |
| `is_heterogeneous` | `TRUE` or `FALSE`. |
| `supporting_genes` | The key supporting genes, separated by `;`. |
| `cell_type_abbrv_label` | The abbreviated label shown on the embedding. |
| `color_hex` | The cluster's color in the report. |
| `author_label` | Your label for the cluster. |
| `author_label_similarity_score` | The Match score from 0 to 100; empty when there is none. |
| `l1`, `l2`, `l3`, `l4` | The annotation hierarchy levels. For an identity outside the hierarchy, `l1` holds the identity and `l2` to `l4` are empty. |
| `qualifier` | The distinguishing feature added by label refinement, if any. |

The CSV leaves out the evidence text, the reviewer panel, the heterogeneity type and level L5. The Python client returns these (see [Results in AnnData and Seurat](results-in-your-data.md)).

## See also

- [Interpret and review results](interpreting-results.md): how to read the labels, confidence and heterogeneity.
- [Results in AnnData and Seurat](results-in-your-data.md): the same results in your object.
- [Visibility and sharing](../platform/sharing.md), [Dashboard](../platform/dashboard.md) and [Chat](../platform/chat.md).
