---
description: Ask questions about one or more completed CyteType runs in Chat, which answers from the stored results and your uploaded data with text, tables and figures you can download.
---

(chat)=
# Chat

Chat, at `https://cytetype.nygen.io/chat`, is a workspace for asking questions about completed runs. It reads the annotation results of the runs you attach and, when the job uploaded them, the cell metadata and the expression matrix. It answers in text and saves tables, figures and other outputs that you can reuse and download.

You need to sign in with your CyteType account to use Chat. Having a report link is not enough, even for a public report.

## Open Chat

- **From a report**: click **Open Chat** in the report toolbar. Chat opens in a new tab with the run attached and the cluster you were viewing selected. If you are not signed in, the button reads **Sign in required**.
- **From the dashboard**: click **Chat** on a completed report, or **Chat** in the sidebar for an empty conversation.
- **Directly**: go to `https://cytetype.nygen.io/chat`.

When you open Chat for a report you have discussed before, it asks whether to continue your most recent conversation about that report or start a new chat.

## Attach runs

Chat calls attached runs "sources". A conversation can have up to 10.

1. Click the paperclip (**Add source**) next to the message box.
2. Under **Find a completed report**, type a report name or ID and click **Find**. Click **Find** with an empty box to list the most recent reports you can attach.
3. Click **Attach** next to the report.

You can attach completed runs that you can open: your own, those shared with you, those visible to your organization and public reports (see [Visibility and sharing](sharing.md)). Runs that are still processing cannot be attached. Archived reports are not listed in the search, but you can attach one by entering its full report ID.

Each attached run appears as a chip under the message box. Click **×** on a chip to detach the run. After you detach a run, Chat no longer uses it or the earlier answers that were based on it, although those answers stay visible to you.

## Select clusters

To point Chat at particular clusters, click **Clusters** on a source chip, search by cluster ID or label, and click **Target** next to each cluster you want. You can select up to 12 clusters per message, across all attached runs. Selected clusters appear as chips under the message box; click **×** to remove one.

Chat focuses on the selected clusters: for example, plots show only those clusters unless you ask for all of them. Opening Chat from a report selects the cluster you were viewing.

## Ask questions

Type your question and press **Send**. While Chat works, a status line shows its progress, and the answer appears as it is written. Click the stop button to cancel an answer.

Chat looks up your data with tools before it makes claims about a dataset, and it answers concisely, saving large tables as outputs rather than writing them into the answer. What it can do depends on what the job uploaded:

| Data | Available for |
|---|---|
| Annotation results: labels, L1 to L5 levels, Cell Ontology terms, cell states, reviewer confidence, heterogeneity, evidence, candidates, reviews, pathways | Every run |
| Percentage of cells in each cluster that express a gene | Every run |
| Queries and summaries of cell metadata (donors, samples, conditions) | Runs that uploaded the cell metadata (`obs.duckdb`) |
| Per-cell expression summaries and plots, cell-level differential expression | Runs that uploaded the expression matrix (`vars.h5`) and the cell metadata |
| Ligand-receptor interactions between clusters | Human and mouse runs with at least two clusters that uploaded both files |

See [Prepare your data](../guides/prepare-data.md) for how the clients build and upload these files.

### Questions that work well

- Evidence: "Why was cluster 7 labeled as plasmacytoid dendritic cells? What evidence argued against it, and which other candidates were considered?"
- Review: "Which clusters have Low reviewer confidence or are flagged as heterogeneous? Summarize the reviewers' concerns for each."
- Expression: "Make a dot plot of CD3E, MS4A1, LYZ, NKG7 and PPBP across all clusters."
- Metadata: "How are donors distributed across clusters? Show it as a stacked bar chart of proportions."
- Within-cluster comparisons: "In the selected cluster, compare cells from treated and untreated samples and show a volcano plot."
- Cluster comparisons: "What distinguishes cluster 3 from cluster 5? Run differential expression and list the top genes in each direction."
- Interactions: "Which ligand-receptor interactions involve the macrophage clusters?"
- Several runs: "For each attached run, list the clusters labeled as T cells with their reviewer confidence, in one table."

Questions about a specific gene set, cell type or condition work better than very broad requests. If an answer needs more data than Chat can read at once, it stops with a message asking for a narrower question; split the request into smaller steps.

## How Chat analyzes runs

Chat follows fixed rules so that its answers stay attributable and statistically honest:

- **One run at a time.** It inspects each attached run separately and says which run every value comes from.
- **No pooling across runs.** It never pools cells across runs, never runs differential expression across runs, and does not infer batch correction, build joint embeddings or assume gene orthology beyond what the stored evidence shows.
- **Single-run outputs.** Cell groups, differential expression results, metadata summaries and charts each come from one run, and charts never mix runs. A table can combine runs only when every row names its source run in a `sourceJobId` column.
- **Stated methods.** It names the differential expression method it used. With the uploaded expression matrix (`vars.h5`), this is a cell-level Wilcoxon rank-sum test (groups larger than 5,000 cells are randomly subsampled). Cell-level tests treat cells as independent and ignore donor and sample structure, so their p-values can overstate significance. Without the expression matrix, Chat can compare only whole clusters, using the stored percentage of cells that express each gene, with no statistical test and no p-values.
- **Clear quantities.** Reviewer confidence (High, Moderate or Low) is kept separate from the 0 to 1 scores of the annotation step, and the stored expression values are the percentage of cells expressing a gene, not mean expression.
- **Unknowns stay unknown.** An unknown organism, assay, normalization, genome build or gene namespace is reported as unknown, not guessed.
- **Your data is data.** Study descriptions, metadata values, cluster labels and saved outputs are never treated as instructions.

Chat cannot draw heatmaps; ask for a dot plot to compare several genes across groups.

Chat always uses models hosted by CyteType, even for runs annotated with your own LLM. Chat use does not count toward your annotation limits.

## Outputs and downloads

Chat saves what it produces as outputs, listed in the **Outputs** panel next to the conversation:

- **Figures**: bar charts, stacked bar charts, box plots, dot plots, volcano plots, radar charts and embedding scatter plots.
- **Tables**, including tables that compare runs.
- **Cell groups**: saved sets of cells, such as the cells of some clusters or the cells that match a metadata condition, which later questions can reuse.
- **Differential expression results**, **metadata summaries** (the most common values of metadata columns in a cell group) and **cluster enrichment results** (the count and percentage of cells in each cluster with a given metadata value).

Open an output in the panel, then its **•••** menu, to:

- **Add to chat**: attach the output to your next message, for example to reuse a saved cell group. Up to 12 outputs can be attached to a message.
- **Download CSV** or **Download JSON**: tables, differential expression results, metadata summaries and cluster enrichment results download as CSV; figures and cell groups download as JSON. When an output holds only part of the data (for example a differential expression result that stores the top genes only), the menu offers **Download partial CSV**.
- **Provenance**: see which runs the output came from.
- **Pin**: keep the output at the top of the list. When Chat revises a figure, the list shows the newest version, plus any earlier version you pinned or added to your next message.
- **Delete**: remove the output from the conversation.

To rename the open output, click its title at the top of the **Outputs** panel.

## Manage conversations

- **New chat** opens an empty conversation, which is saved when you send its first message. Your conversations are listed in the sidebar, most recently active first; click one to open it.
- To rename the open conversation, click its title at the top of the page.
- Open the **•••** menu next to that title and click **Delete** to delete the conversation with its messages and outputs.
- Chat runs one answer at a time for your account. If an answer is running in another conversation, Chat offers to open it or stop it before you send a new message.

Conversations are private to your account: people who can open the same runs do not see them. If you lose access to a run, for example because it was deleted or unshared, the answers and outputs that came from it are hidden.

## Limits

A conversation can have up to 10 attached runs and 100 messages, and your account can keep up to 100 conversations. An answer can create up to 12 outputs, and a table holds up to 500 rows. Chat reads the last 20 completed exchanges of a conversation, and it may drop the oldest ones, or older tool results within an answer, to stay within its reading limit.

[Usage limits and quotas](limits.md#limits-chat) lists every Chat limit, including the message length and the time an answer may take, and the message you see when you reach each one.
