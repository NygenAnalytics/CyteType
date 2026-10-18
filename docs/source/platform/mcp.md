---
description: Connect Claude, Claude Code, Cursor or another MCP client to CyteType to explore and compare your annotation jobs from your own AI assistant.
---

(mcp)=
# Connect AI assistants (MCP)

CyteType runs a [Model Context Protocol](https://modelcontextprotocol.io/) (MCP) server. Add it to an AI assistant that supports MCP, such as Claude, Claude Code or Cursor, and the assistant can open your CyteType jobs, read their annotations and the evidence behind them, query their data and compare jobs with each other. All tools are read-only: they never change your jobs.

**Server URL:** `https://mcp.cytetype.nygen.io`

The main use is synthesis across jobs. Connect several reports and ask the assistant which findings agree, which differ, and which differences are supported by marker genes, metadata, expression, differential expression or cell-cell interactions. Typical uses:

- comparing disease and control cohorts;
- checking whether an annotation reproduces across repeated runs;
- comparing similar tissues, donors, batches or time points;
- finding shared and study-specific cell populations.

[Chat](chat.md) is the alternative inside CyteType: it uses CyteType-hosted models and saves figures and tables. The MCP connector brings the same data to the assistant you already use.

## Before you start

- A CyteType account. You sign in with your institutional email address, as on the website (see [Sign in and API keys](../getting-started/sign-in.md)).
- An MCP client that supports remote servers over HTTP with OAuth sign-in.

You do not need an API key: the connector uses its own sign-in.

## Set up your client

### Claude (web, desktop and mobile)

1. In Claude, open **Customize**, then **Connectors** (`https://claude.ai/customize/connectors`).
2. Click **Add custom connector**.
3. Enter the name `CyteType` and the server URL `https://mcp.cytetype.nygen.io`, and continue through the dialog. If it asks how to register the OAuth client, choose **Register automatically**: the CyteType server supports automatic (dynamic) registration.
4. Click **Add**, then sign in (see [Sign in](#mcp-sign-in)). If the connector shows **Connect** instead, click it to sign in.
5. In a conversation, click **+** next to the message box, choose **Connectors** and turn CyteType on.

On Team and Enterprise plans, an organization owner first adds the connector under **Organization settings**, **Connectors**. Each member then finds it under **Customize**, **Connectors**, clicks **Connect** and signs in with their own CyteType account. A connector added in Claude is available in the web app, Claude Desktop and the mobile apps.

### Claude Code

Add the server from your terminal:

```bash
claude mcp add --transport http --scope user cytetype https://mcp.cytetype.nygen.io
```

`--scope user` makes the connector available in all your projects; without it, it is added to the current project only. Then start Claude Code, run `/mcp`, select `cytetype` and choose **Authenticate**. Your browser opens the CyteType sign-in page.

`claude mcp list` shows whether the server is connected or still needs authentication. `claude mcp remove cytetype` removes it, together with its stored sign-in.

### Cursor

Add the server to `~/.cursor/mcp.json` (all projects) or to `.cursor/mcp.json` in a project folder:

```json
{
  "mcpServers": {
    "cytetype": {
      "url": "https://mcp.cytetype.nygen.io"
    }
  }
}
```

Cursor supports OAuth sign-in for remote servers: when it connects to the server, it sends you to the CyteType sign-in page. You can enable or disable the server in Cursor's MCP settings.

### Other clients

Any MCP client can connect if it supports:

- remote servers over Streamable HTTP;
- OAuth 2.0 sign-in with dynamic client registration and PKCE;
- a sign-in redirect to a local address on your computer (`http://localhost` or `http://127.0.0.1`, any port). Redirects to Claude and Cursor are also accepted; a client that registers any other redirect address is refused.

A client that only runs local (stdio) servers can connect through the `mcp-remote` bridge, which needs Node.js and handles the sign-in through your browser. In the client's MCP configuration:

```json
{
  "mcpServers": {
    "cytetype": {
      "command": "npx",
      "args": ["mcp-remote", "https://mcp.cytetype.nygen.io"]
    }
  }
}
```

(mcp-sign-in)=
## Sign in

The first time your client connects, it opens the CyteType sign-in page ("Sign in to connect"):

1. Enter your email address and click **Continue**.
2. Enter the 8-digit code from the email. It is checked as soon as you have entered all eight digits, or when you click **Verify**.
3. The page confirms and sends you back to your client.

Signing in follows the website's rules: you need an institutional email address, a first sign-in creates your account, and reports already shared with that address become available.

Finish signing in within 30 minutes; after that, start the connection again from your client. Your client then stays signed in while you use it: access renews automatically and expires after 6 days without use, after which you sign in again. To stop a client from using your account, remove or disconnect the connector in that client. CyteType cannot end a client's sign-in from its side: the access a client already holds lapses on its own, within an hour, and its renewal within 6 days.

You must sign in to use the connector, also for public reports. Once signed in, you can open:

- public reports;
- your own reports;
- reports shared with you;
- reports visible to your organization, while you are a member.

These are the same rules as for the report page and Chat (see [Visibility and sharing](sharing.md)). Unlike Chat, the connector also opens jobs that are still running.

## Connect jobs

Ask the assistant to connect a job by its report ID or report URL (`https://cytetype.nygen.io/report/<job_id>`). Without an ID, it can list your jobs or, in clients that support it, let you choose one. The response gives the job's status:

| Status | Meaning |
|---|---|
| `completed` | The job is finished |
| `processing`, `pending` | The job is still running; its annotation results appear when it reaches its final step, so treat what you see as partial |
| `failed` | The job failed |
| `incomplete` | An older job with missing results |

It also gives the number of clusters and annotated clusters, the start of the study context and a summary per cluster.

- Up to 10 jobs can be connected at once.
- Connected jobs belong to your account and your client, so every conversation in that client sees the same list. They are forgotten after 12 hours in which no tool call uses them.
- With one job connected, the assistant does not need to name it. With several, it passes the job ID to each tool.
- Access is checked when a job is connected and again at least once a minute while it is used. If a job is deleted or your access ends, it is disconnected and the next tool call says why.

## Tools

The tools fall into six groups. Your assistant chooses them from your request; you do not need to name them.

### Session

| Tool | What it does |
|---|---|
| `connect_job` | Connects a job by ID or report URL, or lets you choose from your jobs; reports its status and whether results are partial |
| `list_my_jobs` | Lists up to 100 recent completed jobs you can open: your own, those shared with you and those visible to your organization (not other people's public jobs, and not archived jobs) |
| `list_connected_jobs` | Lists the connected jobs with status, cluster counts and labels, and any that were disconnected or could not be loaded |
| `disconnect_job` | Removes a job from the connected jobs; the job itself is not changed |

### Job overview and comparison

| Tool | What it does |
|---|---|
| `get_job_overview` | One row per cluster: label, classification, L1 to L5 levels, Cell Ontology term, cell state, reviewer confidence, selection score, heterogeneity, status and your original label |
| `list_clusters` | Clusters with label, status, top 10 marker genes, reviewer confidence, selection score, heterogeneity, original label and label-comparison score |
| `get_cluster_detail` | Everything about one cluster, including marker genes and the context narrative |
| `get_study_context` | The study context used during annotation, with the submitted study description, organism, keywords and header tags |
| `compare_connected_jobs` | Lines up clusters across connected jobs, with each job's organism, keywords and study summary |

### Expression and metadata

| Tool | What it does |
|---|---|
| `lookup_gene_expression` | Percentage of cells in each cluster that express each gene (not mean expression), for up to 30 genes |
| `summarize_gene_expression` | Per-cell expression statistics for groups of cells defined by metadata columns, for up to 10 genes and 40 groups; needs the uploaded files |
| `describe_metadata` | The cell metadata columns, their types, value counts and likely roles (cluster, sample, donor, disease, treatment) |
| `query_metadata` | One read-only SQL `SELECT` over the cell metadata table `obs`; at most 200 rows are returned |
| `cluster_metadata_enrichment` | Per-cluster count and percentage of cells with a given metadata value |
| `run_differential_expression` | Differential expression between two non-overlapping groups of clusters, stating the method used |
| `get_cell_interactions` | Ligand-receptor interactions for the whole dataset or one cluster |

### Annotation evidence

| Tool | What it does |
|---|---|
| `get_annotation_evidence` | Supporting, missing and unexpected genes behind a label, with reasoning |
| `get_annotation_context_fit` | How the label fits the tissue, study and cluster context |
| `get_cell_identity` | Cell state with its gene sets, the Cell Ontology term and why it was chosen, and the granular cell type |
| `get_v2_annotation_trace` | How the final label was reached: the L1 to L5 levels, classification, why no deeper level was supported, and the decision record |
| `get_annotation_candidates` | The candidate cell types considered, with reasoning and reasons for rejection |
| `get_candidate_gene_groups` | Supporting and refuting gene sets for one candidate |
| `get_embedding_coordinates` | Subsampled 2D coordinates with cluster centroids: 2,000 points by default, at most 5,000 |

### Quality and review

| Tool | What it does |
|---|---|
| `get_quality_control` | Reviewer confidence with reasoning, selection score, heterogeneity, classifier accuracy, marker validation against CellGuide and SingleCellBase, and the predictions of an independent cell-type model compared with the label |
| `get_expert_reviews` | Specialist reviewer assessments: the evidence examined, strengths, weaknesses, interpretation, alternatives and conclusion |
| `get_neighbor_degs` | Pre-computed differential genes against the most similar clusters |
| `get_co_expression_evidence` | Gene pairs tested for co-expression within a cluster's cells, with statistics and rationale |
| `get_pathway_details` | Pathway enrichment per gene-set library, with normalized enrichment scores |

### Targets and labels

| Tool | What it does |
|---|---|
| `get_target_discovery` | The job's latest [target analysis](target-analysis.md) run: status, options, design, populations in scope and the top candidates of each ranking (10 by default, at most 25) |
| `compare_labels_cyteonto` | Scores the job's original cluster labels against CyteType's labels from 0 to 100 with CyteOnto, and gives the reason for each skipped cluster |

### Reading the results

- `reviewerConfidence` (High, Moderate or Low) comes from the independent review step. `selectionScore` is the annotation step's 0 to 1 score for its final choice and is on a different scale.
- `classification` is `catalogue` (the label matches an entry in the reference cell-type hierarchy), `outOfCatalog` (a coherent identity without a reference entry) or `unresolved` (no identity could be assigned).
- A cluster without a result has `annotated: false`. Your own label is kept in `originalLabel` and is never presented as CyteType's label.
- `compare_connected_jobs` matches clusters by Cell Ontology term first, then by identical L1 to L5 levels, and finally by label text alone. Rows matched by text alone are flagged `textOnlyMatch`; treat them as weak evidence. Clusters without an annotation or with an unresolved identity are listed as excluded.
- `run_differential_expression` names its method. With the job's uploaded expression data, it runs a cell-level Wilcoxon rank-sum test on up to 5,000 cells per group, with Benjamini-Hochberg adjusted p-values; cells are treated as independent observations, so the p-values can overstate significance. Without uploaded data, it compares the percentage of cells expressing each gene, with no statistical test and no p-values.
- `query_metadata` accepts one `SELECT` from `obs`, with `WHERE`, `GROUP BY`, `HAVING`, `ORDER BY`, `LIMIT` and window functions; joins, subqueries, `WITH`, multiple statements and file or network functions are rejected.

## Example prompts

- "Connect to job `<job_id>` and summarize its clusters with their reviewer confidence."
- "What share of cells express CD3E, MS4A1 and LYZ in each cluster?"
- "Compare clusters 0 and 3 against clusters 1 and 5 with differential expression, and tell me which method was used."
- "Show the evidence for cluster 2's label, the other candidates that were considered and how the final label was chosen."
- "How many cells per condition are in each cluster?"
- "How does IL7R expression differ between conditions within each cluster?"
- "Which clusters have Low reviewer confidence, and what did the expert reviewers flag?"
- "Connect jobs `<job_id_1>` and `<job_id_2>` and compare their cell populations. Check each job's organism and tissue first, and treat label-only matches as tentative."

## Limits

Besides the limits on connected jobs and query sizes described above, two tools have hourly limits per user: cell-level differential expression with `run_differential_expression` (20 runs per hour) and `compare_labels_cyteonto` (10 requests per hour). When an hourly limit is reached, the tool says when to try again. [Usage limits and quotas](limits.md#limits-mcp) lists every connector limit.

## Privacy

- **Read-only.** The tools cannot change your jobs, their visibility or their sharing. `disconnect_job` only changes the list of connected jobs.
- **Your access, nothing more.** The assistant can open only the jobs your account can open, and access is checked again at least once a minute while a job is in use.
- **What the assistant sees.** Everything a tool returns is sent to your AI client and the model provider behind it: annotations, study descriptions and metadata tags and, for jobs with uploaded data, rows of the cell metadata table (up to 200 per query) and expression statistics. How that data is stored and used is governed by your AI provider's terms and settings, not by CyteType. Check that your data-use agreements allow it before you connect jobs with sensitive metadata.
- **What CyteType sees.** CyteType receives the tool calls your assistant makes, with their arguments (such as job IDs, gene lists and SQL queries), not the rest of your conversation.

See [Security and privacy](../support/security-and-privacy.md).
