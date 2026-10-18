---
description: What the CyteType clients upload, how CyteType protects and uses your data, who can see your reports, how deletion and retention work, and what to check before you upload sensitive data.
---

(security-and-privacy)=
# Security and privacy

CyteType prepares your data on your computer with the Python or R client and runs the annotation on a hosted service. This page describes what leaves your computer, how it is protected and used, who can see your results, and what to check before you upload sensitive data. If your institution has contractual, regulatory or model-provider requirements, contact contact@nygen.io before you upload any data.

## At a glance

- New reports are private. The owner decides whether to share a report with specific people, with the organization or publicly.
- Data is encrypted in transit (HTTPS with TLS) and encrypted at rest by the infrastructure providers that store it.
- Your scientific data is not used to train or benchmark models.
- You sign in with a one-time code sent to your email. The clients use API keys that you can disable at any time.
- You can run the annotation with your own model provider.
- Owners can delete reports. Deletion removes all access but keeps the underlying data (see [Deletion and retention](#security-and-privacy-deletion)).
- A data processing agreement (DPA) is available with a commercial license.

(security-and-privacy-uploads)=
## What leaves your computer

A job sends a request and uploads two files.

### The job request

- Your study context and the metadata tags you pass as `metadata`.
- Your cluster labels and the top marker genes of each cluster.
- For every gene, the percentage of cells in each cluster that express it.
- Per-cluster summaries of your categorical cell metadata: the values found in more than `min_percentage` of a cluster's cells. Turn them off with `aggregate_metadata=False` (R: `aggregate_metadata = FALSE`).
- Embedding coordinates of up to `max_cells_per_group` sampled cells per cluster.
- Settings such as `n_parallel_clusters` and the report visibility, and the name and version of the client.
- If you use your own model: your `llm_configs`, including the provider API key.

The client also saves the request in `query.json` in your working directory. The Python client never writes `llm_configs` to that file. CyteTypeR does, including your provider keys, so pass `save_query = FALSE` when you use `llm_configs` in R.

### The two files

- **`vars.h5`:** the normalized expression matrix for all cells and genes (`adata.X`; in R, the normalized data of the `RNA` assay), the gene metadata table (`adata.var`; in R, the assay's feature metadata), and raw counts when the client finds them (Python: whole-number counts in `adata.layers["counts"]`, `adata.raw.X` or `adata.X`; R: the `counts` layer of the `RNA` assay, when it holds whole numbers).
- **`obs.duckdb`:** the complete cell metadata table, every column of `adata.obs` (R: `meta.data`), and the 2D embedding coordinates of all cells.

### What is not uploaded

FASTQ or BAM files, and your AnnData or Seurat object as a whole, are not uploaded. Neither are other layers or assays, other embeddings, neighbor graphs, or unstructured data such as `adata.uns`, apart from the marker genes the client takes from your marker results.

The CyteType files still contain your expression data and all of your cell metadata, so treat them with the same care as the source data. Some settings do less than their names suggest:

- `aggregate_metadata=False` removes the metadata summaries from the request, but `obs.duckdb` still contains every metadata column.
- `require_artifacts=False` only lets a run continue when a file fails to build or upload. It is not a way to avoid uploading: when the files build, they are uploaded.
- There is no setting that sends only the summaries. To limit what is uploaded, remove columns from a copy of your object before you prepare it.

The files stay on your computer after the run until you delete them: `annotator.cleanup()` (R: `CleanUpArtifacts(prepped_data)`) removes `vars.h5` and `obs.duckdb`. Delete `query.json` yourself if you do not want to keep it.

### Other connections the clients make

- Browser sign-in opens CyteType pages in your browser and uses a temporary listener on `127.0.0.1` on your own computer.
- When you import the Python package or run the `cytetype` command, it asks PyPI whether a newer version exists. CyteTypeR asks GitHub when it is loaded. These checks send no information about your data or analysis.

## How your data is used

CyteType uses your data to run your jobs, to host your reports, and for the features you use on them, such as [Chat](../platform/chat.md), the [MCP connector](../platform/mcp.md) and [target analysis](../platform/target-analysis.md). Your scientific data is not used to train or benchmark models.

CyteType also collects operational telemetry to run, secure, support and improve the service: service health, performance, errors and how the product is used. When a job fails, diagnostic records can include information about that job.

Authorized Nygen staff may access service data when this is needed to operate, secure, support or troubleshoot CyteType.

(security-and-privacy-services)=
## Models and other services

The annotation workflow sends focused context to language models, not your files: the study context, marker genes, expression percentages, metadata summaries, candidate annotations and the reasoning produced along the way. Tools in the workflow can look up values in your uploaded files, and those values can then be part of a prompt. The complete `vars.h5` and `obs.duckdb` files are not sent to model providers.

Unless you bring your own model, CyteType runs the models it hosts through commercial model providers. Which provider serves a step can vary with the workflow step and with availability.

Other services that process data for CyteType:

- **Infrastructure:** cloud providers for sign-in and account records, computing and storage, and file transfer. They encrypt the data they store.
- **Scientific resources:** workflow tools query public gene-set, ontology and cell-marker databases (and, for target analysis, a public drug-target knowledge base) with focused terms such as gene names, library names or ontology terms, never with your files.
- **Operations:** separate services for product analytics, application monitoring, email delivery and checking email addresses at sign-in. They receive only the account, usage or diagnostic information they need.

(security-and-privacy-own-llm)=
## Your own model provider

When you pass `llm_configs`, the annotation agents of that job use only the models you configure, under your own agreement with that provider. Keep in mind:

- Your provider key and the prompts pass through CyteType's service, which calls your provider for you. Use a key made for CyteType, with a spending limit if your provider offers one, and disable it when you no longer need it.
- Your files are still uploaded to CyteType, and CyteType still hosts the report.
- [Chat](../platform/chat.md) always uses CyteType's hosted models, also for jobs annotated with your own model.
- How long your provider keeps prompts, and whether it trains on them, depends on your agreement and settings with that provider.
- A ChatGPT Enterprise workspace is separate from OpenAI API access. You need an API key for the provider's API.

See [Use your own LLM](../guides/llm-configuration.md).

(security-and-privacy-access)=
## Access controls

- **Sign-in:** the dashboard, reports and Chat use a one-time code sent to your email; there are no passwords. Browser sessions use secure, HTTP-only cookies that page scripts cannot read.
- **Command-line sign-in:** `cytetype setup` and `SetupCyteTypeR()` ask you to confirm in the browser every time (**Authorize CLI**), even when you are already signed in. The key is handed to a temporary listener on your own computer with a one-time code that only that listener can exchange, and it never appears in a URL or in the terminal.
- **API keys:** personal keys start with `cyt_p_` and organization keys with `cyt_`. A key created in the dashboard is shown once, when it is created; browser sign-in (`cytetype setup`, `SetupCyteTypeR()`) gives each computer you sign in from the same "CyteType CLI" key (see [Sign in and API keys](../getting-started/sign-in.md)). Keys do not expire, but you can disable a key at once in the [dashboard](../platform/dashboard.md). `cytetype logout` and `LogoutCyteTypeR()` only delete the copy saved on your computer, where the file is readable only by you on Linux and macOS.
- **Report access:** when you open a report, the server checks ownership, sharing and visibility, then gives your browser a signed, short-lived access token that works only for that report.
- **On the server:** every request for results is checked against the job's owner, shares and visibility. Database connections use TLS, and the service uses rate limits and monitoring.

(security-and-privacy-visibility)=
## Who can see your results

| Visibility | Who has access |
|---|---|
| Private (the default) | You, and the people you share the report with by email. A recipient signs in with that email address to open it. |
| Organization | Also all members of your organization. Only for reports submitted with an organization API key. |
| Public | Anyone with the link can open the report page, without signing in. In addition, any signed-in CyteType user can find the report in Chat's report search, use it in [Chat](../platform/chat.md) and the [MCP connector](../platform/mcp.md), and fetch its results with their own API key. |

Keep in mind:

- People who have access to a job in Chat or through the MCP connector can query its uploaded cell-level data, including individual rows of the cell metadata table with every column you uploaded. For a public report, that means any signed-in user.
- Anyone who can open a report can also read the results of the job's target analysis, although the report page does not display them.
- The MCP connector always requires sign-in. When you open results in an AI assistant through it, they are sent to that assistant under your agreement with its provider.
- Organization admins have no extra access to other members' private reports.
- The owner can change the visibility in the dashboard at any time, in both directions. That cannot recall copies that people already made, such as CSV downloads.

Before you make a report public, check that you may publish everything it exposes. [Visibility and sharing](../platform/sharing.md) lists exactly what a public link exposes, and which API keys can fetch a job's results.

(security-and-privacy-deletion)=
## Deletion and retention

Only the owner can delete a report, from the [dashboard](../platform/dashboard.md). Deletion cannot be undone from the dashboard.

After deletion:

- nobody can open the report, including you: the link shows that the owner deleted it;
- all shares and pending invitations are removed;
- requests for its results with an API key fail with the error code `REPORT_DELETED`;
- Chat and the MCP connector can no longer use it.

Deletion is a soft delete: the job's underlying data is kept on CyteType's storage and is not erased. Deleted jobs still count toward your usage limits or your organization's quota. If you delete a job that is still running, it disappears from your views but runs to the end.

Archiving is different: it can be reversed and does not change who can open the report. It only hides the report from your dashboard list, and from the report search in Chat and the job list of the MCP connector.

**Retention.** CyteType does not delete reports or job data automatically after a period of time. Uploaded files that never become part of a job expire after 24 hours and are then deleted. If you need data erased from CyteType's storage, or your institution requires a specific deletion procedure, contact contact@nygen.io before you upload data.

(security-and-privacy-sensitive)=
## Sensitive and regulated data

:::{important}
CyteType does not detect or remove personal identifiers or protected health information. Complete your institution's review before you upload confidential human data.
:::

Before you upload confidential human data:

1. Work on a reviewed copy of the AnnData or Seurat object.
2. Remove direct identifiers, and any cell metadata columns that you do not need, from that copy. Every remaining column is uploaded in `obs.duckdb`.
3. Keep the study context, metadata tags and cluster labels free of names, record numbers and other identifying text.
4. Confirm the transfer with your institution's data protection, ethics or security team.
5. Contact contact@nygen.io before you start if you need a DPA, a specific model provider, or another access or deployment arrangement.

## Security reviews and agreements

Customers with a commercial license can request a data processing agreement. Nygen can also answer security questionnaires and discuss requirements for model providers, access and deployment. Contact contact@nygen.io, ideally before you upload data that falls under a formal institutional review.
