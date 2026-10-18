---
description: Control who can open a CyteType report with private, organization or public visibility and email sharing, and understand what each level exposes.
---

(sharing)=
# Visibility and sharing

Every report has a visibility level, and its owner can also share it with individual people by email. Together they decide who can open the report page, use the run in [Chat](chat.md) and connect it from an AI assistant through the [MCP connector](mcp.md). Fetching results with an API key follows slightly different rules, described below.

## Visibility levels

| Visibility | Who can open the report |
|---|---|
| **Private** (default) | You, and the people you share it with |
| **Organization** | Every member of the report's organization, and the people you share it with |
| **Public** | Anyone with the link, without signing in |

New reports are private unless you choose otherwise. Organization visibility applies only to reports in an organization workspace, which are the jobs submitted with an organization API key (see [Dashboard](dashboard.md)).

## Choose the visibility when you submit

Pass `visibility` to `run()` (cytetype 0.21.0 or later) or to `CyteTypeR()` (CyteTypeR 0.6.0 or later): `"private"`, `"org"` or `"public"`. Leaving it out (`None` in Python, `NULL` in R) uses the server default, private. The requested value is also saved with the job details, in `adata.uns["cytetype_jobDetails"]["visibility"]` or `obj@misc$cytetype_jobDetails$visibility` (with the default `results_prefix`).

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
adata = annotator.run(
    study_context="Human PBMC from a healthy donor",
    visibility="public",
)
```
:::
:::{tab-item} R
:sync: r
```r
obj <- CyteTypeR(
  obj, prepped_data,
  study_context = "Human PBMC from a healthy donor",
  visibility = "public"
)
```
:::
::::

:::{warning}
Organization visibility needs an organization API key (`cyt_...`). A personal key is refused even if your account belongs to an organization. The server checks this only when the job is submitted, which happens after your files are uploaded: with a personal key, the run stops with an `APIError` ("Org visibility requires an organization license"; in R, `Job submission failed: [HTTP_ERROR] Org visibility requires an organization license`), no job is created and the upload is discarded. Run it again with an organization key, or submit it as private and change the visibility in the dashboard.
:::

## Change the visibility later

In the dashboard, open **My Reports**, click the edit icon next to the report's visibility, choose the new level and click **Save**. You can change it in either direction at any time; **Organization** is offered only in the organization workspace.

The change applies at once:

- The report page checks access each time it loads data, so people who no longer qualify lose access straight away.
- Chat checks access before each answer, and hides earlier answers and outputs that came from runs the person can no longer open.
- The MCP connector checks access at least once a minute while a job is in use, and disconnects jobs that are no longer accessible.

Changing the visibility cannot recall copies that people already made, such as CSV downloads, Chat downloads or results an AI assistant has already received.

## Share with specific people

The owner can share a report with people by email, from the dashboard or with the **Share** button on the report page. The steps are in [Dashboard](dashboard.md).

- Recipients get view access: they can open the report and use the run in Chat and over MCP. They cannot rename it, change its visibility, share it further, transfer it or delete it.
- Each recipient gets an email with the link and the address to sign in with.
- A recipient without a CyteType account is listed as **Pending** until they sign in for the first time with that address; the invite then turns into access.
- Recipients need an institutional email address, the same rule as for signing in.
- Shares work at every visibility level, including for people outside your organization.
- You can remove a share or a pending invite at any time in the dashboard. Deleting a report removes all of its shares and invites.

Recipients find personal reports shared with them under **Shared with Me** in their dashboard, and reports of their own organization under **Workspace Reports**. A report from an organization they do not belong to is not listed in their dashboard; they open it from the link in the email.

## Who can open what

| Access | Private | Organization | Public |
|---|---|---|---|
| Report page | The owner and the people it is shared with | Also all members of the organization | Anyone with the link, signed in or not |
| Chat (sign-in required) | The owner and the people it is shared with | Also all members of the organization | Any signed-in user |
| MCP connector (sign-in required) | The owner and the people it is shared with | Also all members of the organization | Any signed-in user |
| Results with an API key (`get_results()`, `GetResults()`) | The owner's API keys | Also organization API keys of the same organization | Any API key |

Points to note:

- Email shares give access to the report page, Chat and MCP, but not to results through an API key. A recipient who calls `get_results()` or `GetResults()` with their own key gets an `ACCESS_DENIED` error (an `AuthorizationError` in Python). To give a collaborator the results, they can download the CSV from the report page, or you can make the report Organization (for members of your organization who use an organization key) or Public.
- For organization reports, API-key access works only with the organization's keys (`cyt_...`), not with members' personal keys.
- Chat opens completed runs only. The MCP connector also opens jobs that are still running and marks their results as partial.
- Organization admins have no extra access to reports (see [Organization roles](dashboard.md)).
- Nobody can open a deleted report.

## What a public link exposes

A public report publishes more than the page you see. Anyone with the link, without signing in, can:

- view the whole report: labels, Cell Ontology terms, cell states, confidence and heterogeneity calls, every evidence and review section, the study summary and the marker dot plot;
- see the inputs the report is built from: the study description and metadata tags you submitted (including any links), your own cluster labels if you supplied them, the per-cluster summaries of categorical cell metadata, the marker genes, the 2D embedding coordinates sent with the job (a subsample of cells), and the percentage of cells in each cluster that express each gene, for every gene in the dataset;
- download the annotation table as CSV.

The job's latest [target analysis](target-analysis.md) run and its results belong to the report's data, so they are public too, even though the report page does not display them.

In addition, anyone with a CyteType account can:

- fetch the results with their own API key (`get_results()`, `GetResults()`);
- find the report in Chat's report search, which lists public reports alongside the user's own, so a public report can be found by name and not only through its link;
- attach the run in Chat and analyze the uploaded cell-level data: queries over the cell metadata table, per-cell expression summaries and plots, and differential expression between groups of cells;
- connect the job through the MCP connector, whose tools include read-only SQL queries over the cell metadata (up to 200 rows per query) and per-cell expression statistics.

The uploaded files cannot be downloaded as files, but Chat and MCP return results computed from them, including individual rows of the cell metadata table. That table holds every cell metadata column you uploaded (`adata.obs` or the Seurat `meta.data`), not only the columns shown in the report.

A public link does not let anyone change the report, its visibility or its sharing, and it does not let them start a target analysis.

:::{important}
Before making a report public, check that you could publish everything listed above, including every cell metadata column (for example donor or sample identifiers), under your data-use agreements and institutional policies. See [Security and privacy](../support/security-and-privacy.md).
:::
