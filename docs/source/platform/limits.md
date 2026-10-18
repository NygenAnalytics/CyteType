---
description: CyteType usage limits for personal accounts and organizations, the data, job, Chat and MCP limits, the errors you get when you reach them, and how to get more capacity.
---

(limits)=
# Usage limits and quotas

CyteType is free for academic and non-commercial research within the limits on this page. Personal accounts have rolling job limits; organizations have a cluster quota set by their license. A few technical limits apply to everyone. The [Account page](dashboard.md) shows how much you have left.

## Personal accounts

Jobs submitted with a personal API key (`cyt_p_...`) count toward three rolling windows:

| Window | Jobs |
|---|---|
| Last 24 hours | 10 |
| Last 7 days | 50 |
| Last 30 days | 100 |

The windows roll: a job stops counting 24 hours, 7 days or 30 days after you submitted it. When you reach a limit, the error tells you which window is full and roughly when it frees up.

What counts:

- Jobs that are running or completed, including jobs you later archive or delete.
- Jobs annotated with your own LLM. Passing `llm_configs` does not bypass the limits (see [Use your own LLM](../guides/llm-configuration.md)).

What does not count:

- Failed jobs.
- Requests rejected at submission, for example by an invalid key or a failed model check.
- Jobs submitted with an organization key, which count toward the organization's quota instead.
- Viewing reports, fetching results, Chat, the MCP connector and target analysis.

## Organizations

Jobs submitted with an organization API key (`cyt_...`) count toward the organization's cluster quota:

- The quota is a number of clusters per license period, set by the license. For a license without an end date, all of the organization's jobs count.
- The clusters of running and completed jobs count, including jobs that were later deleted. Failed jobs do not count.
- A job is accepted only if all of its clusters fit in the remaining quota.
- When the license expires, uploads and new jobs with the organization's keys are refused.

The Account page in the organization workspace shows the clusters annotated in the current period and the clusters remaining.

## Per-key quotas

Any API key, personal or organization, can have its own **Quota limit**: the maximum number of clusters that jobs submitted with the key can contain, over the key's lifetime. Running and completed jobs count, including deleted ones; failed jobs do not. A key quota applies in addition to your personal limits or your organization's quota.

You set quotas on your personal keys yourself. Organization admins set them on organization keys, for example to give each member or pipeline a share of the organization's quota. See [Dashboard](dashboard.md).

## When limits are checked

The clients upload your data files before they submit the job:

1. When an upload starts, the server checks that your API key is valid and active and that its organization's license has not expired.
2. When the job is submitted, after the upload, the server checks your personal limits, the key's quota, the organization's quota and the requested visibility.

A job refused at step 2 does not run, and its upload is lost: you upload again when you resubmit. Before submitting a large dataset, check the Account page. Organization visibility needs a key with an organization license; when you request it with a key that looks like a personal key (`cyt_p_...`), both clients warn before the upload.

## Data and upload limits

| File | Maximum size |
|---|---|
| `obs.duckdb` (cell metadata) | 2 GiB (2,147,483,648 bytes) |
| `vars.h5` (expression matrix) | 50 GiB (53,687,091,200 bytes) |

Both clients check the file sizes before uploading. A file that is too large stops the run with an error that names the file ("exceeds upload limit"); in Python it is a `ValueError`. With `require_artifacts=False`, the job continues without that file (CyteTypeR without either file), and features that need it (cell-level queries in Chat and MCP, cell-level differential expression, target analysis) are unavailable for the job. [Large datasets](../guides/large-datasets.md) explains how to reduce the file sizes.

An uploaded file expires if no job uses it within 24 hours. The clients submit the job right after uploading and upload again on every run, so you do not need to act on this.

## Job limits

| Limit | Value |
|---|---|
| Marker genes used per cluster | 250 (any more are ignored) |
| `n_parallel_clusters` | 1 to 50; used only with your own LLM |
| Server time limit per job | 4 hours; the job is then stopped and reported as failed |
| Client wait (`timeout_seconds`) | 7,200 seconds (2 hours) by default |
| Target analysis run | 3 hours |

When the client stops waiting, Python raises `CyteTypeTimeoutError` and CyteTypeR stops with an error of class `cytetype_timeout_error` that names the job. The job keeps running on the server. Fetch the results later with `get_results()` in Python or `GetResults(job_id = ...)` in R (see [Run, monitor and rerun jobs](../guides/running-jobs.md)).

(limits-chat)=
## Chat limits

| Limit | Value |
|---|---|
| Runs attached to a conversation | 10 |
| Clusters selected per message | 12 |
| Outputs attached to a message | 12 |
| Message length | 16,000 characters |
| Answers running at the same time | 1 per account, across all conversations |
| Time per answer | 10 minutes waiting to start, 30 minutes running |
| Outputs created per answer | 12 |
| Rows per table | 500 |
| Conversations | 100 per account |
| Messages per conversation | 100, including stopped and failed ones |

When you reach one of these limits, Chat tells you what to do:

| Situation | Message |
|---|---|
| Another answer is running | "Another Chat turn is still running. Open that conversation or wait for it to finish." |
| 100 conversations | "You have reached the saved-conversation limit. Delete an older conversation before starting another." |
| 100 messages in a conversation | "This conversation has reached its turn limit. Start a new conversation to continue." |
| An answer waited or ran too long | "The turn exceeded its execution window. Please retry.", "The chat turn exceeded its execution deadline." or "The Chat turn exceeded its time limit." |
| An answer gathered too much data | "This turn gathered more tool output than the model can read at once. Ask a narrower question or split it into smaller steps." |

See [Chat](chat.md).

(limits-mcp)=
## MCP connector limits

| Limit | Value |
|---|---|
| Connected jobs | 10 per signed-in user and AI client |
| Connected jobs forgotten after | 12 hours without a tool call that uses them |
| Cell-level differential expression | 20 runs per hour per user |
| Label comparison (`compare_labels_cyteonto`) | 10 requests per hour per user |
| Rows per metadata query | 200 |
| Genes per expression lookup | 30 (`lookup_gene_expression`), 10 (`summarize_gene_expression`) |
| Groups per expression summary | The 40 largest groups, with up to 5,000 sampled cells per group |
| Cells per group in cell-level differential expression | 5,000 (larger groups are subsampled) |
| Genes per direction in differential expression | 1 to 100 |
| Embedding points | 2,000 by default, at most 5,000 |
| Target candidates per ranking | 10 by default, at most 25 |
| Sign-in | Renewed automatically while in use; sign in again after 6 days without use |

When an hourly limit is reached, the tool returns a message such as "Rate limit reached: at most 20 differential expression runs on cell-level data per hour. Try again in about 12 minute(s)." See [Connect AI assistants (MCP)](mcp.md).

## Errors when you reach a limit

| Situation | Error code | Python exception | Message |
|---|---|---|---|
| Personal limit reached | `RATE_LIMIT_EXCEEDED` | `RateLimitError` | "Daily limit reached (10/10). Resets in about 3 hours. For a free academic license with higher limits, email contact@nygen.io" ("Weekly" and "Monthly" for the other windows) |
| Organization quota too small for the job | `QUOTA_EXCEEDED` | `QuotaExceededError` | "Quota exceeded: This job requires 24 clusters but you only have 10 remaining (990/1000 used). Please contact your organization admin" |
| Key quota too small for the job | `QUOTA_EXCEEDED` | `QuotaExceededError` | The same message, without the last sentence |
| Organization license expired | `LICENSE_EXPIRED` | `AuthorizationError` | "License has expired. Please contact your organization admin(s): ..." |
| Key disabled | `TOKEN_INACTIVE` | `AuthenticationError` | "API token has expired or been deactivated" |
| Organization visibility with a key without an organization license, such as a personal key | `HTTP_ERROR` | `APIError` | "Org visibility requires an organization license" |
| Upload file too large | (checked by the client) | `ValueError` | "... exceeds upload limit: ..." |
| Job stopped after 4 hours | (job status) | `JobFailedError` | "Job ... failed" |
| Client stopped waiting | (checked by the client) | `CyteTypeTimeoutError` | "Job ... did not complete within ...s" |

The numbers in the messages are examples. The CyteType exceptions have `.error_code` and `.message` attributes and can be imported from `cytetype`:

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
from cytetype import QuotaExceededError, RateLimitError

try:
    adata = annotator.run(study_context="Human PBMC from a healthy donor")
except RateLimitError as err:
    print(err.error_code)  # RATE_LIMIT_EXCEEDED
    print(err.message)     # which window is full and when it frees up
except QuotaExceededError as err:
    print(err.message)     # clusters required and clusters remaining
```
:::
:::{tab-item} R
:sync: r
```r
obj <- tryCatch(
  CyteTypeR(obj, prepped_data, study_context = "Human PBMC from a healthy donor"),
  error = function(e) {
    message(conditionMessage(e))  # the server's message
    obj
  }
)
```
:::
::::

In R, a refused submission stops with a message that starts with "Job submission failed:" and includes the error code in square brackets, for example `[QUOTA_EXCEEDED]`. A personal-limit error shows the server's message directly, followed by a reminder that the limits apply to every job, including jobs that use your own LLM. [Errors and error codes](../reference/errors.md) lists every error.

## Get more capacity

- **Academic and non-commercial research**: email contact@nygen.io about a free academic license with higher limits.
- **Commercial use** requires a license. Email contact@nygen.io.
- **Organizations**: admins can see the organization's usage on the Account page and rebalance per-key quotas in the dashboard. To raise the quota or renew the license, contact contact@nygen.io.
