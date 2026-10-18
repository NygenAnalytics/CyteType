---
description: The exceptions raised by the CyteType Python client and the error conditions signalled by CyteTypeR, the server error codes behind them, and what to do about each.
---

(errors)=
# Errors and error codes

When something goes wrong, the Python client raises an exception and the R client stops with an error condition. Errors that come from the CyteType server carry a machine-readable code, such as `RATE_LIMIT_EXCEEDED`, and a message that explains the problem. This page lists the exceptions, the codes and what to do about each one. For problems described by their symptoms, see [Troubleshooting](../support/troubleshooting.md).

(errors-python)=
## Python exceptions

All CyteType exceptions can be imported from `cytetype` (and from `cytetype.api`); the [Python API reference](../reference/api/exceptions.html) documents each class:

```python
from cytetype import (
    CyteTypeError,
    APIError,
    AuthenticationError,
    AuthorizationError,
    RateLimitError,
    QuotaExceededError,
    JobNotFoundError,
    JobFailedError,
    LLMValidationError,
    ClientDisconnectedError,
    NetworkError,
    CyteTypeTimeoutError,
)
```

They form this hierarchy:

```text
CyteTypeError                      base class: .message, .error_code
├── APIError                       errors reported by the server (any code not listed below)
│   ├── AuthenticationError        INVALID_TOKEN, TOKEN_INACTIVE, AUTHENTICATION_REQUIRED
│   ├── AuthorizationError         ACCESS_DENIED, LICENSE_EXPIRED, VPN_BLOCKED
│   ├── RateLimitError             RATE_LIMIT_EXCEEDED
│   ├── QuotaExceededError         QUOTA_EXCEEDED
│   ├── JobNotFoundError           JOB_NOT_FOUND
│   ├── JobFailedError             JOB_FAILED, or the job's status became "failed"
│   ├── LLMValidationError         LLM_VALIDATION_FAILED
│   └── ClientDisconnectedError    CLIENT_DISCONNECTED
├── NetworkError                   the server could not be reached, or its reply could not be read
└── TimeoutError                   a request, or the wait for a job, took too long
                                   (import it as CyteTypeTimeoutError)
```

Every CyteType exception has two attributes:

- `message`: the human-readable message, also returned by `str(err)`;
- `error_code`: the server's code, such as `"QUOTA_EXCEEDED"`, or a code set by the client. It is `None` for some client-side errors, for example a job that failed while you were waiting for it or a job that did not finish within `timeout_seconds`.

A few points to keep in mind:

- **`CyteTypeTimeoutError` is not Python's built-in `TimeoutError`.** An `except TimeoutError:` clause does not catch it. Catch `CyteTypeTimeoutError` or `CyteTypeError`. The class is also available as `cytetype.api.TimeoutError`.
- **Codes without their own class** (for example `INTERNAL_ERROR`, `HTTP_ERROR` and `REPORT_DELETED`) raise a plain `APIError`. Check `err.error_code` to tell them apart.
- **`JobNotFoundError` is rarely raised.** When the results of a job cannot be found, the client raises a plain `APIError` with `error_code == "JOB_NOT_FOUND"`. Check the code rather than the class.
- **`get_results()` does not raise server errors.** If fetching fails, it logs the reason (for example "Failed to retrieve results for job ...") and returns `None`. It raises only `AuthenticationError` when no API key is available for the job's server, and `ValueError` when the server address stored with the job is not valid.

(errors-python-local)=
### Errors raised before anything is sent

Problems with your data or arguments are found on your computer and raise standard Python exceptions:

| Exception | Raised by | Typical causes |
|---|---|---|
| `KeyError` | `CyteType()` | `group_key` is not a column of `adata.obs`; the `rank_genes_groups` results are not in `adata.uns` under `rank_key`; `gene_symbols_column` is not a column of `adata.var`. |
| `ValueError` | `CyteType()` | Missing cluster labels (`NaN`) without `label_na=True`; gene IDs instead of gene symbols ("Gene Symbol Detection Error"); `rank_genes_groups` computed with another `groupby`; marker genes that do not match `adata.var_names`; an invalid `api_url`. |
| `ValueError` | `run()` | Results with the same `results_prefix` already exist; `visibility` is not `"private"`, `"org"`, `"public"` or `None`; `upload_timeout_seconds` is not positive; `cleanup()` deleted the files and `require_artifacts=True`; an invalid `api_url`. |
| `pydantic.ValidationError` | `run()` | `metadata` values that are not strings; an `llm_configs` entry with an unknown provider or without `apiKey` (or the AWS fields); `n_parallel_clusters` outside 1 to 50; a `study_context` that is not a string. |

`pydantic.ValidationError` is a subclass of `ValueError`, so `except ValueError` catches it too. These checks run before any file is uploaded. For the fixes, see [Troubleshooting](../support/troubleshooting.md#troubleshooting-data).

With `require_artifacts=True` (the default), a failure while building or uploading `vars.h5` or `obs.duckdb` re-raises the original exception, for example an `OSError` when the disk is full or a `ValueError` when a file is larger than the upload limit. See [Artifacts and uploads](../support/troubleshooting.md#troubleshooting-artifacts).

(errors-when)=
### When each error can occur in `run()`

`run()` works through these steps in order, so the step tells you how far a job got:

1. **Checks on the arguments and the `CyteType` object:** `ValueError`.
2. **Finding your API key:** `AuthenticationError` with `AUTHENTICATION_REQUIRED` if no key is saved for the server and none was passed.
3. **Validating the request:** `pydantic.ValidationError`.
4. **Uploading `vars.h5` and `obs.duckdb`:** a rejected key raises `AuthenticationError` and a refused account (for example an expired license) raises `AuthorizationError`, even with `require_artifacts=False`. Other upload failures raise only with `require_artifacts=True`.
5. **Submitting the job:** `RateLimitError`, `QuotaExceededError`, `LLMValidationError`, or `APIError` (for example when `visibility="org"` is used with a personal API key instead of an organization key). All job limits and quotas are checked here, after the upload, so fixing these errors means uploading again.
6. **Waiting for the job:** `JobFailedError` if the job fails, `CyteTypeTimeoutError` if it does not finish within `timeout_seconds`. Brief network problems while waiting are retried. Once the job has been submitted, its ID is stored in `adata.uns` and you can fetch the results later with `get_results()`.

(errors-codes)=
## Server error codes

The table summarizes the codes; the sections below give the typical causes and fixes.

| Code | Python exception | In short |
|---|---|---|
| `INVALID_TOKEN` | `AuthenticationError` | The API key is not recognized. |
| `TOKEN_INACTIVE` | `AuthenticationError` | The API key was disabled. |
| `AUTHENTICATION_REQUIRED` | `AuthenticationError` | No API key was available. |
| `ACCESS_DENIED` | `AuthorizationError` | You may not access this job or resource. |
| `LICENSE_EXPIRED` | `AuthorizationError` | Your organization's license has expired. |
| `VPN_BLOCKED` | `AuthorizationError` | The request was refused because of the network it came from. |
| `RATE_LIMIT_EXCEEDED` | `RateLimitError` | Your personal account reached a job limit. |
| `QUOTA_EXCEEDED` | `QuotaExceededError` | Not enough clusters left in the quota. |
| `JOB_NOT_FOUND` | `APIError` or `JobNotFoundError` | No job with this ID on this server. |
| `JOB_FAILED` | `JobFailedError` | The job failed on the server. |
| `LLM_VALIDATION_FAILED` | `LLMValidationError` | A model in `llm_configs` failed the check before the job started. |
| `CLIENT_DISCONNECTED` | `ClientDisconnectedError` | The connection closed before the server had the whole request. |
| Other codes | `APIError` | See [Other codes](#errors-other-codes). |

(errors-invalid-token)=
### `INVALID_TOKEN`

The server does not recognize the API key: "Invalid API token".

- **Typical causes:** a mistyped or truncated key; a key from another CyteType server than the one in `api_url`.
- **What to do:** sign in again with `cytetype setup --force` (R: `SetupCyteTypeR(force = TRUE)`), or save a correct key with `cytetype login` (R: `LoginCyteTypeR()`). If you pass `auth_token`, check its value and that `api_url` points to the server that issued it.

(errors-token-inactive)=
### `TOKEN_INACTIVE`

The key exists but has been disabled: "API token has expired or been deactivated".

- **Typical causes:** you, or an admin of your organization, disabled the key in the dashboard. Disabling a key cannot be undone. Because browser sign-in reuses your "CyteType CLI" key, disabling it affects every computer that signed in with `cytetype setup`.
- **What to do:** run `cytetype setup --force` (R: `SetupCyteTypeR(force = TRUE)`) to get a new key, or save another active key with `cytetype login`. Jobs that were already running still finish; fetch their results with a valid key.

(errors-authentication-required)=
### `AUTHENTICATION_REQUIRED`

No API key was available.

- **Typical causes:** you have not signed in on this computer; the saved sign-in belongs to a different server than the one in `api_url` or `CYTETYPE_API_URL`; the saved credentials file is damaged or cannot be read (for example in a container without a home directory).
- **Messages:** "CyteType sign-in is required. Run `cytetype setup` first." (R: "Run SetupCyteTypeR() first."), or "Invalid CyteType credentials file: ...". An outdated client may show the server's message, "Authentication is required. Upgrade to the latest CyteType or CyteTypeR and run setup, or provide a valid API token."
- **What to do:** run `cytetype setup` (R: `SetupCyteTypeR()`) for the same server you use in your code, or pass `auth_token`. For a damaged file, run `cytetype logout` (R: `LogoutCyteTypeR()`) and sign in again. See [Sign in and API keys](../getting-started/sign-in.md).

(errors-access-denied)=
### `ACCESS_DENIED`

You may not access this job or resource: "You do not have permission to access this job".

- **Typical causes:** fetching the results of a job that you do not own and that is neither public nor an organization report of your organization. Organization reports can be fetched only with an organization API key of the same organization, not with a member's personal key. Sharing a report by email lets the recipient open it in the browser but does not give their API key access to its results. When you fetch results, an unknown key is also reported this way ("Invalid token or you do not have permission to access this job").
- **What to do:** fetch the results with the owner's key or an organization key, ask the owner to change the [visibility](../platform/sharing.md), or open the report in your browser if it was shared with you and use **Download as CSV**.

(errors-license-expired)=
### `LICENSE_EXPIRED`

The organization license that the key belongs to has expired: "License has expired. Please contact your organization admin(s): ...".

- **Typical causes:** you use an organization key (starting with `cyt_`) after the license period ended. An expired license blocks uploads and new jobs.
- **What to do:** contact your organization's admins, or contact@nygen.io to renew the license. In Python, this error stops `run()` during the upload with "Artifact upload was not authorized".

(errors-vpn-blocked)=
### `VPN_BLOCKED`

The server refused the request because of the network it came from.

- **Typical causes:** the request was sent through a VPN, proxy or similar service.
- **What to do:** try again from your institution's network without the VPN or proxy. If you think the request was refused by mistake, contact contact@nygen.io.

(errors-rate-limit)=
### `RATE_LIMIT_EXCEEDED`

Your personal account reached one of its job limits. The message names the limit and when it resets, for example "Daily limit reached (10/10). Resets in about 5 hours."

- **Typical causes:** personal accounts can run 10 jobs in any rolling 24 hours, 50 in 7 days and 100 in 30 days. Running and completed jobs count, including deleted ones; failed jobs do not. The limits also apply when you use your own LLM.
- **What to do:** wait until the time given in the message, or email contact@nygen.io about a free academic license with higher limits. The limit is checked when the job is submitted, after the files are uploaded. See [Usage limits and quotas](../platform/limits.md).

(errors-quota-exceeded)=
### `QUOTA_EXCEEDED`

The job has more clusters than remain in the quota: "Quota exceeded: This job requires 24 clusters but you only have 10 remaining (990/1000 used)".

- **Typical causes:** your organization's license quota is used up, or the key has its own cluster limit (set when the key was created or by an admin). Clusters of running and completed jobs count, including deleted jobs.
- **What to do:** ask your organization's admins to raise the key's limit or the license quota, or use another key. The quota is checked when the job is submitted, after the files are uploaded. See [Usage limits and quotas](../platform/limits.md).

(errors-job-not-found)=
### `JOB_NOT_FOUND`

The server has no job with this ID: "Job ... not found" (Python: "Results not found for job ...").

- **Typical causes:** a mistyped job ID; asking a different server (`api_url`) than the one that ran the job.
- **What to do:** copy the job ID from the dashboard or from the report link, and use the server the job ran on. In Python, this is usually a plain `APIError` with `error_code == "JOB_NOT_FOUND"`.

(errors-job-failed)=
### `JOB_FAILED`

The job failed on the server. In Python, `run()` raises `JobFailedError` with "Job ... failed" as soon as the job's status becomes failed, and prints the report link.

- **Typical causes:** a problem on the server, or a job that ran longer than the server allows (4 hours).
- **What to do:** open the report to see how far the job got, then run it again, with `override_existing_results=True` or another `results_prefix` if you reuse the same `CyteType` object. Failed jobs do not count toward the personal job limits. If it fails again, contact support with the job ID (see [Get help](../support/troubleshooting.md#troubleshooting-help)).

(errors-llm-validation)=
### `LLM_VALIDATION_FAILED`

A model in `llm_configs` did not pass the check that runs before the job starts: "All provided models failed validation: ... Ensure the model is reachable and supports tool calling."

- **Typical causes:** a wrong model name or API key; a `baseUrl` that CyteType's servers cannot reach (for example a local address); a model without tool calling. If any of the models fails, the whole request is rejected.
- **What to do:** check each configuration, make sure the endpoint is reachable from the internet and that the model supports tool calling, then run again. The check happens at submission, after the files are uploaded. See [Use your own LLM](../guides/llm-configuration.md).

(errors-client-disconnected)=
### `CLIENT_DISCONNECTED`

The connection closed before the server had received the whole request: "Client closed the connection before the request completed".

- **Typical causes:** an unstable network, or a proxy or firewall that cuts long requests, usually during an upload.
- **What to do:** check your connection and run again.

(errors-other-codes)=
### Other codes

These codes raise a plain `APIError` (or `NetworkError`, where noted):

| Code | Meaning | What to do |
|---|---|---|
| `INTERNAL_ERROR` | An unexpected error on the server ("An unexpected error occurred: ..."). | Try again later. The Python client retries each part of a file upload up to three times on its own. If the error persists, contact support with the message. |
| `HTTP_ERROR` | A request was refused for another reason; the message says why, for example "Org visibility requires an organization license". It is a `NetworkError` when the reply could not be read at all, for example an error page from a proxy. | Act on the message. |
| `REPORT_DELETED` | The owner deleted the report ("This report has been deleted by the owner"). | Deletion cannot be undone. Ask the owner if they have another copy of the results. |
| `JOB_PROCESSING`, `JOB_NOT_COMPLETED` | The job has not finished yet. | Wait and fetch the results again. |
| `INVALID_RESPONSE` | Set by the client: the server's reply was not in the expected form. | Try again; if it persists, update the client and contact support. |
| `PRESIGNED_URL_REJECTED`, `MISSING_ETAG` | Set by the client: the file storage service rejected part of an upload (`MISSING_ETAG` is a `NetworkError`). | Run again. If it keeps happening on a slow connection, see [Artifacts and uploads](../support/troubleshooting.md#troubleshooting-artifacts). |

(errors-python-example)=
## Catch and handle errors in Python

Catch the specific classes you can act on, then `CyteTypeError` for the rest:

```python
from cytetype import (
    AuthenticationError,
    CyteType,
    CyteTypeError,
    CyteTypeTimeoutError,
    JobFailedError,
    LLMValidationError,
    QuotaExceededError,
    RateLimitError,
)

annotator = CyteType(adata, group_key="leiden", rank_key="rank_genes_leiden")

try:
    adata = annotator.run(study_context="Human PBMC from a healthy donor")
except AuthenticationError as err:
    print(f"API key problem ({err.error_code}): {err.message}")
    print("Run `cytetype setup --force`, or pass a valid auth_token.")
except RateLimitError as err:
    print(f"Job limit reached: {err.message}")
except QuotaExceededError as err:
    print(f"Quota exceeded: {err.message}")
except LLMValidationError as err:
    print(f"A model in llm_configs failed validation: {err.message}")
except JobFailedError as err:
    print(f"The job failed on the server: {err.message}")
except CyteTypeTimeoutError as err:
    # A submitted job keeps running on the server.
    print(f"Timed out: {err.message}")
    print("Fetch the results later with annotator.get_results().")
except CyteTypeError as err:
    print(f"CyteType error ({err.error_code}): {err.message}")
    raise
```

`run()` returns the annotated AnnData object, so assign its return value (here back to `adata`). With `label_na=True`, it can be a new object.

To retry automatically after a rate limit, read the reset time from the message and run again later. Do not retry `AuthenticationError`, `AuthorizationError` or `QuotaExceededError` in a loop: they need a change on your side.

(errors-r)=
## R errors

CyteTypeR signals errors in these classes:

```text
error
├── cytetype_api_error         a problem with the CyteType API: e$message, e$error_code, e$call ("network" or "api")
│   └── cytetype_auth_error    sign-in and API key problems: e$message, e$error_code
├── cytetype_timeout_error     CyteTypeR() stopped waiting at timeout_seconds: e$message, e$job_id
└── simple R errors            everything else, including most refusals from the server
```

**`cytetype_auth_error`** conditions (which also inherit `cytetype_api_error`) carry an `error_code`:

| `error_code` | When |
|---|---|
| `AUTHENTICATION_REQUIRED` | No key is saved for the server ("CyteType sign-in is required. Run SetupCyteTypeR() first."), the credentials file is damaged ("Invalid CyteType credentials file: ... Run LogoutCyteTypeR(), then SetupCyteTypeR()."), or `GetResults()` was given `auth_token` for a job stored in the object without `api_url`. |
| `API_ORIGIN_MISMATCH` | `GetResults()` was given an `api_url` that differs from the server stored with the job. |
| `INVALID_TOKEN`, `TOKEN_INACTIVE`, `LICENSE_EXPIRED` | `SetupCyteTypeR()` or `LoginCyteTypeR()` was refused by the server. `INVALID_TOKEN` is also used when a request for a job's results is refused with "Authentication failed: Invalid or expired auth token". |
| `ACCESS_DENIED` | The server refused to return a job's results ("Authorization failed: Access denied"): the key is unknown or disabled, or the job is neither yours, nor an organization report of the key's organization, nor public. Also used when you select **Cancel** on the authorization page during setup. |
| `TIMEOUT`, `NETWORK_ERROR`, `CALLBACK_SERVER_ERROR`, `INVALID_RESPONSE`, `AUTHENTICATION_FAILED` | Browser sign-in or key validation failed: it was not finished within 5 minutes, the server could not be reached, the local listener could not start, the server's reply was not as expected, or another refusal. |

**`cytetype_timeout_error`** means that `CyteTypeR()` stopped waiting after `timeout_seconds` (default 7200). The job keeps running on the server. `e$job_id` holds its ID, and the message names the call that fetches its results: "Job <job_id> did not finish within 7200 seconds (timeout_seconds). The job keeps running on the server; retrieve its results later with GetResults(job_id = "<job_id>")."

**`cytetype_api_error`** is the parent class of `cytetype_auth_error`, so a handler for it also catches API key problems. CyteTypeR signals it when a job's status or results cannot be fetched: by `GetResults()`, and by `CyteTypeR()` while it waits (after five consecutive errors, or when a completed job's results cannot be fetched or are not in the expected form). `e$error_code` holds the server's code when it sent one, such as `REPORT_DELETED` or `RATE_LIMIT_EXCEEDED`. An uncaught one prints as `Error: <message>`; `print(e)` shows `CyteType API Error: <message>`.

**Simple R errors** cover the rest. When the server refuses an upload or a job submission, the message names the step and gives the server's explanation, with its code in brackets when the server sent one:

```text
Upload initiate failed: [TOKEN_INACTIVE] API token has expired or been deactivated
Job submission failed: [QUOTA_EXCEEDED] Quota exceeded: This job requires 24 clusters but you only have 10 remaining (990/1000 used)
Job submission failed: [LLM_VALIDATION_FAILED] All provided models failed validation: ...
```

Other messages you may see:

- **Job limit:** the server's message, such as "Daily limit reached (10/10). Resets in about 5 hours. ...", followed by "The limits apply to every job, including jobs that use your own LLM through llm_configs. Wait for the limit to reset before retrying."
- **Job failed:** "Server error: Job failed" from `CyteTypeR()`, "Job failed" from `GetResults()`.
- **Polling stopped:** a `cytetype_api_error` "Stopping after 5 consecutive errors: ..." with "The job may still be running; retrieve its results later with GetResults(job_id = "<job_id>").", or the simple error "Job '...' not found after 5 attempts. Verify the job_id is correct." When the server cannot be reached, `CyteTypeR()` keeps retrying until `timeout_seconds` runs out instead.
- **Results not fetched:** `GetResults()` signals a `cytetype_api_error` with the reason, for example "Error while checking job status: HTTP 410 [REPORT_DELETED] This report has been deleted by the owner" (with `e$error_code` set to `"REPORT_DELETED"`); a failed results request reads "Error while fetching job results: ...". While a job is still running, and for a job ID the server does not know, it returns `NULL` with a message instead.
- **Upload failures:** for example "Upload chunk 3 failed: ..." or "Upload complete failed: ...".
- **Checks before submission:** for example "All infoTags values must be character strings", "n_parallel_clusters must be between 1 and 50", "visibility must be "private", "org", "public" or NULL", "Results with prefix 'cytetype' already exist in this object ...", "Artifact build did not succeed. Set require_artifacts = FALSE to continue without artifacts." and the `llm_configs` checks, such as "provider must be one of: ..." or "Either apiKey or all AWS credentials must be provided". See [Troubleshooting](../support/troubleshooting.md).

(errors-r-example)=
## Catch and handle errors in R

Handle `cytetype_auth_error`, `cytetype_timeout_error` and `cytetype_api_error` first (their fields name the problem), then match server codes in the message of other errors:

```r
annotated <- tryCatch(
  CyteTypeR(
    obj = pbmc,
    prepped_data = prepped_data,
    study_context = "Human PBMC from a healthy donor"
  ),
  cytetype_auth_error = function(e) {
    message("API key problem (", e$error_code, "): ", conditionMessage(e))
    message("Run SetupCyteTypeR(force = TRUE) or LoginCyteTypeR(), then try again.")
    NULL
  },
  cytetype_timeout_error = function(e) {
    # The job keeps running on the server.
    message("Timed out: ", conditionMessage(e))
    message("Fetch the results later with GetResults(job_id = \"", e$job_id, "\").")
    NULL
  },
  cytetype_api_error = function(e) {
    # The job's status or results could not be fetched; the job may still be running.
    code <- if (is.null(e$error_code)) "no code" else e$error_code
    message("CyteType API error (", code, "): ", conditionMessage(e))
    NULL
  },
  error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("INVALID_TOKEN|TOKEN_INACTIVE", msg)) {
      message("The API key was rejected. Run SetupCyteTypeR(force = TRUE).")
    } else if (grepl("limit reached", msg, fixed = TRUE)) {
      message("Job limit reached: ", msg)
    } else if (grepl("QUOTA_EXCEEDED", msg, fixed = TRUE)) {
      message("Quota exceeded: ", msg)
    } else {
      stop(e)
    }
    NULL
  }
)

if (!is.null(annotated)) pbmc <- annotated
```

`CyteTypeR()` returns the annotated Seurat object, so keep its return value. If it stops with an error, your original object is unchanged.
