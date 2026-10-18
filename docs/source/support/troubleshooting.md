---
description: Fixes for common CyteType problems with signing in, preparing data, building and uploading files, running jobs and installing the clients, from the message you see to the cause and the fix.
---

(troubleshooting)=
# Troubleshooting

Find the message or behavior you see, then read why it happens and how to fix it. For the meaning of each exception class and server error code, see [Errors and error codes](../reference/errors.md). If your problem is not listed, see [Get help](#troubleshooting-help).

(troubleshooting-sign-in)=
## Signing in and API keys

### The browser does not open during setup

**You see:** `Was not able to launch web browser` (R: `Could not launch a web browser.`), or nothing opens.

**Why:** the command could not start a browser, for example on a computer without a desktop, or the browser started somewhere you cannot see it.

**Fix:** keep the command running and copy the link printed under "If it does not open automatically, use this URL:" into a browser on the same computer. If the computer has no browser, see [Sign in on a remote machine](#troubleshooting-remote).

### Signing in from WSL

On Windows Subsystem for Linux, the clients open your Windows default browser. The browser then returns the key to the command through `127.0.0.1` inside WSL, which Windows reaches through WSL's localhost forwarding (on by default).

- If no browser opens, copy the printed link into a browser in Windows.
- If the browser cannot connect to `127.0.0.1` after you select **Authorize CLI**, localhost forwarding is probably turned off in your WSL configuration. Turn it back on, or save a key from the dashboard with `cytetype login` (R: `LoginCyteTypeR()`).

(troubleshooting-remote)=
### Sign in on a remote machine (SSH, HPC, remote notebook)

**Why:** browser sign-in hands the key back to a temporary listener at `127.0.0.1` on the machine that runs the command. A browser on your laptop cannot reach that address on a remote server.

**Fix:** create an API key on the **API keys** page of the dashboard, then save it on the remote machine with `cytetype login` (R: `LoginCyteTypeR()`), or pass it to `run()` or `CyteTypeR()` as `auth_token`. [Sign in and API keys](../getting-started/sign-in.md#sign-in-remote) has the steps.

- `cytetype login` asks for the key in a terminal. A notebook has no terminal to type into, so open one (in JupyterLab: **File**, **New**, **Terminal**), or pass the key in code, for example `auth_token=getpass.getpass("API key: ")`.
- If you prefer the browser sign-in over SSH, you can forward the listener's port. Run `cytetype setup` on the server and note the port: the number after `127.0.0.1%3A` in the printed link. In a second terminal on your own computer, run `ssh -N -L <port>:127.0.0.1:<port> <you>@<server>`, then open the printed link in your local browser within 5 minutes.

### Setup waits and stops with "CyteType setup timed out"

**Why:** the browser step was not finished within 5 minutes, or the browser could not deliver the result to the listener. Security software, a firewall or a browser policy that blocks `http://127.0.0.1` has this effect, and so does a browser running on another computer.

**Fix:** run the command again, keep it running, and finish in a browser on the same computer within 5 minutes. If loopback addresses are blocked on your computer, ask your IT team to allow `127.0.0.1`, or save a key from the dashboard with `cytetype login` (R: `LoginCyteTypeR()`).

### "Please use your university/institute/company email"

**Why:** CyteType accounts need an institutional email address. Addresses from personal email providers (such as Gmail or Outlook) are not accepted.

**Fix:** sign in with your university, institute or company address. If your institutional address is refused, contact contact@nygen.io.

### "CyteType sign-in is required" or "No saved CyteType sign-in found"

**Why:** no API key is saved for the server the client is using. Either you have not signed in on this computer as this user, or the saved sign-in belongs to a different server: it is tied to one server, and your `api_url` or `CYTETYPE_API_URL` points to another. A batch job or container that runs as another user, or with another home directory or `XDG_CONFIG_HOME`, does not see your saved key either.

**Fix:**

- Sign in for the server you use in your code: `cytetype setup`, or `cytetype setup --api-url <address>` for another server (R: `SetupCyteTypeR()`).
- If you set `CYTETYPE_API_URL` by mistake, unset it. In Python, an empty `CYTETYPE_API_URL` is an error too.
- In batch jobs and containers, pass the key as `auth_token` instead.

In Python, the note printed by `CyteType()` names the server it looked for, and `annotator.api_url` shows the server in use.

### "Invalid CyteType credentials file"

**Why:** the saved `credentials.json` is damaged or incomplete, for example after a manual edit or a full disk.

**Fix:** delete it and sign in again.

::::{tab-set}
:::{tab-item} Python
:sync: python
```bash
cytetype logout
cytetype setup
```

`cytetype setup --force` also works: it does not read the old file.
:::
:::{tab-item} R
:sync: r
```r
LogoutCyteTypeR()
SetupCyteTypeR()
```

`SetupCyteTypeR(force = TRUE)` also works: it does not read the old file.
:::
::::

### "Saved CyteType credentials could not be read"

**Why (Python):** the client could not reach its configuration folder, for example in a container or batch job without a home directory, or with a home directory you cannot read.

**Fix:** pass the key as `auth_token`, or set `XDG_CONFIG_HOME` to a readable folder that contains `cytetype/credentials.json`.

### The saved key is rejected

**You see:** `Invalid API token` or `API token has expired or been deactivated`.

**Why:** the key was disabled in the dashboard, by you or, for an organization key, by an admin. All computers that signed in with browser sign-in share your "CyteType CLI" key, so disabling it on one signs out all of them. A key from another CyteType server is also rejected.

::::{tab-set}
:::{tab-item} Python
:sync: python
`cytetype setup` stops with "... Run `cytetype setup --force` to re-authenticate with a new key". `run()` raises `AuthenticationError` with "Artifact upload was rejected: ...". Get a new key:

```bash
cytetype setup --force
```

Or save another active key with `cytetype login`. `run()` reads the saved key every time, so you can call it again on the same `CyteType` object without preparing the data again.
:::
:::{tab-item} R
:sync: r
`SetupCyteTypeR()` stops with "... Run SetupCyteTypeR(force = TRUE) or `cytetyper setup --force` to re-authenticate with a new key". `CyteTypeR()` stops with "Upload initiate failed: [TOKEN_INACTIVE] API token has expired or been deactivated". Get a new key:

```r
SetupCyteTypeR(force = TRUE)
```

Or save another active key with `LoginCyteTypeR()`. Then run `CyteTypeR()` again with the same `prepped_data`.
:::
::::

### "CyteType credentials directory is not owned by the current user"

**Why:** the credentials folder was created by another user, for example by running `sudo cytetype setup`, or `XDG_CONFIG_HOME` points to a shared folder.

**Fix:** give the folder back to your user (for example `sudo chown -R "$USER" ~/.config/cytetype`), or set `XDG_CONFIG_HOME` to a folder of your own.

### The `cytetype` or `cytetyper` command is not found

- **Python:** activate the environment where `cytetype` is installed, or run the commands through Python, for example `python -m cytetype.cli setup`. If `python -m cytetype.cli` fails too, you may have an old release; see [An old version is installed](#troubleshooting-old-version).
- **R:** install the launcher with `InstallCyteTypeRCli()`. It goes into the folder of `Rscript`, which must be on your `PATH`. The launcher starts R with `--vanilla`, so it ignores `.Rprofile` and `.Renviron`; if CyteTypeR is installed only in a project library such as renv, it stops with `Error: Unable to load CyteTypeR`. You can always use the R functions instead (`SetupCyteTypeR()` and so on). See [the `cytetyper` command](../getting-started/installation.md#installation-cytetyper).

### The server cannot be reached

**You see:** `Could not validate the saved API key`, `NetworkError: Connection failed: ...`, or uploads that fail with connection errors.

**Why:** the machine has no connection to `https://cytetype.nygen.io`. Typical cases are compute nodes without internet access, firewalls, and networks that require a proxy.

**Fix:** run CyteType on a machine with outbound HTTPS access. If your network uses a proxy, set the standard `HTTPS_PROXY` environment variable; both clients use it.

(troubleshooting-data)=
## Preparing your data

### "Gene Symbol Detection Error"

**Why (Python):** about half or more of the gene names that `CyteType()` selected look like database IDs (Ensembl, RefSeq, Entrez or similar accession numbers). CyteType needs readable gene symbols, such as `CD3E`.

**Fix:** pass `gene_symbols_column` with the name of the `adata.var` column that holds symbols, or convert the IDs to symbols first. Without it, CyteType looks for symbols in the columns `feature_name`, `gene_symbols`, `gene_symbol`, `gene_short_name`, `gene_name` and `symbol`, then in `adata.var_names`, then in the other columns. Values that combine a symbol and an ID, such as `TSPAN6_ENSG00000000003`, are split automatically. A `KeyError` "Column '...' not found in `adata.var`" means the name you passed does not exist.

**In R:** CyteTypeR only warns and continues, which makes the annotation less reliable: `More than half of the gene values (...%) appear to be gene IDs` for the gene names of the `RNA` assay, or `...% of the genes in the marker table look like gene IDs`. Use gene symbols as the row names of the Seurat object before you compute markers. See [Prepare your data](../guides/prepare-data.md#prepare-data-gene-symbols).

### "`rank_genes_groups` run with groupby='...', expected '...'"

**Why:** the marker genes stored under `rank_key` were computed for another cluster column than `group_key`.

**Fix:** compute markers for the same column and pass their key:

```python
sc.tl.rank_genes_groups(adata, groupby="leiden", key_added="rank_genes_leiden")
annotator = CyteType(adata, group_key="leiden", rank_key="rank_genes_leiden")
```

A `KeyError` "'rank_genes_groups' not found in `adata.uns`" means no markers are stored under that key: run `sc.tl.rank_genes_groups` first, or pass the right `rank_key`.

### Marker genes are lost

**You see:** warnings "Number of lost genes (...) for group ..." or the error "All marker gene lists are empty. Gene names in rank_genes_groups could not be matched to adata.var_names."

**Why:** the marker genes do not match `adata.var_names`. Usually `var_names` were changed after the markers were computed (for example from Ensembl IDs to symbols), or the markers were ranked on `adata.raw` (Scanpy's default when `adata.raw` is set), which holds genes that are no longer in `adata.var`, for example after subsetting to highly variable genes.

**Fix:** compute the markers again on the same AnnData object you pass to `CyteType()`, after any change to its genes (with `use_raw=False` if `adata.raw` holds other genes).

### "Internal inconsistency: Group name '...' from rank_genes_groups results was not found"

**Why:** the clusters changed after the markers were computed: you clustered again into the same column, renamed clusters, or removed cells so that a cluster no longer exists.

**Fix:** compute the markers again after changing the clusters.

### Cells without a cluster label

**You see (Python):** `... cells (...%) have NaN values in 'leiden'. Either fix the data or set label_na=True ...`

**You see (R):** `... cells (...%) have NA values in 'seurat_clusters'. Assign these cells to a cluster, or remove them (for example with subset()), before calling PrepareCyteTypeR().`

**Why:** some cells have no cluster label, for example after merging objects or transferring labels.

::::{tab-set}
:::{tab-item} Python
:sync: python
Give those cells a label, remove them, or let CyteType put them in a cluster named "Unknown" (or "Unknown 2" if "Unknown" exists):

```python
annotator = CyteType(adata, group_key="leiden", label_na=True)
adata = annotator.run(study_context="...")
```

With `label_na=True`, `run()` returns a new AnnData object, and only that object holds the results. Always assign the return value.
:::
:::{tab-item} R
:sync: r
CyteTypeR has no `label_na` option. Remove the cells without a label, or give them a label, before `FindAllMarkers()` and `PrepareCyteTypeR()`:

```r
pbmc <- subset(pbmc, cells = colnames(pbmc)[!is.na(pbmc$seurat_clusters)])
```
:::
::::

### No embedding in the report

**You see (Python):** `No suitable 2D coordinates found in adata.obsm ... Visualization will be disabled.` (R: `No suitable 2D coordinates found in the object's reductions. ... Visualization will be disabled.`)

**Why:** the client found no 2D embedding. Python tries `coordinates_key`, then `X_umap`, `X_tsne` and `X_pca` in `adata.obsm` (in any capitalization) and needs one row per cell and at least 2 columns. R tries `coordinates_key`, then `umap`, `tsne` and `pca` among the object's reductions, in the same way.

**Fix:** the job still runs, but the report has no embedding plot. Compute a UMAP before preparing the data (`sc.tl.umap(adata)`, R: `RunUMAP()`), or pass the name of the embedding you have as `coordinates_key` (R: one of `names(pbmc@reductions)`).

### "No integer raw counts found"

**You see (Python):** `No integer raw counts found in adata.layers['counts'], adata.raw.X, or adata.X. Skipping raw counts in vars.h5.` (R: `The raw counts are not integer-valued. Skipping raw counts in vars.h5.`)

**Why:** `CyteType()` includes raw counts in `vars.h5` when it finds a matrix of whole numbers in `adata.layers["counts"]`, then `adata.raw.X` (only if it contains all of `adata.var_names`), then `adata.X`. It found none. CyteTypeR includes the `counts` layer of the `RNA` assay only when every value in it is a whole number.

**Fix:** annotation works without counts, but features that need them, such as [target analysis](../platform/target-analysis.md), are not available for the job. To include counts, keep a copy before you normalize:

```python
adata.layers["counts"] = adata.X.copy()  # before sc.pp.normalize_total and sc.pp.log1p
```

In R, create the Seurat object from the raw counts with `CreateSeuratObject(counts = ...)`; `NormalizeData()` keeps them in the `counts` layer.

(troubleshooting-artifacts)=
## Building and uploading files

`CyteType()` (R: `PrepareCyteTypeR()`) writes two files, `vars.h5` (expression) and `obs.duckdb` (cell metadata and coordinates), and `run()` (R: `CyteTypeR()`) uploads them before it submits the job. [How CyteType works](../getting-started/how-cytetype-works.md) explains what they are used for.

### A file could not be built or uploaded

**You see (Python):** a warning such as `vars.h5 artifact failed during build: ...` when you create `CyteType`, or `obs.duckdb upload failed: ...`, then `Artifact build/upload failed for: ... Rerun with require_artifacts=False to skip this error.` and the original error from `run()`.

**You see (R):** `Error building artifacts: ...` from `PrepareCyteTypeR()` and `Artifact build did not succeed. Set require_artifacts = FALSE to continue without artifacts.` from `CyteTypeR()`, or `Uploading artifacts failed: ...`.

**Why:** with the default `require_artifacts=True` (R: `TRUE`), any problem building or uploading the files stops the run, so that you do not get a report with features missing without noticing.

**Fix:** read the original error at the end of the message; the sections below cover the common causes. If you need the annotation now, run with `require_artifacts=False` (R: `require_artifacts = FALSE`). The job then runs without the failed file; in Python, a file that was built and uploaded is still used, while CyteTypeR runs without either file when one of them fails. Without the files, the evidence that CyteType computes from your cells is missing from the report, Chat has less data to work with, and target analysis is not available (see [`require_artifacts`](../guides/running-jobs.md#running-jobs-options)). In Python, a rejected API key or a refused account stops the run even with `require_artifacts=False`. If a failure looks like a bug, please report it on GitHub (see [Get help](#troubleshooting-help)).

### "Failed to cast value: Could not convert string ..."

**Why (Python):** a column of `adata.obs` with the generic `object` type mixes numbers and text, for example donor ages with one `"unknown"`. The database file infers the column type (such as `INT32` or `DOUBLE`) from the first rows and fails when the text appears, so `obs.duckdb` is not built.

**Fix:** convert such columns to text before you create `CyteType`:

```python
adata.obs["donor_age"] = adata.obs["donor_age"].astype(str)
```

### Full disk, read-only folder or no temporary space

**You see:** `[Errno 28] No space left on device`, `Permission denied` or `Read-only file system` while the files are built, or when `query.json` is saved.

**Why:** the files are written to the current working directory by default, and `vars.h5` can be large: it holds the whole normalized matrix and, when found, the raw counts. Building `obs.duckdb` can also write temporary files. `run()` saves `query.json` in the working directory after the upload and before the job is submitted.

::::{tab-set}
:::{tab-item} Python
:sync: python
Write the files to a disk with space, and remove them after the run:

```python
annotator = CyteType(
    adata,
    group_key="leiden",
    vars_h5_path="/scratch/me/vars.h5",
    obs_duckdb_path="/scratch/me/obs.duckdb",
)
adata = annotator.run(
    study_context="...",
    query_filename="/scratch/me/query.json",  # or save_query=False
)
annotator.cleanup()  # deletes vars.h5 and obs.duckdb
```

Temporary files for `obs.duckdb` go to the system temporary folder. Set the `TMPDIR` environment variable before you start Python to use another folder.
:::
:::{tab-item} R
:sync: r
Write the files to a disk with space, and remove them after the run:

```r
prepped_data <- PrepareCyteTypeR(
  pbmc,
  pbmc.markers,
  vars_h5_path = "/scratch/me/vars.h5",
  obs_duckdb_path = "/scratch/me/obs.duckdb"
)
pbmc <- CyteTypeR(
  obj = pbmc,
  prepped_data = prepped_data,
  study_context = "...",
  query_filename = "/scratch/me/query.json"  # or save_query = FALSE
)
CleanUpArtifacts(prepped_data)  # deletes vars.h5 and obs.duckdb
```

Temporary files for `obs.duckdb` go to R's temporary folder (`tempdir()`). Set the `TMPDIR` environment variable before you start R to use another folder.
:::
::::

### Running out of memory with large data

**Why:** when `adata.X` is held in memory, building `vars.h5` needs the whole matrix plus working copies in memory. `CyteType()` prints `For large datasets, use AnnData backed mode ...` whenever `adata.X` is in memory; for small datasets you can ignore this.

**Fix (Python):** open the file in backed mode, `sc.read_h5ad("data.h5ad", backed="r")`, compute markers with `cytetype.rank_genes_groups_backed()`, and lower `pcent_batch_size` if needed. CyteTypeR has no backed mode: give R enough memory for the object. See [Large datasets](../guides/large-datasets.md).

### A file is larger than the upload limit

**You see:** `obs_duckdb exceeds upload limit: ...` or `vars_h5 exceeds upload limit: ...`.

**Why:** `obs.duckdb` can be up to 2 GiB and `vars.h5` up to 50 GiB (see [Usage limits and quotas](../platform/limits.md)).

**Fix:** `obs.duckdb` holds every column of `adata.obs` (R: `meta.data`), so drop columns you do not need, especially long per-cell text. For `vars.h5`, reduce the data, for example by subsampling cells per cluster; see [Large datasets](../guides/large-datasets.md).

### Uploads are slow or fail part-way

Files are uploaded in parts, several at a time: up to `upload_max_workers` parts (Python default 4, R default 6). Both clients retry a part that fails with a network error, a timeout or a server error up to three times, after 1, 5 and 20 seconds; Python logs each retry as `Chunk ... upload failed (attempt .../4), retrying in ...`.

- On a slow or unstable connection, lower `upload_max_workers` and raise `upload_timeout_seconds`, the time allowed for each request (default 3600 seconds).
- In R, parts are uploaded in parallel only when the furrr package is installed, and also progressr unless you pass `show_progress = FALSE`; otherwise one at a time.
- `Presigned URL upload rejected (HTTP 403)` (R: `Presigned upload rejected with HTTP 403`) means the file storage service rejected a part, for example when an upload has been running for many hours. Run again, ideally from a faster connection.

### "The artifact files were deleted by cleanup()"

**Why (Python):** `cleanup()` deleted `vars.h5` and `obs.duckdb`, and `run()` needs them.

**Fix:** create a new `CyteType` object to build them again, or pass `require_artifacts=False` to run without them.

(troubleshooting-jobs)=
## Submitting and running jobs

### `metadata` values are rejected

::::{tab-set}
:::{tab-item} Python
:sync: python
**You see:** `pydantic.ValidationError` ... `infoTags.<key> Input should be a valid string`.

Every value in `metadata` must be a string:

```python
metadata = {"n_donors": 4, "GEO": "GSE123456"}
adata = annotator.run(
    study_context="...",
    metadata={key: str(value) for key, value in metadata.items()},
)
```
:::
:::{tab-item} R
:sync: r
**You see:** `All infoTags values must be character strings`.

Every value in `metadata` must be a character string:

```r
metadata <- list(n_donors = 4, GEO = "GSE123456")
pbmc <- CyteTypeR(
  obj = pbmc,
  prepped_data = prepped_data,
  study_context = "...",
  metadata = lapply(metadata, as.character)
)
```
:::
::::

### "Results with prefix 'cytetype' already exist"

**Why:** the object already holds a CyteType job under this `results_prefix`. In Python, the job ID is stored as soon as the job is submitted, so a run that later timed out, failed or was interrupted counts too.

**Fix:** to get the results of that job, call `annotator.get_results()` (R: `GetResults(pbmc)`). To keep both runs, pass a new prefix such as `results_prefix="cytetype_v2"`. To replace the old results, pass `override_existing_results=True` (R: `override_existing_results = TRUE`).

### `llm_configs` is rejected

- **Before upload (Python):** a `pydantic.ValidationError` such as "Either apiKey or all AWS credentials must be provided", or a provider name that is not in the list. **R:** "provider must be one of: ..." or "Either apiKey or all AWS credentials must be provided".
- **At submission:** `LLMValidationError` (R: "Job submission failed: [LLM_VALIDATION_FAILED] ..."). A model could not be reached or does not support tool calling. See [LLM_VALIDATION_FAILED](../reference/errors.md#errors-llm-validation) and [Use your own LLM](../guides/llm-configuration.md).

### Job limit reached

**You see:** `RateLimitError` with a message such as `Daily limit reached (10/10). Resets in about 5 hours. ...` (R: the same message, then "The limits apply to every job, including jobs that use your own LLM through llm_configs. Wait for the limit to reset before retrying.").

**Why:** personal accounts can run 10 jobs in any rolling 24 hours, 50 in 7 days and 100 in 30 days. The limit is checked when the job is submitted, after the files are uploaded. The limits also apply when you use your own LLM.

**Fix:** wait until the time in the message, or email contact@nygen.io about a free academic license with higher limits. See [Usage limits and quotas](../platform/limits.md).

### Quota exceeded or license expired

`QuotaExceededError` and `AuthorizationError` with `LICENSE_EXPIRED` concern organization licenses and per-key limits. An expired license stops the run as soon as the files start to upload; a used-up quota is detected when the job is submitted, after the upload. See [QUOTA_EXCEEDED](../reference/errors.md#errors-quota-exceeded), [LICENSE_EXPIRED](../reference/errors.md#errors-license-expired) and [Usage limits and quotas](../platform/limits.md).

### "Org visibility requires an organization license"

**You see (R):** `Job submission failed: [HTTP_ERROR] Org visibility requires an organization license`

**Why:** you passed `visibility="org"` (R: `visibility = "org"`) with an API key that does not belong to an organization license. Only organization API keys (starting with `cyt_`, created by an organization admin) can submit organization reports. A personal key, including the one that `cytetype setup` or `SetupCyteTypeR()` saves, is refused even if your account belongs to an organization. The server checks this when the job is submitted, after the files are uploaded. Python raises an `APIError` with the code `HTTP_ERROR`.

**Fix:** run `run()` (R: `CyteTypeR()`) again with an organization key as `auth_token`, or without `visibility` (the report is then private). The files are uploaded again. See [Visibility and sharing](../platform/sharing.md).

### The client stopped waiting, or your session ended

The job keeps running on the server until it finishes (at most 4 hours) whatever happens to your session. Its report link, printed when the job starts, shows the progress.

::::{tab-set}
:::{tab-item} Python
:sync: python
**You see:** `CyteTypeTimeoutError: Job ... did not complete within 7200s`.

`run()` waits `timeout_seconds` (default 7200, 2 hours). The job ID is already stored in `adata.uns`, so in the same session you can fetch the results later from the same object:

```python
results = annotator.get_results()  # None while the job is still running
```

When the job has finished, `get_results()` adds the same columns to `adata.obs` as `run()`. It makes one check and does not wait, so call it again later if it returns `None`. Next time, pass a larger `timeout_seconds`.

If the Python session ended but you saved the AnnData object after the job started (for example with `annotator.adata.write_h5ad(...)`), load it, create `CyteType` again with the same `group_key` and `rank_key` (and `label_na`, if you used it), and call `get_results()`; creating `CyteType` prepares the data again, which takes some time. You can always open the report from the dashboard or with `cytetype view <JOB_ID>`. See [Run, monitor and rerun jobs](../guides/running-jobs.md).
:::
:::{tab-item} R
:sync: r
**You see:** `Job <JOB_ID> did not finish within 7200 seconds (timeout_seconds). The job keeps running on the server; retrieve its results later with GetResults(job_id = "<JOB_ID>").`

`CyteTypeR()` waits `timeout_seconds` (default 7200, 2 hours), then stops with an error of class `cytetype_timeout_error`. As when you interrupt it or the R session ends, the annotated object is not returned, and your original object holds no job details. Take the job ID from the message, from the line `CyteType job (id: ...) submitted.` or from the dashboard, then:

- open the report with `ViewCyteTypeJob("<JOB_ID>")`; **Download as CSV** in the report gives a table whose `author_label` column holds your original cluster labels, so you can join it to `pbmc@meta.data`;
- or fetch the results table with `GetResults(job_id = "<JOB_ID>")`. It returns `NULL` with a message while the job is still running. Once the job has finished, it returns the table and also writes the raw results to `cytetypeR_results_<JOB_ID>.json` in the working directory.

Next time, pass a larger `timeout_seconds`. See [Run, monitor and rerun jobs](../guides/running-jobs.md).
:::
::::

(troubleshooting-job-failed)=
### The job failed

**You see:** `JobFailedError: Job ... failed` (R: `Server error: Job failed`).

**Fix:** open the report link to see how far the job got, then run the job again; in Python, pass `override_existing_results=True` (or another `results_prefix`) to reuse the same object. Failed jobs do not count toward the personal job limits. If it fails again, contact support with the job ID.

### `get_results()` returns `None`

The log line before it says why:

- `Job ... is still processing. Results not yet available.` (or `still pending`): try again later.
- `Job ... results not yet available (404).`: the server does not know this job ID. Check the ID and that the client uses the server the job ran on.
- `Job ... failed`: see [The job failed](#troubleshooting-job-failed).
- `No results found locally and no job details available for retrieval.`: this object has no job under this `results_prefix`. Check the prefix.
- `Failed to retrieve results for job ...: <reason>`: the reason is a server message, for example from [ACCESS_DENIED](../reference/errors.md#errors-access-denied) or `REPORT_DELETED` (the owner deleted the report).

**In R,** `GetResults()` returns `NULL` with a message when the job is still pending or processing (`Job ... is still processing; its results are not available yet. ...`): try again later. It also returns `NULL` when the server does not know the job ID (`Job ... was not found on ...`): check the ID and the `api_url`. For a failed job or a deleted report, it stops with an error instead.

(troubleshooting-r)=
## CyteTypeR

### The marker table is rejected

| Message | Why | Fix |
|---|---|---|
| `Please check if cluster labels are consistent between marker table and seurat obj!` | The `cluster` values of the marker table are not the same set as the values of `group_key`. Common causes: markers computed for another clustering; a filter that removed every marker of a cluster; cells removed after `FindAllMarkers()`. | Compute markers with `FindAllMarkers()` on the same object and clustering, and filter them so that every cluster keeps some. |
| `Invalid marker genes, some clusters have fewer than 5 markers` | After `n_top_genes`, some clusters have fewer than 5 markers, usually because of a strict filter such as `avg_log2FC > 1`. | Loosen the filter, or check those clusters. |
| An error about a missing `cluster`, `gene` or `avg_log2FC` column | The table lacks a required column. | Use the `FindAllMarkers()` output, or rename your columns to `cluster`, `gene` and `avg_log2FC`. |

### Problems in CyteTypeR 0.5.1 and earlier

In CyteTypeR 0.5.1 and earlier, you may see these messages and problems. Upgrade to CyteTypeR 0.6.0, which fixes them (see [Installation](../getting-started/installation.md)):

- `Package 'duckdb' is required to build obs.duckdb`: duckdb is not installed with the package.
- `GetAssayData doesn't work for multiple layers in v5 assay`: the `RNA` assay of a Seurat v5 object is split into layers.
- `feature_names length (...) does not match matrix columns (...)`: the default assay is not `RNA`.
- `duplicate 'row.names' are not allowed` at the end of `CyteTypeR()`: two or more clusters have no result.
- `Error: API key is required` from `cytetyper login`.
- `CyteTypeR()` does not stop at `timeout_seconds`, the evidence columns of `cytetype_results` are empty, there is no reviewer confidence column, and the job limit message wrongly suggests that your own LLM avoids the limits.

(troubleshooting-install)=
## Installation and startup

(troubleshooting-old-version)=
### An old version is installed

**You see:** `cytetype --version` or `cytetype.__version__` shows an old release, the `cytetype` command does not exist, or a job is refused with "Authentication is required. Upgrade to the latest CyteType or CyteTypeR and run setup, or provide a valid API token."

**Why:** pip only offers releases that support the Python you run it with, and it does not warn you. On Python 3.11 it installs 0.12.0, the last release for that version, which predates the `cytetype` command and the saved sign-in. Older CyteTypeR versions likewise lack `SetupCyteTypeR()` and the other sign-in functions.

**Fix:** check `python --version`, then install the current release into an environment with Python 3.12 or newer (`pip install -U cytetype`). For R, update CyteTypeR from GitHub. See [Installation](../getting-started/installation.md) and its [Python notes](../getting-started/installation.md#installation-python-notes).

### Installing CyteTypeR fails

The most common causes are an R version older than 4.2, which CyteTypeR requires, and building the Bioconductor package rhdf5 from source. See the [R notes](../getting-started/installation.md#installation-r-notes) on the installation page.

### Importing the package pauses for a moment

**Why:** `import cytetype` (and each `cytetype` command) asks PyPI whether a newer version exists, and `library(CyteTypeR)` asks GitHub. Each check waits at most about 2 seconds for an answer. On a machine without internet access, the import pauses until the check gives up, then continues silently.

**Fix:** nothing is needed; the check does not affect your analysis. If a newer version exists, you see `A newer CyteType version is available on PyPI` (R: `A newer version of CyteTypeR is available`); see [Installation](../getting-started/installation.md) to update.

(troubleshooting-help)=
## Get help

If your problem is not listed here:

- Check [Errors and error codes](../reference/errors.md) and the [FAQ](faq.md).
- Open an issue on GitHub: [CyteType](https://github.com/NygenAnalytics/CyteType/issues) for the Python client, [CyteTypeR](https://github.com/NygenAnalytics/CyteTypeR/issues) for the R package.
- Email contact@nygen.io (the address behind **Contact support** in the dashboard).
- Ask the community on [Discord](https://discord.gg/V6QFM4AN).

Include:

- the client and its version (`cytetype --version`, or `packageVersion("CyteTypeR")` in R) and your Python or R version and operating system;
- the full error message and traceback;
- the job ID or report link, if a job was submitted.

Never include your API key, your `credentials.json` file, or a CyteTypeR `query.json` that contains LLM API keys.
