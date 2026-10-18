---
description: What run() and CyteTypeR() do from upload to results, how to follow a job, what happens when the client stops waiting, and how to fetch results later and rerun jobs.
---

(running-jobs)=
# Run, monitor and rerun jobs

`run()` in Python and `CyteTypeR()` in R submit an annotation job and wait for its results. This page explains each step, how to follow a job while it runs, what happens when the client stops waiting, and how to fetch the results later or annotate the same data again.

## What a run does

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
adata = annotator.run(
    study_context="Peripheral blood mononuclear cells from a healthy human donor, 10x Genomics 3' scRNA-seq.",
    metadata={"title": "PBMC 3k"},
)
```
:::
:::{tab-item} R
:sync: r
```r
pbmc <- CyteTypeR(
  pbmc, prepped_data,
  study_context = "Peripheral blood mononuclear cells from a healthy human donor, 10x Genomics 3' scRNA-seq.",
  metadata = list(title = "PBMC 3k")
)
```
:::
::::

1. **Checks.** The client looks for results already stored under the same `results_prefix`, finds your API key, and validates the study context, the metadata tags, the visibility and any `llm_configs`. A problem here stops the run before anything is uploaded.
2. **Upload.** It uploads `vars.h5` and `obs.duckdb` in chunks, several at a time (see [Large datasets](large-datasets.md#large-datasets-upload)).
3. **Copy of the request.** With `save_query=True` (the default), it writes the request to `query.json`.
4. **Submission.** The server checks your API key, your [usage limits](../platform/limits.md) and the requested visibility, and tests your own models if you passed any. If everything passes, it creates the job and returns its ID.
5. **Job details.** The client records the job ID and the report link. Python stores them in `adata.uns["cytetype_jobDetails"]` straight away; CyteTypeR stores them in the object it returns at the end.
6. **Report link.** The client prints the link to the job's report.
7. **Waiting.** It checks the job's status every `poll_interval_seconds` (10 by default) and shows the progress.
8. **Results.** When the job completes, it downloads the results, adds them to your object and returns the object.

Always keep the returned object: `adata = annotator.run(...)` in Python, `obj <- CyteTypeR(obj, ...)` in R (see [Keep the returned object](../results/results-in-your-data.md#results-return-value)).

## Follow a job

The job runs on CyteType's servers. You can follow it in several places:

- **The report.** Its link, `https://cytetype.nygen.io/report/<job_id>`, is printed as soon as the job is submitted. While the job runs, the report shows the current stage (Preparing, Understanding context, Annotating clusters, Refining annotations, Reviewing, Finalizing); the labels and the evidence appear when the job reaches its final step. Reports are private by default: open the link in a browser where you are signed in to CyteType, or run `cytetype view <job_id>` (R: `ViewCyteTypeJob("<job_id>")`), which opens it through sign-in.
- **The dashboard** lists all your jobs, under the `title` from your [metadata tags](study-context.md#study-context-metadata-tags).
- **Email.** CyteType emails you when the job completes or fails.

In your session, the client shows its progress while it waits (`show_progress=True`, the default):

::::{tab-set}
:::{tab-item} Python
:sync: python
```text
CyteType job submitted.
If your session disconnects, results can still be fetched later with:
`results = annotator.get_results()`

[TRACK PROGRESS]
https://cytetype.nygen.io/report/<job_id>
⠹ CyteType job running... 04:12 elapsed
```

In a terminal, the last line updates in place. In a Jupyter notebook, a single output line updates. When the output goes to a file, as in a batch job, a new line is printed only when the status changes. At the end, the client prints `[DONE] CyteType job completed.` (or `[FAILED] CyteType job failed.`).
:::
:::{tab-item} R
:sync: r
CyteTypeR prints the job ID, the report link and the call that fetches the results if you get disconnected, `GetResults(job_id = "<job_id>")`. Then it prints a line with one symbol per cluster, completed (✔), running (◉), waiting (◯) or failed (✖), and the number of completed clusters. The line updates in place. The clusters are reported as completed together, when the job reaches its final step.
:::
::::

`show_progress=False` turns this display off; messages such as the report link are still printed. In R, it also hides the upload progress bars.

## Waiting and timeouts

The job does not depend on your session. The client only waits for it, and the job keeps running when the client stops waiting for any reason: a timeout, an interrupted call, a restarted notebook kernel, a computer that went to sleep or a lost connection.

- **`timeout_seconds`** (default 7200, so 2 hours) is how long the client waits. When the time runs out, both clients print the report link. Python raises `CyteTypeTimeoutError` (importable from `cytetype`); CyteTypeR stops with an error of class `cytetype_timeout_error`, whose message names the job and the `GetResults(job_id = ...)` call that fetches its results.
- **The server stops a job** that runs for more than 4 hours, and reports it as failed.
- **Connection problems.** While waiting, both clients retry when the server cannot be reached, until `timeout_seconds` runs out. Python stops at once on an error reported by the server; CyteTypeR stops after five such errors in a row, or at once for a problem with your API key.

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
from cytetype import CyteTypeTimeoutError

try:
    adata = annotator.run(study_context="...", timeout_seconds=3600)
except CyteTypeTimeoutError:
    # The job is still running. Its ID is already stored in annotator.adata.uns.
    results = annotator.get_results()  # None until the job has completed
```
:::
:::{tab-item} R
:sync: r
```r
result <- tryCatch(
  CyteTypeR(pbmc, prepped_data, study_context = "...", timeout_seconds = 3600),
  cytetype_timeout_error = function(e) e
)
if (inherits(result, "cytetype_timeout_error")) {
  # The job is still running; the error holds its ID.
  res <- GetResults(job_id = result$job_id)  # NULL until the job has completed
} else {
  pbmc <- result
}
```
:::
::::

### If the job fails

Python raises `JobFailedError` and prints the report link; CyteTypeR stops with `Server error: Job failed`. Failed jobs do not count toward your usage limits. Submit the job again (in Python, from the same object, with `override_existing_results=True` or another `results_prefix`, because the failed job's details are already stored); if it keeps failing, contact support with the job ID, which is in the report link and the email.

A job can also complete with some clusters that have no result. Both clients still add the results to your object, mark the cells of those clusters as `Unknown` and log a warning that names the clusters. In R, these clusters have no row in the results table.

(running-jobs-fetch-later)=
## Fetch results later

Each job's results stay on the server, so you can fetch them after the client stopped waiting, in the same session or a new one. [Results in AnnData and Seurat](../results/results-in-your-data.md#results-per-cluster) describes what is returned.

::::{tab-set}
:::{tab-item} Python
:sync: python
`get_results()` first returns any results stored in the object. If there are none but the object holds job details, it checks the job once: when the job has completed, it downloads the results, adds the columns as `run()` would (to `annotator.adata`) and returns them; while the job is still running, it returns `None`.

**Same session.** Call it on the same `annotator`:

```python
results = annotator.get_results()
if results is None:
    print("The job is still running; try again later.")
```

**New session.** The job details are part of the AnnData, so when the client stops waiting, save the object with `annotator.adata.write_h5ad(...)`. Later, reload it and create a `CyteType` object with the same `group_key` (and `label_na`, if you used it) so that the cluster IDs line up:

```python
import anndata as ad
from cytetype import CyteType

adata = ad.read_h5ad("pbmc_submitted.h5ad")
annotator = CyteType(adata, group_key="leiden", rank_key="rank_genes_leiden")
results = annotator.get_results()
annotator.cleanup()
```

Creating `CyteType` prepares the data again, including the upload files, which `cleanup()` then deletes.

**Only the job ID.** If you lost the object, take the job ID from the report link, the dashboard or the email, and add it to the AnnData before you create `CyteType`:

```python
adata.uns["cytetype_jobDetails"] = {"job_id": "<job_id>"}
annotator = CyteType(adata, group_key="leiden", rank_key="rank_genes_leiden")
results = annotator.get_results()
```

With another `results_prefix`, use the key `f"{results_prefix}_jobDetails"` and call `get_results(results_prefix=...)`.
:::
:::{tab-item} R
:sync: r
`CyteTypeR()` returns the object only when the job completes, so when the call stops early, your object holds no job details. Use the job ID instead: CyteTypeR prints it at submission, together with the `GetResults(job_id = "<job_id>")` call to use, and repeats it when the call times out or is interrupted. It is also in the report link, the dashboard and the email.

**A saved object.** If you saved the object that `CyteTypeR()` returned, the results are in it, in any session:

```r
pbmc <- readRDS("pbmc_cytetype.rds")
GetResults(pbmc)
```

**Only the job ID.**

```r
res <- GetResults(job_id = "<job_id>")
```

- It checks the job once and does not wait. While the job is still running, or when the server does not know the job ID, it returns `NULL` with a message; call it again later. For a failed job, it stops with an error.
- When the job has completed, it returns a data frame with one row per annotated cluster and the same columns as `obj@misc$cytetype_results`; the attributes `"job_id"` and `"result"` hold the job ID and the full results. It also saves the server's response as `cytetypeR_results_<job_id>.json` in the working directory, unless you pass `save_json = FALSE`.
- It does not add columns to a Seurat object, and its `clusterId` column holds CyteType's cluster IDs (1 to N in the sorted order of your labels). To add the labels to your object, map the IDs back to your clusters:

  ```r
  labels <- as.character(sort(unique(pbmc$seurat_clusters)))  # the order CyteTypeR numbered them in
  res$cluster <- labels[as.integer(res$clusterId)]
  pbmc$cytetype_annotation_seurat_clusters <-
    res$annotation[match(as.character(pbmc$seurat_clusters), res$cluster)]
  ```

  This works only while the cluster column is unchanged since the run.
:::
::::

(running-jobs-rerun)=
## Rerun a job

You may want to annotate the same data again, for example with a more detailed study context or with [your own LLM](llm-configuration.md). Runs are not deterministic, so a rerun can also differ in wording or, for ambiguous clusters, in the call itself.

- **Keep both results** by giving the new run another `results_prefix`, such as `"cytetype_v2"`. Its columns are then named `cytetype_v2_annotation_leiden` and so on, next to the first run's.
- **Replace the earlier results** with `override_existing_results=True`.
- **Without either,** `run()` and `CyteTypeR()` stop with an error that results with that prefix already exist. In Python the check looks at the job details, which are stored at submission, so this also applies after a run that stopped waiting.
- **Reuse the prepared object.** The same `annotator` (Python) or `prepped_data` (R) can submit any number of jobs without preparing the data again. Each job uploads the files again.
- **Clean up after the last run.** `cleanup()` and `CleanUpArtifacts()` delete the upload files. After that, both `run()` and `CyteTypeR()` stop with an error because the files are missing, unless you pass `require_artifacts=False` to run without them. Create a new `CyteType` object, or call `PrepareCyteTypeR()` again, to rebuild them.
- **In CyteTypeR,** the results table `obj@misc$cytetype_results` belongs to the default prefix. A run with another prefix, such as `"cytetype_v2"`, stores its results in `obj@misc$cytetype_v2_results` and leaves the table as it is (see [Results in AnnData and Seurat](../results/results-in-your-data.md#results-per-cluster)).

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
adata = annotator.run(study_context="Human PBMC", results_prefix="cytetype")
adata = annotator.run(
    study_context="Peripheral blood mononuclear cells from a healthy human donor, 10x Genomics 3' scRNA-seq.",
    results_prefix="cytetype_v2",
)
annotator.cleanup()  # after the last run
```
:::
:::{tab-item} R
:sync: r
```r
pbmc <- CyteTypeR(pbmc, prepped_data, study_context = "Human PBMC", results_prefix = "cytetype")
pbmc <- CyteTypeR(
  pbmc, prepped_data,
  study_context = "Peripheral blood mononuclear cells from a healthy human donor, 10x Genomics 3' scRNA-seq.",
  results_prefix = "cytetype_v2"
)
CleanUpArtifacts(prepped_data)  # after the last run
```
:::
::::

(running-jobs-options)=
## Other run options

Every argument is described in the [Python API reference](../reference/api/cytetype.html#cytetype.CyteType.run), and in R with `?CyteTypeR`. These are the ones that change how a run behaves.

### `require_artifacts`

By default (`True`), a failure to build or upload `vars.h5` or `obs.duckdb` stops the run. With `require_artifacts=False`, the job runs without the missing files: CyteType still annotates from the marker genes and cluster summaries, but evidence computed from your cells is missing, Chat has less data to work with, and [target analysis](../platform/target-analysis.md) is not available. In Python, a file that was uploaded is still used; CyteTypeR uses the files only when both uploads succeed. Errors with your API key or account (authentication and authorization errors) always stop the run.

Use the default unless an upload keeps failing for a reason you cannot fix, and report the error (see [Troubleshooting](../support/troubleshooting.md)).

### `save_query` and `query_filename`

With `save_query=True` (the default), each run writes the request to `query_filename` (default `query.json` in the working directory), overwriting an earlier file of that name. The Python client never writes LLM API keys to it. In CyteTypeR the file also contains your `llm_configs`, including their keys, so pass `save_query = FALSE` whenever you use your own models (see [Use your own LLM](llm-configuration.md#llm-configuration-keys)).

### `visibility`

`visibility` sets who can open the report: `"private"` (the default), `"org"` or `"public"`. Organization visibility needs an organization API key, issued under your organization's license. The personal key that `cytetype setup` or `SetupCyteTypeR()` saves is refused, even if your account belongs to an organization, and the server checks this only at submission, after the upload. You can also change the visibility in the dashboard at any time. See [Visibility and sharing](../platform/sharing.md).

### `n_parallel_clusters`

The number of clusters annotated at the same time, from 1 to 50 (default 2). It applies only when you use your own models; see [Use your own LLM](llm-configuration.md#llm-configuration-parallel).

## When a submission is refused

| Situation | Python | R |
|---|---|---|
| No API key, or a rejected key | `AuthenticationError` | `CyteType sign-in is required. Run SetupCyteTypeR() first.`, or the server's message |
| Personal usage limit reached | `RateLimitError` | Error with the server's message |
| Organization quota used up | `QuotaExceededError` | Error with the server's message |
| One of your own models failed the test | `LLMValidationError` | `Job submission failed: [LLM_VALIDATION_FAILED] ...` |
| `visibility="org"` with a personal API key | `APIError` (code `HTTP_ERROR`) | `Job submission failed: [HTTP_ERROR] Org visibility requires an organization license` |

A missing or rejected API key is caught before the files are uploaded. The other checks run at submission, so when a job is refused for them, the upload time has already been spent, and the files are uploaded again when you resubmit. The message for a personal limit says which limit you reached and when it resets. See [Usage limits and quotas](../platform/limits.md) for the limits, and [Errors and error codes](../reference/errors.md) for every error.

## Next steps

- [Results in AnnData and Seurat](../results/results-in-your-data.md): what the client adds to your object.
- [The CyteType report](../results/report.md): read the evidence for each cluster.
- [Automation, HPC and CI](automation.md): run jobs from scripts, batch systems and CI.
