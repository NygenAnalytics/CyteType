---
description: Run CyteType from scripts, HPC batch jobs and CI, with an API key from a secret store, writable file paths, wall-time limits, saved results and many datasets.
---

(automation)=
# Automation, HPC and CI

CyteType works the same way in a script as in a notebook. This page covers what is different when no one is at the keyboard: how the job gets an API key, where files go, how to fit a job into a wall-time limit, how to save the results, and how to run many datasets.

Before you start, check that:

- **The job has an API key** (see [below](#automation-api-keys)).
- **The machine can reach CyteType.** The process that calls `run()` or `CyteTypeR()` uploads the files and then checks the job until it finishes, so it needs outbound HTTPS access the whole time, to `https://cytetype.nygen.io` and to the storage addresses it hands out for uploads. On clusters whose compute nodes have no internet access, ask your administrators which nodes or partitions can reach external services.
- **The files have somewhere to go:** a writable directory with enough space (see [Where files are written](#automation-files)).
- **The wall time covers the wait** (see [Long jobs and wall-time limits](#automation-long-jobs)).

(automation-api-keys)=
## API keys in non-interactive environments

Browser sign-in does not work in batch jobs or CI. Instead, create a separate key for each CI system or cluster in the [dashboard](../platform/dashboard.md) (**API keys**, then **Create API key**; the key is shown only once). Do not reuse the key that browser sign-in saved: browser sign-in gives every computer the same "CyteType CLI" key, so disabling it signs them all out (see [Sign in and API keys](../getting-started/sign-in.md#sign-in-remote)). Then use the key in one of two ways.

**Pass the key in code.** Store the key in your CI system's secrets or another secret store, expose it to the job as an environment variable, and pass it as `auth_token`:

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
import os

from cytetype import CyteType

annotator = CyteType(
    adata,
    group_key="leiden",
    rank_key="rank_genes_leiden",
    auth_token=os.environ["CYTETYPE_API_TOKEN"],
)
```

Passing the key to `CyteType()` rather than to `run()` also avoids the note about a missing saved sign-in.
:::
:::{tab-item} R
:sync: r
```r
obj <- CyteTypeR(
  obj, prepped_data,
  study_context = "...",
  auth_token = Sys.getenv("CYTETYPE_API_TOKEN")
)
```

If the variable is not set, `Sys.getenv()` returns an empty string, and CyteTypeR stops with `auth_token must be a non-empty character string`.
:::
::::

`CYTETYPE_API_TOKEN` is a name you choose; the clients never read it on their own.

**Save the key once (HPC).** On a login node, run `cytetype login` (R: `LoginCyteTypeR()`) and paste a key from the dashboard. The key is saved in `~/.config/cytetype/credentials.json`, or under `$XDG_CONFIG_HOME/cytetype/` when that variable is set. Batch jobs on compute nodes find it there when they share your home directory and see the same `XDG_CONFIG_HOME`; they then need no `auth_token`. See [Sign in and API keys](../getting-started/sign-in.md#sign-in-remote) for details and for [the order in which the clients choose a key](../getting-started/sign-in.md#sign-in-key-order).

## A minimal batch script

These scripts annotate one dataset and save the result. They take an input and an output path, write their working files next to the output, and name each file after the output, so several copies can run at the same time.

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
"""Annotate one dataset: python annotate.py input.h5ad output.h5ad"""
import os
import sys

import scanpy as sc
from cytetype import CyteType, CyteTypeTimeoutError

input_path, output_path = sys.argv[1], sys.argv[2]
name = os.path.splitext(os.path.basename(output_path))[0]
work_dir = os.path.dirname(os.path.abspath(output_path))

# Log-normalized expression in adata.X, clusters in adata.obs["leiden"].
adata = sc.read_h5ad(input_path)
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", key_added="rank_genes_leiden")

annotator = CyteType(
    adata,
    group_key="leiden",
    rank_key="rank_genes_leiden",
    vars_h5_path=os.path.join(work_dir, f"{name}_vars.h5"),
    obs_duckdb_path=os.path.join(work_dir, f"{name}_obs.duckdb"),
    auth_token=os.environ["CYTETYPE_API_TOKEN"],  # leave out to use a saved key
)

try:
    annotator.run(
        study_context="Peripheral blood mononuclear cells from a healthy human donor, 10x Genomics 3' scRNA-seq.",
        metadata={"title": name},
        query_filename=os.path.join(work_dir, f"{name}_query.json"),
        timeout_seconds=3 * 60 * 60,
    )
except CyteTypeTimeoutError:
    print("The job is still running; fetch the results later with get_results().", file=sys.stderr)

annotator.cleanup()
# The saved file holds the results, or the job details needed to fetch them later.
annotator.adata.write_h5ad(output_path)
```
:::
:::{tab-item} R
:sync: r
```r
#!/usr/bin/env Rscript
# Annotate one dataset: Rscript annotate.R input.rds output.rds
library(Seurat)
library(CyteTypeR)

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]
name <- tools::file_path_sans_ext(basename(output_path))
work_dir <- dirname(normalizePath(output_path, mustWork = FALSE))

# Normalized RNA assay, clusters in meta.data$seurat_clusters.
obj <- readRDS(input_path)
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(obj, only.pos = TRUE)

prepped_data <- PrepareCyteTypeR(
  obj, markers,
  group_key = "seurat_clusters",
  vars_h5_path = file.path(work_dir, paste0(name, "_vars.h5")),
  obs_duckdb_path = file.path(work_dir, paste0(name, "_obs.duckdb"))
)

annotated <- tryCatch(
  CyteTypeR(
    obj, prepped_data,
    study_context = "Peripheral blood mononuclear cells from a healthy human donor, 10x Genomics 3' scRNA-seq.",
    metadata = list(title = name),
    auth_token = Sys.getenv("CYTETYPE_API_TOKEN"),  # leave out to use a saved key
    query_filename = file.path(work_dir, paste0(name, "_query.json")),
    timeout_seconds = 3 * 60 * 60,
    show_progress = FALSE
  ),
  cytetype_timeout_error = function(e) {
    # The job is still running: log the message and save the job ID for later.
    message(conditionMessage(e))
    writeLines(e$job_id, file.path(work_dir, paste0(name, "_job_id.txt")))
    NULL
  }
)

CleanUpArtifacts(prepped_data)
if (!is.null(annotated)) saveRDS(annotated, output_path)
```

`show_progress = FALSE` keeps the log readable: otherwise CyteTypeR rewrites its progress line twice a second. When the job is still running after `timeout_seconds`, the script saves its ID next to the output instead of the annotated object (see [Long jobs and wall-time limits](#automation-long-jobs)).
:::
::::

Submit the script with your scheduler, for example with `sbatch` on Slurm, and give it enough memory for the preparation (see [Large datasets](large-datasets.md)) and enough wall time for the wait.

(automation-files)=
## Where files are written

By default, the clients write into the working directory. In a shared or read-only directory, set the paths explicitly:

| File | Default location | Change it with |
|---|---|---|
| `vars.h5`, `obs.duckdb` | Working directory | `vars_h5_path` and `obs_duckdb_path` in `CyteType()` or `PrepareCyteTypeR()` |
| `query.json` | Working directory | `query_filename`, or `save_query=False` to skip it |
| `cytetypeR_results_<job_id>.json`, written by `GetResults(job_id = ...)` in R | Working directory | `save_json = FALSE` to skip it; to move it, change the working directory with `setwd()` |
| Temporary DuckDB files while `obs.duckdb` is written | The system temporary directory (in R, under the session's `tempdir()`) | `TMPDIR`, set before Python or R starts |
| The saved API key, which the clients read during a run | `~/.config/cytetype/credentials.json` | `XDG_CONFIG_HOME` |

- **Give every run its own paths.** The default file names are fixed, so runs that start in the same directory at the same time, such as the tasks of an array job, overwrite each other's files. Include the dataset name in each path, as the scripts above do.
- **Clean up.** `cleanup()` and `CleanUpArtifacts()` delete the upload files. A run that fails leaves them on disk; delete them yourself when you no longer need them.
- **Mind the size.** The upload files can be large; [Large datasets](large-datasets.md#large-datasets-disk) explains how to estimate them.

(automation-long-jobs)=
## Long jobs and wall-time limits

The annotation runs on CyteType's servers, and your batch job only waits for it. Its wall time must cover the preparation, the upload and the wait. The server stops a job after 4 hours, and CyteType emails you when a job completes or fails.

::::{tab-set}
:::{tab-item} Python
:sync: python
Set `timeout_seconds` (default 2 hours) to fit within your wall time. When it runs out, `run()` raises `CyteTypeTimeoutError`; the job keeps running, and the object already holds its job details. The script above saves the object in that case, so a later job can fetch the results.

With short wall times, submit and save right away, then fetch the results in a later job:

```python
# Job 1: submit, wait at most 5 minutes, save.
try:
    annotator.run(study_context="...", timeout_seconds=300)
except CyteTypeTimeoutError:
    pass  # still running; the job details are stored in annotator.adata.uns
annotator.cleanup()
annotator.adata.write_h5ad("sample1_submitted.h5ad")
```

```python
# Job 2, later: fetch the results.
import anndata as ad
from cytetype import CyteType

adata = ad.read_h5ad("sample1_submitted.h5ad")
annotator = CyteType(adata, group_key="leiden", rank_key="rank_genes_leiden")
results = annotator.get_results()  # None while the job is still running
annotator.cleanup()
if results is not None:
    adata.write_h5ad("sample1_cytetype.h5ad")
```

Creating `CyteType` in the second job prepares the data again, so give that job the same memory as the first. See [Fetch results later](running-jobs.md#running-jobs-fetch-later).
:::
:::{tab-item} R
:sync: r
Set `timeout_seconds` (default 2 hours) to fit within your wall time. When it runs out, `CyteTypeR()` stops with an error of class `cytetype_timeout_error`, and the job keeps running. R returns no object in that case, so keep the job ID, which the error holds as `job_id`; the script above saves it to a file. If the scheduler stops the batch job first, take the job ID from the log (CyteTypeR prints it at submission) or from the dashboard.

With short wall times, submit and save the job ID right away, then fetch the results in a later job:

```r
# Job 1: submit, wait at most 5 minutes, save the job ID.
job_id <- tryCatch(
  {
    obj <- CyteTypeR(obj, prepped_data, study_context = "...", timeout_seconds = 300)
    saveRDS(obj, "sample1_cytetype.rds")  # the job finished in time
    obj@misc$cytetype_jobDetails$job_id
  },
  cytetype_timeout_error = function(e) e$job_id  # still running
)
CleanUpArtifacts(prepped_data)
writeLines(job_id, "sample1_job_id.txt")
```

```r
# Job 2, later: fetch the results.
library(CyteTypeR)

res <- GetResults(job_id = readLines("sample1_job_id.txt"))  # NULL while the job is still running
if (!is.null(res)) write.csv(res, "sample1_cytetype_clusters.csv", row.names = FALSE)
```

`GetResults(job_id = ...)` needs no preparation, but it returns a table with one row per cluster rather than an annotated object; [Fetch results later](running-jobs.md#running-jobs-fetch-later) shows how to add the labels to your object.
:::
::::

## Save the annotated object

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
annotator.cleanup()  # also removes the temporary __cytetype_gene_symbols column
adata.write_h5ad("sample1_cytetype.h5ad")
```

The h5ad file keeps the result columns in `adata.obs` and the entries `cytetype_results` and `cytetype_jobDetails` in `adata.uns` (the results are stored as a JSON string, so they survive the round trip). On the reloaded object, `get_results()` returns them without contacting the server.
:::
:::{tab-item} R
:sync: r
```r
CleanUpArtifacts(prepped_data)
saveRDS(obj, "sample1_cytetype.rds")
```

The RDS file keeps the result columns in `obj@meta.data` and the results and job details in `obj@misc`.
:::
::::

[Results in AnnData and Seurat](../results/results-in-your-data.md) lists everything that is stored.

## Run many datasets

- **Usage limits.** Personal accounts can submit 10 jobs in any 24 hours, with further weekly and monthly limits; running and completed jobs count, and so do jobs you delete, but failed jobs do not. Organization licenses have a cluster quota instead, to which deleted jobs also count. See [Usage limits and quotas](../platform/limits.md), and check your current usage on the **Account** page of the [dashboard](../platform/dashboard.md).
- **Limits are checked at submission,** after the files are uploaded, so a refused job has already spent its upload time. Python raises `RateLimitError` (personal limits) or `QuotaExceededError` (organization quota). Waiting a few minutes does not help, so stop the loop rather than retrying straight away; for a personal limit, the message says when it resets. CyteTypeR stops with the server's message.
- **Separate paths and names.** Give each dataset its own file paths, and a `title` tag so that you can tell the runs apart in the dashboard.
- **Keep the job IDs.** Each object stores its report link, in `adata.uns["cytetype_jobDetails"]["report_url"]` or `obj@misc$cytetype_jobDetails$report_url`.

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
import scanpy as sc
from cytetype import CyteType, QuotaExceededError, RateLimitError

for name, path in datasets.items():  # for example {"sample1": "data/sample1.h5ad", ...}
    adata = sc.read_h5ad(path)
    annotator = CyteType(
        adata,
        group_key="leiden",
        rank_key="rank_genes_leiden",
        vars_h5_path=f"scratch/{name}_vars.h5",
        obs_duckdb_path=f"scratch/{name}_obs.duckdb",
    )
    try:
        adata = annotator.run(
            study_context=study_contexts[name],
            metadata={"title": name},
            query_filename=f"scratch/{name}_query.json",
        )
    except (RateLimitError, QuotaExceededError) as error:
        print(f"Stopped before {name}: {error.message}")
        annotator.cleanup()
        break
    annotator.cleanup()
    adata.write_h5ad(f"results/{name}_cytetype.h5ad")
```
:::
:::{tab-item} R
:sync: r
```r
for (name in names(datasets)) {  # for example c(sample1 = "data/sample1.rds", ...)
  obj <- readRDS(datasets[[name]])
  Idents(obj) <- "seurat_clusters"
  markers <- FindAllMarkers(obj, only.pos = TRUE)
  prepped_data <- PrepareCyteTypeR(
    obj, markers,
    group_key = "seurat_clusters",
    vars_h5_path = file.path("scratch", paste0(name, "_vars.h5")),
    obs_duckdb_path = file.path("scratch", paste0(name, "_obs.duckdb"))
  )
  result <- tryCatch(
    CyteTypeR(
      obj, prepped_data,
      study_context = study_contexts[[name]],
      metadata = list(title = name),
      query_filename = file.path("scratch", paste0(name, "_query.json")),
      show_progress = FALSE
    ),
    error = function(e) {
      message("Stopped at ", name, ": ", conditionMessage(e))
      NULL
    }
  )
  CleanUpArtifacts(prepped_data)
  if (is.null(result)) break
  saveRDS(result, file.path("results", paste0(name, "_cytetype.rds")))
}
```

This loop stops at the first error of any kind, including a usage limit or a timeout.
:::
::::

## Logs and progress output

- **Python** logs to standard output. When the output goes to a file, the progress display prints a new line only when the job's status changes, and the preparation and upload steps print progress bars to standard error. `show_progress=False` turns off the job status display.
- **R** logs through the logger package. Use `show_progress = FALSE` in batch jobs, as in the scripts above.
- **Keep the logs.** They contain the job ID and the report link, which you need to fetch results from a run that stopped early.

## See also

- [Run, monitor and rerun jobs](running-jobs.md): what `run()` and `CyteTypeR()` do, and how to fetch results later.
- [Large datasets](large-datasets.md): memory, disk space and upload size for big jobs.
- [Usage limits and quotas](../platform/limits.md): the limits that apply when you submit many jobs.
