---
description: The exported functions of the CyteTypeR R package, with their usage, arguments, return values and examples, generated from the package's help pages.
---

(r-api)=
# R function reference

This page is generated from the help pages of CyteTypeR 0.6.0. To regenerate it, run `Rscript docs/render_r_reference.R ../CyteTypeR` from the root of the CyteType repository. In R, the same help is available with `?` and the function name, for example `?PrepareCyteTypeR`.

| Function | Section | Purpose |
| --- | --- | --- |
| [`PrepareCyteTypeR()`](#r-api-preparecytetyper) | Prepare and run | Prepare Seurat Object for CyteType Analysis |
| [`CyteTypeR()`](#r-api-cytetyper) | Prepare and run | Submit CyteType Analysis Job |
| [`CleanUpArtifacts()`](#r-api-cleanupartifacts) | Prepare and run | Clean Up Artifacts |
| [`GetResults()`](#r-api-getresults) | Results | Retrieve CyteType Analysis Results |
| [`SetupCyteTypeR()`](#r-api-setupcytetyper) | Sign-in and the cytetyper command | Set up CyteType authentication |
| [`LoginCyteTypeR()`](#r-api-logincytetyper) | Sign-in and the cytetyper command | Save an existing CyteType API key |
| [`LogoutCyteTypeR()`](#r-api-logoutcytetyper) | Sign-in and the cytetyper command | Remove local CyteType credentials |
| [`OpenCyteTypeDashboard()`](#r-api-opencytetypedashboard) | Sign-in and the cytetyper command | Open the CyteType dashboard |
| [`ViewCyteTypeJob()`](#r-api-viewcytetypejob) | Sign-in and the cytetyper command | Open a CyteType job |
| [`InstallCyteTypeRCli()`](#r-api-installcytetypercli) | Sign-in and the cytetyper command | Install the cytetyper terminal command |
| [`cytetype_api_error()`](#r-api-cytetype-api-error) | Errors | Custom API Error |

(r-api-prepare-and-run)=
## Prepare and run

Prepare a Seurat object, submit the annotation job and wait for the results, then delete the uploaded files.

(r-api-preparecytetyper)=
### `PrepareCyteTypeR()`

*Prepare Seurat Object for CyteType Analysis*

Prepares a Seurat object and a gene marker table for CyteType API by extracting and formatting cluster labels, metadata, marker genes, visualization coordinates, and expression data into the required data structure. It also writes the two files that [`CyteTypeR()`](#r-api-cytetyper) uploads with the job: `vars.h5` (expression) and `obs.duckdb` (cell metadata and coordinates).

**Usage**

```r
PrepareCyteTypeR(
  obj,
  marker_table,
  group_key = "seurat_clusters",
  gene_symbols = "gene_symbols",
  n_top_genes = 50,
  aggregate_metadata = TRUE,
  min_percentage = 10,
  max_metadata_categories = 500L,
  pcent_batch_size = 5000,
  coordinates_key = "umap",
  max_cells_per_group = 1000,
  vars_h5_path = "vars.h5",
  obs_duckdb_path = "obs.duckdb"
)
```

**Arguments**

| Argument | Description |
| --- | --- |
| `obj` | A Seurat object with normalized data and clustering results. It needs an assay named `"RNA"` with a normalized `data` layer. |
| `marker_table` | Data frame containing marker genes with the columns `cluster`, `gene` and `avg_log2FC`, such as the output of Seurat's `FindAllMarkers()`. Its clusters must be exactly the clusters in `group_key`, each with at least five marker genes. |
| `group_key` | Character string specifying the metadata column containing cluster assignments. Every cell needs a cluster: `NA` values stop the preparation with an error. Default is `"seurat_clusters"`. |
| `gene_symbols` | Character string naming the column of the RNA assay's feature metadata that holds gene symbols. It is used only to check that the genes are symbols rather than database identifiers; when the column does not exist, the row names of the RNA assay are checked instead. Default is `"gene_symbols"`. |
| `n_top_genes` | Integer specifying the maximum number of top marker genes per cluster to include after ordering by `avg_log2FC`. Filter `marker_table` before calling this function when a fold-change threshold is required. Default is 50. |
| `aggregate_metadata` | Logical indicating whether to aggregate metadata across cells within each cluster. Default is `TRUE`. |
| `min_percentage` | Numeric threshold for minimum percentage. Default is 10. |
| `max_metadata_categories` | Integer. Maximum number of unique values a categorical metadata column may have to be included in aggregation. Columns exceeding this limit are silently skipped. Default is 500. |
| `pcent_batch_size` | Integer specifying batch size for expression percentage calculations. Default is 5000. |
| `coordinates_key` | Character string naming the dimensional reduction used for the report's visualization, such as `"umap"` or `"tsne"`. This reduction is tried first, then `"umap"`, `"tsne"` and `"pca"`, matched case-insensitively; the first one with a row for every cell and at least two dimensions is used. When none qualifies, the job runs without visualization. Default is `"umap"`. |
| `max_cells_per_group` | Integer specifying maximum cells per cluster for visualization subsampling. Default is 1000. |
| `vars_h5_path` | Character string specifying the local file path for the generated vars.h5 artifact (feature expression). Default is `"vars.h5"`. |
| `obs_duckdb_path` | Character string specifying the local file path for the generated obs.duckdb artifact (cell metadata). Default is `"obs.duckdb"`. |

**Value**

A named list to pass to [`CyteTypeR()`](#r-api-cytetyper):

- `clusterLabels`: List mapping the CyteType cluster IDs (`"1"`, `"2"`, ...) to your cluster labels
- `clusterMetadata`: Aggregated metadata by cluster (if enabled)
- `markerGenes`: List of top marker genes per cluster
- `visualizationData`: Sampled coordinates and cluster assignments for plotting, or `NULL` when no coordinates were found
- `expressionData`: Expression percentages by cluster
- `group_key`: The metadata column name used for cluster assignments
- `build_succeeded`: Logical indicating whether artifact files were built successfully
- `vars_h5_path`: Path to the generated vars.h5 file
- `obs_duckdb_path`: Path to the generated obs.duckdb file
- `coordinates_key`: The name of the reduction whose coordinates were used for the visualization, or `NULL` when none was used

**Details**

The function performs the following steps:

1. Validates the Seurat object structure and the cluster labels
2. Creates a mapping from original cluster labels to sequential numeric IDs
3. Aggregates metadata across cells within clusters (optional)
4. Filters and selects top marker genes per cluster
5. Extracts dimensional reduction coordinates for visualization
6. Calculates expression percentages across clusters
7. Builds a vars.h5 artifact from the normalized expression matrix, with the raw counts when they are whole numbers
8. Builds an obs.duckdb artifact from the cell metadata

Clusters are renumbered sequentially (1, 2, 3, ...) in the sorted order of your labels (level order for a factor) to ensure consistent formatting regardless of original cluster naming.

Expression always comes from the `"RNA"` assay, whatever the default assay is, and Seurat v5 layers split by sample are joined on a copy (your object is not changed). The gene names come from the same assay. When the default assay is not `"RNA"`, a warning says that `"RNA"` is used.

Errors while building the two files are logged and recorded as `build_succeeded = FALSE`; [`CyteTypeR()`](#r-api-cytetyper) then stops unless `require_artifacts = FALSE`.

**Examples**

```r
# Prepare Seurat object with UMAP coordinates
prepped_data <- PrepareCyteTypeR(
  obj = seurat_obj,
  marker_table = markers_df,
  group_key = "seurat_clusters",
  n_top_genes = 30,
  coordinates_key = "umap"
)
```

**See also**

[`CyteTypeR()`](#r-api-cytetyper) to submit the job, [`CleanUpArtifacts()`](#r-api-cleanupartifacts) to delete the files afterwards.

(r-api-cytetyper)=
### `CyteTypeR()`

*Submit CyteType Analysis Job*

Submits a prepared dataset to CyteType for automated cell type annotation using large language models. Handles job submission, polling for results, and integration of annotations back into the Seurat object.

**Usage**

```r
CyteTypeR(
  obj,
  prepped_data,
  study_context = NULL,
  llm_configs = NULL,
  metadata = NULL,
  n_parallel_clusters = 2L,
  results_prefix = "cytetype",
  poll_interval_seconds = NULL,
  timeout_seconds = NULL,
  api_url = NULL,
  auth_token = NULL,
  save_query = TRUE,
  query_filename = "query.json",
  upload_timeout_seconds = 3600L,
  upload_max_workers = 6L,
  require_artifacts = TRUE,
  show_progress = TRUE,
  override_existing_results = FALSE,
  visibility = NULL
)
```

**Arguments**

| Argument | Description |
| --- | --- |
| `obj` | A Seurat object (will be returned with annotations added). |
| `prepped_data` | Named list containing prepared data from `PrepareCyteTypeR()`. |
| `study_context` | Optional character. Biological context for the experimental setup (e.g. organisms, tissues, diseases, developmental stages, single-cell methods, experimental conditions). Default is `NULL`. |
| `llm_configs` | Optional list of LLM configs. Each element must match the LLMModelConfig schema with required `provider` and `name`; either `apiKey` or all AWS credentials (`awsAccessKeyId`, `awsSecretAccessKey`, `awsDefaultRegion`) must be provided. `provider` is one of `"anthropic"`, `"bedrock"`, `"deepseek"`, `"fireworks"`, `"google"`, `"groq"`, `"huggingface"`, `"mistral"`, `"nebius"`, `"openai"`, `"openrouter"`, `"together"`, `"vertex"` or `"xai"`. Default is `NULL` (API default model). |
| `metadata` | Optional named list. Custom metadata tags for the report header; URL-like values are made clickable in the report. Default is `NULL`. |
| `n_parallel_clusters` | Integer. Number of parallel requests to the model (max 50). High values can trigger rate limits. Default is 2. |
| `results_prefix` | Character. Prefix for the results added to `obj`. The factor columns `<prefix>_annotation_<group_key>`, `<prefix>_cellOntologyTerm_<group_key>`, `<prefix>_cellOntologyTermID_<group_key>`, `<prefix>_cellState_<group_key>` and `<prefix>_reviewerConfidence_<group_key>` are added to `obj@meta.data`, with three copies under older names: `<prefix>_<group_key>` (label), `<prefix>_ontologyTerm_<group_key>` and `<prefix>_ontologyID_<group_key>`. The job details are stored in `obj@misc[[paste0(prefix, "_jobDetails")]]` and the results as described under Value. Default is `"cytetype"`. |
| `poll_interval_seconds` | Integer. How often (seconds) to poll the API for results. Default from options (10). |
| `timeout_seconds` | Integer. Maximum time (seconds) to wait for results. When it is reached, `CyteTypeR()` stops with an error of class `cytetype_timeout_error` that names the job ID; the job keeps running on the server and [`GetResults()`](#r-api-getresults) retrieves its results later. Default from options (7200). |
| `api_url` | Optional character. CyteType API server origin. If `NULL`, uses `CYTETYPE_API_URL`, the configured R option, or `https://cytetype.nygen.io`. |
| `auth_token` | Optional character. Bearer token for API authentication. If `NULL`, uses credentials saved by [`SetupCyteTypeR()`](#r-api-setupcytetyper) or [`LoginCyteTypeR()`](#r-api-logincytetyper). |
| `save_query` | Logical. Whether to save the request payload to a JSON file. Default is `TRUE`. Set this to `FALSE` when `llm_configs` contains provider credentials because the saved payload includes the complete LLM configuration. |
| `query_filename` | Character. Filename for the saved query when `save_query` is `TRUE`. Default is `"query.json"`. |
| `upload_timeout_seconds` | Integer. Socket read timeout (seconds) for each artifact upload. Default is 3600. |
| `upload_max_workers` | Integer. Maximum number of parallel upload workers. Default is 6. |
| `require_artifacts` | Logical. If `TRUE`, an error during artifact build or upload stops the run; if `FALSE`, failures are skipped and annotation continues without artifacts. Default is `TRUE`. |
| `show_progress` | Logical. Whether to show progress (spinner and cluster status). Set `FALSE` to disable. Default is `TRUE`. |
| `override_existing_results` | Logical. If `TRUE`, allow overwriting existing results with the same `results_prefix`. If `FALSE` and results exist, an error is raised. Default is `FALSE`. |
| `visibility` | Optional character. Who can view the report: `"private"`, `"org"` (your organization; needs an organization API key) or `"public"` (anyone with the link). Default is `NULL`, which uses the server default (`"private"`). With `"org"` and an API key that looks like a personal key (it starts with `cyt_p_`), a warning is logged before the files are uploaded and the run continues; the server checks the key when the job is submitted and rejects personal keys. |

**Value**

The Seurat object `obj` with the results added. Assign it, for example `obj <- CyteTypeR(obj, prepped_data)`; the object you pass in is not changed.

- `obj@meta.data` gains the factor columns listed under `results_prefix`. Cells in clusters without a result get `"Unknown"` (`""` for the cell state) in the first five columns and `NA` in the three older ones. The reviewer confidence column holds `"High"`, `"Moderate"`, `"Low"`, or `"Unknown"` when it is missing; its levels follow that order, and only the levels that occur are kept.
- `obj@misc[[paste0(prefix, "_jobDetails")]]` holds `job_id`, `report_url`, `api_url`, `group_key`, `cluster_labels` and, when set, `visibility`.
- With the default prefix, `obj@misc$cytetype_results` is a data.frame with one row per annotated cluster and the columns `clusterId` (your cluster label), `annotation`, `ontologyTerm`, `ontologyTermID`, `granularAnnotation`, `cellState`, `justification`, `supportingMarkers`, `conflictingMarkers`, `missingExpression`, `unexpectedExpression`, `reviewerConfidence`, `isHeterogeneous` and `heterogeneityType`. Its attributes `"job_id"` and `"result"` hold the job ID and the full results in the form returned by [`GetResults()`](#r-api-getresults).
- With another prefix, `obj@misc[[paste0(prefix, "_results")]]` is a list with `job_id` and `result` (the full results), and `obj@misc$cytetype_results` is left as it is.

**Details**

The function performs the following workflow:

1. Constructs the analysis query from prepared data
2. Optionally saves the query to a JSON file
3. Uploads the artifact files and submits the job
4. Polls for job completion with progress updates
5. Retrieves results and integrates annotations into the Seurat object

The evidence columns of the results table (`justification`, `supportingMarkers`, `conflictingMarkers`, `missingExpression` and `unexpectedExpression`) are filled from the structured output of current jobs when the stored fields are empty. `supportingMarkers` and `conflictingMarkers` list genes separated by `"; "`.

While waiting, network errors are retried until `timeout_seconds` is reached. Other errors while checking the job stop the run after five in a row. When the results of a completed job cannot be fetched, the rate limit and server errors (HTTP 500, 502, 503 and 504) are retried in the same way, and other errors, such as a deleted report (HTTP 410), stop the run at once. These errors signal a [`cytetype_api_error()`](#r-api-cytetype-api-error) condition. If the run stops or is interrupted, the job keeps running on the server: retrieve its results with `GetResults(job_id = "<job_id>")`, using the job ID printed when the job was submitted.

**Examples**

```r
# Basic usage with prepared data
annotated_seurat <- CyteTypeR(
  obj = seurat_obj,
  prepped_data = prepped_data,
  study_context = "Single-cell analysis of human brain tissue"
)

# With custom settings and metadata
annotated_seurat <- CyteTypeR(
  obj = seurat_obj,
  prepped_data = prepped_data,
  study_context = "Alzheimer's disease brain samples",
  metadata = list(
    "DOI" = "10.1038/example",
    "GEO_Accession" = "GSE123456"
  ),
  timeout_seconds = 1800,
  poll_interval_seconds = 60,
  visibility = "public"
)
```

**See also**

[`PrepareCyteTypeR()`](#r-api-preparecytetyper) for data preparation, [`GetResults()`](#r-api-getresults) for manual result retrieval

(r-api-cleanupartifacts)=
### `CleanUpArtifacts()`

*Clean Up Artifacts*

Deletes the artifact files (`vars.h5` and `obs.duckdb`) that [`PrepareCyteTypeR()`](#r-api-preparecytetyper) wrote, once [`CyteTypeR()`](#r-api-cytetyper) has uploaded them.

**Usage**

```r
CleanUpArtifacts(prepped_data)
```

**Arguments**

| Argument | Description |
| --- | --- |
| `prepped_data` | Named list containing prepared data from `PrepareCyteTypeR()`. |

**Value**

Invisibly, `TRUE` when the artifact files had been built and are now removed, and `FALSE` when they were never built (nothing is deleted).

**Examples**

```r
prepped_data <- PrepareCyteTypeR(pbmc, pbmc.markers)
pbmc <- CyteTypeR(pbmc, prepped_data, study_context = "Human PBMCs")
CleanUpArtifacts(prepped_data)
```

(r-api-results)=
## Results

Read the results stored in a Seurat object, or fetch the results of a job by its ID.

(r-api-getresults)=
### `GetResults()`

*Retrieve CyteType Analysis Results*

If `obj` and `results_prefix` are given, tries to load from `obj@misc` first. If no local result exists but job details are available, fetches the result from the API. If only `job_id` is given, fetches the result from the API.

**Usage**

```r
GetResults(
  obj = NULL,
  job_id = NULL,
  results_prefix = "cytetype",
  auth_token = NULL,
  api_url = NULL,
  save_json = TRUE
)
```

**Arguments**

| Argument | Description |
| --- | --- |
| `obj` | Seurat object that may contain stored results (optional). |
| `job_id` | Job ID from submission (optional if `obj` + `results_prefix` have stored job details). |
| `results_prefix` | Prefix used when storing results. Default `"cytetype"`. |
| `auth_token` | Optional bearer token for API fetch. If `NULL`, uses saved CyteType credentials. When retrieving a job stored in `obj`, an explicit token also requires a matching `api_url`. |
| `api_url` | Optional CyteType server origin for standalone or legacy job retrieval. When `NULL`, uses `CYTETYPE_API_URL` or `https://cytetype.nygen.io`. |
| `save_json` | Logical. When the results are fetched by `job_id`, whether to save the server's response as `cytetypeR_results_<job_id>.json` in the working directory. Default is `TRUE`. |

**Value**

Depends on where the results come from:

- `GetResults(obj)` after a run with the default prefix: the data.frame in `obj@misc$cytetype_results` (see [`CyteTypeR()`](#r-api-cytetyper)).
- Results stored under another prefix, or fetched for the job stored in `obj`: a list with `annotations` (one list per annotated cluster with `clusterId`, the CyteType cluster ID, `annotation`, `ontologyTerm`, `ontologyTermID`, `granularAnnotation`, `cellState`, `justification`, `supportingMarkers`, `conflictingMarkers`, `missingExpression`, `unexpectedExpression`, `reviewerConfidence`, `isHeterogeneous` and `heterogeneityType`), `raw_annotations` (the server's per-cluster results, unchanged) and the server's `studyContext`, `summary`, `semanticOrder`, `clusterCategories` and `displayLevel`.
- `GetResults(job_id = ...)`: a data.frame with the same columns as `obj@misc$cytetype_results`, where `clusterId` is the CyteType cluster ID. Its attributes `"job_id"` and `"result"` hold the job ID and the list described above.
- `NULL`, with a message, when the job has not finished yet or is not found.

**Details**

Fetching asks the server for the job's status and, when the job has completed, for its results: two requests, without waiting. A failed job stops with an error. When the status or the results cannot be fetched, `GetResults()` signals a [`cytetype_api_error()`](#r-api-cytetype-api-error) condition whose `error_code` holds the server's code when it sent one, such as `"REPORT_DELETED"`. Problems with the API key signal its subclass `cytetype_auth_error`.

`reviewerConfidence` is `"High"`, `"Moderate"` or `"Low"` (`NULL` in the list and `NA` in the table when missing), `isHeterogeneous` is `TRUE` or `FALSE` (or missing), and `heterogeneityType` names the kind of heterogeneity the reviewer found. In the list, `supportingMarkers` and `conflictingMarkers` are lists of gene names.

Fetched results are not added to `obj`; R does not modify the object you pass in.

**Examples**

```r
# Results stored in the object by CyteTypeR()
results <- GetResults(annotated_seurat)

# Results of a job by its ID, without saving the JSON file
results <- GetResults(job_id = "<JOB_ID>", save_json = FALSE)
```

**See also**

[`CyteTypeR()`](#r-api-cytetyper) for job submission

(r-api-sign-in)=
## Sign-in and the cytetyper command

Save an API key, open the dashboard or a job's report, and install the `cytetyper` terminal command. Each `cytetyper` command calls one of these functions (see the [Command-line reference](cli.md)).

(r-api-setupcytetyper)=
### `SetupCyteTypeR()`

*Set up CyteType authentication*

Validates matching saved credentials. When none are available, opens CyteType sign-in in a browser and stores the resulting personal API key.

**Usage**

```r
SetupCyteTypeR(api_url = NULL, force = FALSE)
```

**Arguments**

| Argument | Description |
| --- | --- |
| `api_url` | Optional CyteType server origin. When `NULL`, uses `CYTETYPE_API_URL` or `https://cytetype.nygen.io`. |
| `force` | Whether to skip saved credential validation and authenticate with a new key. |

**Value**

The saved credential metadata, invisibly.

**Details**

The browser flow verifies the user's email with a one-time code and returns to a local callback server. CyteTypeR and the Python client use the same credentials file. For development, pass the development server origin as `api_url`. If a saved key is invalid or inactive, use `force = TRUE`.

**Examples**

```r
SetupCyteTypeR()
SetupCyteTypeR(force = TRUE)
SetupCyteTypeR(api_url = "<URL>")
```

(r-api-logincytetyper)=
### `LoginCyteTypeR()`

*Save an existing CyteType API key*

Validates an API key and stores it for later CyteTypeR requests.

**Usage**

```r
LoginCyteTypeR(api_token = NULL, api_url = NULL)
```

**Arguments**

| Argument | Description |
| --- | --- |
| `api_token` | Optional API key. When omitted, a hidden prompt is shown in an interactive session. Otherwise the key is read from the first line of standard input, so `cytetyper login` and scripts can pipe it in. Typed into a terminal, the key is hidden on Linux and macOS but shown on Windows. |
| `api_url` | Optional CyteType server origin. When `NULL`, uses `CYTETYPE_API_URL` or `https://cytetype.nygen.io`. |

**Value**

The saved credential metadata, invisibly.

**Examples**

```r
LoginCyteTypeR()
LoginCyteTypeR(api_url = "<URL>")
```

(r-api-logoutcytetyper)=
### `LogoutCyteTypeR()`

Remove local CyteType credentials

**Usage**

```r
LogoutCyteTypeR()
```

**Value**

`TRUE` when credentials were removed, otherwise `FALSE`, invisibly.

**Details**

This removes only the local credentials file. Revoke an API key from the dashboard when it must also be invalidated on the server.

**Examples**

```r
LogoutCyteTypeR()
```

(r-api-opencytetypedashboard)=
### `OpenCyteTypeDashboard()`

Open the CyteType dashboard

**Usage**

```r
OpenCyteTypeDashboard()
```

**Value**

The dashboard URL, invisibly.

**Details**

Uses the dashboard associated with saved credentials. When credentials are absent, opens the production dashboard at `https://cytetype.nygen.io/dashboard`.

**Examples**

```r
OpenCyteTypeDashboard()
```

(r-api-viewcytetypejob)=
### `ViewCyteTypeJob()`

Open a CyteType job

**Usage**

```r
ViewCyteTypeJob(job_id, api_url = NULL)
```

**Arguments**

| Argument | Description |
| --- | --- |
| `job_id` | CyteType job identifier. |
| `api_url` | Optional CyteType server origin. When `NULL`, uses the origin from saved credentials or the configured default. |

**Value**

The sign-in URL for the job, invisibly.

**Examples**

```r
ViewCyteTypeJob("<JOB_ID>")
ViewCyteTypeJob("<JOB_ID>", api_url = "<URL>")
```

(r-api-installcytetypercli)=
### `InstallCyteTypeRCli()`

Install the cytetyper terminal command

**Usage**

```r
InstallCyteTypeRCli(destination = NULL, overwrite = FALSE)
```

**Arguments**

| Argument | Description |
| --- | --- |
| `destination` | Directory in which to install the command. Defaults to the directory containing the active `Rscript` executable. The default is recommended because the launcher uses the adjacent `Rscript` executable. |
| `overwrite` | Whether to replace an existing command. |

**Value**

The installed command path, invisibly.

**Details**

After installation, run `cytetyper setup` in a terminal to authenticate against `https://cytetype.nygen.io`. During development, use `cytetyper setup --api-url "<URL>"`.

**Examples**

```r
InstallCyteTypeRCli()
```

(r-api-errors)=
## Errors

The condition that CyteTypeR signals when the status or the results of a job cannot be fetched (see [Errors and error codes](errors.md)).

(r-api-cytetype-api-error)=
### `cytetype_api_error()`

*Custom API Error*

Creates the condition that CyteTypeR signals when the status or the results of a job cannot be fetched: by [`GetResults()`](#r-api-getresults), and by [`CyteTypeR()`](#r-api-cytetyper) when it stops waiting for a job because of such errors. Errors about the API key are signalled as `cytetype_auth_error`, a subclass.

**Usage**

```r
cytetype_api_error(message, call = NULL, error_code = NULL, status_code = NULL)
```

**Arguments**

| Argument | Description |
| --- | --- |
| `message` | Error message. |
| `call` | Where the error came from: `"network"` when the server could not be reached, `"api"` for other errors. Default is `NULL`. |
| `error_code` | The error code the server sent, such as `"REPORT_DELETED"`, or `NULL` when there is none. Default is `NULL`. |
| `status_code` | The HTTP status of the server's response, such as `503`, or `NULL` when there is none, for example when the server could not be reached. Default is `NULL`. |

**Value**

An error condition of class `cytetype_api_error`, with the elements `message`, `call`, `error_code` and `status_code`. Printing it shows `CyteType API Error: <message>`.

**Examples**

```r
error <- cytetype_api_error(
  paste(
    "Error while checking job status: HTTP 410 [REPORT_DELETED]",
    "This report has been deleted by the owner"
  ),
  call = "api",
  error_code = "REPORT_DELETED",
  status_code = 410L
)
print(error)

# Handle the condition by its class and read the server's code.
tryCatch(
  stop(error),
  cytetype_api_error = function(e) e$error_code
)
```

