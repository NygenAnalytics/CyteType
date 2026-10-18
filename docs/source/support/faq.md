---
description: Short answers to common questions about CyteType, its data requirements, analysis, privacy, access and pricing, and support.
---

(faq)=
# FAQ

Short answers to common questions, with links to the pages that cover each topic in depth. For problems and error messages, see [Troubleshooting](troubleshooting.md).

(faq-general)=
## General

### What is CyteType?

CyteType annotates the clusters of single-cell RNA-seq data. You prepare your clustered data with the Python client (AnnData and Scanpy) or the R client (Seurat), and CyteType's servers run a structured workflow: they read your study context, weigh competing cell type candidates for each cluster against the marker and expression evidence, link the chosen identity to the Cell Ontology, compare related clusters, and have an independent reviewer assess each call. The results come back into your data object and into an interactive [report](../results/report.md). CyteType is developed by Nygen Analytics AB in Sweden. See [How CyteType works](../getting-started/how-cytetype-works.md).

The [2025 preprint](citing.md) describes an earlier version of the workflow. The hosted service has changed since, so record the client version and run date when you publish (see [Cite CyteType and report versions](citing.md)).

### Which environments are supported?

- **Python:** the `cytetype` package (Python 3.12 or newer) for AnnData objects, usually prepared with Scanpy.
- **R:** the `CyteTypeR` package for Seurat objects.

Both clients share one sign-in. Start with [Installation](../getting-started/installation.md), then the [Python quickstart](../getting-started/quickstart.ipynb) or the [R quickstart](../getting-started/quickstart-r.md).

### Does CyteType run on my computer?

Partly. The client prepares your data on your computer: it extracts the marker genes, summarizes expression and metadata per cluster, and writes two files with your expression data and cell metadata. The annotation itself runs on CyteType's servers, so you need internet access and an API key. The clients cannot annotate offline.

### Do I need an account?

Yes. You sign in with an institutional email address (university, institute or company); addresses from personal email providers are not accepted. Signing in from the client saves an API key that it then uses for you. See [Sign in and API keys](../getting-started/sign-in.md).

(faq-data)=
## Data

### What do I need before I run CyteType?

- Log-normalized expression: `adata.X` in Python, the `RNA` assay after `NormalizeData()` in R.
- A cluster column with a label for every cell.
- Marker genes for those clusters: `sc.tl.rank_genes_groups()` with `groupby` set to the cluster column in Python, `FindAllMarkers()` in R.
- Gene symbols, such as `CD3E`, rather than database IDs.

A 2D embedding such as a UMAP, raw counts and cell metadata are optional but add to the report and the evidence. See [Prepare your data](../guides/prepare-data.md).

### Does CyteType compute marker genes?

No. You provide them, because the right marker method often depends on the study. CyteType takes the top markers of each cluster (`n_top_genes`, 50 by default) as the starting point. It adds the share of cells that express each gene in each cluster, computed for all genes, and uses the uploaded files to look up cell-level expression while it checks a call.

### Which species are supported?

Human and mouse. CyteType has not been validated on other organisms, and the clients do not map genes between species.

### Can I use Ensembl IDs instead of gene symbols?

No, CyteType needs gene symbols. The Python client looks for them in common `adata.var` columns (such as `gene_symbols` or `feature_name`) and in `adata.var_names`, and splits values like `TSPAN6_ENSG00000000003`. It stops with an error if about half or more of the names still look like IDs; then pass `gene_symbols_column` or convert the IDs first. CyteTypeR only warns, so make sure the row names of your Seurat object are symbols.

### Does CyteType work on disease, tumor or perturbed samples?

Yes. CyteType uses the disease, treatment and experimental conditions you describe in the [study context](../guides/study-context.md) instead of assuming healthy tissue. Check unusual or transitional populations against the evidence in the report.

### Does CyteType handle spatial, multimodal or CITE-seq data?

The clients take an expression matrix with clusters: an AnnData object, or the `RNA` assay of a Seurat object. There is no dedicated workflow for spatial, multimodal or protein data.

### Does CyteType detect doublets or low-quality cells?

It is not a quality-control tool. The reviewer flags heterogeneous clusters, including clusters that look like doublet artifacts, but you should remove doublets and low-quality cells during preprocessing.

(faq-analysis)=
## Analysis

### Are annotations made per cell or per cluster?

Per cluster. Every cell of a cluster gets that cluster's labels. To annotate at a finer level, cluster at a higher resolution (or subcluster a population) and run CyteType on those clusters.

### Why do repeated runs give slightly different results?

CyteType uses language models, so runs are not deterministic: wording, reasoning, confidence and sometimes the label can change between runs on the same data. Treat agreement across runs as supporting evidence, and record the job ID of each run you report.

### How should I read the confidence label?

The reviewer rates each annotation `High`, `Moderate` or `Low`. It is a categorical review signal, not a calibrated probability that the label is correct. Heterogeneity is assessed separately. Use both to decide which clusters to review first. See [Interpret and review results](../results/interpreting-results.md).

### How long does a job take?

It depends on the number of clusters, the size of the data and the load on the service. The report link printed when the job starts shows its progress live. Both clients wait up to 2 hours by default (`timeout_seconds`); if they stop waiting, the job keeps running. The server stops a job after 4 hours.

### Does the cell sampling in the report change the annotation?

No. For the report's embedding plot, the client samples up to `max_cells_per_group` cells per cluster (1000 by default). The expression summaries sent with the job are computed from all cells, and the uploaded files contain every cell.

### Can I tell CyteType which tissue or cell types to expect?

Describe the organism, tissue, disease and experiment in `study_context`. CyteType uses it to weigh the candidates, but it does not restrict the answer to a list you provide, so that plausible alternatives are not excluded too early. See [Study context and metadata](../guides/study-context.md).

### Can I use my own language model?

Yes. Pass `llm_configs` with a supported provider (for example `openai`, `anthropic`, `bedrock` or `openrouter`), a model name and your provider key; the model must support tool calling. Every model is checked before the job starts. The usage limits still apply, and [Chat](../platform/chat.md) always uses CyteType's hosted models. See [Use your own LLM](../guides/llm-configuration.md).

### What if I disagree with an annotation?

Open the cluster in the report and check its evidence: the supporting, missing and unexpected genes, the candidates that were considered and the reviewer's assessment. You can ask follow-up questions in [Chat](../platform/chat.md). Record your own conclusion in your analysis. To try again, for example with a more specific study context, run a new job with a different `results_prefix` so that both results are kept.

(faq-privacy)=
## Privacy

### What data leaves my computer?

A request with your study context, report metadata, cluster labels, marker genes, per-cluster expression and metadata summaries and sampled embedding coordinates, plus two files: `vars.h5` with the normalized expression matrix (and raw counts when found) and `obs.duckdb` with all cell metadata and the 2D coordinates. FASTQ or BAM files and your AnnData or Seurat object as a whole are not uploaded. See [Security and privacy](security-and-privacy.md).

### Is my data used to train models?

No. Your scientific data is used to run your jobs and the features you use on them, not to train or benchmark models. CyteType does collect operational telemetry to run and improve the service. See [Security and privacy](security-and-privacy.md).

### Who can see my report?

Only you, by default: new reports are private. You can share a report with specific people by email, make it visible to your organization (for reports submitted with an organization API key), or make it public. A public report can be opened by anyone with the link, and any signed-in CyteType user can also find it through Chat's report search and analyze it in Chat or through the MCP connector. See [Visibility and sharing](../platform/sharing.md), which lists what a public report exposes.

### Can I delete a report?

Yes, in the [dashboard](../platform/dashboard.md). Deletion cannot be undone from the dashboard: everyone, including you, loses access to the report, and its shares are removed. The underlying job data is kept on CyteType's storage, and deleted jobs still count toward usage limits. See [Security and privacy](security-and-privacy.md#security-and-privacy-deletion).

### Can I use CyteType with confidential human data?

CyteType does not detect or remove personal identifiers or protected health information. Remove identifiers and unneeded columns from your cell metadata (all of it is uploaded in `obs.duckdb`), keep the study context and metadata free of identifying text, and complete your institution's review first. There is no mode that uploads only summaries. A data processing agreement is available with a commercial license. See [Security and privacy](security-and-privacy.md).

(faq-access)=
## Access and pricing

### Is CyteType free?

CyteType is free for academic and non-commercial research within the usage limits. Personal accounts can run 10 jobs in any rolling 24 hours, 50 in 7 days and 100 in 30 days; failed jobs do not count. Commercial use requires a license: contact contact@nygen.io. See [Usage limits and quotas](../platform/limits.md).

### How can I run more jobs?

Email contact@nygen.io about a free academic license with higher limits. Organizations can get a license with a cluster quota shared by their members. Using your own language model does not raise the limits.

### Under which licenses are the clients released?

Both the Python package and the R package are released under the BSD 3-Clause License. Use of the hosted service is covered by separate terms. See [Licenses](citing.md#citing-licenses).

(faq-support)=
## Support

### My job is taking a long time. Is it stuck?

Open the report link printed when the job started (or the job in the [dashboard](../platform/dashboard.md)) to see its progress. The job keeps running on the server even if your session disconnects or the client stops waiting. Fetch the results later with `annotator.get_results()` in Python, or `GetResults(job_id = "<JOB_ID>")` in R. A job that runs for more than 4 hours is stopped and marked as failed. See [Run, monitor and rerun jobs](../guides/running-jobs.md).

### Where can I get help?

Start with [Troubleshooting](troubleshooting.md). You can also open an issue on GitHub for [CyteType](https://github.com/NygenAnalytics/CyteType/issues) or [CyteTypeR](https://github.com/NygenAnalytics/CyteTypeR/issues), email contact@nygen.io, or ask on [Discord](https://discord.gg/V6QFM4AN).

### How do I cite CyteType?

Cite the bioRxiv preprint by Ahuja et al. (2025) and report the client version, the run date and the job ID. [Cite CyteType and report versions](citing.md) has the reference, BibTeX and a checklist.
