---
description: How to cite CyteType, what to record so that others can inspect and reproduce your annotations, and the licenses that apply to the clients and the service.
---

(citing)=
# Cite CyteType and report versions

This page gives the reference to cite, the details to record for each run so that others can inspect and reproduce your annotations, a methods template, and the licenses of the clients and the service.

(citing-preprint)=
## Cite the preprint

If you use CyteType in your work, please cite the preprint:

> Ahuja G, Antill A, Su Y, Dall'Olio GM, Basnayake S, Karlsson G, Dhapola P. Multi-agent AI enables evidence-based cell annotation in single-cell transcriptomics. *bioRxiv* 2025. doi:[10.1101/2025.11.06.686964](https://doi.org/10.1101/2025.11.06.686964)

BibTeX:

```bibtex
@article{cytetype2025,
  title   = {Multi-agent {AI} enables evidence-based cell annotation in single-cell transcriptomics},
  author  = {Ahuja, Gautam and Antill, Alex and Su, Yi and Dall'Olio, Giovanni Marco and Basnayake, Sukhitha and Karlsson, G{\"o}ran and Dhapola, Parashar},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.11.06.686964},
  url     = {https://www.biorxiv.org/content/10.1101/2025.11.06.686964v1}
}
```

The **Citation** link at the bottom of the report's sidebar also opens the preprint.

The preprint describes an earlier version of the CyteType workflow, and the service has changed since. Report the client version and the date of your run as well, as described below.

(citing-record)=
## Record what you ran

The annotation runs on CyteType's servers, and the workflow and the models it uses are updated over time: a client you installed earlier sends its jobs to the current workflow. Runs are also not deterministic, so the same data can give slightly different results. Record these items for every run you report:

| Item | Where to find it |
|---|---|
| Job ID | Python: `adata.uns["cytetype_jobDetails"]["job_id"]`. R: `pbmc@misc$cytetype_jobDetails$job_id`. It is also part of the report link and listed in the [dashboard](../platform/dashboard.md). |
| Report link | Python: `adata.uns["cytetype_jobDetails"]["report_url"]`. R: `pbmc@misc$cytetype_jobDetails$report_url`. |
| Client and version | Python: `cytetype.__version__` or `cytetype --version`. R: `packageVersion("CyteTypeR")`. The report's study information shows them too. |
| Date of the run | The **Created** column of the report list in the dashboard. |
| Request | `query.json`, saved in your working directory unless you pass `save_query=False` (R: `save_query = FALSE`). It holds the study context, metadata tags, cluster labels, marker genes, expression summaries and the client version. |
| Study context | The `study_context` text you passed. |
| Your settings | The cluster column (`group_key`), how you computed the marker genes (method and parameters), `n_top_genes`, `results_prefix` and, if you used your own LLM, the provider and model names. Never record API keys. |
| Results | Python: the `cytetype_*` columns of `adata.obs` and `adata.uns["cytetype_results"]`. R: the `cytetype_*` columns of `pbmc@meta.data` and `pbmc@misc$cytetype_results`. In the report: **Download as CSV**. See [Results in AnnData and Seurat](../results/results-in-your-data.md). |

If you passed another `results_prefix`, replace `cytetype` in these names with your prefix.

To print the main items:

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
import cytetype

details = adata.uns["cytetype_jobDetails"]
print("cytetype", cytetype.__version__)
print("Job ID:", details["job_id"])
print("Report:", details["report_url"])
```
:::
:::{tab-item} R
:sync: r
```r
details <- pbmc@misc$cytetype_jobDetails
packageVersion("CyteTypeR")
details$job_id
details$report_url
```
:::
::::

`query.json` records what you sent, but not the server version or the models that served the job. If you need that level of provenance for a publication, contact contact@nygen.io with the job ID.

:::{warning}
In CyteTypeR, `query.json` also contains your `llm_configs`, including provider API keys. Remove the keys before you share the file, or pass `save_query = FALSE`. The Python client never writes `llm_configs` to `query.json`.
:::

(citing-methods)=
## Describe CyteType in your methods

Adapt this template to your analysis:

::::{tab-set}
:::{tab-item} Python
:sync: python
```text
Clusters were annotated with CyteType (Ahuja et al., 2025) using the cytetype
Python client, version <version>, on <date> (job <job ID>). Marker genes were the
top <n_top_genes> genes per cluster from Scanpy's rank_genes_groups (<method>).
The study context was: "<study context>". The annotation report is available at
<report link>.
```
:::
:::{tab-item} R
:sync: r
```text
Clusters were annotated with CyteType (Ahuja et al., 2025) using the CyteTypeR
package, version <version>, on <date> (job <job ID>). Marker genes were the top
<n_top_genes> genes per cluster by average log2 fold change from Seurat's
FindAllMarkers (<parameters>). The study context was: "<study context>". The
annotation report is available at <report link>.
```
:::
::::

If you used your own LLM, also name the provider and model.

(citing-share)=
## Share the report with readers

- **Make it public.** New reports are private. To let reviewers and readers open the report, set its visibility to **Public** in the dashboard, then cite the link: anyone with the link can view a public report. A public report is also listed in Chat's report search for every signed-in CyteType user, who can then query its uploaded cell-level data, including all of your cell metadata columns. Check that you may publish all of it first; [Visibility and sharing](../platform/sharing.md) lists what a public report exposes.
- **Keep it.** Deleting a report breaks its link for everyone and cannot be undone from the dashboard. Keep the reports you cite.
- **Add a table.** **Download as CSV** in the report's toolbar gives one row per cluster, which works well as a supplementary table. See [Download as CSV](../results/report.md#report-csv) for its columns.

(citing-licenses)=
## Licenses

| Component | License or terms |
|---|---|
| `cytetype`, the Python client | BSD 3-Clause |
| `CyteTypeR`, the R client | BSD 3-Clause |
| The CyteType service (annotation, report, dashboard, Chat and MCP connector) | Free for academic and non-commercial research within the [usage limits](../platform/limits.md). Commercial use requires a license from Nygen Analytics AB. |

The client licenses cover the code of the two packages. The service terms cover the annotations and other data that the service generates: data generated under free academic use is licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/), and a commercial license includes commercial rights to the generated data. For commercial licensing, contact contact@nygen.io.
