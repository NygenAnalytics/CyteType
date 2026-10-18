---
description: Write a study context that tells CyteType what your sample is, label runs with metadata tags, and choose cell metadata columns that add useful evidence.
---

(study-context)=
# Study context and metadata

The study context is a short description of your experiment that you pass to `run()` or `CyteTypeR()` as `study_context`. It is your main way to tell CyteType what the sample is. This page shows what to write and how CyteType uses it, and how metadata tags and your cell metadata add to it.

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
adata = annotator.run(
    study_context=(
        "Peripheral blood mononuclear cells from three healthy adult human donors, "
        "profiled with 10x Genomics 3' v3 single-cell RNA-seq."
    ),
    metadata={
        "title": "PBMC healthy donors",
        "GEO": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<accession>",
    },
)
```
:::
:::{tab-item} R
:sync: r
```r
obj <- CyteTypeR(
  obj, prepped_data,
  study_context = paste(
    "Peripheral blood mononuclear cells from three healthy adult human donors,",
    "profiled with 10x Genomics 3' v3 single-cell RNA-seq."
  ),
  metadata = list(
    title = "PBMC healthy donors",
    GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<accession>"
  )
)
```
:::
::::

## How CyteType uses the study context

At the start of every job, CyteType builds the background of your study:

1. **Organism.** It identifies the organism from your study context and the marker genes.
2. **Tissue.** It infers the likely tissue from the expression data.
3. **Study summary.** It writes a short summary of the biological context, with keywords, from your study context, the organism and the inferred tissue.
4. **Pathway libraries.** It chooses the pathway libraries that fit the study. They are used to describe the pathway activity of each cluster.

The summary then goes with every cluster through annotation and review, including the comparison of CyteType's labels with yours. In the report, the Context Fit section of each cluster shows how its label fits the study (see [The CyteType report](../results/report.md)).

The study context guides the reasoning but does not restrict it. CyteType is not limited to the tissue or the cell types you mention: when the evidence supports it, it can call cell types that you did not mention or that are unexpected for the stated tissue.

In Python, `study_context` is a required argument. In R it is optional, but without it CyteType has only the genes to go on, so always write one.

## What to include

Describe the sample in a sentence or two, or a short paragraph, in plain language. Cover the items that apply to your study:

| Item | Examples |
|---|---|
| Organism | human; mouse (C57BL/6J) |
| Tissue or organ, with the region when it matters | lung, distal parenchyma; hippocampus; bone marrow aspirate |
| Disease or condition | healthy donors; idiopathic pulmonary fibrosis; tumor with matched adjacent normal tissue |
| Assay | 10x Genomics 3' single-cell RNA-seq; single-nucleus RNA-seq; Smart-seq2 |
| Enrichment or sorting | CD45+ cells sorted by FACS; dead-cell removal; epithelial enrichment |
| Notable perturbations | LPS stimulation for 6 hours; anti-PD-1 treatment; knockout of a named gene; organoids after 30 days of differentiation |
| Developmental stage or age | embryonic day 14.5; fetal, 18 weeks post-conception; aged mice, 24 months |

Some guidance:

- **Be specific about the tissue.** A region or compartment narrows the plausible cell types more than an organ name.
- **Mention sorting and enrichment.** They tell CyteType which populations to expect and which were depleted.
- **Say whether you profiled nuclei or whole cells.** Single-nucleus and single-cell data look different.
- **Describe the sample, not the expected result.** Leave out the cell types you expect to find: whatever you write ends up in the study summary that every cluster is read against. Mention cell types only when they describe how the sample was made, such as sorted CD4+ T cells.
- **Leave out personal or identifying details,** such as patient names or record numbers. The study context is sent with the job as you write it.

## Examples

**Blood**

- Weak: `"PBMC"`
- Better: `"Peripheral blood mononuclear cells from three healthy adult human donors, profiled with 10x Genomics 3' v3 single-cell RNA-seq."`

The better version adds the organism, the condition and the assay.

**Brain**

- Weak: `"Mouse brain, disease model"`
- Better: `"Single-nucleus RNA-seq of hippocampus from 6-month-old male 5xFAD mice, an amyloid model of Alzheimer's disease, and wild-type littermates."`

The better version names the region, the model and its controls, the age, and that nuclei were profiled.

**Tumor**

- Weak: `"Tumor microenvironment"`
- Better: `"Human colorectal adenocarcinoma: treatment-naive primary tumors and matched adjacent normal mucosa from eight patients. CD45+ immune cells were enriched by FACS before 10x Genomics 5' single-cell RNA-seq."`

The better version adds the organism, the disease, the sample types, the immune-cell enrichment and the assay. Without the enrichment, the near absence of epithelial and stromal cells would be hard to interpret.

(study-context-metadata-tags)=
## Metadata tags

`metadata` labels a run: a dictionary in Python, a named list in R. Useful tags include a title, the dataset accession, the publication, the sample or batch, and an analysis version.

- **Values must be strings.** Write numbers and dates as text, for example `"2025"`. Other values stop the run before anything is uploaded: Python raises a `pydantic.ValidationError`, and CyteTypeR stops with `All infoTags values must be character strings`.
- **`title` names the job.** A `title` tag (in any capitalization) becomes the job's name in the [dashboard](../platform/dashboard.md). Without one, the job is listed as "Untitled"; you can rename it in the dashboard.
- **Links.** Any `http://` or `https://` address in a value becomes a clickable link in the report. Write a DOI as a full address, such as `https://doi.org/10.1000/xyz123`.
- **Where tags appear.** The report's Study Information panel lists them under "User-Provided Metadata", next to the client and its version. [Chat](../platform/chat.md) and the [MCP connector](../platform/mcp.md) also see them as part of the run's details.
- **Tags are not evidence.** CyteType does not use them to annotate clusters. Put anything that should inform the labels into the study context.

(study-context-cell-metadata)=
## Cell metadata columns

Besides the study context, CyteType reads the categorical columns of your cell metadata (`adata.obs` or `obj@meta.data`). For each cluster, the client records which values make up more than `min_percentage` of its cells (see [Cell metadata](prepare-data.md#prepare-data-metadata)). CyteType then uses them in two ways:

- **A description of the cluster's sample composition,** for example the share of cells from each condition and donor. CyteType writes it to describe the sample rather than cell identity, and it becomes part of the cluster's context.
- **A comparison with the rest of the study.** CyteType notes when a cluster is clearly enriched in or depleted of a condition, tissue region or treatment compared with the other cells, and when most of its cells come from one donor. With the uploaded `obs.duckdb`, this comparison uses exact proportions; without it, it works from the per-cluster summaries.

The most useful columns describe where the cells came from:

- **Condition or disease status**, such as `disease` with the values `healthy` and `IPF`.
- **Tissue, region or sample site**, such as `site` with `tumor` and `normal`.
- **Treatment and time point.**
- **Donor, patient or sample**, so that CyteType can tell whether a cluster comes mostly from one donor.

Some guidance:

- **Use readable names and values,** such as `condition` with `healthy` and `COVID-19` rather than `grp` with `A` and `B`. CyteType reads them as text.
- **Keep the number of values small.** The comparison works best for columns with a handful of values. Columns with more than `max_metadata_categories` values, such as barcodes, are left out of the summaries.
- **Per-cell labels are not evidence.** Columns that hold an earlier cell-type annotation or cluster assignments are left out of the comparison, and the composition description is written to avoid cell identities. Such columns are still part of the per-cluster summaries, so remove them before you prepare the job if they must not influence the labels. To compare CyteType's labels with your own, use your labels as the cluster column (see [Clusters](prepare-data.md#prepare-data-clusters)).
- **Everything is uploaded.** `obs.duckdb` contains all your metadata columns, whether or not they are summarized. Remove columns that must not leave your computer before you prepare the job.

## Next steps

- [Run, monitor and rerun jobs](running-jobs.md): submit the job with your study context and metadata tags.
- [The CyteType report](../results/report.md): where the study summary, your metadata tags and each cluster's Context Fit appear.
