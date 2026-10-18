---
description: Rank cell-type-resolved target hypotheses for a completed CyteType job from its uploaded counts and metadata and public target knowledge, and read the results and their limits.
---

(target-analysis)=
# Target analysis

Target analysis turns a completed annotation job into cell-type-resolved target hypotheses, each backed by evidence you can audit. It answers two kinds of question:

- **Discovery**: which genes, in which annotated populations, are the strongest candidates across the dataset?
- **Target assessment**: you name the targets, and it profiles each one in every population in scope: expression, disease association, consistency across donors, specificity and public facts.

The output is a ranked, annotated list of hypotheses for further study. It does not establish that any target causes disease, works or is safe; see [Limits](#target-limits).

:::{note}
Target analysis is available on licensed plans. Contact contact@nygen.io for access.
:::

## What CyteType adds

- **Reviewed cell identities.** Population labels are anchored to the Cell Ontology, lineage runs from L1 to L5, and cell state is kept apart from identity. The reviewer's confidence and heterogeneity calls are used to flag weak or mixed populations, because a target claim for a cell type depends on that cell type being called correctly.
- **Cell-level data the job already holds.** The raw counts and the donor and condition metadata uploaded with the job.
- **Provenance.** Every number traces to a computation, and every routing and design choice is recorded with who made it: you, a rule, a model or a default.

## Requirements

- **A licensed plan** that includes target analysis.
- **Permission to start a run.** You own the job, or the job belongs to your organization. Anyone who can open the report can view its results (see [Who can see the results](#target-access)).
- **A completed job.** Every cluster must have finished annotation.
- **Uploaded raw counts and cell metadata.** The job needs its cell metadata file (`obs.duckdb`) and an expression file (`vars.h5`) that contains raw integer counts and gene symbols. The clients upload both by default; see [Raw counts](../guides/prepare-data.md#prepare-data-raw-counts) for where they look for the counts. Without them, the run fails with "Target analysis needs uploaded cell metadata (obs.duckdb) and raw counts (vars.h5 with a raw layer and gene symbols)".
- **For disease claims, a comparison in the metadata.** Disease association needs columns that identify donors (or samples) and the disease and reference conditions, with at least 2 donors in each group; a donor counts for a population only with at least 10 of its cells there. Without such a comparison, or when disease status is confounded with a technical, anatomical or treatment column, or when the groups share no sampling site, the run continues with a descriptive design (see [Claims ladder](#target-claims)).

## Start an analysis

1. Open the [dashboard](dashboard.md) and go to **Reports**.
2. Find the completed report and click **Target analysis**.
3. The dialog shows the status of the job's latest run. If you can start a run, it also shows an **Analysis options (JSON)** box, filled with the options of the latest run or with `{}`.
4. Keep `{}` or edit the options (see [Analysis options](#target-options)), then click **Run analysis**. After a completed run the button reads **Run again**.
5. Click **Refresh** to update the status.

If you cannot start a run, the dialog says why: your plan does not include target analysis, your organization's license has expired, or you are neither the owner nor a member of the job's organization. You can still view existing results.

A job has one current run at a time:

- You cannot start a run while another one for the same job is queued or running.
- A new run replaces the previous one as the job's latest run. The dialog and the MCP connector show results only for the latest run, once it completes, so save the results JSON before you start another run.
- A run that has not reported progress for 20 minutes is treated as stopped, and the dialog lets you start another.
- A run is stopped after 3 hours.

(target-options)=
## Analysis options

The options are a JSON object. `{}` lets CyteType infer everything from the study. Each key you add narrows the analysis; unknown keys are rejected.

| Key | Value | Effect |
|---|---|---|
| `disease` | Text, for example `"ulcerative colitis"` | Names the indication instead of inferring it from the study description and metadata. It is matched to a disease term in the Open Targets knowledge base; if no term matches, the disease stays unresolved. |
| `diseaseArea` | `"oncology"`, `"immunology"`, `"neurology"`, `"metabolic"` or `"general"` | Overrides the area inferred from the disease. The area sets the scoring weights, which evidence counts as genetic support, whether a known cancer driver counts as a safety liability (it does not in oncology), and the area-specific interpretation questions. |
| `targets` | A list of up to 50 objects `{"gene": "...", "direction": "..."}`; `direction` is optional: `"inhibit"`, `"activate"` or `"deplete"` | Switches from discovery to target assessment of these genes. Each gene may appear once. `deplete` means removing the cells that carry the target. |
| `modalities` | A list from `"smallMolecule"`, `"antibody"`, `"antibodyDrugConjugate"`, `"cellTherapy"`, `"oligonucleotide"`, `"proteinDegrader"`, `"geneTherapy"` | Restricts candidates to targets that these modalities can act on, and adds a question about modality fit. |
| `cellTypes` | A list of up to 20 cell types | Restricts the populations in scope. Each entry is matched to the CyteType label, granular label, lineage levels or Cell Ontology term of the populations (a Cell Ontology ID such as `"CL:0000236"` also works); without an exact match, a model chooses among the populations. A cell type that matches nothing stops the run with an error rather than widening the scope. |
| `comparison` | `{"conditionColumn": "...", "diseaseLabels": [...], "referenceLabels": [...], "donorColumn": "..."}`; `donorColumn` is optional | Fixes the condition column, the disease and reference labels and, optionally, the donor column instead of inferring them. The column and labels are cell metadata column names and values; disease and reference labels must not overlap, and a column or label that does not exist stops the run. The other design checks still apply, so the design can still become descriptive. |

Discovery with an explicit disease and cell types:

```json
{
  "disease": "ulcerative colitis",
  "cellTypes": ["fibroblast", "CL:0000236"]
}
```

Target assessment for antibody targets:

```json
{
  "targets": [
    {"gene": "TNF", "direction": "inhibit"},
    {"gene": "IL23A"}
  ],
  "modalities": ["antibody"]
}
```

An explicit comparison:

```json
{
  "disease": "ulcerative colitis",
  "comparison": {
    "conditionColumn": "disease",
    "diseaseLabels": ["ulcerative colitis"],
    "referenceLabels": ["normal"],
    "donorColumn": "donor_id"
  }
}
```

When you leave `disease` out and the study profiles more than one disease in separate patient groups, the run stops and asks you to name one: "The study profiles several diseases: ... Name the one to analyse in the request's disease field; each needs its own run."

## How a run works

The status line shows the current stage:

| Stage | What it does |
|---|---|
| load | Reads the job's annotated populations and uploaded files |
| route | Decides the objective, disease, disease area, modalities and populations in scope |
| design | Plans the disease-versus-reference comparison from the metadata, or explains why there is none |
| measurements | Computes expression, specificity, disease association and changes in population abundance from the raw counts |
| knowledge | Collects public facts from Open Targets for every gene that can become a candidate |
| scoring | Scores and ranks candidates under two lenses and checks how stable the rankings are |
| interpretation | Asks a model short questions about the top candidates |

Only the route, design and interpretation stages use a language model. In the route and design stages, the model's choices are checked against the options offered. Every number comes from computation, and model interpretations never change scores.

## Results

When a run is complete, click **Open results JSON** in the dialog. It opens the run record and its result as JSON in a new tab, where you can save it. The report page does not display target analysis results.

The result holds:

| Key | Content |
|---|---|
| `request` | The options the run used |
| `route` | The objective, disease, disease areas, populations in scope and the routing decisions, each with its source |
| `design` | `comparative` or `descriptive`; the donor and condition columns, disease and reference labels, donors per group, covariates, compared sites, checks and decisions |
| `populations` | Per population: cell counts, the disease association and abundance-shift summaries, the top 5 candidates under each lens, and any conflict between the annotation and where the population is found |
| `rankings` | One ranking per lens (`novelty` and `validation`), with the criterion weights and the top 25 candidates in discovery, or every candidate in target assessment |
| `interpretations` | The model's answers about the interpreted candidates (at most 30) |
| `genes` | Public facts for the ranked genes |
| `knowledge` | The Open Targets release used and what could not be resolved |
| `exclusions`, `unmeasuredTargets`, `limitations` | How many candidates were excluded and why, requested targets that could not be measured, and the limitations computed for this run |
| `versions`, `stageSeconds` | Versions of the result schema, playbooks, thresholds and scoring, and the time each stage took |

### Which genes become candidates

In discovery, a gene becomes a candidate in a population when at least 5% of the population's cells, and at least 5 cells, express it (in a comparative design, of the compared disease and reference cells). Non-protein-coding genes, gene symbols shared by more than one gene in the data, genes that none of the requested modalities can act on, and genes characteristic of another cell type (most likely doublets or contamination) are dropped; `exclusions` counts them. In target assessment, every requested gene is kept in every population in scope, and these conditions become flags instead.

### Read a ranked candidate

Each candidate is a gene in a population (in discovery, each gene appears once, at its best population, with the other populations that score nearly as well listed). It has:

- **score**: the weighted mean of its criteria under the lens.
- **criteria**: the criteria in the table below, each scaled from 0 to 1, with the measured value and the value used in the score. A criterion that cannot be measured for any candidate in the run, such as disease association in a descriptive design, is left out and the other weights are rescaled.
- **coverage**: the share of the applicable weight backed by measurements. Unknown values take a fixed prior, the base rate for a typical gene (zero for most evidence criteria), and a criterion that cannot apply (such as disease association for a population found only in disease) is left out.
- **topShare** and **rankInterval**: how stable the rank is when the criterion weights are redrawn 200 times. `topShare` is how often the candidate stays in the top list; `rankInterval` is its 5th to 95th percentile rank.
- **direction**: the intervention (inhibit, activate or deplete) you requested, or else the one the data imply, if any.
- **evidence**: the fraction of cells expressing, counts per million, the donor fraction, the specificity ratio and the runner-up population, the disease-versus-reference log2 fold change (raw and shrunk) with adjusted p-values, donor consistency and donor counts, the abundance shift, and, when the signal may come from another population, which one.
- **flags**: cautions (see the second table below). A flag does not change the score by itself, but the conditions behind three flags already weaken the disease association: it is halved for `exploratoryDonorCount` and treated as unknown for `ambientLikeShift` and `otherCellTypeGene`.

| Criterion | Meaning |
|---|---|
| Disease association | Donor-level expression change between disease and reference within the population |
| Population shift | Donor-level change in the population's abundance with disease |
| Population specificity | Expression in the population relative to other cell types in the dataset, scaled by the share of donors who express it |
| Tissue specificity | Body-wide tissue specificity from public expression atlases |
| Genetic support | Human genetic evidence linking the gene to the disease |
| Tractability | Evidence that a requested (or any feasible) modality can act on the target |
| Safety | Absence of known on-target safety liabilities |
| Clinical precedent | The most advanced clinical program against the target in any indication |
| Indication activity | Clinical programs against the target in this disease or its subtypes |

The two lenses weight the criteria differently. **novelty** favors strong data signals on targets with little clinical precedent or indication activity, and in most disease areas gives population specificity more weight; **validation** favors strong data signals on targets with clinical precedent. Each ranking in the result lists its weights.

| Flag | Meaning |
|---|---|
| `notExpressed` | Expressed by fewer than 5% of the population's cells, or by fewer than 5 cells |
| `ambiguousSymbol` | The gene symbol is shared by more than one gene in the uploaded data |
| `notProteinCoding` | Not a protein-coding gene |
| `notInKnowledgeBase` | No single matching entry in the knowledge base |
| `modalityIncompatible` | None of the requested modalities can act on it |
| `directionConflict` | The requested direction is the opposite of the one the data imply |
| `exploratoryDonorCount` | The comparison has fewer than 3 donors in a group |
| `lowAnnotationConfidence` | The population's reviewer confidence is Low |
| `heterogeneousPopulation` | The population was flagged as heterogeneous |
| `ambientLikeShift` | The disease shift may come from another population's RNA (ambient RNA or doublets) |
| `otherCellTypeGene` | The gene is characteristic of another cell type and most likely comes from doublets or contamination |
| `annotationContextConflict` | The population's annotation conflicts with where it is found, for example a non-malignant cell type found only in disease |

(target-claims)=
## Claims ladder

What a result can support depends on the data the job holds:

| Data available | Claims allowed | Claims not allowed |
|---|---|---|
| Annotation only, no uploaded files | None: the run fails | |
| Raw counts and cell metadata | Expression level, fraction of cells expressing, specificity in the profiled tissue | Disease association; specificity outside the tissue |
| Donors with disease and reference labels, at least 3 donors per group | Donor-level disease difference within a population, donor consistency, change in abundance | Causality; benefit of an intervention |

Comparisons with 2 donors per group run, but count as exploratory: their disease-association score is halved and they get the `exploratoryDonorCount` flag. Small cohorts detect only large changes, so a gene without a significant change may still differ.

(target-limits)=
## Limits

Target analysis never claims:

1. **Causality.** An expression change can follow from disease rather than cause it.
2. **Efficacy or benefit** of an intervention.
3. **Safety clearance.** RNA from one tissue cannot clear safety.
4. **Protein presence, surface accessibility or activity.** Surface location comes from public annotations, never from the RNA data.
5. **Novelty, freedom to operate, patent position, licensing or market size.** The novelty lens reflects only how little public clinical activity a target has.
6. **Generalization beyond the cohort** without replication.
7. **Absence from missing evidence.** An untested gene is not a null result.

No output approves a target for development. Also keep in mind:

- Public facts describe human genes and are matched by gene symbol. For a non-human study, such as mouse, orthology is assumed rather than checked, a gene whose human ortholog has a different symbol gets no public facts, and the result lists this limitation.
- Abundance changes are relative: one population expanding makes the others shrink.
- The criterion weights are research defaults, not fitted values. Use `topShare` and `rankInterval` to see how much a ranking depends on them.
- Interpretations are model-generated text based on the computed evidence and the model's general knowledge. Treat them as prompts for follow-up, not as findings.

(target-access)=
## Who can see the results

Target analysis results belong to the report's data, so everyone who can open the report can read them: the people you share it with, members of your organization for an Organization report, and anyone with the link for a public report (see [Visibility and sharing](sharing.md)). In the dashboard, the **Target analysis** button of a personal report shared with someone shows them the results without letting them start a run. Members of the job's organization can start runs on its reports if their plan includes target analysis.

AI assistants connected through the [MCP connector](mcp.md) can read the latest run with the `get_target_discovery` tool: its status and options, the disease and comparison design, the populations in scope with their disease-associated gene counts and abundance shifts, the top candidates of each lens (10 by default, at most 25) with their scores, coverage, direction, flags and key evidence, and the run's limitations. While the latest run is queued, running or failed, the tool reports that status instead of results. The connector cannot start a run.
