<h1 align="left">CyteType</h1>
<h3 align="left">Agentic, Evidence-Based Cell Type Annotation for Single-Cell RNA-seq</h3>

<p align="left">
  <a href="https://github.com/NygenAnalytics/cytetype/actions/workflows/publish.yml">
    <img src="https://github.com/NygenAnalytics/cytetype/actions/workflows/publish.yml/badge.svg" alt="CI Status">
  </a>
  <img src="https://img.shields.io/badge/python-≥3.12-blue.svg" alt="Python Version">
  <a href="https://pypi.org/project/cytetype/">
    <img src="https://img.shields.io/pypi/v/cytetype.svg" alt="PyPI version">
  </a>
  <a href="https://raw.githubusercontent.com/NygenAnalytics/CyteType/refs/heads/master/LICENSE.md">
    <img src="https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg" alt="License: CC BY-NC-SA 4.0">
  </a>
  <a href="https://pypi.org/project/cytetype/">
    <img src="https://img.shields.io/pypi/dm/cytetype" alt="PyPI downloads">
  </a>
</p>

**CyteType** performs **automated cell type annotation** in **single-cell RNA sequencing (scRNA-seq)** data. It uses a multi-agent AI architecture to deliver transparent, evidence-based annotations with Cell Ontology mapping.

Integrates with **Scanpy** and **Seurat** workflows.

---

> **Preprint published:** Nov. 7, 2025: [bioRxiv link](https://www.biorxiv.org/content/10.1101/2025.11.06.686964v1) - Dive into benchmarking results

---

## Why CyteType?

Cell type annotation is one of the most time-consuming steps in single-cell analysis. It typically requires weeks of expert curation, and the results often vary between annotators. When annotations do get done, the reasoning is rarely documented; this makes it difficult to reproduce or audit later.

CyteType addresses this with a novel agentic architecture: specialized AI agents collaborate on marker gene analysis, literature evidence retrieval, and ontology mapping. The result is consistent, reproducible annotations with a full evidence trail for every decision.

<img width="800" alt="CyteType multi-agent AI architecture for single-cell RNA-seq cell type annotation" src="https://github.com/user-attachments/assets/c4cc4f67-9c63-4590-9717-c2391b3e5faf" />

---

## Key Features

| Feature | Description |
|---------|-------------|
| **Cell Ontology Integration** | Automatic CL ID assignment for standardized terminology and cross-study comparison |
| **Confidence Scores** | Numeric certainty values (0–1) for cell type, subtype, and activation state — useful for flagging ambiguous clusters |
| **Linked Literature** | Each annotation includes supporting publications and condition-specific references — see exactly why a call was made |
| **Annotation QC via Match Scores** | Compare CyteType results against your existing annotations to quickly identify discrepancies and validate previous work |
| **Embedded Chat Interface** | Explore results interactively; chat is connected to your expression data for on-the-fly queries |

Also included: interactive HTML reports, Scanpy/Seurat compatibility (R wrapper via [CyteTypeR](https://github.com/NygenAnalytics/CyteTypeR)), and no API keys required out of the box.

📹 [Watch CyteType intro video](https://vimeo.com/nygen/cytetype)

---

## Quick Start

### Installation

```bash
pip install cytetype
```

### Basic Usage with Scanpy

```python
import scanpy as sc
from cytetype import CyteType

# Assumes preprocessed AnnData with clusters and marker genes
group_key = 'clusters'
annotator = CyteType(
    adata, 
    group_key=group_key, 
    rank_key='rank_genes_' + group_key, 
    n_top_genes=100
)
adata = annotator.run(study_context="Human PBMC from healthy donor")
sc.pl.umap(adata, color='cytetype_annotation_clusters')
```
🚀 [Try it in Google Colab](https://colab.research.google.com/drive/1aRLsI3mx8JR8u5BKHs48YUbLsqRsh2N7?usp=sharing)

> **Note:** No API keys required for default configuration. See [custom LLM configuration](docs/configuration.md#llm-configuration) for advanced options.
>
> `run()` now handles artifact packaging and upload automatically (`vars.h5` + `obs.duckdb`) before annotation.
> Generated artifact files are kept on disk by default; use `cleanup_artifacts=True` to remove them after run completion/failure.

**Using R/Seurat?** → [CyteTypeR](https://github.com/NygenAnalytics/CyteTypeR)

---

## Documentation

| Resource | Description |
|----------|-------------|
| [Configuration](docs/configuration.md) | LLM settings, parameters, and customization |
| [Output Columns](docs/results.md) | Understanding annotation results and metadata |
| [Troubleshooting](docs/troubleshooting.md) | Common issues and solutions |
| [Development](docs/development.md) | Contributing and local setup |
| [Discord](https://discord.gg/V6QFM4AN) | Community support |

---

## Output Reports

Each analysis generates an HTML report documenting annotation decisions, reviewer comments and an embedded chat interface for further exploration.

<img width="1000" alt="CyteType HTML report showing cell type annotations marker genes and confidence scores" src="https://github.com/user-attachments/assets/9f0f4b36-2dd7-4cb8-93e3-ecda9c97a930" />

[View example report](https://prod.cytetype.nygen.io/report/6420a807-8bf3-4c33-8731-7617edfc2ad0?v=251124)

---

## Benchmarks

Validated across PBMC, bone marrow, tumor microenvironment, and cross-species datasets. CyteType's agentic architecture consistently outperforms existing annotation methods:

| Comparison | Improvement |
|------------|-------------|
| vs GPTCellType | +388% |
| vs CellTypist | +268% |
| vs SingleR | +101% |

<img width="500" alt="CyteType benchmark comparison against GPTCellType CellTypist SingleR" src="https://github.com/user-attachments/assets/a63cadc1-d8c5-4ac0-bba7-af36f9b3c46d" />

[Browse CyteType results on atlas scale datasets](docs/examples.md)

---

## Citation

If you use CyteType in your research, please cite our preprint:

> Ahuja G, Antill A, Su Y, Dall'Olio GM, Basnayake S, Karlsson G, Dhapola P. Multi-agent AI enables evidence-based cell annotation in single-cell transcriptomics. *bioRxiv* 2025. doi: [10.1101/2025.11.06.686964](https://www.biorxiv.org/content/10.1101/2025.11.06.686964v1)

```bibtex
@article{cytetype2025,
  title={Multi-agent AI enables evidence-based cell annotation in single-cell transcriptomics},
  author={Gautam Ahuja, Alex Antill, Yi Su, Giovanni Marco Dall'Olio, Sukhitha Basnayake, Göran Karlsson, Parashar Dhapola},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/2025.11.06.686964},
  url={https://www.biorxiv.org/content/10.1101/2025.11.06.686964v1}
}
```

---

## License

CyteType is free for academic and non-commercial research under [CC BY-NC-SA 4.0](LICENSE.md).

For commercial licensing, contact [contact@nygen.io](mailto:contact@nygen.io).

---
