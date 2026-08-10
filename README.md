
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

**CyteType** is an end-to-end cell type annotation system for **single-cell RNA sequencing (scRNA-seq)**, designed for repeatable analysis pipelines rather than one-off prompting. It combines cluster-level marker genes, expression context, study metadata, literature retrieval, ontology mapping, and a dedicated review step in a structured workflow that operates directly on AnnData.

For Seurat workflows, use [CyteTypeR](https://github.com/NygenAnalytics/CyteTypeR).

> [!IMPORTANT]
> CyteType requires an API key. Use is free for academic and non-commercial research. Commercial use requires a [license](#license).

## Quick Start

### 1. Install

```bash
pip install cytetype
```

### 2. Set up your API key

```bash
cytetype setup
```

This opens passwordless CyteType sign-in in your browser and saves the API key locally for automatic use from Python. You can also [create or manage API keys in the dashboard](https://cytetype.nygen.io/dashboard).

Already have an API key? Save and validate it locally once:

```bash
cytetype login
```

### 3. Annotate with Scanpy

```python
import scanpy as sc
from cytetype import CyteType

# Assumes preprocessed AnnData with clusters and marker genes
group_key = "clusters"
annotator = CyteType(
    adata,
    group_key=group_key,
    rank_key=f"rank_genes_{group_key}",
    n_top_genes=100,
)
adata = annotator.run(study_context="Human PBMC from a healthy donor")
sc.pl.umap(adata, color="cytetype_annotation_clusters")
```

[Try CyteType in Google Colab](https://colab.research.google.com/drive/1aRLsI3mx8JR8u5BKHs48YUbLsqRsh2N7?usp=sharing).

## What You Get

- **Annotations:** Cell type, subtype, and activation state for every cluster
- **Cell Ontology mapping:** Standardized CL IDs for comparison across studies
- **Confidence and quality control:** Confidence values, plus match scores against your existing labels
- **Supporting evidence:** Publications and condition-specific references behind each call

## Example Report

Each analysis generates an HTML report with annotation decisions, reviewer comments, supporting evidence, and an embedded chat interface connected to your expression data.

<img width="1000" alt="CyteType HTML report showing cell type annotations marker genes" src="https://github.com/user-attachments/assets/e5373fdd-7173-42db-b863-76a1e8ecfe01" />

[View example report](https://cytetype.nygen.io/report/e70e2883-7713-4121-94f2-5b57eabd1468?v=260303)

## Benchmarks

Across PBMC, bone marrow, tumor microenvironment, and cross-species datasets, the multi-agent approach outperforms existing annotation methods:

| Compared with | Improvement |
|---------------|-------------|
| GPTCellType | +388% |
| CellTypist | +268% |
| SingleR | +101% |

Methods and full results are in the [preprint](https://www.biorxiv.org/content/10.1101/2025.11.06.686964v1). You can also [browse results on atlas-scale datasets](docs/examples.md).

## Resources

- 🔐 [CLI and Authentication](docs/cli.md): Set up API keys and manage saved credentials.
- ⚙️ [Configuration](docs/configuration.md): Customize annotation settings, LLM providers, and artifacts.
- 📋 [Output Columns](docs/results.md): Understand annotations and metadata added to AnnData.
- 🛠️ [Troubleshooting](docs/troubleshooting.md): Resolve authentication, API, artifact, and LLM issues.
- 🧑‍💻 [Development](docs/development.md): Configure a local environment and contribute.
- 🎥 [Introduction video](https://vimeo.com/nygen/cytetype): Watch a quick overview of CyteType.
- 💬 [Discord community](https://discord.gg/V6QFM4AN): Ask questions and get support.

## Citation

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

## License

CyteType is free for academic and non-commercial research under [CC BY-NC-SA 4.0](LICENSE.md).

For commercial licensing, contact [contact@nygen.io](mailto:contact@nygen.io).
