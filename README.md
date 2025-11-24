<h1 align="left">CyteType</h1>
<h3 align="left">Automated, evidence-based cell type annotation for single-cell transcriptomics</h3>

<p align="left">
  <a href="https://github.com/NygenAnalytics/cytetype/actions/workflows/publish.yml">
    <img src="https://github.com/NygenAnalytics/cytetype/actions/workflows/publish.yml/badge.svg" alt="CI Status">
  </a>
  <img src="https://img.shields.io/badge/python-≥3.11-blue.svg" alt="Python Version">
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

---

> 🎉 **NEW:** [Preprint published November 7, 2025](https://www.biorxiv.org/content/10.1101/2025.11.06.686964v1) on bioRxiv  
> 📅 **FREE Webinar:** [Register now](https://attendee.gotowebinar.com/register/1731194703386732893) — Learn CyteType from the developers

---

## Why CyteType?

Manual cell type annotation takes weeks and varies between experts. CyteType delivers consistent, expert-level annotations in minutes using a multi-agent AI system where specialized agents collaborate on marker analysis, literature evidence, and Cell Ontology mapping.

<img width="800" alt="CyteType Overview" src="https://github.com/user-attachments/assets/c4cc4f67-9c63-4590-9717-c2391b3e5faf" />

- **Save weeks of manual curation** — Annotate entire datasets at expert level in minutes, not days
- **Drop-in integration** — 3 lines of code, works with existing Scanpy/Seurat workflows
- **No setup friction** — No API keys required; built-in LLM with optional custom configuration
- **Standards-compliant output** — Automatic Cell Ontology term mapping (CL IDs)
- **Comprehensive annotations** — Cell types, subtypes, activation states, confidence scores, and lineage
- **Transparent and auditable** — Interactive HTML reports show evidence, reasoning, and confidence for every annotation

**🚀 [Try it now in Colab](https://colab.research.google.com/drive/1aRLsI3mx8JR8u5BKHs48YUbLsqRsh2N7?usp=sharing)** • **[See example report](https://prod.cytetype.nygen.io/report/6420a807-8bf3-4c33-8731-7617edfc2ad0?v=251124)**

---

## Quick Start

```bash
pip install cytetype
```

```python
import scanpy as sc
from cytetype import CyteType

# Assumes preprocessed AnnData with clusters and marker genes
annotator = CyteType(adata, group_key='clusters')
adata = annotator.run(study_context="Human PBMC from healthy donor")
sc.pl.umap(adata, color='cytetype_annotation_clusters')
```

> **Note:** No API keys required for default configuration. See [custom LLM configuration](docs/configuration.md#llm-configuration) for advanced options.

**Using R/Seurat?** → [CyteTypeR](https://github.com/NygenAnalytics/CyteTypeR)

---

## Output Reports

Each analysis generates an HTML report documenting annotation decisions, marker genes, confidence scores, and Cell Ontology mappings:

<img width="1000" alt="CyteType Report Example" src="https://github.com/user-attachments/assets/9f0f4b36-2dd7-4cb8-93e3-ecda9c97a930" />

**[View example report with embedded chat interface](https://prod.cytetype.nygen.io/report/6420a807-8bf3-4c33-8731-7617edfc2ad0?v=251124)**

---

## Benchmarks

Validated across multiple datasets, tissues, and organisms. CyteType's agentic architecture consistently outperforms other methods across multiple LLMs:

**📊 Performance:** 388% improvement over GPTCellType, 268% over CellTypist, 101% over SingleR

<img width="500" alt="CyteType Benchmark Results" src="https://github.com/user-attachments/assets/a63cadc1-d8c5-4ac0-bba7-af36f9b3c46d" />

**[Browse results from single-cell atlases →](docs/examples.md)**

## Need Help?

📖 [Configuration options](docs/configuration.md) • [Understanding output columns](docs/results.md) • [Common issues & fixes](docs/troubleshooting.md) • [Development guide](docs/development.md)  
💬 [Join Discord](https://discord.gg/V6QFM4AN) for support

---

## Citation

If you use CyteType in your research, please cite our preprint:

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

CyteType is free for academic and non-commercial research use under CC BY‑NC‑SA 4.0 — see [LICENSE.md](LICENSE.md)

For commercial use, please contact us at [contact@nygen.io](mailto:contact@nygen.io)
