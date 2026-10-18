
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
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-BSD%203--Clause-blue.svg" alt="License: BSD 3-Clause">
  </a>
  <a href="https://pypi.org/project/cytetype/">
    <img src="https://img.shields.io/pypi/dm/cytetype" alt="PyPI downloads">
  </a>
</p>

**CyteType** is an end-to-end cell type annotation system for **single-cell RNA sequencing (scRNA-seq)**, designed for repeatable analysis pipelines rather than one-off prompting. It combines cluster-level marker genes, expression context, study metadata, ontology mapping, and a dedicated review step in a structured workflow that operates directly on AnnData.

For Seurat workflows, use [CyteTypeR](https://github.com/NygenAnalytics/CyteTypeR).

📚 **Documentation:** [docs.nygen.io/CyteType](https://docs.nygen.io/CyteType/) covers the Python and R clients, the report, the dashboard, Chat and the MCP connector.

> [!IMPORTANT]
> Already installed CyteType? Open a terminal and run `pip install -U cytetype` to update to the latest version before you start.
>
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

# Assumes log-normalized AnnData with clusters in adata.obs["clusters"]
group_key = "clusters"
sc.tl.rank_genes_groups(adata, groupby=group_key, key_added=f"rank_genes_{group_key}")
annotator = CyteType(
    adata,
    group_key=group_key,
    rank_key=f"rank_genes_{group_key}",
    n_top_genes=100,
)
adata = annotator.run(study_context="Human PBMC from a healthy donor")
sc.pl.umap(adata, color="cytetype_annotation_clusters")
```

Follow the [full quickstart](https://docs.nygen.io/CyteType/getting-started/quickstart/), or [try CyteType in Google Colab](https://colab.research.google.com/drive/1aRLsI3mx8JR8u5BKHs48YUbLsqRsh2N7?usp=sharing).

## What You Get

- **Annotations:** Cell type, subtype, and activation state for every cluster
- **Cell Ontology mapping:** Standardized CL IDs for comparison across studies
- **Confidence and quality control:** Reviewer confidence (High, Moderate, Low) and heterogeneity flags, plus match scores against your existing labels
- **Supporting evidence:** Supporting, missing and unexpected marker genes, pathways, and the reasoning behind each call

## Example Report

Each analysis generates an interactive report with annotation decisions, reviewer assessments and supporting evidence. From the report, open [Chat](https://docs.nygen.io/CyteType/platform/chat/) to ask follow-up questions about your data.

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

- 📚 [Documentation](https://docs.nygen.io/CyteType/): Guides for the Python and R clients, the report, the dashboard, Chat and the MCP connector.
- 🔐 [Sign in and API keys](https://docs.nygen.io/CyteType/getting-started/sign-in/): Set up API keys and manage saved credentials.
- 🧬 [Prepare your data](https://docs.nygen.io/CyteType/guides/prepare-data/): Clusters, marker genes, gene symbols and other inputs.
- ⚙️ [Run, monitor and rerun jobs](https://docs.nygen.io/CyteType/guides/running-jobs/): Annotation settings, timeouts, reruns and upload files.
- 🤖 [Use your own LLM](https://docs.nygen.io/CyteType/guides/llm-configuration/): LLM providers, including local models.
- 📋 [Results in AnnData](https://docs.nygen.io/CyteType/results/results-in-your-data/): Columns and keys added to AnnData.
- 🛠️ [Troubleshooting](https://docs.nygen.io/CyteType/support/troubleshooting/): Resolve sign-in, data, upload and job issues.
- 🐍 [Python API reference](https://docs.nygen.io/CyteType/api/index.html): Every public class and function.
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

The CyteType Python client is open source software released under the [BSD 3-Clause License](LICENSE) and maintained by [Nygen](https://nygen.io).

### CyteType service terms

Annotation runs on the CyteType service, which requires an API key.

- **Academic and non-commercial use** is free within the [usage limits](https://docs.nygen.io/CyteType/platform/limits/). Data generated under academic use is licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).
- **Commercial use and higher limits** require a commercial license from Nygen Analytics AB, which includes commercial rights to the generated data.

For licensing, contact [contact@nygen.io](mailto:contact@nygen.io).
