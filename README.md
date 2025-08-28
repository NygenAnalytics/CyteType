<h1 align="left">CyteType</h1>

<p align="left">
  <a href="https://github.com/NygenAnalytics/cytetype/actions/workflows/publish.yml">
    <img src="https://github.com/NygenAnalytics/cytetype/actions/workflows/publish.yml/badge.svg" alt="CI Status">
  </a>
  <img src="https://img.shields.io/badge/python-≥3.11-blue.svg" alt="Python Version">
  <a href="https://pypi.org/project/cytetype/">
    <img src="https://img.shields.io/pypi/v/cytetype.svg" alt="PyPI version">
  </a>
  <a href="https://github.com/NygenAnalytics/cytetype/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg" alt="License: CC BY-NC-SA 4.0">
  </a>
  <a href="https://pypi.org/project/cytetype/">
    <img src="https://img.shields.io/pypi/dm/cytetype" alt="PyPI downloads">
  </a>
  <a href="https://discord.com/channels/1339594966300622908/1398003605545422998">
    <img src="https://img.shields.io/discord/1339594966300622908" alt="PyPI version">
  </a>
</p>

---

CyteType is a Python client for single‑cell RNA‑seq cluster annnotation using multi-agent workflow.

Try <a href="https://colab.research.google.com/drive/1aRLsI3mx8JR8u5BKHs48YUbLsqRsh2N7?usp=sharing" target="_blank">Colab Notebook</a> or browse this <a href="https://nygen-labs-prod--cytetype-api.modal.run/report/5b4eb3e1-fde7-4609-8be0-2bea015c241d?v=250722" target="_blank">example HTML report</a>

**Atlas scale examples**: <a href="docs/examples.md">docs/examples.md</a>

---
<img width="2063" height="1857" alt="CyteType architecture" src="https://github.com/user-attachments/assets/c55f00a2-c4d1-420a-88c2-cdb507898383" />

---


## Installation
`pip install cytetype`

## Quick Start

```python
import anndata
import scanpy as sc
from cytetype import CyteType

# ------ Example Scanpy Pipeline ------
#  Skip this step if you already have clusters and marker genes in an AnnData object. 
adata = anndata.read_h5ad("path/to/your/data.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1000)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, key_added="clusters")
sc.tl.rank_genes_groups(adata, groupby="clusters", method="t-test")
# ------ Example Scanpy Pipeline ------

# ------ CyteType ------
annotator = CyteType(adata, group_key="clusters")
adata = annotator.run(
    study_context="Brief study description (e.g., Human brain tissue ...)",
)

# View results
print(adata.obs.cytetype_annotation_clusters)
print(adata.obs.cytetype_cellOntologyTerm_clusters)
```

## Documentation
- Configuration (LLMs, auth, advanced): <a href="docs/configuration.md">docs/configuration.md</a>
- Results: <a href="docs/results.md">docs/results.md</a>
- Troubleshooting: <a href="docs/troubleshooting.md">docs/troubleshooting.md</a>
- Development: <a href="docs/development.md">docs/development.md</a>
- Server Overview (high‑level): <a href="docs/server-overview.md">docs/server-overview.md</a>
- Ollama Integration: <a href="docs/ollama.md">docs/ollama.md</a>

## License
Licensed under CC BY‑NC‑SA 4.0 — see <a href="LICENSE.md">LICENSE.md</a>.
