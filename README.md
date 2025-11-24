<h1 align="left">CyteType</h1>
<h4 align="left">Multi-agent AI driven cell annotation in single-cell transcriptomics</h4>

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
  <a href="https://discord.gg/V6QFM4AN">
    <img src="https://img.shields.io/discord/1339594966300622908" alt="PyPI version">
  </a>
</p>

---
<a href="https://colab.research.google.com/drive/1aRLsI3mx8JR8u5BKHs48YUbLsqRsh2N7?usp=sharing" target="_blank">Example Notebook</a> | 
<a href="https://prod.cytetype.nygen.io/report/6420a807-8bf3-4c33-8731-7617edfc2ad0?v=251124" target="_blank">Example output</a> | 
<a href="docs/examples.md">Atlas scale results</a>

Switch to R/Seurat package: <a href="https://github.com/NygenAnalytics/CyteTypeR">CyteTypeR</a>

Read Preprint: <a href="https://www.biorxiv.org/content/10.1101/2025.11.06.686964v1">BioRxiv</a>

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

group_key = 'clusters'  # Wherever you want to store or already have clusters in adata.obs

sc.tl.leiden(adata, key_added=group_key)
sc.tl.umap(adata)
sc.tl.rank_genes_groups(adata, groupby=group_key, method="t-test")
# ------ Example Scanpy Pipeline ------

# ------ CyteType ------
annotator = CyteType(adata, group_key=group_key)
adata = annotator.run(
    study_context="Brief study description (e.g., Human brain tissue ...)",
)

# Visualize results
sc.pl.embedding(adata, basis='umap', color=f'cytetype_annotation_{group_key}')
sc.pl.embedding(adata, basis='umap', color=f'cytetype_cellOntologyTerm_{group_key}')
sc.pl.embedding(adata, basis='umap', color=f'cytetype_cellOntologyTermID_{group_key}')
sc.pl.embedding(adata, basis='umap', color=f'cytetype_cellState_{group_key}')
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

---
<img width="2063" height="1857" alt="CyteType architecture" src="https://github.com/user-attachments/assets/c55f00a2-c4d1-420a-88c2-cdb507898383" />
