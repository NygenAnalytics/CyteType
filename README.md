<h1 align="left">CyteType</h1>

<p align="left">
  <!-- GitHub Actions CI Badge -->
  <a href="https://github.com/NygenAnalytics/cytetype/actions/workflows/publish.yml">
    <img src="https://github.com/NygenAnalytics/cytetype/actions/workflows/publish.yml/badge.svg" alt="CI Status">
  </a>
  <a href="https://github.com/NygenAnalytics/cytetype/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg" alt="License: CC BY-NC-SA 4.0">
  </a>
  <img src="https://img.shields.io/badge/python-≥3.12-blue.svg" alt="Python Version">
</p>

---

**CyteType** is a Python package for automated cell type annotation of single-cell RNA-seq clusters using large language models.

## Quick Start

```python
import anndata
import scanpy as sc
import cytetype

# Load and preprocess your data
adata = anndata.read_h5ad("path/to/your/data.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.leiden(adata, key_added='leiden')
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')

# Initialize CyteType (performs data preparation)
annotator = cytetype.CyteType(adata, group_key='leiden')

# Run annotation
adata = annotator.run(
    study_context="Human brain tissue from Alzheimer's disease patients"
)

# View results
print(adata.obs.cytetype_leiden)
```

## Installation

```bash
pip install cytetype
# or
uv add cytetype
```

## How It Works

CyteType uses a two-step process:

1. **Data Preparation** (during `__init__`): Validates data, calculates expression percentages, and extracts marker genes
2. **API Annotation** (during `run()`): Sends data to CyteType API for LLM-powered cell type annotation

This design allows efficient reuse when making multiple annotation requests with different parameters.

## Basic Usage

### Required Preprocessing

Your `AnnData` object must have:
- Log-normalized expression data in `adata.X`
- Cluster labels in `adata.obs` 
- Differential expression results from `sc.tl.rank_genes_groups`

```python
import scanpy as sc

# Standard preprocessing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Clustering
sc.tl.leiden(adata, key_added='clusters')

# Differential expression (required)
sc.tl.rank_genes_groups(adata, groupby='clusters', method='t-test')
```

### Annotation

```python
from cytetype import CyteType

# Initialize (data preparation happens here)
annotator = CyteType(adata, group_key='clusters')

# Run annotation
adata = annotator.run(
    study_context="Mouse brain cortex, postnatal development"
)

# Results are stored in:
# - adata.obs.cytetype_clusters (cell type annotations)
# - adata.uns['cytetype_results'] (full API response)
```

## Configuration Options

### Initialization Parameters

```python
annotator = CyteType(
    adata,
    group_key='leiden',                    # Required: cluster column name
    rank_key='rank_genes_groups',          # DE results key (default)
    gene_symbols_column='gene_symbols',    # Gene symbols column (default)
    n_top_genes=50,                        # Top marker genes per cluster
    results_prefix='cytetype'              # Prefix for result columns
)
```

### Custom LLM Configuration

```python
# Use your own API key
adata = annotator.run(
    study_context="Human PBMC from COVID-19 patients",
    model_config=[{
        'provider': 'openai',
        'name': 'gpt-4o',
        'apiKey': 'your-api-key'
    }]
)
```

Supported providers: `openai`, `anthropic`, `google`, `xai`, `groq`

## Example Report

View a sample annotation report: [CyteType Report](https://cytetype.nygen.io/report/97ba2a69-ccfa-4b57-8614-746ce2024333)

The report includes:
- Detailed cell type annotations with confidence scores
- Marker gene analysis and supporting evidence
- Alternative annotations and biological justifications

## Advanced Usage

### Multiple Annotations

```python
# Initialize once, run multiple times
annotator = CyteType(adata, group_key='leiden')

# Different contexts
adata1 = annotator.run(study_context="Healthy tissue")
adata2 = annotator.run(study_context="Disease tissue")
```

### Custom Gene Symbols

```python
# If gene symbols are in a different column
annotator = CyteType(
    adata, 
    group_key='leiden',
    gene_symbols_column='gene_name'  # instead of default 'gene_symbols'
)
```

## Development

### Setup

```bash
git clone https://github.com/NygenAnalytics/cytetype.git
cd cytetype
uv sync --all-extras
uv run pip install -e .
```

### Testing

```bash
uv run pytest              # Run tests
uv run ruff check .        # Linting
uv run ruff format .       # Formatting
uv run mypy .              # Type checking
```

## License

Licensed under CC BY-NC-SA 4.0 - see [LICENSE](LICENSE) for details.
