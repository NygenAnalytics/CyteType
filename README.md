<p align="center">
  <img src="logo.svg" alt="CyteType Logo" width="200"/>
</p>

<h1 align="center">CyteType Python Client</h1>

<p align="center">
  <!-- GitHub Actions CI Badge -->
  <a href="https://github.com/paras5/cytetype/actions/workflows/ci.yml">
    <img src="https://github.com/paras5/cytetype/actions/workflows/ci.yml/badge.svg" alt="CI Status">
  </a>
  <a href="https://pypi.org/project/cytetype/">
    <img src="https://img.shields.io/pypi/v/cytetype.svg" alt="PyPI version">
  </a>
  <a href="https://github.com/paras5/cytetype/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg" alt="License: CC BY-NC-SA 4.0">
  </a>
</p>

---

**CyteType** is a Python package providing a convenient interface to the CyteType API for automated cell type annotation of single-cell RNA-seq data. It works directly with `AnnData` objects commonly used in bioinformatics workflows.

## Key Features

*   Seamless integration with `AnnData` objects.
*   Submits marker genes derived from `scanpy.tl.rank_genes_groups`.
*   Handles API communication, job submission, and results polling.
*   Adds annotation results directly back into your `AnnData` object (`adata.obs` and `adata.uns`).
*   Configurable API endpoint, polling interval, and timeout.

## Installation

You can install CyteType using `pip` or `uv`:

```bash
pip install cytetype
```

or

```bash
uv pip install cytetype
```

## Basic Usage

Here's a minimal example demonstrating how to use CyteType after running standard Scanpy preprocessing and marker gene identification:

```python
import anndata
import scanpy as sc
from cytetype.main import annotate_anndata
from cytetype.config import logger

# --- Preprocessing ---
# 1. Load your data
adata = anndata.read_h5ad("path/to/your/data.h5ad")

# 2. Perform standard preprocessing (filtering, normalization, etc.)
#    Example steps (adapt to your data):
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# ... other steps like HVG selection, scaling, PCA, neighbors ...

# 3. Perform clustering (e.g., Leiden)
sc.tl.leiden(adata, key_added='leiden_clusters')

# 4. Find marker genes using rank_genes_groups
#    Ensure you use the same key as your clustering ('leiden_clusters' here)
sc.tl.rank_genes_groups(adata, groupby='leiden_clusters', method='t-test', key_added='rank_genes_leiden')

# --- Annotation ---
adata = annotate_anndata(
    adata=adata,
    cell_group_key='leiden_clusters',    # Key in adata.obs with cluster labels
    rank_genes_key='rank_genes_leiden',  # Key in adata.uns with rank_genes_groups results
    results_key_added='CyteType',        # Prefix for keys added by CyteType
    n_top_genes=50                      # Number of top marker genes per cluster to submit
)

```

## Development

To set up for development:

1.  Clone the repository:
    ```bash
    git clone https://github.com/paras5/cytetype.git
    cd cytetype
    ```
2.  Create and activate a virtual environment (recommended):
    ```bash
    python -m venv .venv
    source .venv/bin/activate # or .venv\Scripts\activate on Windows
    ```
3.  Install dependencies using `uv` (includes development tools):
    ```bash
    pip install uv # Install uv if you don't have it
    uv pip sync --all-extras
    ```
4.  Install the package in editable mode:
    ```bash
    pip install -e .
    ```

### Running Checks and Tests

*   **Mypy (Type Checking):** `uv run mypy .`
*   **Ruff (Linting & Formatting):** `uv run ruff check .` and `uv run ruff format .`
*   **Pytest (Unit Tests):** `uv run pytest`


## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
