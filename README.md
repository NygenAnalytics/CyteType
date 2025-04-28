<p align="center">
  <img src="logo.svg" alt="CyteType Logo" width="200"/>
</p>

<h1 align="left">CyteType</h1>

<p align="left">
  <!-- GitHub Actions CI Badge -->
  <a href="https://github.com/NygenAnalytics/CyteType/actions/workflows/publish.yml">
    <img src="https://github.com/NygenAnalytics/CyteType/actions/workflows/publish.yml/badge.svg" alt="CI Status">
  </a>
  <a href="https://pypi.org/project/cytetype/">
    <img src="https://img.shields.io/pypi/v/cytetype.svg" alt="PyPI version">
  </a>
  <a href="https://github.com/NygenAnalytics/CyteType/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg" alt="License: CC BY-NC-SA 4.0">
  </a>
  <img src="https://img.shields.io/badge/python-≥3.12-blue.svg" alt="Python Version">
</p>

---

**CyteType** is a Python package for automated cell type annotation of single-cell RNA-seq data. It integrates directly with `AnnData` object.

## Key Features

*   Seamless integration with `AnnData` objects.
*   Submits marker genes derived from `scanpy.tl.rank_genes_groups`.
*   Adds annotation results directly back into your `AnnData` object (`adata.obs` and `adata.uns`).

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
from cytetype import annotate_anndata

# --- Preprocessing ---
adata = anndata.read_h5ad("path/to/your/data.h5ad")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# ... other steps like HVG selection, scaling, PCA, neighbors ...

sc.tl.leiden(adata, key_added='leiden_clusters')

sc.tl.rank_genes_groups(adata, groupby='leiden_clusters', method='t-test', key_added='rank_genes_leiden')

# --- Annotation ---
adata = annotate_anndata(
    adata=adata,
    cell_group_key='leiden_clusters',    # Key in adata.obs with cluster labels
    rank_genes_key='rank_genes_leiden',  # Key in adata.uns with rank_genes_groups results
    results_key_added='CyteType',        # Prefix for keys added by CyteType
    n_top_genes=50                      # Number of top marker genes per cluster to submit
)

# Access the cell type annotations that were added to the AnnData object
print(adata.obs.CyteType_leiden_clusters)

# Get detailed information about the annotation results
print (adata.uns['CyteType_leiden_clusters'])

```

## Development

To set up for development:

1.  Clone the repository:
    ```bash
    git clone https://github.com/NygenAnalytics/CyteType.git
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

This project is licensed under the CC BY-NC-SA 4.0 License - see the [LICENSE](LICENSE) file for details.
