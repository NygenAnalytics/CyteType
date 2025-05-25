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

**CyteType** is a Python package for automated cell type annotation of single-cell RNA-seq data.

## Quick Start

```python
import anndata
import scanpy as sc
import cytetype

# Load and preprocess your data
adata = anndata.read_h5ad("path/to/your/data.h5ad")
sc.tl.rank_genes_groups(adata, groupby='leiden_clusters', method='t-test', key_added='rank_genes_leiden')

# Annotate cell types
cytetype_annotator = cytetype.CyteType(adata, cell_group_key='leiden_clusters', rank_genes_key='rank_genes_leiden')
adata = cytetype_annotator.run(
    bio_context={
        'organisms': ['Homo sapiens'],
        'tissues': ['Brain'],
        'diseases': ['Alzheimer disease']
    }
)

# View results
print(adata.obs.CyteType_leiden_clusters)
```

## Example Report

CyteType generates comprehensive annotation reports with detailed justifications for each cell type assignment. You can see an example of the report structure and analysis depth at: [Hosted Report](https://nygen-labs--cell-annotation-agent-fastapi-app.modal.run/report/97ba2a69-ccfa-4b57-8614-746ce2024333)

The report includes:
- **Detailed cluster annotations** with confidence scores
- **Marker gene analysis** and supporting evidence
- **Alternative annotations** and conflicting markers
- **Biological justifications** for each cell type assignment

## Key Features

*   **Seamless integration** with `AnnData` objects
*   **Type-safe configuration** using Pydantic models
*   **Configurable LLM support** for annotation (optional)

## Installation

You can install CyteType using `pip` or `uv`:

```bash
pip install cytetype
```

## Usage

The `CyteType` class separates expensive data preparation from API calls, making it efficient for multiple annotation requests:

```python
import anndata
import scanpy as sc
import cytetype

# --- Preprocessing ---
adata = anndata.read_h5ad("path/to/your/data.h5ad")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# ... other steps like HVG selection, scaling, PCA, neighbors ...

sc.tl.leiden(adata, key_added='leiden_clusters')
sc.tl.rank_genes_groups(adata, groupby='leiden_clusters', method='t-test', key_added='rank_genes_leiden')

# --- Initialize CyteType ---
cytetype_annotator = cytetype.CyteType(
    adata,
    cell_group_key='leiden_clusters',
    rank_genes_key='rank_genes_leiden'
)

# --- Run annotation ---
adata = cytetype_annotator.run(
    bio_context={
        'organisms': ['Homo sapiens'],
        'tissues': ['Brain', 'Nervous system'],
        'diseases': ['Alzheimer disease']
    },
    n_top_genes=50
)

# Access annotations
print(adata.obs.CyteType_leiden_clusters)
print(adata.uns['CyteType_results'])
```

## Advanced: Configuring the Annotation Model

You can specify different LLM providers and models using the `model_config` parameter. **Important:** Make sure to use an LLM that supports tool use, as CyteType relies on function calling capabilities for accurate annotations.

The currently supported providers are: `google`, `openai`, `xai`, `anthropic`, and `groq`.

```python
# Initialize once
cytetype_annotator = cytetype.CyteType(adata, cell_group_key='leiden_clusters')

# Run annotation with custom model
adata = cytetype_annotator.run(
    bio_context={
        'organisms': ['Homo sapiens'],
        'tissues': ['Brain']
    },
    model_config=[{
        'provider': 'openai',
        'name': 'gpt-4o',
        'apiKey': 'YOUR_API_KEY'
    }]
)
```

**Important:** Handle API keys securely using environment variables or secrets management tools.

## Development

To set up for development:

1.  Clone the repository:
    ```bash
    git clone https://github.com/NygenAnalytics/CyteType.git
    cd cytetype
    ```
2.  Install dependencies using `uv` (includes development tools):
    ```bash
    pip install uv # Install uv if you don't have it
    uv pip sync --all-extras
    ```
3.  Install the package in editable mode:
    ```bash
    uv run pip install -e .
    ```

### Running Checks and Tests

*   **Mypy (Type Checking):** `uv run mypy .`
*   **Ruff (Linting & Formatting):** `uv run ruff check .` and `uv run ruff format .`
*   **Pytest (Unit Tests):** `uv run pytest`


## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## License

This project is licensed under the CC BY-NC-SA 4.0 License - see the [LICENSE](LICENSE) file for details.
