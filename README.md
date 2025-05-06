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
*   Allows configuration of the underlying Large Language Model (LLM) used for annotation (optional).

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
    n_top_genes=50,                      # Number of top marker genes per cluster to submit
)

# Access the cell type annotations that were added to the AnnData object
print(adata.obs.CyteType_leiden_clusters)

# Get detailed information about the annotation results
print (adata.uns['CyteType_leiden_clusters'])

```

## Advanced: Configuring the Annotation Model

By default, CyteType uses the default model configured on the backend API service. However, you can specify a different LLM provider and model using the `model_config` parameter in `annotate_anndata`. This is useful if you want to experiment with different models or use a specific provider you have access to.

The currently supported providers are: `google`, `openai`, `xai`, `anthropic`, and `groq`.

The `model_config` parameter expects a list of dictionaries, where each dictionary specifies a model configuration. The structure follows these fields:

*   **Provider (`provider`):** Must be one of the supported providers listed above.
*   **Model Name (`name`, optional):** The specific model you want to use (e.g., `"gpt-4o"`, `"claude-3-opus-20240229"`). If omitted, the default model for the chosen provider will be used by the backend service.
*   **API Key (`apiKey`, optional):** Your API key for the chosen provider. If omitted, you can still access the default models for each provider, but you may be subject to stricter rate limits imposed by the backend service.
*   **Base URL (`baseUrl`, optional):** Use this to specify a custom API endpoint, such as a self-hosted model or a proxy service like OpenRouter.

```python
# Example using a different provider (adjust keys as needed per provider)
adata = annotate_anndata(
    adata=adata,
    cell_group_key='leiden_clusters',
    rank_genes_key='rank_genes_leiden',
    model_config=[
        {
            "provider": "groq", # Or "google", "openai", etc.
            "name": "meta-llama/llama-4-maverick-17b-128e-instruct",
            "apiKey": "YOUR_GROQ_API_KEY", # Optional, but recommended to avoid rate limits
            # "baseUrl": "https://openrouter.ai/api/v1" # Example using OpenRouter
        }
    ]
)
```

**Important:** Ensure you handle API keys securely. Avoid hardcoding them directly in your scripts. Use environment variables or a secrets management tool.

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
