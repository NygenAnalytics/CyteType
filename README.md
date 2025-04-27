# CyteType Python Client

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) <!-- TODO: Update license badge if needed -->

A Python package to interact with the CyteType API for cell type annotation of single-cell RNA-seq data stored in AnnData objects.

## Installation

Currently, you can install the package directly from the source code:

```bash
# Clone the repository (if you haven't already)
# git clone https://github.com/yourusername/cytetype.git
# cd cytetype

pip install .
```

Alternatively, if the package is published to PyPI:

```bash
pip install cytetype
```

**Dependencies:**

The package requires the following Python libraries:

*   `anndata`
*   `numpy`
*   `pandas`
*   `scanpy`
*   `requests`
*   `tqdm`

These will be installed automatically when using `pip`.

## Usage

Here's a basic example of how to use the package:

```python
import scanpy as sc
import anndata
from cytetype.annotator import annotate_anndata

# Load your AnnData object
# Ensure it has raw counts in adata.raw.X
# Ensure it has log1p normalized data in adata.X
# Ensure you have run clustering (e.g., Leiden) and stored in adata.obs['leiden']
adata = anndata.read_h5ad("path/to/your/data.h5ad")

# --- Preprocessing (Example - adapt as needed) ---
# Make sure raw counts are stored if not already
if adata.raw is None:
    adata.raw = adata.copy()

# Basic filtering, normalization, PCA, neighbors, clustering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, key_added='leiden') # Use 'leiden' as cell_group_key
# -------------------------------------------------

# Define biological context (optional but recommended)
org = "Homo sapiens"
tissue_list = ["Peripheral blood mononuclear cell", "Blood"]
disease_list = ["Healthy"]

# Run CyteType annotation
try:
    adata = annotate_anndata(
        adata=adata,
        cell_group_key='leiden',  # Key in adata.obs containing cluster labels
        organism=org,
        tissues=tissue_list,
        diseases=disease_list,
        # Optional: specify a different API endpoint
        # api_url="http://your-custom-api-url/"
    )

    # Results are added to:
    # adata.obs['CyteType_leiden'] (annotations)
    # adata.uns['CyteType_results'] (raw API response and status)

    print("Annotation successful!")
    print(adata.obs['CyteType_leiden'].value_counts())

except Exception as e:
    print(f"Annotation failed: {e}")
    # Check adata.uns['CyteType_results'] for error details if available

```

## Development & Testing

To set up for development:

```bash
pip install -r requirements-dev.txt
pip install -e .
```

Run tests using pytest:

```bash
pytest
```

<!-- TODO: Add contribution guidelines if applicable -->
<!-- TODO: Add link to API documentation if available -->
