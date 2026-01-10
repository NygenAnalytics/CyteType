"""Shared test fixtures for CyteType tests."""

import pytest
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from typing import Any


@pytest.fixture
def mock_adata() -> anndata.AnnData:
    """Create minimal but realistic AnnData for testing.

    Returns:
        AnnData with 100 cells, 50 genes, 3 clusters, with all required fields.
    """
    n_cells, n_genes = 100, 50
    rng = np.random.default_rng(42)  # Fixed seed for reproducibility

    # Expression matrix (log1p normalized)
    X = rng.lognormal(0, 0.5, (n_cells, n_genes))
    X = np.log1p(X)

    # Observations (cells)
    obs = pd.DataFrame(
        {
            "leiden": rng.choice(["0", "1", "2"], n_cells),
            "batch": rng.choice(["batch1", "batch2"], n_cells),
            "donor": rng.choice(["donor1", "donor2", "donor3"], n_cells),
        }
    )

    # Variables (genes) - with gene symbols
    var = pd.DataFrame(
        {"gene_symbols": [f"GENE{i}" for i in range(n_genes)]},
        index=[f"ENSG{i:011d}" for i in range(n_genes)],
    )

    # Create AnnData
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    # Add UMAP coordinates
    adata.obsm["X_umap"] = rng.standard_normal((n_cells, 2))

    # Run rank_genes_groups (required by CyteType)
    sc.tl.rank_genes_groups(adata, "leiden", method="t-test")

    return adata


@pytest.fixture
def mock_api_response() -> dict[str, Any]:
    """Standard successful API response for mocking.

    Returns:
        Dictionary matching API response format with 3 cluster annotations.
    """
    return {
        "annotations": [
            {
                "clusterId": "1",
                "annotation": "T cells",
                "ontologyTerm": "T cell",
                "ontologyTermID": "CL:0000084",
                "granularAnnotation": "CD8+ T cells",
                "cellState": "activated",
                "justification": "High expression of CD8A, CD8B",
                "supportingMarkers": ["CD8A", "CD8B", "CD3D"],
                "conflictingMarkers": [],
                "missingExpression": "",
                "unexpectedExpression": "",
            },
            {
                "clusterId": "2",
                "annotation": "B cells",
                "ontologyTerm": "B cell",
                "ontologyTermID": "CL:0000236",
                "granularAnnotation": "Naive B cells",
                "cellState": "",
                "justification": "High expression of CD19, MS4A1",
                "supportingMarkers": ["CD19", "MS4A1", "CD79A"],
                "conflictingMarkers": [],
                "missingExpression": "",
                "unexpectedExpression": "",
            },
            {
                "clusterId": "3",
                "annotation": "NK cells",
                "ontologyTerm": "natural killer cell",
                "ontologyTermID": "CL:0000623",
                "granularAnnotation": "",
                "cellState": "",
                "justification": "High expression of NCAM1, NKG7",
                "supportingMarkers": ["NCAM1", "NKG7", "KLRD1"],
                "conflictingMarkers": [],
                "missingExpression": "",
                "unexpectedExpression": "",
            },
        ],
        "summary": {"totalClusters": 3, "annotatedClusters": 3},
        "clusterCategories": [],
        "studyContext": "Test study",
    }
