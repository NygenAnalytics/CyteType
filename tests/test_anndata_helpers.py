import pytest
import anndata
import numpy as np
import pandas as pd

# Import helpers to test
from cytetype.anndata_helpers import (
    _validate_adata,
    _calculate_pcent,
    _get_markers,
)

# --- Fixtures ---


# TODO: Consider if this fixture should be shared via conftest.py
# if it's also needed in test_main.py or other files.
@pytest.fixture
def mock_adata() -> anndata.AnnData:
    """Creates a basic AnnData object suitable for testing helpers."""
    n_obs, n_vars = 100, 50
    rng = np.random.default_rng(0)
    # Simulate log1p normalized data directly in X
    X = rng.poisson(1, size=(n_obs, n_vars)).astype(np.float32)
    X = np.log1p(X)

    obs = pd.DataFrame(
        {
            "cell_type": [f"type_{i % 3}" for i in range(n_obs)],
            "leiden": [f"{i % 3}" for i in range(n_obs)],  # Clusters as strings
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    # Ensure var_names match X dimensions
    var_index = [f"gene_{i}" for i in range(n_vars)]
    var = pd.DataFrame(index=var_index)
    var.index.name = "gene_id"  # Use 'gene_id' for index name for clarity
    var["gene_symbols"] = var.index  # Add the required gene_symbols column

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    # Simulate rank_genes_groups results for cluster '0', '1', '2'
    group_names = ["0", "1", "2"]
    n_genes_ranked = 20
    names_list = [
        [f"gene_{i * 3 + j}" for i in range(n_genes_ranked)]
        for j in range(len(group_names))
    ]

    dtype = [(name, "U20") for name in group_names]
    names_arr = np.array(list(zip(*names_list)), dtype=dtype)

    # Default rank_genes_key
    default_rank_key = "rank_genes_groups"
    adata.uns[default_rank_key] = {
        "params": {"groupby": "leiden", "method": "t-test"},
        "names": names_arr,
    }
    # Add another key for testing the rank_genes_key parameter
    custom_rank_key = "custom_rank_genes"
    adata.uns[custom_rank_key] = {
        "params": {"groupby": "leiden", "method": "t-test"},
        "names": names_arr,  # Use same data for simplicity
    }

    return adata


# --- Test Helper Functions ---


def test_validate_adata_success(mock_adata: anndata.AnnData) -> None:
    """Test validation passes with a correctly formatted AnnData object."""
    _validate_adata(
        mock_adata, "leiden", "rank_genes_groups", gene_symbols_col="gene_symbols"
    )  # Should not raise


def test_validate_adata_missing_group(mock_adata: anndata.AnnData) -> None:
    """Test validation fails if cell_group_key is missing."""
    with pytest.raises(KeyError, match="not found in `adata.obs`"):
        _validate_adata(
            mock_adata,
            "unknown_key",
            "rank_genes_groups",
            gene_symbols_col="gene_symbols",
        )


def test_validate_adata_missing_x(mock_adata: anndata.AnnData) -> None:
    """Test validation fails if adata.X is missing."""
    mock_adata.X = None
    with pytest.raises(ValueError, match="`adata.X` is required"):
        _validate_adata(
            mock_adata, "leiden", "rank_genes_groups", gene_symbols_col="gene_symbols"
        )


def test_validate_adata_rank_key_missing(mock_adata: anndata.AnnData) -> None:
    """Test validation fails if rank_genes_key is missing in uns."""
    with pytest.raises(KeyError, match="not found in `adata.uns`"):
        _validate_adata(
            mock_adata,
            "leiden",
            "nonexistent_rank_key",
            gene_symbols_col="gene_symbols",
        )


def test_calculate_pcent(mock_adata: anndata.AnnData) -> None:
    """Test percentage calculation using adata.X."""
    # Map original cluster labels ('0', '1', '2') to API cluster IDs ('1', '2', '3') as strings
    ct_map = {
        str(x): str(n + 1)
        for n, x in enumerate(sorted(mock_adata.obs["leiden"].unique().tolist()))
    }
    clusters_str = [ct_map[str(x)] for x in mock_adata.obs["leiden"].values.tolist()]

    pcent = _calculate_pcent(
        mock_adata,
        clusters_str,  # Pass list of strings
        gene_names=mock_adata.var_names.to_list(),
        batch_size=10,
    )
    assert isinstance(pcent, dict)
    assert len(pcent) == mock_adata.n_vars  # Should have entry for each gene
    # Check a specific gene and cluster (values depend on mock data & log1p)
    assert "gene_0" in pcent
    assert "1" in pcent["gene_0"]  # Cluster IDs are strings '1', '2', '3'
    # Since input is log1p(counts+1), (X > 0) should be equivalent to (raw > 0)
    # for typical count data, so percentage should still be reasonable.
    assert 0 <= pcent["gene_0"]["1"] <= 100


def test_get_markers(mock_adata: anndata.AnnData) -> None:
    """Test marker gene extraction."""
    # Map original cluster labels ('0', '1', '2') to API cluster IDs ('1', '2', '3') as strings
    ct_map = {"0": "1", "1": "2", "2": "3"}
    n_top = 5
    rank_key = "rank_genes_groups"
    markers = _get_markers(
        mock_adata,
        "leiden",
        rank_key,
        ct_map,
        gene_symbols_col="gene_symbols",
        n_top_genes=n_top,
    )
    assert isinstance(markers, dict)
    assert list(markers.keys()) == ["1", "2", "3"]  # API cluster IDs as strings
    assert len(markers["1"]) == n_top
    assert markers["1"][0] == "gene_0"  # Based on mock rank_genes_groups
    assert markers["2"][0] == "gene_1"
    assert markers["3"][0] == "gene_2"


# Add a test for validation failure due to rank_genes groupby mismatch
def test_validate_adata_groupby_mismatch(mock_adata: anndata.AnnData) -> None:
    """Test validation fails if rank_genes_groups groupby mismatches cell_group_key."""
    # Modify the mock adata to have a mismatch
    mock_adata.uns["rank_genes_groups"]["params"]["groupby"] = "different_group"
    # Update the regex to be more specific to the expected error message format
    expected_error_msg = r"`rank_genes_groups` run with groupby=\'different_group\', expected \'leiden\'."
    with pytest.raises(ValueError, match=expected_error_msg):
        _validate_adata(
            mock_adata, "leiden", "rank_genes_groups", gene_symbols_col="gene_symbols"
        )


# Add a test for _get_markers failure due to ct_map mismatch
def test_get_markers_ct_map_mismatch(mock_adata: anndata.AnnData) -> None:
    """Test _get_markers fails if rank_genes group is not in ct_map."""
    # Use a ct_map that is missing a mapping for one of the groups ('2')
    ct_map = {"0": "1", "1": "2"}  # Missing mapping for group '2'
    n_top = 5
    rank_key = "rank_genes_groups"
    with pytest.raises(ValueError, match="Internal inconsistency"):
        _get_markers(
            mock_adata,
            "leiden",
            rank_key,
            ct_map,
            gene_symbols_col="gene_symbols",
            n_top_genes=n_top,
        )
