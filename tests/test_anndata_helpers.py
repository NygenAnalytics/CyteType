import pytest
import anndata
import numpy as np
import pandas as pd
from io import StringIO
from loguru import logger

# Import helpers to test
from cytetype.anndata_helpers import (
    _validate_adata,
    _calculate_pcent,
    _get_markers,
    _is_gene_id_like,
    _validate_gene_symbols_column,
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
        mock_adata,
        "leiden",
        "rank_genes_groups",
        gene_symbols_col="gene_symbols",
        coordinates_key="X_umap",
    )  # Should not raise


def test_validate_adata_missing_group(mock_adata: anndata.AnnData) -> None:
    """Test validation fails if cell_group_key is missing."""
    with pytest.raises(KeyError, match="not found in `adata.obs`"):
        _validate_adata(
            mock_adata,
            "unknown_key",
            "rank_genes_groups",
            gene_symbols_col="gene_symbols",
            coordinates_key="X_umap",
        )


def test_validate_adata_missing_x(mock_adata: anndata.AnnData) -> None:
    """Test validation fails if adata.X is missing."""
    mock_adata.X = None
    with pytest.raises(ValueError, match="`adata.X` is required"):
        _validate_adata(
            mock_adata,
            "leiden",
            "rank_genes_groups",
            gene_symbols_col="gene_symbols",
            coordinates_key="X_umap",
        )


def test_validate_adata_rank_key_missing(mock_adata: anndata.AnnData) -> None:
    """Test validation fails if rank_genes_key is missing in uns."""
    with pytest.raises(KeyError, match="not found in `adata.uns`"):
        _validate_adata(
            mock_adata,
            "leiden",
            "nonexistent_rank_key",
            gene_symbols_col="gene_symbols",
            coordinates_key="X_umap",
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
            mock_adata,
            "leiden",
            "rank_genes_groups",
            gene_symbols_col="gene_symbols",
            coordinates_key="X_umap",
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


def test_get_markers_zero_n_top_genes_error(mock_adata: anndata.AnnData) -> None:
    """Test _get_markers raises error when n_top_genes is 0 for all groups."""
    ct_map = {"0": "1", "1": "2", "2": "3"}

    # Capture loguru logs
    log_capture = StringIO()
    handler_id = logger.add(log_capture, level="WARNING", format="{message}")

    try:
        with pytest.raises(ValueError, match="No marker genes found for any group"):
            _get_markers(
                mock_adata,
                "leiden",
                "rank_genes_groups",
                ct_map,
                gene_symbols_col="gene_symbols",
                n_top_genes=0,  # This will result in empty top_genes for all groups
            )

        # Check that warnings were logged for all groups
        log_output = log_capture.getvalue()
        assert "No top genes found for group '0' (cluster '1')" in log_output
        assert "No top genes found for group '1' (cluster '2')" in log_output
        assert "No top genes found for group '2' (cluster '3')" in log_output
    finally:
        logger.remove(handler_id)


def test_get_markers_empty_mdf_error(mock_adata: anndata.AnnData) -> None:
    """Test _get_markers raises error when rank_genes dataframe is empty."""
    # Create an empty dataframe to simulate no ranked genes
    empty_df = pd.DataFrame()
    mock_adata.uns["rank_genes_groups"]["names"] = empty_df

    ct_map = {"0": "1", "1": "2", "2": "3"}

    # Capture loguru logs
    log_capture = StringIO()
    handler_id = logger.add(log_capture, level="WARNING", format="{message}")

    try:
        with pytest.raises(ValueError, match="No marker genes found for any group"):
            _get_markers(
                mock_adata,
                "leiden",
                "rank_genes_groups",
                ct_map,
                gene_symbols_col="gene_symbols",
                n_top_genes=5,
            )
    finally:
        logger.remove(handler_id)


# --- Test Gene Symbols Validation Functions ---


def test_is_gene_id_like_ensembl_ids() -> None:
    """Test detection of Ensembl gene IDs."""
    # Human Ensembl IDs
    assert _is_gene_id_like("ENSG00000000003")
    assert _is_gene_id_like("ENSG99999999999")

    # Mouse Ensembl IDs
    assert _is_gene_id_like("ENSMUSG00000000001")

    # Other species Ensembl IDs
    assert _is_gene_id_like("ENSDARG00000000001")  # zebrafish

    # Non-Ensembl IDs should return False
    assert not _is_gene_id_like("ENSG0000000000")  # too short
    assert not _is_gene_id_like("ENSG000000000003")  # too long
    assert not _is_gene_id_like("ENSG0000000000A")  # letter in number part


def test_is_gene_id_like_refseq_ids() -> None:
    """Test detection of RefSeq gene IDs."""
    # RefSeq mRNA IDs
    assert _is_gene_id_like("NM_000001")
    assert _is_gene_id_like("NM_123456")
    assert _is_gene_id_like("NM_999999")

    # RefSeq protein IDs
    assert _is_gene_id_like("XM_000001")
    assert _is_gene_id_like("XR_000001")
    assert _is_gene_id_like("NR_000001")

    # Non-RefSeq patterns should return False
    assert not _is_gene_id_like("NM_ABC")  # letters in number part
    assert not _is_gene_id_like("AM_000001")  # wrong prefix
    assert not _is_gene_id_like("NM000001")  # missing underscore


def test_is_gene_id_like_numeric_ids() -> None:
    """Test detection of purely numeric gene IDs."""
    assert _is_gene_id_like("12345")
    assert _is_gene_id_like("999999")
    assert _is_gene_id_like("1")

    # Non-numeric should return False
    assert not _is_gene_id_like("12345A")
    assert not _is_gene_id_like("A12345")


def test_is_gene_id_like_other_patterns() -> None:
    """Test detection of other database-style gene IDs."""
    # Long alphanumeric with dots/underscores
    assert _is_gene_id_like("A123456789.123456")
    assert _is_gene_id_like("ABC123_DEF456_GHI789")

    # Should return False for shorter patterns
    assert not _is_gene_id_like("A123.456")
    assert not _is_gene_id_like("ABC_DEF")


def test_is_gene_id_like_gene_symbols() -> None:
    """Test that common gene symbols are not detected as gene IDs."""
    # Common human gene symbols
    assert not _is_gene_id_like("TSPAN6")
    assert not _is_gene_id_like("DPM1")
    assert not _is_gene_id_like("SCYL3")
    assert not _is_gene_id_like("TP53")
    assert not _is_gene_id_like("BRCA1")
    assert not _is_gene_id_like("GAPDH")

    # Gene symbols with numbers
    assert not _is_gene_id_like("H3F3A")
    assert not _is_gene_id_like("HIST1H1A")


def test_is_gene_id_like_edge_cases() -> None:
    """Test edge cases for gene ID detection."""
    # Empty/None/whitespace
    assert not _is_gene_id_like("")
    assert not _is_gene_id_like("   ")
    assert not _is_gene_id_like("None")

    # Mixed case (should still work)
    assert _is_gene_id_like("ensg00000000003")
    assert not _is_gene_id_like("tspan6")


def test_validate_gene_symbols_column_valid_symbols(
    mock_adata: anndata.AnnData,
) -> None:
    """Test validation passes with valid gene symbols."""
    # Use default mock_adata which has gene symbols like "gene_0", "gene_1", etc.
    # These should not be detected as gene IDs

    # Capture loguru logs
    log_capture = StringIO()
    handler_id = logger.add(log_capture, level="WARNING", format="{message}")

    try:
        _validate_gene_symbols_column(mock_adata, "gene_symbols")

        # Should not raise any exceptions or warnings
        log_output = log_capture.getvalue()
        assert log_output == ""
    finally:
        logger.remove(handler_id)


def test_validate_gene_symbols_column_high_gene_id_percentage(
    mock_adata: anndata.AnnData,
) -> None:
    """Test validation warns when >50% of values look like gene IDs."""
    # Replace gene symbols with Ensembl IDs for most genes
    ensembl_ids = [f"ENSG{i:011d}" for i in range(len(mock_adata.var))]
    mock_adata.var["gene_symbols"] = ensembl_ids

    # Capture loguru logs
    log_capture = StringIO()
    handler_id = logger.add(log_capture, level="WARNING", format="{message}")

    try:
        _validate_gene_symbols_column(mock_adata, "gene_symbols")

        # Should warn about gene IDs
        log_output = log_capture.getvalue()
        assert "appears to contain gene IDs rather than gene symbols" in log_output
        assert "The annotation might not be accurate" in log_output
        assert "100.0% of values look like gene IDs" in log_output
    finally:
        logger.remove(handler_id)


def test_validate_gene_symbols_column_mixed_content(
    mock_adata: anndata.AnnData,
) -> None:
    """Test validation with mixed gene symbols and IDs."""
    # Mix of valid gene symbols and Ensembl IDs (30% gene IDs)
    gene_symbols = []
    for i in range(len(mock_adata.var)):
        if i % 10 < 3:  # 30% gene IDs
            gene_symbols.append(f"ENSG{i:011d}")
        else:  # 70% gene symbols
            gene_symbols.append(f"GENE{i}")

    mock_adata.var["gene_symbols"] = gene_symbols

    # Capture loguru logs
    log_capture = StringIO()
    handler_id = logger.add(log_capture, level="WARNING", format="{message}")

    try:
        _validate_gene_symbols_column(mock_adata, "gene_symbols")

        # Should warn but not raise error (30% > 20% but < 50%)
        log_output = log_capture.getvalue()
        assert "contains 30.0% values that look like gene IDs" in log_output
        assert "Please verify this column contains gene symbols" in log_output
    finally:
        logger.remove(handler_id)


def test_validate_gene_symbols_column_empty_column(mock_adata: anndata.AnnData) -> None:
    """Test validation with empty gene symbols column."""
    # Create empty column
    mock_adata.var["empty_gene_symbols"] = [None] * len(mock_adata.var)

    # Capture loguru logs
    log_capture = StringIO()
    handler_id = logger.add(log_capture, level="WARNING", format="{message}")

    try:
        _validate_gene_symbols_column(mock_adata, "empty_gene_symbols")

        # Should warn about empty column
        log_output = log_capture.getvalue()
        assert "is empty or contains only NaN values" in log_output
    finally:
        logger.remove(handler_id)


def test_validate_gene_symbols_column_missing_column(
    mock_adata: anndata.AnnData,
) -> None:
    """Test validation with missing gene symbols column."""
    # This should be caught by the KeyError check in _validate_adata
    # but test the function directly
    with pytest.raises(KeyError):
        _validate_gene_symbols_column(mock_adata, "nonexistent_column")


def test_validate_adata_with_gene_symbols_validation(
    mock_adata: anndata.AnnData,
) -> None:
    """Test that _validate_adata calls gene symbols validation."""
    # Replace gene symbols with mostly Ensembl IDs to trigger warning
    ensembl_ids = [f"ENSG{i:011d}" for i in range(len(mock_adata.var))]
    mock_adata.var["gene_symbols"] = ensembl_ids

    # Capture loguru logs
    log_capture = StringIO()
    handler_id = logger.add(log_capture, level="WARNING", format="{message}")

    try:
        _validate_adata(
            mock_adata,
            "leiden",
            "rank_genes_groups",
            gene_symbols_col="gene_symbols",
            coordinates_key="X_umap",
        )

        # Should include gene symbols validation warning
        log_output = log_capture.getvalue()
        assert "appears to contain gene IDs rather than gene symbols" in log_output
    finally:
        logger.remove(handler_id)


def test_validate_adata_missing_gene_symbols_column(
    mock_adata: anndata.AnnData,
) -> None:
    """Test validation fails if gene_symbols_col is missing."""
    with pytest.raises(
        KeyError, match="Column 'nonexistent_column' not found in `adata.var`"
    ):
        _validate_adata(
            mock_adata,
            "leiden",
            "rank_genes_groups",
            gene_symbols_col="nonexistent_column",
            coordinates_key="X_umap",
        )
