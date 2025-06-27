"""Tests for missing functionality in cytetype.anndata_helpers module."""

import pytest
import anndata
import numpy as np
import pandas as pd
from unittest.mock import patch, MagicMock

from cytetype.anndata_helpers import _extract_sampled_coordinates, _aggregate_metadata


# --- Fixtures ---


@pytest.fixture
def mock_adata_with_coords() -> anndata.AnnData:
    """Creates AnnData object with coordinates for testing coordinate extraction."""
    n_obs, n_vars = 200, 30  # Larger for sampling tests
    rng = np.random.default_rng(42)

    X = rng.poisson(1, size=(n_obs, n_vars)).astype(np.float32)
    X = np.log1p(X)

    # Create groups with different sizes to test sampling
    obs = pd.DataFrame(
        {
            "leiden": [f"{i % 3}" for i in range(n_obs)],  # 3 groups
            "cell_type": [f"type_{i % 2}" for i in range(n_obs)],  # 2 types
            "batch": [f"batch_{i % 4}" for i in range(n_obs)],  # 4 batches
            "treatment": ["treated" if i < 100 else "control" for i in range(n_obs)],
            "numeric_col": rng.random(n_obs),  # Should be ignored in metadata
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )

    var = pd.DataFrame(
        {"gene_symbols": [f"gene_{i}" for i in range(n_vars)]},
        index=[f"gene_{i}" for i in range(n_vars)],
    )

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    # Add 2D coordinates
    adata.obsm["X_umap"] = rng.random((n_obs, 2)) * 10
    # Add higher-dimensional coordinates to test dimension reduction
    adata.obsm["X_pca"] = rng.random((n_obs, 50))

    return adata


# --- Test _extract_sampled_coordinates ---


def test_extract_sampled_coordinates_basic(
    mock_adata_with_coords: anndata.AnnData,
) -> None:
    """Test basic coordinate extraction with sampling."""
    cluster_map = {"0": "1", "1": "2", "2": "3"}
    max_cells = 50  # Small value to force sampling

    coords, labels = _extract_sampled_coordinates(
        adata=mock_adata_with_coords,
        coordinates_key="X_umap",
        group_key="leiden",
        cluster_map=cluster_map,
        max_cells_per_group=max_cells,
    )

    assert coords is not None
    assert len(coords) == len(labels)
    assert len(coords) <= 3 * max_cells  # At most max_cells per group

    # Check that coordinates are 2D
    for coord in coords[:5]:  # Check first few
        assert len(coord) == 2
        assert isinstance(coord[0], float)
        assert isinstance(coord[1], float)

    # Check that labels are mapped correctly
    unique_labels = set(labels)
    assert unique_labels.issubset({"1", "2", "3"})


def test_extract_sampled_coordinates_no_coordinates_key(
    mock_adata_with_coords: anndata.AnnData,
) -> None:
    """Test coordinate extraction when coordinates_key is None."""
    cluster_map = {"0": "1", "1": "2", "2": "3"}

    coords, labels = _extract_sampled_coordinates(
        adata=mock_adata_with_coords,
        coordinates_key=None,
        group_key="leiden",
        cluster_map=cluster_map,
    )

    assert coords is None
    assert labels == []


def test_extract_sampled_coordinates_high_dimensions(
    mock_adata_with_coords: anndata.AnnData,
) -> None:
    """Test coordinate extraction with high-dimensional coordinates (should use first 2)."""
    cluster_map = {"0": "1", "1": "2", "2": "3"}

    coords, labels = _extract_sampled_coordinates(
        adata=mock_adata_with_coords,
        coordinates_key="X_pca",  # 50 dimensions
        group_key="leiden",
        cluster_map=cluster_map,
        max_cells_per_group=1000,  # Don't sample
    )

    assert coords is not None
    assert len(coords) == len(labels)

    # Should still be 2D even though source was 50D
    for coord in coords[:5]:
        assert len(coord) == 2


def test_extract_sampled_coordinates_reproducible_sampling(
    mock_adata_with_coords: anndata.AnnData,
) -> None:
    """Test that sampling is reproducible with same random_state."""
    cluster_map = {"0": "1", "1": "2", "2": "3"}
    max_cells = 30
    random_state = 123

    coords1, labels1 = _extract_sampled_coordinates(
        adata=mock_adata_with_coords,
        coordinates_key="X_umap",
        group_key="leiden",
        cluster_map=cluster_map,
        max_cells_per_group=max_cells,
        random_state=random_state,
    )

    coords2, labels2 = _extract_sampled_coordinates(
        adata=mock_adata_with_coords,
        coordinates_key="X_umap",
        group_key="leiden",
        cluster_map=cluster_map,
        max_cells_per_group=max_cells,
        random_state=random_state,
    )

    assert coords1 == coords2
    assert labels1 == labels2


def test_extract_sampled_coordinates_no_sampling_needed(
    mock_adata_with_coords: anndata.AnnData,
) -> None:
    """Test coordinate extraction when no sampling is needed."""
    cluster_map = {"0": "1", "1": "2", "2": "3"}
    max_cells = 1000  # Much larger than any group

    coords, labels = _extract_sampled_coordinates(
        adata=mock_adata_with_coords,
        coordinates_key="X_umap",
        group_key="leiden",
        cluster_map=cluster_map,
        max_cells_per_group=max_cells,
    )

    assert coords is not None
    assert len(coords) == mock_adata_with_coords.n_obs  # All cells included
    assert len(labels) == mock_adata_with_coords.n_obs


@patch("cytetype.anndata_helpers.logger")
def test_extract_sampled_coordinates_logging(
    mock_logger: MagicMock, mock_adata_with_coords: anndata.AnnData
) -> None:
    """Test that appropriate logging occurs during sampling."""
    cluster_map = {"0": "1", "1": "2", "2": "3"}
    max_cells = 20  # Small to force sampling and logging

    _extract_sampled_coordinates(
        adata=mock_adata_with_coords,
        coordinates_key="X_umap",
        group_key="leiden",
        cluster_map=cluster_map,
        max_cells_per_group=max_cells,
    )

    # Check that info logging was called
    assert mock_logger.info.call_count >= 1
    # Should have logged about extracted coordinates and possibly about sampling


# --- Test _aggregate_metadata ---


def test_aggregate_metadata_basic(mock_adata_with_coords: anndata.AnnData) -> None:
    """Test basic metadata aggregation functionality."""
    result = _aggregate_metadata(
        adata=mock_adata_with_coords, group_key="leiden", min_percentage=10
    )

    assert isinstance(result, dict)

    # Should have entries for each group
    groups = mock_adata_with_coords.obs["leiden"].unique()
    for group in groups:
        assert str(group) in result

    # Should have analyzed categorical columns (excluding group_key)
    for group_name in result:
        assert isinstance(result[group_name], dict)
        # Check that leiden (group_key) is not included
        assert "leiden" not in result[group_name]


def test_aggregate_metadata_min_percentage_filtering(
    mock_adata_with_coords: anndata.AnnData,
) -> None:
    """Test that min_percentage filtering works correctly."""
    # Test with high threshold
    result_high = _aggregate_metadata(
        adata=mock_adata_with_coords,
        group_key="leiden",
        min_percentage=90,  # Very high threshold
    )

    # Test with low threshold
    result_low = _aggregate_metadata(
        adata=mock_adata_with_coords,
        group_key="leiden",
        min_percentage=1,  # Very low threshold
    )

    # High threshold should have fewer or equal entries
    total_entries_high = sum(
        len(group_data.get("treatment", {})) for group_data in result_high.values()
    )
    total_entries_low = sum(
        len(group_data.get("treatment", {})) for group_data in result_low.values()
    )

    assert total_entries_high <= total_entries_low


def test_aggregate_metadata_column_type_filtering(
    mock_adata_with_coords: anndata.AnnData,
) -> None:
    """Test that only categorical/string columns are processed."""
    result = _aggregate_metadata(
        adata=mock_adata_with_coords, group_key="leiden", min_percentage=1
    )

    # Check that numeric_col (float) is not included
    for group_name in result:
        assert "numeric_col" not in result[group_name]

    # Check that categorical/string columns are included
    for group_name in result:
        # At least one of the categorical columns should be present
        categorical_cols = ["cell_type", "batch", "treatment"]
        # Note: might not be present if percentages are too low
        # So we'll just check the structure is correct if present
        for col in categorical_cols:
            if col in result[group_name]:
                assert isinstance(result[group_name][col], dict)
                for value, percentage in result[group_name][col].items():
                    assert isinstance(value, str)
                    assert isinstance(percentage, int)
                    assert 0 <= percentage <= 100


def test_aggregate_metadata_empty_adata() -> None:
    """Test metadata aggregation with minimal AnnData."""
    # Create minimal AnnData
    n_obs = 10
    obs = pd.DataFrame(
        {
            "group": ["A"] * 5 + ["B"] * 5,  # First 5 are A, next 5 are B
            "category": ["X", "Y", "X", "Y", "X"]
            + ["X", "Y", "X", "Y", "X"],  # Each group has 3 X and 2 Y
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )

    var = pd.DataFrame(index=[f"gene_{i}" for i in range(5)])
    X = np.random.random((n_obs, 5))

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    result = _aggregate_metadata(
        adata=adata,
        group_key="group",
        min_percentage=30,  # Should include both X (60%) and Y (40%)
    )

    assert "A" in result
    assert "B" in result
    # Each group should have category with X at 60% and Y at 40%
    assert "category" in result["A"]
    assert "X" in result["A"]["category"]
    assert result["A"]["category"]["X"] == 60
    assert "Y" in result["A"]["category"]
    assert result["A"]["category"]["Y"] == 40


def test_aggregate_metadata_no_categorical_columns() -> None:
    """Test metadata aggregation when there are no categorical columns besides group_key."""
    n_obs = 20
    obs = pd.DataFrame(
        {
            "group": ["A", "B"] * 10,
            "numeric1": np.random.random(n_obs),
            "numeric2": np.random.random(n_obs),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )

    var = pd.DataFrame(index=[f"gene_{i}" for i in range(5)])
    X = np.random.random((n_obs, 5))

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    result = _aggregate_metadata(adata=adata, group_key="group", min_percentage=10)

    # Should have entries for groups but no metadata
    assert "A" in result
    assert "B" in result
    assert result["A"] == {}
    assert result["B"] == {}


def test_aggregate_metadata_percentage_calculations() -> None:
    """Test that percentage calculations are correct."""
    # Create controlled data for precise percentage testing
    obs = pd.DataFrame(
        {
            "group": ["A"] * 100,  # Single group for simplicity
            "treatment": ["control"] * 80
            + ["treated"] * 20,  # 80% control, 20% treated
        },
        index=[f"cell_{i}" for i in range(100)],
    )

    var = pd.DataFrame(index=["gene_1"])
    X = np.random.random((100, 1))

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    result = _aggregate_metadata(
        adata=adata,
        group_key="group",
        min_percentage=15,  # Should include both values
    )

    assert result["A"]["treatment"]["control"] == 80
    assert result["A"]["treatment"]["treated"] == 20
