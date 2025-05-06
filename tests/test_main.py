import pytest
import anndata
import numpy as np
import pandas as pd
from unittest.mock import patch, MagicMock
from typing import Any

# Import the main function to test
from cytetype.main import annotate_anndata

# Import helpers (though testing them here is not ideal long-term)

# Import config and exceptions
from cytetype.config import DEFAULT_API_URL, DEFAULT_POLL_INTERVAL, DEFAULT_TIMEOUT
from cytetype.exceptions import (
    CyteTypeAPIError,
    CyteTypeTimeoutError,
)

# --- Fixtures ---


@pytest.fixture
def mock_adata() -> anndata.AnnData:
    """Creates a basic AnnData object suitable for testing."""
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
    # No need for adata.raw based on current anndata_helpers

    # Simulate rank_genes_groups results for cluster '0', '1', '2'
    # Assume this was run beforehand by the user
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


# --- Test Main Function --- #


@patch("cytetype.main.submit_annotation_job")
@patch("cytetype.main.poll_for_results")
def test_annotate_anndata_success(
    mock_poll: MagicMock, mock_submit: MagicMock, mock_adata: anndata.AnnData
) -> None:
    """Test the main annotation function end-to-end with mocks."""
    job_id = "mock_job_success"
    mock_submit.return_value = job_id
    mock_result: dict[str, list[dict[str, str]]] = {
        "annotations": [
            {"clusterId": "1", "annotation": "Cell Type A"},  # Corresponds to '0'
            {"clusterId": "2", "annotation": "Cell Type B"},  # Corresponds to '1'
            {"clusterId": "3", "annotation": "Cell Type C"},  # Corresponds to '2'
        ]
    }
    mock_poll.return_value = mock_result

    group_key = "leiden"
    result_prefix = "TestCyteType"
    rank_key = "rank_genes_groups"  # Use default key

    adata_result = annotate_anndata(
        mock_adata,
        cell_group_key=group_key,
        rank_genes_key=rank_key,
        results_key_added=result_prefix,
        n_top_genes=5,
    )

    # Check mocks called correctly
    mock_submit.assert_called_once()
    query_arg, url_arg = mock_submit.call_args[0]
    assert url_arg == DEFAULT_API_URL
    assert list(query_arg["markerGenes"].keys()) == ["1", "2", "3"]

    mock_poll.assert_called_once()
    job_id_arg, url_arg, interval_arg, timeout_arg = mock_poll.call_args[0]
    assert job_id_arg == job_id
    assert url_arg == DEFAULT_API_URL
    assert interval_arg == DEFAULT_POLL_INTERVAL
    assert timeout_arg == DEFAULT_TIMEOUT

    # Check results added to AnnData
    assert f"{result_prefix}_results" in adata_result.uns
    assert "job_id" in adata_result.uns[f"{result_prefix}_results"]
    assert adata_result.uns[f"{result_prefix}_results"]["result"] == mock_result

    obs_key = f"{result_prefix}_{group_key}"
    assert obs_key in adata_result.obs
    assert pd.api.types.is_categorical_dtype(adata_result.obs[obs_key])

    # Check annotation mapping
    expected_annotations = []
    ct_map = {"0": 1, "1": 2, "2": 3}
    anno_map = {1: "Cell Type A", 2: "Cell Type B", 3: "Cell Type C"}
    for cluster_label in mock_adata.obs[group_key]:
        cluster_int = ct_map[cluster_label]
        expected_annotations.append(anno_map[cluster_int])

    pd.testing.assert_series_equal(
        adata_result.obs[obs_key],
        pd.Series(expected_annotations, index=mock_adata.obs.index, dtype="category"),
        check_names=False,
    )


@patch("cytetype.main.submit_annotation_job")
def test_annotate_anndata_submit_fails(
    mock_submit: MagicMock, mock_adata: anndata.AnnData
) -> None:
    """Test main function when submission fails (exception propagates)."""
    mock_submit.side_effect = CyteTypeAPIError("Submit failed")

    with pytest.raises(CyteTypeAPIError, match="Submit failed"):
        annotate_anndata(mock_adata, cell_group_key="leiden")


@patch("cytetype.main.submit_annotation_job", return_value="mock_job_poll_fail")
@patch("cytetype.main.poll_for_results")
def test_annotate_anndata_poll_fails(
    mock_poll: MagicMock, mock_submit: MagicMock, mock_adata: anndata.AnnData
) -> None:
    """Test main function when polling fails (exception propagates)."""
    mock_poll.side_effect = CyteTypeTimeoutError("Poll timed out")

    with pytest.raises(CyteTypeTimeoutError, match="Poll timed out"):
        annotate_anndata(mock_adata, cell_group_key="leiden")


@patch("cytetype.main.submit_annotation_job", return_value="mock_job_custom_url")
@patch("cytetype.main.poll_for_results")
def test_annotate_anndata_custom_url_and_params(
    mock_poll: MagicMock, mock_submit: MagicMock, mock_adata: anndata.AnnData
) -> None:
    """Test using a custom API URL and non-default poll/timeout params."""
    custom_url = "http://my-custom-api.com"
    custom_poll = 5
    custom_timeout = 60
    job_id = "mock_job_custom_url"
    mock_result: dict[str, list[Any]] = {"annotations": []}  # Added type hint
    mock_poll.return_value = mock_result

    annotate_anndata(
        mock_adata,
        cell_group_key="leiden",
        api_url=custom_url,
        poll_interval_seconds=custom_poll,
        timeout_seconds=custom_timeout,
    )

    # Check submit was called with custom url
    mock_submit.assert_called_once()
    _, url_arg = mock_submit.call_args[0]
    assert url_arg == custom_url

    # Check poll was called with custom url and params
    mock_poll.assert_called_once()
    job_id_arg, url_arg, interval_arg, timeout_arg = mock_poll.call_args[0]
    assert job_id_arg == job_id
    assert url_arg == custom_url
    assert interval_arg == custom_poll
    assert timeout_arg == custom_timeout


@patch("cytetype.main.submit_annotation_job")
@patch("cytetype.main.poll_for_results")
def test_annotate_anndata_custom_rank_key(
    mock_poll: MagicMock, mock_submit: MagicMock, mock_adata: anndata.AnnData
) -> None:
    """Test using a custom rank_genes_key."""
    job_id = "mock_job_custom_rank"
    mock_submit.return_value = job_id
    mock_result: dict[str, list[Any]] = {"annotations": []}  # Added type hint
    mock_poll.return_value = mock_result
    custom_rank_key = "custom_rank_genes"

    # Run annotation using the custom key
    annotate_anndata(
        mock_adata, cell_group_key="leiden", rank_genes_key=custom_rank_key
    )

    # Check that _get_markers (called internally) would have used the correct key
    # We can check the markerGenes part of the query submitted
    mock_submit.assert_called_once()
    query_arg, _ = mock_submit.call_args[0]
    # The markers should be the same as before because we used the same mock data
    # in adata.uns["custom_rank_genes"]
    assert "markerGenes" in query_arg
    assert list(query_arg["markerGenes"].keys()) == ["1", "2", "3"]
    assert query_arg["markerGenes"]["1"][0] == "gene_0"


# --- TODO ---
# - Add tests specifically for cytetype/anndata_helpers.py
# - Add tests specifically for cytetype/client.py (e.g., more nuanced API responses)
# - Test cases where _get_markers raises ValueError (e.g., group mismatch)
# - Test case where API result format is invalid (should raise CyteTypeAPIError from poll_for_results)
# - Test case for annotation processing error (e.g., non-integer clusterId in response)
