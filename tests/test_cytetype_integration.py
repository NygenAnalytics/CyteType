"""Integration tests for CyteType class."""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
from pydantic import ValidationError
from typing import Any
import anndata

from cytetype import CyteType
from cytetype.api.exceptions import RateLimitError, AuthenticationError
from cytetype.api import UploadResponse


def test_fixture_works(mock_adata: anndata.AnnData) -> None:
    """Test that our mock_adata fixture is properly structured."""
    # Verify required fields exist
    assert "leiden" in mock_adata.obs.columns
    assert "gene_symbols" in mock_adata.var.columns
    assert "X_umap" in mock_adata.obsm
    assert "rank_genes_groups" in mock_adata.uns
    assert mock_adata.X is not None
    assert mock_adata.shape == (100, 50)


def test_cytetype_initialization(mock_adata: anndata.AnnData) -> None:
    """Test CyteType initializes successfully with valid data."""
    ct = CyteType(mock_adata, group_key="leiden")

    # Verify instance attributes set correctly
    assert ct.adata is mock_adata
    assert ct.group_key == "leiden"
    assert ct.rank_key == "rank_genes_groups"
    assert ct.gene_symbols_column == "gene_symbols"

    # Verify data preparation completed
    assert len(ct.clusters) == len(mock_adata)
    assert len(ct.cluster_map) == 3  # 3 unique leiden clusters
    assert ct.marker_genes is not None
    assert len(ct.marker_genes) == 3
    assert ct.expression_percentages is not None
    assert ct.visualization_data is not None
    assert ct.visualization_data["coordinates"] is not None


@pytest.fixture(autouse=True)
def mock_internal_artifact_pipeline(monkeypatch: pytest.MonkeyPatch) -> None:
    """Avoid file and network work for run() in tests by mocking internals."""

    def _save_vars(*args: Any, **kwargs: Any) -> None:
        return None

    def _save_obs(*args: Any, **kwargs: Any) -> None:
        return None

    def _upload_obs(*args: Any, **kwargs: Any) -> UploadResponse:
        return UploadResponse(
            upload_id="obs_upload_123",
            file_kind="obs_duckdb",
            file_name="obs.duckdb",
            size_bytes=100,
            expires_at="2099-01-01T00:00:00Z",
        )

    def _upload_vars(*args: Any, **kwargs: Any) -> UploadResponse:
        return UploadResponse(
            upload_id="vars_upload_456",
            file_kind="vars_h5",
            file_name="vars.h5",
            size_bytes=100,
            expires_at="2099-01-01T00:00:00Z",
        )

    monkeypatch.setattr("cytetype.main.save_features_matrix", _save_vars)
    monkeypatch.setattr("cytetype.main.save_obs_duckdb_file", _save_obs)
    monkeypatch.setattr("cytetype.main.upload_obs_duckdb_file", _upload_obs)
    monkeypatch.setattr("cytetype.main.upload_vars_h5_file", _upload_vars)


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_success(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
) -> None:
    """Test CyteType.run() end-to-end with mocked API."""
    # Setup mocks
    mock_submit.return_value = "test_job_123"
    mock_wait.return_value = mock_api_response

    # Initialize and run
    ct = CyteType(mock_adata, group_key="leiden")
    result_adata = ct.run(study_context="Test human PBMC study")

    # Verify API calls made
    mock_submit.assert_called_once()
    mock_wait.assert_called_once()
    payload = mock_submit.call_args.args[2]
    assert payload["uploaded_files"] == {
        "obs_duckdb": "obs_upload_123",
        "vars_h5": "vars_upload_456",
    }

    # Verify result is the modified adata
    assert result_adata is mock_adata

    # Verify annotation columns created
    assert "cytetype_annotation_leiden" in result_adata.obs.columns
    assert "cytetype_cellOntologyTerm_leiden" in result_adata.obs.columns
    assert "cytetype_cellOntologyTermID_leiden" in result_adata.obs.columns
    assert "cytetype_cellState_leiden" in result_adata.obs.columns

    # Verify annotations applied correctly
    annotations = result_adata.obs["cytetype_annotation_leiden"].unique()
    assert (
        "T cells" in annotations
        or "B cells" in annotations
        or "NK cells" in annotations
    )

    # Verify results stored in adata.uns
    assert "cytetype_results" in result_adata.uns
    assert "cytetype_jobDetails" in result_adata.uns
    assert result_adata.uns["cytetype_jobDetails"]["job_id"] == "test_job_123"


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_auto_uploads_artifacts(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
) -> None:
    """Test run() auto-creates upload references in /annotate payload."""
    mock_submit.return_value = "job_with_uploads"
    mock_wait.return_value = mock_api_response

    ct = CyteType(mock_adata, group_key="leiden")
    ct.run(study_context="Test")

    call_args = mock_submit.call_args.args
    payload = call_args[2]
    assert payload["uploaded_files"] == {
        "obs_duckdb": "obs_upload_123",
        "vars_h5": "vars_upload_456",
    }


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_artifact_failure_raises_by_default(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test run() raises when artifact build fails and require_artifacts=True (default)."""
    mock_submit.return_value = "job_no_artifacts"
    mock_wait.return_value = mock_api_response

    monkeypatch.setattr(
        "cytetype.main.save_features_matrix",
        MagicMock(side_effect=RuntimeError("disk full")),
    )

    ct = CyteType(mock_adata, group_key="leiden")
    with pytest.raises(RuntimeError, match="disk full"):
        ct.run(study_context="Test")


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_artifact_failure_continues_when_not_required(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test run() proceeds with partial uploads when require_artifacts=False.

    vars.h5 save fails but obs.duckdb still succeeds independently.
    """
    mock_submit.return_value = "job_no_artifacts"
    mock_wait.return_value = mock_api_response

    monkeypatch.setattr(
        "cytetype.main.save_features_matrix",
        MagicMock(side_effect=RuntimeError("disk full")),
    )

    ct = CyteType(mock_adata, group_key="leiden")
    result = ct.run(study_context="Test", require_artifacts=False)

    # Job should still complete successfully
    assert result is not None
    assert mock_submit.called

    # obs.duckdb should have succeeded independently
    payload = mock_submit.call_args.args[2]
    assert payload["uploaded_files"] == {"obs_duckdb": "obs_upload_123"}


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_cleanup_artifacts(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test run() can cleanup generated artifact files when requested."""
    mock_submit.return_value = "job_cleanup"
    mock_wait.return_value = mock_api_response

    vars_path = tmp_path / "vars.h5"
    obs_path = tmp_path / "obs.duckdb"

    def _save_vars(*args: Any, **kwargs: Any) -> None:
        Path(kwargs["out_file"]).write_bytes(b"vars")

    def _save_obs(*args: Any, **kwargs: Any) -> None:
        Path(kwargs["out_file"]).write_bytes(b"obs")

    monkeypatch.setattr("cytetype.main.save_features_matrix", _save_vars)
    monkeypatch.setattr("cytetype.main.save_obs_duckdb_file", _save_obs)

    ct = CyteType(mock_adata, group_key="leiden")
    ct.run(
        study_context="Test",
        vars_h5_path=str(vars_path),
        obs_duckdb_path=str(obs_path),
        cleanup_artifacts=True,
    )

    assert not vars_path.exists()
    assert not obs_path.exists()


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_override_existing_results(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
) -> None:
    """Test override_existing_results parameter prevents/allows overwriting."""
    mock_submit.return_value = "job_first"
    mock_wait.return_value = mock_api_response

    ct = CyteType(mock_adata, group_key="leiden")

    # First run - should succeed
    ct.run(study_context="First run")
    assert "cytetype_jobDetails" in ct.adata.uns

    # Second run without override - should raise error
    mock_submit.return_value = "job_second"
    with pytest.raises(ValueError, match="already exist"):
        ct.run(study_context="Second run", override_existing_results=False)

    # Second run with override - should succeed
    ct.run(study_context="Second run", override_existing_results=True)
    assert ct.adata.uns["cytetype_jobDetails"]["job_id"] == "job_second"


def test_cytetype_get_results_local(
    mock_adata: anndata.AnnData, mock_api_response: dict[str, Any]
) -> None:
    """Test get_results() retrieves from local storage."""
    import json

    ct = CyteType(mock_adata, group_key="leiden")

    # Manually store results (simulating previous run)
    ct.adata.uns["cytetype_results"] = {
        "job_id": "stored_job",
        "result": json.dumps(mock_api_response),
    }

    # Retrieve results
    results = ct.get_results()

    assert results is not None
    assert results["annotations"] == mock_api_response["annotations"]
    assert len(results["annotations"]) == 3


@patch("cytetype.main.fetch_remote_results")
def test_cytetype_get_results_remote(
    mock_fetch: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
) -> None:
    """Test get_results() fetches from API when not local."""
    mock_fetch.return_value = mock_api_response

    ct = CyteType(mock_adata, group_key="leiden")

    # Store job details (simulating previous run)
    ct.adata.uns["cytetype_jobDetails"] = {
        "job_id": "remote_job",
        "api_url": "https://api.test",
    }

    # Retrieve results (should fetch from API)
    results = ct.get_results()

    assert results is not None
    assert results == mock_api_response
    mock_fetch.assert_called_once()


def test_cytetype_initialization_with_auth_token(mock_adata: anndata.AnnData) -> None:
    """Test CyteType can store auth_token in __init__."""
    ct = CyteType(mock_adata, group_key="leiden", auth_token="test_token_123")

    assert ct.auth_token == "test_token_123"
    assert ct.api_url == "https://prod.cytetype.nygen.io"


def test_cytetype_no_coordinates(mock_adata: anndata.AnnData) -> None:
    """Test CyteType works without coordinates."""
    # Remove coordinates
    del mock_adata.obsm["X_umap"]

    ct = CyteType(mock_adata, group_key="leiden")

    # Should still initialize
    assert ct.coordinates_key is None
    assert ct.visualization_data["coordinates"] is None
    assert ct.visualization_data["clusters"] == []


def test_cytetype_no_metadata(mock_adata: anndata.AnnData) -> None:
    """Test CyteType works with aggregate_metadata=False."""
    ct = CyteType(mock_adata, group_key="leiden", aggregate_metadata=False)

    assert ct.group_metadata == {}


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_with_custom_prefix(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
) -> None:
    """Test run() with custom results_prefix."""
    mock_submit.return_value = "job_custom"
    mock_wait.return_value = mock_api_response

    ct = CyteType(mock_adata, group_key="leiden")
    result = ct.run(study_context="Test", results_prefix="custom_prefix")

    # Verify custom prefix used
    assert "custom_prefix_annotation_leiden" in result.obs.columns
    assert "custom_prefix_results" in result.uns
    assert "custom_prefix_jobDetails" in result.uns


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_with_api_url_override(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
) -> None:
    """Test run() with api_url override."""
    mock_submit.return_value = "job_override"
    mock_wait.return_value = mock_api_response

    ct = CyteType(mock_adata, group_key="leiden", api_url="https://original.api")

    # Run with different API URL
    ct.run(study_context="Test", api_url="https://override.api")

    # Verify API URL was updated
    assert ct.api_url == "https://override.api"


def test_cytetype_get_results_no_job_details(mock_adata: anndata.AnnData) -> None:
    """Test get_results() returns None when no job details exist."""
    ct = CyteType(mock_adata, group_key="leiden")

    # No results stored, no job details
    results = ct.get_results()

    assert results is None


def test_cytetype_get_results_missing_job_id(mock_adata: anndata.AnnData) -> None:
    """Test get_results() handles missing job_id in job details."""
    ct = CyteType(mock_adata, group_key="leiden")

    # Store job details without job_id
    ct.adata.uns["cytetype_jobDetails"] = {
        "api_url": "https://api.test",
    }

    results = ct.get_results()

    assert results is None


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_with_auth_token_override(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
) -> None:
    """Test run() with auth_token override."""
    mock_submit.return_value = "job_auth"
    mock_wait.return_value = mock_api_response

    ct = CyteType(mock_adata, group_key="leiden", auth_token="token_init")

    # Run with different auth token
    ct.run(study_context="Test", auth_token="token_override")

    # Verify auth token was updated
    assert ct.auth_token == "token_override"


@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_handles_rate_limit_error(
    mock_submit: MagicMock, mock_adata: anndata.AnnData
) -> None:
    """Test run() catches and re-raises RateLimitError with clean message."""
    # Mock API to raise RateLimitError
    mock_submit.side_effect = RateLimitError(
        "Too many requests. Rate limit is 5 jobs per 24 hours. Try again after 3600 seconds.",
        error_code="RATE_LIMIT_EXCEEDED",
    )

    ct = CyteType(mock_adata, group_key="leiden")

    # Should raise RateLimitError (after logging clean message)
    with pytest.raises(RateLimitError) as exc_info:
        ct.run(study_context="Test")

    # Verify exception has correct attributes (server's message preserved)
    assert exc_info.value.error_code == "RATE_LIMIT_EXCEEDED"
    assert "Too many requests" in exc_info.value.message
    assert "3600 seconds" in exc_info.value.message


@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_handles_auth_error(
    mock_submit: MagicMock, mock_adata: anndata.AnnData
) -> None:
    """Test run() catches and re-raises AuthenticationError with clean message."""
    # Mock API to raise AuthenticationError
    mock_submit.side_effect = AuthenticationError(
        "API token is invalid or not found. Re-authenticate with valid token.",
        error_code="INVALID_TOKEN",
    )

    ct = CyteType(mock_adata, group_key="leiden")

    # Should raise AuthenticationError (after logging clean message)
    with pytest.raises(AuthenticationError) as exc_info:
        ct.run(study_context="Test")

    # Verify exception has correct attributes (server's message preserved)
    assert exc_info.value.error_code == "INVALID_TOKEN"
    assert "invalid or not found" in exc_info.value.message


def test_cytetype_initialization_missing_required_key(
    mock_adata: anndata.AnnData,
) -> None:
    """Test CyteType raises clear error when required key is missing."""
    # Remove required key
    del mock_adata.obs["leiden"]

    # Should raise KeyError with helpful message
    with pytest.raises(KeyError, match="not found in `adata.obs`"):
        CyteType(mock_adata, group_key="leiden")


def test_cytetype_initialization_missing_rank_genes(
    mock_adata: anndata.AnnData,
) -> None:
    """Test CyteType raises clear error when rank_genes_groups is missing."""
    # Remove rank_genes_groups
    del mock_adata.uns["rank_genes_groups"]

    # Should raise KeyError with helpful message
    with pytest.raises(KeyError, match="Run `sc.tl.rank_genes_groups` first"):
        CyteType(mock_adata, group_key="leiden")


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_validation_error_llm_config(
    mock_submit: MagicMock, mock_wait: MagicMock, mock_adata: anndata.AnnData
) -> None:
    """Test helpful error when LLMModelConfig validation fails."""
    ct = CyteType(mock_adata, group_key="leiden")

    # Invalid llm_configs (missing required fields)
    invalid_configs = [
        {"provider": "openai"}  # Missing 'name' and credentials
    ]

    # Should raise ValidationError with helpful message
    with pytest.raises(ValidationError):
        ct.run(study_context="Test", llm_configs=invalid_configs)


def test_cytetype_get_results_backwards_compat_dict_format(
    mock_adata: anndata.AnnData, mock_api_response: dict[str, Any]
) -> None:
    """Test get_results() handles old dict format (backwards compat)."""
    ct = CyteType(mock_adata, group_key="leiden")

    # Store results in old format (dict, not JSON string)
    ct.adata.uns["cytetype_results"] = {
        "job_id": "old_job",
        "result": mock_api_response,  # Dict instead of JSON string
    }

    # Should still load successfully
    results = ct.get_results()

    assert results is not None
    assert results == mock_api_response


def test_cytetype_get_results_malformed_json(mock_adata: anndata.AnnData) -> None:
    """Test get_results() handles corrupted JSON gracefully."""
    ct = CyteType(mock_adata, group_key="leiden")

    # Store malformed JSON string
    ct.adata.uns["cytetype_results"] = {
        "job_id": "bad_job",
        "result": "not valid json{[",
    }

    # Should return None (warning logged but no crash)
    results = ct.get_results()

    assert results is None


@patch("cytetype.core.results.get_job_status")
@patch("cytetype.main.fetch_remote_results")
def test_cytetype_get_results_job_status_failed(
    mock_fetch: MagicMock, mock_status: MagicMock, mock_adata: anndata.AnnData
) -> None:
    """Test get_results() when remote job has failed status."""
    mock_fetch.return_value = None  # Failed job returns None

    ct = CyteType(mock_adata, group_key="leiden")

    # Store job details
    ct.adata.uns["cytetype_jobDetails"] = {
        "job_id": "failed_job",
        "api_url": "https://api.test",
    }

    # Get results (job failed)
    results = ct.get_results()

    assert results is None
