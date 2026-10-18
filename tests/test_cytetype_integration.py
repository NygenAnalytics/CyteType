"""Integration tests for CyteType class."""

import copy
import inspect
import json
import tempfile
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import anndata
import duckdb
import numpy as np
import pytest
import scanpy as sc
from pydantic import ValidationError

import cytetype.main
from cytetype import CyteType
from cytetype.api import UploadResponse
from cytetype.api.exceptions import (
    AuthenticationError,
    AuthorizationError,
    RateLimitError,
)
from cytetype.config import StoredCredentials, load_credentials, save_credentials
from cytetype.core.artifacts import save_obs_duckdb

# Results stored without raw_annotations get the SDK-only keys as None.
NO_REVIEW_FIELDS = dict.fromkeys(
    ("reviewerConfidence", "isHeterogeneous", "heterogeneityType")
)


def _saved_credentials(api_token: str) -> StoredCredentials:
    return StoredCredentials(
        apiUrl="https://cytetype.nygen.io",
        apiToken=api_token,
        tokenId=f"{api_token}-id",
        userId="user-id",
        email="researcher@university.edu",
    )


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
    assert ct.gene_symbols_column is not None
    assert ct.gene_symbols_column.startswith("__cytetype_gene_symbols")
    assert ct.gene_symbols_column in mock_adata.var.columns
    assert (
        mock_adata.var[ct.gene_symbols_column].tolist()
        == mock_adata.var["gene_symbols"].tolist()
    )

    # Verify data preparation completed
    assert len(ct.clusters) == len(mock_adata)
    assert len(ct.cluster_map) == 3  # 3 unique leiden clusters
    assert ct.marker_genes is not None
    assert len(ct.marker_genes) == 3
    assert ct.expression_percentages is not None
    assert ct.visualization_data is not None
    assert ct.visualization_data["coordinates"] is not None


def test_cytetype_initialization_uses_environment_api_url(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv(
        "CYTETYPE_API_URL",
        "https://dev.cytetype.example/",
    )

    ct = CyteType(mock_adata, group_key="leiden")

    assert ct.api_url == "https://dev.cytetype.example"


def test_cytetype_api_url_argument_overrides_environment(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("CYTETYPE_API_URL", "https://environment.example")

    ct = CyteType(
        mock_adata,
        group_key="leiden",
        api_url="https://explicit.example",
    )

    assert ct.api_url == "https://explicit.example"


@pytest.mark.parametrize(
    "api_url",
    [
        "http://example.com",
        "https://user@example.com",
        "https://example.com/path",
        "ftp://example.com",
    ],
)
def test_cytetype_initialization_rejects_unsafe_api_url(
    mock_adata: anndata.AnnData,
    api_url: str,
) -> None:
    with pytest.raises(ValueError, match="API URL|Non-local"):
        CyteType(mock_adata, group_key="leiden", api_url=api_url)


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
    monkeypatch.setattr(
        "cytetype.main.load_credentials",
        lambda api_url: StoredCredentials(
            apiUrl=api_url,
            apiToken="stored_test_token",
            tokenId="stored-token-id",
            userId="stored-user-id",
            email="stored@university.edu",
        ),
    )


def test_cytetype_materializes_canonical_column_from_composite_source(
    mock_adata: anndata.AnnData, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, Any] = {}
    original_var_names = mock_adata.var_names.tolist()
    mock_adata.var["gene_symbols"] = [
        f"GENE{i}|{var_name}" for i, var_name in enumerate(original_var_names)
    ]

    def _save_vars(*args: Any, **kwargs: Any) -> None:
        captured.update(kwargs)

    monkeypatch.setattr("cytetype.main.save_features_matrix", _save_vars)

    ct = CyteType(mock_adata, group_key="leiden")

    assert ct.gene_symbols_column is not None
    assert ct.gene_symbols_column.startswith("__cytetype_gene_symbols")
    assert mock_adata.var[ct.gene_symbols_column].tolist() == [
        f"GENE{i}" for i in range(mock_adata.n_vars)
    ]
    assert captured["gene_symbols_column"] == ct.gene_symbols_column


def test_cytetype_materializes_canonical_column_from_composite_var_names(
    mock_adata: anndata.AnnData, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, Any] = {}
    original_var_names = mock_adata.var_names.tolist()
    expected_gene_names = [f"GENE{i}" for i in range(mock_adata.n_vars)]

    mock_adata.var = mock_adata.var.drop(columns=["gene_symbols"])
    mock_adata.var_names = [
        f"{var_name}|{gene_name}"
        for var_name, gene_name in zip(original_var_names, expected_gene_names)
    ]
    sc.tl.rank_genes_groups(mock_adata, "leiden", method="t-test")

    def _save_vars(*args: Any, **kwargs: Any) -> None:
        captured.update(kwargs)

    monkeypatch.setattr("cytetype.main.save_features_matrix", _save_vars)

    ct = CyteType(mock_adata, group_key="leiden")

    assert ct.gene_symbols_column is not None
    assert ct.gene_symbols_column.startswith("__cytetype_gene_symbols")
    assert mock_adata.var[ct.gene_symbols_column].tolist() == expected_gene_names
    assert captured["gene_symbols_column"] == ct.gene_symbols_column


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
    assert (
        result_adata.uns["cytetype_jobDetails"]["api_url"]
        == "https://cytetype.nygen.io"
    )


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
def test_cytetype_cleanup_deletes_artifact_files(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test cleanup() deletes the artifact files built during initialization."""
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

    ct = CyteType(
        mock_adata,
        group_key="leiden",
        vars_h5_path=str(vars_path),
        obs_duckdb_path=str(obs_path),
    )
    ct.run(study_context="Test")

    assert vars_path.exists()
    assert obs_path.exists()
    assert ct.gene_symbols_column is not None
    assert ct.gene_symbols_column.startswith("__cytetype_gene_symbols")
    assert ct.gene_symbols_column in ct.adata.var.columns

    ct.cleanup()

    assert not vars_path.exists()
    assert not obs_path.exists()
    assert ct._vars_h5_path is None
    assert ct._obs_duckdb_path is None
    assert "__cytetype_gene_symbols" not in ct.adata.var.columns
    assert not any(
        col.startswith("__cytetype_gene_symbols") for col in ct.adata.var.columns
    )
    assert ct.gene_symbols_column == "gene_symbols"


def test_cytetype_init_failure_rolls_back_temporary_gene_symbols_column(
    mock_adata: anndata.AnnData, monkeypatch: pytest.MonkeyPatch
) -> None:
    original_columns = mock_adata.var.columns.tolist()
    monkeypatch.setattr(
        "cytetype.main.aggregate_expression_percentages",
        MagicMock(side_effect=RuntimeError("aggregation failed")),
    )

    with pytest.raises(RuntimeError, match="aggregation failed"):
        CyteType(mock_adata, group_key="leiden")

    assert mock_adata.var.columns.tolist() == original_columns
    assert not any(
        col.startswith("__cytetype_gene_symbols") for col in mock_adata.var.columns
    )


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
    assert results["annotations"] == [
        {**annotation, **NO_REVIEW_FIELDS}
        for annotation in mock_api_response["annotations"]
    ]
    assert len(results["annotations"]) == 3


@patch("cytetype.main.fetch_remote_results")
def test_cytetype_get_results_remote(
    mock_fetch: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test get_results() fetches from API when not local."""
    mock_fetch.return_value = mock_api_response
    credential_loader = MagicMock(
        return_value=StoredCredentials(
            apiUrl="https://api.test",
            apiToken="target_test_token",
            tokenId="target-token-id",
            userId="target-user-id",
            email="target@university.edu",
        )
    )
    monkeypatch.setattr("cytetype.main.load_credentials", credential_loader)

    ct = CyteType(
        mock_adata,
        group_key="leiden",
        api_url="https://trusted.test",
        auth_token="trusted_test_token",
    )

    # Store job details (simulating previous run)
    ct.adata.uns["cytetype_jobDetails"] = {
        "job_id": "remote_job",
        "api_url": "https://api.test/",
    }

    # Retrieve results (should fetch from API)
    results = ct.get_results()

    assert results is not None
    assert results == mock_api_response
    credential_loader.assert_called_once_with("https://api.test")
    mock_fetch.assert_called_once()
    fetch_args = mock_fetch.call_args.args
    assert fetch_args[2] == "https://api.test"
    assert fetch_args[3] == "target_test_token"
    # The saved key is not kept; the token given to CyteType stays bound to its URL.
    assert ct.auth_token == "trusted_test_token"
    assert ct._auth_token_api_url == "https://trusted.test"


@patch("cytetype.main.fetch_remote_results")
def test_cytetype_get_results_legacy_job_uses_instance_api_url(
    mock_fetch: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mock_fetch.return_value = mock_api_response
    credential_loader = MagicMock()
    monkeypatch.setattr("cytetype.main.load_credentials", credential_loader)
    ct = CyteType(
        mock_adata,
        group_key="leiden",
        api_url="https://trusted.test/",
        auth_token="trusted_test_token",
    )
    ct.adata.uns["cytetype_jobDetails"] = {"job_id": "legacy_job"}

    results = ct.get_results()

    assert results == mock_api_response
    credential_loader.assert_not_called()
    fetch_args = mock_fetch.call_args.args
    assert fetch_args[2] == "https://trusted.test"
    assert fetch_args[3] == "trusted_test_token"


def test_cytetype_get_results_rejects_modified_job_url_before_network(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    credential_loader = MagicMock(return_value=None)
    upload_obs = MagicMock()
    upload_vars = MagicMock()
    submit = MagicMock()
    wait = MagicMock()
    get_status = MagicMock()
    fetch_results = MagicMock()
    monkeypatch.setattr("cytetype.main.load_credentials", credential_loader)
    monkeypatch.setattr("cytetype.main.upload_obs_duckdb_file", upload_obs)
    monkeypatch.setattr("cytetype.main.upload_vars_h5_file", upload_vars)
    monkeypatch.setattr("cytetype.main.submit_annotation_job", submit)
    monkeypatch.setattr("cytetype.main.wait_for_completion", wait)
    monkeypatch.setattr("cytetype.core.results.get_job_status", get_status)
    monkeypatch.setattr("cytetype.core.results.fetch_job_results", fetch_results)
    ct = CyteType(
        mock_adata,
        group_key="leiden",
        api_url="https://trusted.test",
        auth_token="trusted_test_token",
    )
    ct.adata.uns["cytetype_jobDetails"] = {
        "job_id": "trusted_job",
        "api_url": "https://modified.test",
    }

    with pytest.raises(AuthenticationError, match="cytetype setup"):
        ct.get_results()

    credential_loader.assert_called_once_with("https://modified.test")
    upload_obs.assert_not_called()
    upload_vars.assert_not_called()
    submit.assert_not_called()
    wait.assert_not_called()
    get_status.assert_not_called()
    fetch_results.assert_not_called()


def test_cytetype_get_results_rejects_invalid_job_url_before_credentials(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    credential_loader = MagicMock()
    get_status = MagicMock()
    fetch_results = MagicMock()
    monkeypatch.setattr("cytetype.main.load_credentials", credential_loader)
    monkeypatch.setattr("cytetype.core.results.get_job_status", get_status)
    monkeypatch.setattr("cytetype.core.results.fetch_job_results", fetch_results)
    ct = CyteType(
        mock_adata,
        group_key="leiden",
        api_url="https://trusted.test",
        auth_token="trusted_test_token",
    )
    ct.adata.uns["cytetype_jobDetails"] = {
        "job_id": "trusted_job",
        "api_url": "https://modified.test/path",
    }

    with pytest.raises(ValueError, match="server origin"):
        ct.get_results()

    credential_loader.assert_not_called()
    get_status.assert_not_called()
    fetch_results.assert_not_called()


def test_cytetype_initialization_with_auth_token(mock_adata: anndata.AnnData) -> None:
    """Test CyteType can store auth_token in __init__."""
    ct = CyteType(mock_adata, group_key="leiden", auth_token="test_token_123")

    assert ct.auth_token == "test_token_123"
    assert ct.api_url == "https://cytetype.nygen.io"
    assert ct._auth_token_api_url == "https://cytetype.nygen.io"


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


def test_cytetype_run_rejects_unsafe_api_url_before_network(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    credential_loader = MagicMock()
    upload_obs = MagicMock()
    upload_vars = MagicMock()
    submit = MagicMock()
    wait = MagicMock()
    monkeypatch.setattr("cytetype.main.load_credentials", credential_loader)
    monkeypatch.setattr("cytetype.main.upload_obs_duckdb_file", upload_obs)
    monkeypatch.setattr("cytetype.main.upload_vars_h5_file", upload_vars)
    monkeypatch.setattr("cytetype.main.submit_annotation_job", submit)
    monkeypatch.setattr("cytetype.main.wait_for_completion", wait)
    ct = CyteType(mock_adata, group_key="leiden")
    # __init__ only checks for saved credentials to warn early.
    credential_loader.assert_called_once_with("https://cytetype.nygen.io")
    credential_loader.reset_mock()

    with pytest.raises(ValueError, match="server origin"):
        ct.run(
            study_context="Test",
            api_url="https://other.example/path",
        )

    credential_loader.assert_not_called()
    upload_obs.assert_not_called()
    upload_vars.assert_not_called()
    submit.assert_not_called()
    wait.assert_not_called()


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
    ct.run(
        study_context="Test",
        api_url="https://override.test/",
        auth_token="token_override",
    )

    # Verify auth token was updated
    assert ct.auth_token == "token_override"
    assert ct.api_url == "https://override.test"
    assert ct._auth_token_api_url == "https://override.test"
    assert mock_submit.call_args.args[1] == "token_override"


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_uses_stored_credentials(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
) -> None:
    mock_submit.return_value = "job_stored_auth"
    mock_wait.return_value = mock_api_response
    ct = CyteType(mock_adata, group_key="leiden")

    ct.run(study_context="Test")

    # Saved credentials are read for each call and not kept on the instance.
    assert ct.auth_token is None
    assert mock_submit.call_args.args[1] == "stored_test_token"
    assert mock_wait.call_args.args[1] == "stored_test_token"


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_stored_credentials_are_not_reused_for_another_server(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def load_for_server(api_url: str) -> StoredCredentials | None:
        if api_url != "https://cytetype.nygen.io":
            return None
        return StoredCredentials(
            apiUrl=api_url,
            apiToken="stored_test_token",
            tokenId="stored-token-id",
            userId="stored-user-id",
            email="stored@university.edu",
        )

    monkeypatch.setattr("cytetype.main.load_credentials", load_for_server)
    mock_submit.return_value = "job_stored_auth"
    mock_wait.return_value = mock_api_response
    ct = CyteType(mock_adata, group_key="leiden")
    ct.run(study_context="Test")

    with pytest.raises(AuthenticationError, match="cytetype setup"):
        ct.run(
            study_context="Test other server",
            api_url="https://other.example",
            results_prefix="other",
        )

    assert mock_submit.call_count == 1


def test_cytetype_run_requires_setup_before_upload(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    upload_obs = MagicMock()
    upload_vars = MagicMock()
    monkeypatch.setattr("cytetype.main.load_credentials", lambda api_url: None)
    monkeypatch.setattr("cytetype.main.upload_obs_duckdb_file", upload_obs)
    monkeypatch.setattr("cytetype.main.upload_vars_h5_file", upload_vars)
    ct = CyteType(mock_adata, group_key="leiden")

    with pytest.raises(AuthenticationError, match="cytetype setup") as error_info:
        ct.run(study_context="Test")

    assert error_info.value.error_code == "AUTHENTICATION_REQUIRED"
    upload_obs.assert_not_called()
    upload_vars.assert_not_called()


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
    before = copy.deepcopy(mock_api_response)

    # Store results in old format (dict, not JSON string)
    ct.adata.uns["cytetype_results"] = {
        "job_id": "old_job",
        "result": mock_api_response,  # Dict instead of JSON string
    }

    # Should still load successfully
    results = ct.get_results()

    assert results is not None
    assert results == {
        **before,
        "annotations": [
            {**annotation, **NO_REVIEW_FIELDS} for annotation in before["annotations"]
        ],
    }
    assert ct.adata.uns["cytetype_results"]["result"] == before


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


@pytest.mark.parametrize("require_artifacts", [True, False])
def test_cytetype_run_raises_upload_authentication_errors_with_setup_hint(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
    log_messages: list[tuple[str, str]],
    require_artifacts: bool,
) -> None:
    upload_error = AuthenticationError(
        "API token is invalid or not found. Re-authenticate with valid token.",
        error_code="INVALID_TOKEN",
    )
    submit = MagicMock()
    monkeypatch.setattr(
        "cytetype.main.upload_vars_h5_file", MagicMock(side_effect=upload_error)
    )
    monkeypatch.setattr("cytetype.main.submit_annotation_job", submit)
    ct = CyteType(mock_adata, group_key="leiden")

    with pytest.raises(
        AuthenticationError, match="cytetype setup --force"
    ) as error_info:
        ct.run(study_context="Test", require_artifacts=require_artifacts)

    assert "cytetype login" in error_info.value.message
    assert error_info.value.error_code == "INVALID_TOKEN"
    assert "API token is invalid or not found" in error_info.value.message
    assert error_info.value.__cause__ is upload_error
    submit.assert_not_called()
    # Not reported as an artifact bug.
    assert not any("github" in message.lower() for _, message in log_messages)


def test_cytetype_run_raises_upload_authorization_errors_directly(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
    log_messages: list[tuple[str, str]],
) -> None:
    upload_error = AuthorizationError("License expired", error_code="LICENSE_EXPIRED")
    submit = MagicMock()
    monkeypatch.setattr(
        "cytetype.main.upload_obs_duckdb_file", MagicMock(side_effect=upload_error)
    )
    monkeypatch.setattr("cytetype.main.submit_annotation_job", submit)
    ct = CyteType(mock_adata, group_key="leiden")

    with pytest.raises(AuthorizationError, match="License expired") as error_info:
        ct.run(study_context="Test")

    assert error_info.value.error_code == "LICENSE_EXPIRED"
    assert "https://cytetype.nygen.io/dashboard" in error_info.value.message
    assert error_info.value.__cause__ is upload_error
    submit.assert_not_called()
    assert not any("github" in message.lower() for _, message in log_messages)


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_sends_visibility(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
) -> None:
    mock_submit.return_value = "job_public"
    mock_wait.return_value = mock_api_response
    ct = CyteType(mock_adata, group_key="leiden")

    ct.run(study_context="Test", visibility="public")

    payload = mock_submit.call_args.args[2]
    assert payload["visibility"] == "public"
    assert "visibility" not in payload["input_data"]
    assert ct.adata.uns["cytetype_jobDetails"]["visibility"] == "public"


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_omits_visibility_by_default(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
    tmp_path: Path,
) -> None:
    mock_submit.return_value = "job_default_visibility"
    mock_wait.return_value = mock_api_response
    ct = CyteType(mock_adata, group_key="leiden")

    ct.run(study_context="Test")

    payload = mock_submit.call_args.args[2]
    assert "visibility" not in payload
    assert "visibility" not in ct.adata.uns["cytetype_jobDetails"]
    ct.adata.write_h5ad(tmp_path / "annotated.h5ad")


def test_cytetype_run_rejects_invalid_visibility_before_upload(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    upload_obs = MagicMock()
    upload_vars = MagicMock()
    submit = MagicMock()
    monkeypatch.setattr("cytetype.main.upload_obs_duckdb_file", upload_obs)
    monkeypatch.setattr("cytetype.main.upload_vars_h5_file", upload_vars)
    monkeypatch.setattr("cytetype.main.submit_annotation_job", submit)
    ct = CyteType(mock_adata, group_key="leiden")

    with pytest.raises(ValueError, match="visibility"):
        ct.run(study_context="Test", visibility="everyone")  # type: ignore[arg-type]

    upload_obs.assert_not_called()
    upload_vars.assert_not_called()
    submit.assert_not_called()
    assert "cytetype_jobDetails" not in ct.adata.uns


def test_visibility_is_the_last_run_parameter() -> None:
    parameters = inspect.signature(CyteType.run).parameters

    assert list(parameters)[-1] == "visibility"
    assert parameters["visibility"].default is None


def test_cytetype_init_warns_when_saved_credentials_file_is_invalid(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
    log_messages: list[tuple[str, str]],
) -> None:
    load_error = ValueError("Invalid CyteType credentials file: /home/u/c.json")
    monkeypatch.setattr(
        "cytetype.main.load_credentials", MagicMock(side_effect=load_error)
    )

    ct = CyteType(mock_adata, group_key="leiden")

    warnings = [message for level, message in log_messages if level == "WARNING"]
    assert any(
        str(load_error) in m
        and "Unless you pass auth_token to run()" in m
        and "`cytetype logout`, then `cytetype setup`" in m
        for m in warnings
    )
    assert ct.auth_token is None


@pytest.mark.parametrize(
    "load_error",
    [
        PermissionError("Permission denied"),
        RuntimeError("Could not determine home directory."),
    ],
)
def test_cytetype_init_warns_when_saved_credentials_cannot_be_read(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
    log_messages: list[tuple[str, str]],
    load_error: Exception,
) -> None:
    monkeypatch.setattr(
        "cytetype.main.load_credentials", MagicMock(side_effect=load_error)
    )

    ct = CyteType(mock_adata, group_key="leiden")

    warnings = [message for level, message in log_messages if level == "WARNING"]
    assert (
        f"Saved CyteType credentials could not be read ({load_error}). Pass "
        "`auth_token`, or make sure XDG_CONFIG_HOME or HOME points to a readable "
        "configuration directory."
    ) in warnings
    assert not any("logout" in message for message in warnings)
    assert ct.auth_token is None


@pytest.mark.parametrize(
    "load_error",
    [
        PermissionError("Permission denied"),
        RuntimeError("Could not determine home directory."),
    ],
)
def test_cytetype_run_reports_unreadable_credentials_as_authentication_error(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
    load_error: Exception,
) -> None:
    upload_vars = MagicMock()
    monkeypatch.setattr(
        "cytetype.main.load_credentials", MagicMock(side_effect=load_error)
    )
    monkeypatch.setattr("cytetype.main.upload_vars_h5_file", upload_vars)
    ct = CyteType(mock_adata, group_key="leiden")

    with pytest.raises(AuthenticationError, match="could not be read") as run_error:
        ct.run(study_context="Test")
    ct.adata.uns["cytetype_jobDetails"] = {"job_id": "job-1"}
    with pytest.raises(AuthenticationError, match="XDG_CONFIG_HOME") as fetch_error:
        ct.get_results()

    for error_info in (run_error, fetch_error):
        assert error_info.value.error_code == "AUTHENTICATION_REQUIRED"
        assert error_info.value.__cause__ is load_error
    upload_vars.assert_not_called()


def test_cytetype_init_reports_missing_sign_in_before_preparation(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
    log_messages: list[tuple[str, str]],
) -> None:
    calls: list[str] = []
    resolve_gene_symbols_column = cytetype.main.resolve_gene_symbols_column

    def _load_credentials(api_url: str) -> None:
        calls.append(f"credentials {api_url}")

    def _resolve(*args: Any, **kwargs: Any) -> str | None:
        calls.append("preparation")
        return resolve_gene_symbols_column(*args, **kwargs)

    monkeypatch.setattr("cytetype.main.load_credentials", _load_credentials)
    monkeypatch.setattr("cytetype.main.resolve_gene_symbols_column", _resolve)

    ct = CyteType(mock_adata, group_key="leiden")

    assert calls == ["credentials https://cytetype.nygen.io", "preparation"]
    # Not a warning: passing auth_token only to run() (as in CI) is fine.
    assert (
        "INFO",
        "No saved CyteType sign-in found for https://cytetype.nygen.io. Unless you "
        "pass auth_token to run(), run `cytetype setup` (or `cytetype login`) before "
        "calling run().",
    ) in log_messages
    assert not any(
        level == "WARNING" and "cytetype setup" in message
        for level, message in log_messages
    )
    assert ct.auth_token is None
    assert ct._auth_token_api_url is None


def test_cytetype_init_with_saved_credentials_does_not_warn_or_keep_token(
    mock_adata: anndata.AnnData,
    log_messages: list[tuple[str, str]],
) -> None:
    ct = CyteType(mock_adata, group_key="leiden")

    assert not any("cytetype setup" in message for _, message in log_messages)
    assert ct.auth_token is None
    assert ct._auth_token_api_url is None


def test_cytetype_init_with_auth_token_skips_credentials_check(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    credential_loader = MagicMock()
    monkeypatch.setattr("cytetype.main.load_credentials", credential_loader)

    CyteType(mock_adata, group_key="leiden", auth_token="explicit_token")

    credential_loader.assert_not_called()


def test_cytetype_get_results_renormalizes_stored_results_without_obs_backfill(
    mock_adata: anndata.AnnData,
) -> None:
    ct = CyteType(mock_adata, group_key="leiden")
    stored = {
        "annotations": [
            {
                "clusterId": "1",
                "annotation": "T Cell",
                "ontologyTerm": "T cell",
                "ontologyTermID": "CL:0000084",
                "granularAnnotation": "",
                "cellState": "",
                "justification": "",
                "supportingMarkers": [],
                "conflictingMarkers": [],
                "missingExpression": "",
                "unexpectedExpression": "",
            }
        ],
        "summary": {},
        "clusterCategories": [],
        "studyContext": "",
        "raw_annotations": {
            "1": {
                "latest": {
                    "annotation": {
                        "clusterId": "1",
                        "annotation": "T Cell",
                        "fullOutput": {
                            "cellType": {
                                "label": "T Cell",
                                "keySupportingGenes": ["CD3D", "CD3E"],
                                "conclusion": "CD3 genes are expressed.",
                            }
                        },
                    },
                    "review": {"confidence": "High", "isHeterogeneous": False},
                },
                "history": [],
            }
        },
    }
    ct.adata.uns["cytetype_results"] = {
        "job_id": "stored_job",
        "result": json.dumps(stored),
    }
    obs_columns = list(ct.adata.obs.columns)

    results = ct.get_results()

    assert results is not None
    annotation = results["annotations"][0]
    assert annotation["supportingMarkers"] == ["CD3D", "CD3E"]
    assert annotation["justification"] == "CD3 genes are expressed."
    assert annotation["reviewerConfidence"] == "High"
    assert annotation["isHeterogeneous"] is False
    assert annotation["heterogeneityType"] is None
    assert json.loads(ct.adata.uns["cytetype_results"]["result"]) == stored
    assert list(ct.adata.obs.columns) == obs_columns


@patch("cytetype.main.wait_for_completion")
def test_cytetype_run_reads_replaced_saved_key_after_rejection(
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    monkeypatch.setattr("cytetype.main.load_credentials", load_credentials)
    mock_wait.return_value = mock_api_response
    submitted_tokens: list[str] = []

    def _upload_vars(api_url: str, token: str, *args: Any, **kwargs: Any) -> Any:
        if token == "revoked_key":
            raise AuthenticationError(
                "API token is invalid or not found.", error_code="INVALID_TOKEN"
            )
        return UploadResponse(
            upload_id="vars_upload_456",
            file_kind="vars_h5",
            file_name="vars.h5",
            size_bytes=100,
            expires_at="2099-01-01T00:00:00Z",
        )

    def _submit(api_url: str, token: str, payload: dict[str, Any]) -> str:
        submitted_tokens.append(token)
        return "job_new_key"

    monkeypatch.setattr("cytetype.main.upload_vars_h5_file", _upload_vars)
    monkeypatch.setattr("cytetype.main.submit_annotation_job", _submit)
    save_credentials(_saved_credentials("revoked_key"))
    ct = CyteType(mock_adata, group_key="leiden")

    with pytest.raises(AuthenticationError, match="cytetype setup --force"):
        ct.run(study_context="Test")
    assert submitted_tokens == []

    # What `cytetype setup --force` does: save a new key for the same server.
    save_credentials(_saved_credentials("new_key"))
    ct.run(study_context="Test")

    assert submitted_tokens == ["new_key"]
    assert mock_wait.call_args.args[1] == "new_key"
    assert ct.auth_token is None


def test_cytetype_get_results_reads_saved_credentials_each_time(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    monkeypatch.setattr("cytetype.main.load_credentials", load_credentials)
    fetch = MagicMock(return_value=None)
    monkeypatch.setattr("cytetype.main.fetch_remote_results", fetch)
    save_credentials(_saved_credentials("first_key"))
    ct = CyteType(mock_adata, group_key="leiden")
    ct.adata.uns["cytetype_jobDetails"] = {
        "job_id": "job-1",
        "api_url": "https://cytetype.nygen.io",
    }

    ct.get_results()
    save_credentials(_saved_credentials("second_key"))
    ct.get_results()

    assert [call.args[3] for call in fetch.call_args_list] == [
        "first_key",
        "second_key",
    ]


def test_cytetype_obs_artifact_builds_when_temp_directory_cannot_be_created(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    def _no_temp_directory(*args: Any, **kwargs: Any) -> Any:
        raise OSError(28, "No space left on device")

    monkeypatch.setattr("cytetype.main.save_obs_duckdb_file", save_obs_duckdb)
    monkeypatch.setattr(tempfile, "TemporaryDirectory", _no_temp_directory)
    obs_path = tmp_path / "obs.duckdb"
    obs_columns = list(mock_adata.obs.columns)

    ct = CyteType(mock_adata, group_key="leiden", obs_duckdb_path=str(obs_path))

    assert ct._artifact_build_errors == []
    assert ct._obs_duckdb_path == str(obs_path)
    assert list(mock_adata.obs.columns) == obs_columns
    with duckdb.connect(str(obs_path), read_only=True) as con:
        columns = [row[0] for row in con.execute("DESCRIBE obs").fetchall()]
    assert "__vis_coordinates_X_umap_1" in columns


def test_cytetype_failed_obs_artifact_leaves_no_columns_or_files(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    def _connect(database: str, *args: Any, **kwargs: Any) -> Any:
        Path(database).write_bytes(b"partial")
        Path(f"{database}.wal").write_bytes(b"partial")
        raise duckdb.IOException("disk full")

    monkeypatch.setattr("cytetype.main.save_obs_duckdb_file", save_obs_duckdb)
    monkeypatch.setattr("cytetype.core.artifacts.duckdb.connect", _connect)
    obs_path = tmp_path / "obs.duckdb"
    obs_columns = list(mock_adata.obs.columns)

    ct = CyteType(mock_adata, group_key="leiden", obs_duckdb_path=str(obs_path))

    assert [name for name, _ in ct._artifact_build_errors] == ["obs_duckdb"]
    assert ct._obs_duckdb_path is None
    assert list(mock_adata.obs.columns) == obs_columns
    assert list(tmp_path.iterdir()) == []


def test_cytetype_failed_vars_artifact_is_deleted(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    def _write_partial_vars(*args: Any, **kwargs: Any) -> None:
        Path(kwargs["out_file"]).write_bytes(b"partial")
        raise RuntimeError("disk full")

    monkeypatch.setattr("cytetype.main.save_features_matrix", _write_partial_vars)
    vars_path = tmp_path / "vars.h5"

    ct = CyteType(mock_adata, group_key="leiden", vars_h5_path=str(vars_path))

    assert [name for name, _ in ct._artifact_build_errors] == ["vars_h5"]
    assert ct._vars_h5_path is None
    assert not vars_path.exists()


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_after_cleanup_requires_a_new_object(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mock_submit.return_value = "job_first"
    mock_wait.return_value = mock_api_response
    ct = CyteType(mock_adata, group_key="leiden")
    ct.run(study_context="First run")
    ct.cleanup()
    credential_loader = MagicMock()
    upload_vars = MagicMock()
    upload_obs = MagicMock()
    monkeypatch.setattr("cytetype.main.load_credentials", credential_loader)
    monkeypatch.setattr("cytetype.main.upload_vars_h5_file", upload_vars)
    monkeypatch.setattr("cytetype.main.upload_obs_duckdb_file", upload_obs)
    mock_submit.reset_mock()

    with pytest.raises(ValueError, match="deleted by cleanup"):
        ct.run(study_context="Second run", override_existing_results=True)

    credential_loader.assert_not_called()
    upload_vars.assert_not_called()
    upload_obs.assert_not_called()
    mock_submit.assert_not_called()


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_after_cleanup_can_continue_without_artifacts(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
    log_messages: list[tuple[str, str]],
) -> None:
    mock_submit.return_value = "job_first"
    mock_wait.return_value = mock_api_response
    ct = CyteType(mock_adata, group_key="leiden")
    ct.run(study_context="First run")
    ct.cleanup()
    mock_submit.return_value = "job_second"

    ct.run(
        study_context="Second run",
        override_existing_results=True,
        require_artifacts=False,
    )

    assert "uploaded_files" not in mock_submit.call_args.args[2]
    assert (
        "WARNING",
        "The artifact files were deleted by cleanup(). Continuing without them.",
    ) in log_messages
    assert ct.adata.uns["cytetype_jobDetails"]["job_id"] == "job_second"


@patch("cytetype.main.wait_for_completion")
@patch("cytetype.main.submit_annotation_job")
def test_cytetype_run_returns_the_relabeled_object_with_label_na(
    mock_submit: MagicMock,
    mock_wait: MagicMock,
    mock_adata: anndata.AnnData,
    mock_api_response: dict[str, Any],
) -> None:
    mock_submit.return_value = "job_label_na"
    mock_wait.return_value = mock_api_response
    mock_adata.obs.loc[mock_adata.obs.index[:5], "leiden"] = np.nan

    ct = CyteType(mock_adata, group_key="leiden", label_na=True)
    result = ct.run(study_context="Test")

    # Results go to the returned object, a new one when NaN labels were relabeled.
    assert result is ct.adata
    assert result is not mock_adata
    assert "cytetype_annotation_leiden" in result.obs.columns
    assert "cytetype_annotation_leiden" not in mock_adata.obs.columns
    assert int(mock_adata.obs["leiden"].isna().sum()) == 5
