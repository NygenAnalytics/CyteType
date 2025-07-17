"""Tests for missing functionality in cytetype.client module."""

import pytest
from typing import Any
from unittest.mock import patch, MagicMock

from cytetype.client import check_job_status
from cytetype.api import make_results_request
from cytetype.exceptions import CyteTypeAPIError
from cytetype.config import DEFAULT_API_URL


# --- Test check_job_status function ---

MOCK_JOB_ID = "test-job-456"


@patch("cytetype.client.make_results_request")
def test_check_job_status_completed(mock_make_results: MagicMock) -> None:
    """Test check_job_status returns completed job status."""
    mock_result = {"annotations": [{"clusterId": "1", "annotation": "Type A"}]}
    mock_response = {
        "status": "completed",
        "result": mock_result,
        "message": "Job completed successfully",
        "raw_response": {},
    }
    mock_make_results.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "completed"
    assert result["result"] == mock_result
    assert result["message"] == "Job completed successfully"
    assert "raw_response" in result

    mock_make_results.assert_called_once_with(MOCK_JOB_ID, DEFAULT_API_URL, None)


@patch("cytetype.client.make_results_request")
def test_check_job_status_processing(mock_make_results: MagicMock) -> None:
    """Test check_job_status returns processing job status."""
    mock_response: dict[str, Any] = {
        "status": "processing",
        "result": None,
        "message": "Job is processing",
        "raw_response": {},
    }
    mock_make_results.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "processing"
    assert result["result"] is None
    assert result["message"] == "Job is processing"


@patch("cytetype.client.make_results_request")
def test_check_job_status_pending(mock_make_results: MagicMock) -> None:
    """Test check_job_status returns pending job status."""
    mock_response: dict[str, Any] = {
        "status": "pending",
        "result": None,
        "message": "Job is pending",
        "raw_response": {},
    }
    mock_make_results.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "pending"
    assert result["result"] is None
    assert result["message"] == "Job is pending"


@patch("cytetype.client.make_results_request")
def test_check_job_status_error(mock_make_results: MagicMock) -> None:
    """Test check_job_status returns error job status."""
    error_msg = "Internal processing error"
    mock_response: dict[str, Any] = {
        "status": "failed",
        "result": None,
        "message": error_msg,
        "raw_response": {},
    }
    mock_make_results.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "failed"
    assert result["result"] is None
    assert result["message"] == error_msg


@patch("cytetype.client.make_results_request")
def test_check_job_status_not_found(mock_make_results: MagicMock) -> None:
    """Test check_job_status handles 404 responses."""
    mock_response = {
        "status": "not_found",
        "result": None,
        "message": "Job not found",
        "raw_response": None,
    }
    mock_make_results.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "not_found"
    assert result["result"] is None
    assert result["message"] == "Job not found"
    assert result["raw_response"] is None


@patch("cytetype.client.make_results_request")
def test_check_job_status_with_auth_token(mock_make_results: MagicMock) -> None:
    """Test check_job_status with authentication token."""
    auth_token = "test-token-123"
    mock_response: dict[str, Any] = {
        "status": "pending",
        "result": None,
        "message": "Job is pending",
        "raw_response": {},
    }
    mock_make_results.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL, auth_token=auth_token)

    assert result["status"] == "pending"
    mock_make_results.assert_called_once_with(MOCK_JOB_ID, DEFAULT_API_URL, auth_token)


@patch("cytetype.client.make_results_request")
def test_check_job_status_network_error(mock_make_results: MagicMock) -> None:
    """Test check_job_status handles network errors."""
    mock_make_results.side_effect = CyteTypeAPIError(
        "Network error while checking job status"
    )

    with pytest.raises(
        CyteTypeAPIError, match="Network error while checking job status"
    ):
        check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)


@patch("cytetype.client.make_results_request")
def test_check_job_status_invalid_json(mock_make_results: MagicMock) -> None:
    """Test check_job_status handles invalid JSON responses."""
    mock_make_results.side_effect = CyteTypeAPIError(
        "Network error while checking job status"
    )

    with pytest.raises(
        CyteTypeAPIError, match="Network error while checking job status"
    ):
        check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)


@patch("cytetype.client.make_results_request")
def test_check_job_status_http_error(mock_make_results: MagicMock) -> None:
    """Test check_job_status handles HTTP errors (non-404)."""
    mock_make_results.side_effect = CyteTypeAPIError(
        "Network error while checking job status"
    )

    with pytest.raises(
        CyteTypeAPIError, match="Network error while checking job status"
    ):
        check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)


@patch("cytetype.client.make_results_request")
def test_check_job_status_unknown_status(mock_make_results: MagicMock) -> None:
    """Test check_job_status handles unknown status values."""
    mock_response: dict[str, Any] = {
        "status": "unknown",
        "result": None,
        "message": "Unknown job status: unknown_status",
        "raw_response": {},
    }
    mock_make_results.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "unknown"
    assert result["result"] is None
    assert "Unknown job status: unknown_status" in result["message"]


@patch("cytetype.client.make_results_request")
def test_check_job_status_invalid_completed_response(
    mock_make_results: MagicMock,
) -> None:
    """Test check_job_status handles completed job with invalid result format."""
    mock_response: dict[str, Any] = {
        "status": "failed",
        "result": None,
        "message": "Invalid response format from API",
        "raw_response": {},
    }
    mock_make_results.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "failed"
    assert result["result"] is None
    assert result["message"] == "Invalid response format from API"


# --- Test make_results_request function directly ---


def test_make_results_request_completed() -> None:
    """Test make_results_request directly for completed job."""
    mock_result = {"annotations": [{"clusterId": "1", "annotation": "Type A"}]}
    with (
        patch("cytetype.api._check_job_status") as mock_check_status,
        patch("cytetype.api._fetch_results") as mock_fetch_results,
        patch("cytetype.api._transform_results") as mock_transform_results,
    ):
        # Mock a completed job
        mock_check_status.return_value = (200, {"jobStatus": "completed"})
        mock_fetch_results.return_value = (200, {"result": mock_result})
        mock_transform_results.return_value = mock_result

        result = make_results_request(MOCK_JOB_ID, DEFAULT_API_URL)

        assert result["status"] == "completed"
        assert result["result"] == mock_result
        assert result["message"] == "Job completed successfully"


def test_make_results_request_with_auth() -> None:
    """Test make_results_request with authentication."""
    with patch("cytetype.api._check_job_status") as mock_check_status:
        auth_token = "test-token-456"

        # Mock a pending job
        mock_check_status.return_value = (200, {"jobStatus": "pending"})

        result = make_results_request(
            MOCK_JOB_ID, DEFAULT_API_URL, auth_token=auth_token
        )

        assert result["status"] == "pending"
        mock_check_status.assert_called_once_with(
            MOCK_JOB_ID, DEFAULT_API_URL, auth_token
        )
