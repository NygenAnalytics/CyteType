"""Tests for missing functionality in cytetype.client module."""

import pytest
import requests
from unittest.mock import patch, MagicMock

from cytetype.client import check_job_status, _make_results_request
from cytetype.exceptions import CyteTypeAPIError
from cytetype.config import DEFAULT_API_URL


# --- Test check_job_status function ---

MOCK_JOB_ID = "test-job-456"


@patch("cytetype.client.requests.get")
def test_check_job_status_completed(mock_get: MagicMock) -> None:
    """Test check_job_status returns completed job status."""
    mock_result = {"annotations": [{"clusterId": "1", "annotation": "Type A"}]}
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {"status": "completed", "result": mock_result}
    mock_get.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "completed"
    assert result["result"] == mock_result
    assert result["message"] == "Job completed successfully"
    assert "raw_response" in result

    expected_url = f"{DEFAULT_API_URL}/results/{MOCK_JOB_ID}"
    mock_get.assert_called_once_with(expected_url, headers={}, timeout=30)


@patch("cytetype.client.requests.get")
def test_check_job_status_processing(mock_get: MagicMock) -> None:
    """Test check_job_status returns processing job status."""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {"status": "processing"}
    mock_get.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "processing"
    assert result["result"] is None
    assert result["message"] == "Job is processing"


@patch("cytetype.client.requests.get")
def test_check_job_status_pending(mock_get: MagicMock) -> None:
    """Test check_job_status returns pending job status."""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {"status": "pending"}
    mock_get.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "pending"
    assert result["result"] is None
    assert result["message"] == "Job is pending"


@patch("cytetype.client.requests.get")
def test_check_job_status_error(mock_get: MagicMock) -> None:
    """Test check_job_status returns error job status."""
    error_msg = "Internal processing error"
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {"status": "error", "message": error_msg}
    mock_get.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "error"
    assert result["result"] is None
    assert result["message"] == error_msg


@patch("cytetype.client.requests.get")
def test_check_job_status_not_found(mock_get: MagicMock) -> None:
    """Test check_job_status handles 404 responses."""
    mock_response = MagicMock()
    mock_response.status_code = 404
    mock_get.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "not_found"
    assert result["result"] is None
    assert result["message"] == "Job results not yet available"
    assert result["raw_response"] is None


@patch("cytetype.client.requests.get")
def test_check_job_status_with_auth_token(mock_get: MagicMock) -> None:
    """Test check_job_status with authentication token."""
    auth_token = "test-token-123"
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {"status": "pending"}
    mock_get.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL, auth_token=auth_token)

    assert result["status"] == "pending"
    expected_headers = {"Authorization": f"Bearer {auth_token}"}
    expected_url = f"{DEFAULT_API_URL}/results/{MOCK_JOB_ID}"
    mock_get.assert_called_once_with(expected_url, headers=expected_headers, timeout=30)


@patch("cytetype.client.requests.get")
def test_check_job_status_network_error(mock_get: MagicMock) -> None:
    """Test check_job_status handles network errors."""
    mock_get.side_effect = requests.exceptions.RequestException("Connection failed")

    with pytest.raises(
        CyteTypeAPIError, match="Network error while checking job status"
    ):
        check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)


@patch("cytetype.client.requests.get")
def test_check_job_status_invalid_json(mock_get: MagicMock) -> None:
    """Test check_job_status handles invalid JSON responses."""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.side_effect = requests.exceptions.JSONDecodeError(
        "Invalid JSON", "", 0
    )
    mock_response.raise_for_status.return_value = None
    mock_get.return_value = mock_response

    with pytest.raises(
        CyteTypeAPIError, match="Network error while checking job status"
    ):
        check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)


@patch("cytetype.client.requests.get")
def test_check_job_status_http_error(mock_get: MagicMock) -> None:
    """Test check_job_status handles HTTP errors (non-404)."""
    mock_response = MagicMock()
    mock_response.status_code = 500
    mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
        "500 Server Error"
    )
    mock_get.return_value = mock_response

    with pytest.raises(
        CyteTypeAPIError, match="Network error while checking job status"
    ):
        check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)


@patch("cytetype.client.requests.get")
def test_check_job_status_unknown_status(mock_get: MagicMock) -> None:
    """Test check_job_status handles unknown status values."""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {"status": "unknown_status"}
    mock_get.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "unknown"
    assert result["result"] is None
    assert "Unknown job status: unknown_status" in result["message"]


@patch("cytetype.client.requests.get")
def test_check_job_status_invalid_completed_response(mock_get: MagicMock) -> None:
    """Test check_job_status handles completed job with invalid result format."""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {
        "status": "completed",
        "result": "invalid_result_format",  # Should be dict with 'annotations'
    }
    mock_get.return_value = mock_response

    result = check_job_status(MOCK_JOB_ID, DEFAULT_API_URL)

    assert result["status"] == "error"
    assert result["result"] is None
    assert result["message"] == "Invalid response format from API"


# --- Test _make_results_request function directly ---


def test_make_results_request_completed() -> None:
    """Test _make_results_request directly for completed job."""
    with patch("cytetype.client.requests.get") as mock_get:
        mock_result = {"annotations": [{"clusterId": "1", "annotation": "Type A"}]}
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"status": "completed", "result": mock_result}
        mock_get.return_value = mock_response

        result = _make_results_request(MOCK_JOB_ID, DEFAULT_API_URL)

        assert result["status"] == "completed"
        assert result["result"] == mock_result
        assert result["message"] == "Job completed successfully"


def test_make_results_request_with_auth() -> None:
    """Test _make_results_request with authentication."""
    with patch("cytetype.client.requests.get") as mock_get:
        auth_token = "test-token-456"
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"status": "pending"}
        mock_get.return_value = mock_response

        result = _make_results_request(
            MOCK_JOB_ID, DEFAULT_API_URL, auth_token=auth_token
        )

        assert result["status"] == "pending"
        expected_headers = {"Authorization": f"Bearer {auth_token}"}
        mock_get.assert_called_once_with(
            f"{DEFAULT_API_URL}/results/{MOCK_JOB_ID}",
            headers=expected_headers,
            timeout=30,
        )
