import pytest
import requests
from unittest.mock import patch, MagicMock, call
from typing import Any

from cytetype.client import submit_job, poll_for_results
from cytetype.exceptions import CyteTypeAPIError, CyteTypeTimeoutError, CyteTypeJobError
from cytetype.config import (
    DEFAULT_API_URL,
)

# --- Test submit_annotation_job ---

MOCK_QUERY: dict[str, Any] = {"bioContext": {}, "markerGenes": {}, "expressionData": {}}
MOCK_JOB_ID = "test-job-123"


@patch("cytetype.client.requests.post")
def test_submit_annotation_job_success(mock_post: MagicMock) -> None:
    """Test successful job submission."""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {"job_id": MOCK_JOB_ID}
    mock_post.return_value = mock_response

    job_id = submit_job(MOCK_QUERY, DEFAULT_API_URL)

    assert job_id == MOCK_JOB_ID
    mock_post.assert_called_once_with(
        f"{DEFAULT_API_URL}/annotate",
        json=MOCK_QUERY,
        headers={"Content-Type": "application/json"},
        timeout=60,  # Added timeout to match the actual call
    )


@patch("cytetype.client.requests.post")
def test_submit_annotation_job_api_error(mock_post: MagicMock) -> None:
    """Test job submission failure due to API error (non-200 status)."""
    # Use spec=requests.Response and mock raise_for_status
    mock_response = MagicMock(spec=requests.Response)
    mock_response.status_code = 500
    mock_response.text = "Internal Server Error"
    mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
        "500 Server Error", response=mock_response
    )
    mock_post.return_value = mock_response

    with pytest.raises(CyteTypeAPIError, match="Network error while submitting job"):
        submit_job(MOCK_QUERY, DEFAULT_API_URL)
    mock_post.assert_called_once()


@patch("cytetype.client.requests.post")
def test_submit_annotation_job_connection_error(mock_post: MagicMock) -> None:
    """Test job submission failure due to connection error."""
    mock_post.side_effect = requests.exceptions.RequestException("Connection failed")

    with pytest.raises(CyteTypeAPIError, match="Network error while submitting job"):
        submit_job(MOCK_QUERY, DEFAULT_API_URL)
    mock_post.assert_called_once()


@patch("cytetype.client.requests.post")
def test_submit_annotation_job_with_model_config(mock_post: MagicMock) -> None:
    """Test job submission with custom model config."""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {"job_id": MOCK_JOB_ID}
    mock_post.return_value = mock_response

    model_config = [
        {
            "provider": "openai",
            "modelName": "gpt-4",
            "apiKey": "sk-testkey",
            "baseURL": "http://custom.openai.api",
        }
    ]

    # Create payload with model config included (matching how main.py does it)
    payload_with_model_config = MOCK_QUERY.copy()
    payload_with_model_config["modelConfig"] = model_config

    job_id = submit_job(payload_with_model_config, DEFAULT_API_URL)

    assert job_id == MOCK_JOB_ID
    mock_post.assert_called_once()
    args, kwargs = mock_post.call_args
    assert args[0] == f"{DEFAULT_API_URL}/annotate"
    # Check that the payload includes modelConfig
    assert kwargs["json"] == payload_with_model_config
    assert kwargs["json"]["modelConfig"] == model_config


# --- Test poll_for_results ---

MOCK_RESULT_PAYLOAD = {"annotations": [{"clusterId": "1", "annotation": "Type X"}]}
# Add a mock log response
MOCK_LOG_PAYLOAD = "Log line 1\nLog line 2"


@patch("cytetype.client.time.sleep", return_value=None)  # Prevent actual sleep
@patch("cytetype.client.make_results_request")
def test_poll_for_results_success(
    mock_make_results: MagicMock, mock_sleep: MagicMock
) -> None:
    """Test successful polling completion."""
    # Mock sequence: pending -> completed
    mock_make_results.side_effect = [
        {
            "status": "pending",
            "result": None,
            "message": "Job is pending",
            "raw_response": {"clusterStatus": {}},
        },
        {
            "status": "completed",
            "result": MOCK_RESULT_PAYLOAD,
            "message": "Job completed successfully",
            "raw_response": {"clusterStatus": {}},
        },
    ]

    result = poll_for_results(MOCK_JOB_ID, DEFAULT_API_URL, poll_interval=1, timeout=5)

    assert result == MOCK_RESULT_PAYLOAD
    # Should be called twice - once for pending, once for completed
    assert mock_make_results.call_count == 2
    # Check that sleep was called for initial delay and poll interval (0.5s x 2 for 1s poll interval)
    expected_sleep_calls = [call(5), call(0.5), call(0.5)]
    mock_sleep.assert_has_calls(expected_sleep_calls)


@patch("cytetype.client.time.sleep", return_value=None)
@patch("cytetype.client.make_results_request")
def test_poll_for_results_error_status(
    mock_make_results: MagicMock, mock_sleep: MagicMock
) -> None:
    """Test polling when API returns an 'error' status."""
    mock_make_results.return_value = {
        "status": "failed",
        "result": None,
        "message": "Annotation failed internally",
        "raw_response": {"clusterStatus": {}},
    }

    with pytest.raises(
        CyteTypeJobError, match="Server error: Annotation failed internally"
    ):
        poll_for_results(MOCK_JOB_ID, DEFAULT_API_URL, poll_interval=1, timeout=5)
    mock_make_results.assert_called_once()
    # Assert only the initial sleep(5) happened
    mock_sleep.assert_called_once_with(5)


@patch("cytetype.client.time.sleep", return_value=None)
@patch("cytetype.client.requests.get")
def test_poll_for_results_timeout(mock_get: MagicMock, mock_sleep: MagicMock) -> None:
    """Test polling timeout."""
    mock_response_pending = MagicMock()
    mock_response_pending.status_code = 200
    mock_response_pending.json.return_value = {"status": "pending"}
    mock_get.return_value = mock_response_pending  # Always pending

    poll_interval = 1
    timeout = 2

    with pytest.raises(CyteTypeTimeoutError, match="Timeout while fetching results"):
        poll_for_results(MOCK_JOB_ID, DEFAULT_API_URL, poll_interval, timeout)

    # Check it called sleep multiple times before timeout
    assert mock_sleep.call_count > 1
    assert mock_get.call_count > 1


@patch("cytetype.client.time.sleep", return_value=None)
@patch("cytetype.client.requests.get")
def test_poll_for_results_api_error(mock_get: MagicMock, mock_sleep: MagicMock) -> None:
    """Test polling failure due to API error (non-200 status)."""
    mock_response = MagicMock(spec=requests.Response)
    mock_response.status_code = 503
    mock_response.text = "Service Unavailable"
    # Mock raise_for_status to raise HTTPError
    mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
        "503 Service Unavailable", response=mock_response
    )
    mock_get.return_value = mock_response

    with pytest.raises(CyteTypeTimeoutError, match="Timeout while fetching results"):
        poll_for_results(MOCK_JOB_ID, DEFAULT_API_URL, poll_interval=1, timeout=2)
    # The function will keep retrying until timeout, so call count will be high
    assert mock_get.call_count > 1
    # Assert initial sleep(10) happened
    assert mock_sleep.call_count >= 1


@patch("cytetype.client.time.sleep", return_value=None)
@patch("cytetype.client.requests.get")
def test_poll_for_results_connection_error(
    mock_get: MagicMock, mock_sleep: MagicMock
) -> None:
    """Test polling failure due to connection error."""
    mock_get.side_effect = requests.exceptions.RequestException("Connection failed")

    with pytest.raises(CyteTypeTimeoutError, match="Timeout while fetching results"):
        poll_for_results(MOCK_JOB_ID, DEFAULT_API_URL, poll_interval=1, timeout=2)
    # The function will keep retrying until timeout, so call count will be high
    assert mock_get.call_count > 1
    # Assert initial sleep(10) happened
    assert mock_sleep.call_count >= 1


@patch("cytetype.client.time.sleep", return_value=None)
@patch("cytetype.client.requests.get")
def test_poll_for_results_invalid_json(
    mock_get: MagicMock, mock_sleep: MagicMock
) -> None:
    """Test polling failure due to invalid JSON response."""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.side_effect = requests.exceptions.JSONDecodeError(
        "msg", "doc", 0
    )
    mock_get.return_value = mock_response

    with pytest.raises(CyteTypeTimeoutError, match="Timeout while fetching results"):
        poll_for_results(MOCK_JOB_ID, DEFAULT_API_URL, poll_interval=1, timeout=2)
    # The function will keep retrying until timeout, so call count will be high
    assert mock_get.call_count > 1
    # Assert initial sleep(10) happened
    assert mock_sleep.call_count >= 1


@patch("cytetype.client.time.sleep", return_value=None)
@patch("cytetype.client.make_results_request")
def test_poll_for_results_missing_keys(
    mock_make_results: MagicMock, mock_sleep: MagicMock
) -> None:
    """Test polling failure due to missing keys in 'completed' response."""
    # Status is failed due to invalid response format
    mock_make_results.return_value = {
        "status": "failed",
        "result": None,
        "message": "Invalid response format from API",
        "raw_response": {"clusterStatus": {}},
    }

    with pytest.raises(
        CyteTypeJobError, match="Server error: Invalid response format from API"
    ):
        poll_for_results(MOCK_JOB_ID, DEFAULT_API_URL, poll_interval=1, timeout=5)
    mock_make_results.assert_called_once()
    # Assert only the initial sleep(5) happened
    mock_sleep.assert_called_once_with(5)


# --- Test auth_token functionality ---


@patch("cytetype.client.requests.post")
def test_submit_job_with_auth_token(mock_post: MagicMock) -> None:
    """Test submit_job with auth_token includes Bearer token in headers."""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {"job_id": MOCK_JOB_ID}
    mock_post.return_value = mock_response

    auth_token = "test-bearer-token-123"
    job_id = submit_job(MOCK_QUERY, DEFAULT_API_URL, auth_token=auth_token)

    assert job_id == MOCK_JOB_ID
    mock_post.assert_called_once()
    args, kwargs = mock_post.call_args
    assert args[0] == f"{DEFAULT_API_URL}/annotate"
    assert kwargs["json"] == MOCK_QUERY

    # Check that Authorization header is included
    expected_headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {auth_token}",
    }
    assert kwargs["headers"] == expected_headers


@patch("cytetype.client.time.sleep", return_value=None)
@patch("cytetype.client.make_results_request")
def test_poll_for_results_with_auth_token(
    mock_make_results: MagicMock, mock_sleep: MagicMock
) -> None:
    """Test poll_for_results with auth_token passes it to make_results_request."""
    mock_make_results.return_value = {
        "status": "completed",
        "result": MOCK_RESULT_PAYLOAD,
        "message": "Job completed successfully",
        "raw_response": {"clusterStatus": {}},
    }

    auth_token = "test-bearer-token-456"
    result = poll_for_results(
        MOCK_JOB_ID, DEFAULT_API_URL, poll_interval=1, timeout=5, auth_token=auth_token
    )

    assert result == MOCK_RESULT_PAYLOAD

    # Check that auth_token is passed to make_results_request
    mock_make_results.assert_called_once_with(MOCK_JOB_ID, DEFAULT_API_URL, auth_token)
