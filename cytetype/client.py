import requests
import time
from typing import Any

from .config import logger
from .exceptions import CyteTypeAPIError, CyteTypeTimeoutError, CyteTypeJobError


def _make_results_request(
    job_id: str,
    api_url: str,
    auth_token: str | None = None,
) -> dict[str, Any]:
    """Make a single request to check job results status.

    This is a shared helper function for both polling and single status checks.

    Args:
        job_id: The job ID to check
        api_url: The API base URL
        auth_token: Bearer token for API authentication

    Returns:
        A dictionary containing:
        - 'status': The job status ('completed', 'processing', 'pending', 'error', 'not_found')
        - 'result': The result data if status is 'completed'
        - 'message': Status message or error message
        - 'raw_response': The raw API response for debugging

    Raises:
        CyteTypeAPIError: For network or API response errors
    """
    results_url = f"{api_url}/results/{job_id}"

    # Prepare headers for authenticated requests
    headers = {}
    if auth_token:
        headers["Authorization"] = f"Bearer {auth_token}"

    try:
        response = requests.get(results_url, headers=headers, timeout=30)

        if response.status_code == 404:
            return {
                "status": "not_found",
                "result": None,
                "message": "Job results not yet available",
                "raw_response": None,
            }

        response.raise_for_status()
        data = response.json()
        status = data.get("status")

        if status == "completed":
            result_data = data.get("result")
            if not isinstance(result_data, dict) or "annotations" not in result_data:
                return {
                    "status": "error",
                    "result": None,
                    "message": "Invalid response format from API",
                    "raw_response": data,
                }
            return {
                "status": "completed",
                "result": result_data,
                "message": "Job completed successfully",
                "raw_response": data,
            }

        elif status == "error":
            error_message = data.get("message", "Unknown error")
            return {
                "status": "error",
                "result": None,
                "message": error_message,
                "raw_response": data,
            }

        elif status in ["processing", "pending"]:
            return {
                "status": status,
                "result": None,
                "message": f"Job is {status}",
                "raw_response": data,
            }

        else:
            return {
                "status": "unknown",
                "result": None,
                "message": f"Unknown job status: {status}",
                "raw_response": data,
            }

    except requests.exceptions.RequestException as e:
        logger.debug(f"Network error during request for job {job_id}: {e}")
        raise CyteTypeAPIError(f"Network error while checking job status: {e}") from e
    except (ValueError, KeyError, requests.exceptions.JSONDecodeError) as e:
        logger.debug(f"Error processing response for job {job_id}: {e}")
        raise CyteTypeAPIError(
            f"Invalid response while checking job status: {e}"
        ) from e


def submit_job(
    payload: dict[str, Any],
    api_url: str,
    auth_token: str | None = None,
) -> str:
    """Submits the job to the API and returns the job ID.

    Args:
        payload: The job payload to submit
        api_url: The API base URL
        auth_token: Bearer token for API authentication

    Returns:
        The job ID returned by the API
    """

    submit_url = f"{api_url}/annotate"
    logger.debug(f"Submitting job to {submit_url}")

    try:
        headers = {"Content-Type": "application/json"}

        # Add bearer token authentication if provided
        if auth_token:
            headers["Authorization"] = f"Bearer {auth_token}"

        response = requests.post(submit_url, json=payload, headers=headers, timeout=60)

        response.raise_for_status()

        job_id = response.json().get("job_id")
        if not job_id:
            raise ValueError("API response did not contain a 'job_id'.")
        logger.debug(f"Job submitted successfully. Job ID: {job_id}")
        return str(job_id)
    except requests.exceptions.Timeout as e:
        raise CyteTypeTimeoutError("Timeout while submitting job") from e
    except requests.exceptions.RequestException as e:
        error_details = ""
        if e.response is not None:
            try:
                error_details = e.response.json()
            except requests.exceptions.JSONDecodeError:
                error_details = e.response.text
        logger.debug(
            f"Network or HTTP error during job submission: {e}. Details: {error_details}"
        )
        raise CyteTypeAPIError("Network error while submitting job") from e
    except (ValueError, KeyError, requests.exceptions.JSONDecodeError) as e:
        logger.debug(f"Error processing submission response: {e}")
        raise CyteTypeAPIError("Invalid response while submitting job") from e


def poll_for_results(
    job_id: str,
    api_url: str,
    poll_interval: int,
    timeout: int,
    auth_token: str | None = None,
) -> dict[str, Any]:
    """Polls the API for results for a given job ID.

    Args:
        job_id: The job ID to poll for results
        api_url: The API base URL
        poll_interval: How often to poll for results (in seconds)
        timeout: Maximum time to wait for results (in seconds)
        auth_token: Bearer token for API authentication

    Returns:
        The result data from the API when the job completes
    """

    logger.info(f"Waiting for results for job ID: {job_id}")

    time.sleep(20)

    logger.info(
        f"View the automatically updating visualization report at: {api_url}/report/{job_id}"
    )

    logs_url = f"{api_url}/logs/{job_id}"
    logger.debug(f"Polling for results for job {job_id}")
    logger.debug(f"Fetching logs for job {job_id} from {logs_url}")
    start_time = time.time()
    last_logs = ""  # Initialize variable to store last fetched logs

    # Prepare headers for authenticated requests (for logs)
    headers = {}
    if auth_token:
        headers["Authorization"] = f"Bearer {auth_token}"

    while True:
        elapsed_time = time.time() - start_time
        if elapsed_time > timeout:
            raise CyteTypeTimeoutError("Timeout while fetching results")

        logger.debug(
            f"Polling attempt for job {job_id}. Elapsed time: {elapsed_time:.1f}s"
        )

        try:
            # Use the shared helper function for the results request
            status_response = _make_results_request(job_id, api_url, auth_token)
            status = status_response["status"]

            if status == "completed":
                logger.info(f"Job {job_id} completed successfully.")
                result = status_response["result"]
                # Ensure we return a proper dict[str, Any] instead of Any
                if not isinstance(result, dict):
                    raise CyteTypeAPIError(
                        f"Expected dict result from API, got {type(result)}"
                    )
                return result

            elif status == "error":
                logger.debug(
                    f"Job {job_id} failed on the server: {status_response['message']}"
                )
                raise CyteTypeJobError(f"Server error: {status_response['message']}")

            elif status in ["processing", "pending"]:
                logger.debug(
                    f"Job {job_id} status: {status}. Checking logs and waiting {poll_interval}s..."
                )
                # Fetch logs (this is specific to polling, not in the shared function)
                try:
                    log_response = requests.get(
                        logs_url, headers=headers, timeout=10
                    )  # Short timeout for logs
                    log_response.raise_for_status()
                    current_logs = log_response.text
                    if current_logs != last_logs:
                        new_log_lines = current_logs[len(last_logs) :].strip()
                        if new_log_lines:
                            for line in new_log_lines.splitlines():
                                logger.info(line)
                        last_logs = current_logs
                except requests.exceptions.RequestException as log_err:
                    logger.warning(f"Could not fetch logs for job {job_id}: {log_err}")

                time.sleep(poll_interval)

            elif status == "not_found":
                logger.debug(
                    f"Results endpoint not ready yet for job {job_id} (404). Waiting {poll_interval}s..."
                )
                time.sleep(poll_interval)

            else:
                logger.warning(
                    f"Job {job_id} has unknown status: '{status}'. Continuing to poll."
                )
                time.sleep(poll_interval)

        except CyteTypeAPIError as e:
            # Handle timeout specifically for polling
            if "Network error" in str(e):
                logger.debug(
                    f"Network error during polling request for {job_id}: {e}. Retrying..."
                )
                time.sleep(min(poll_interval, 5))
            else:
                # Re-raise other API errors
                raise


def check_job_status(
    job_id: str,
    api_url: str,
    auth_token: str | None = None,
) -> dict[str, Any]:
    """Check the status of a job with a single API call (no polling).

    Args:
        job_id: The job ID to check status for
        api_url: The API base URL
        auth_token: Bearer token for API authentication

    Returns:
        A dictionary containing:
        - 'status': The job status ('completed', 'processing', 'pending', 'error', 'not_found')
        - 'result': The result data if status is 'completed'
        - 'message': Error message if status is 'error'
        - 'raw_response': The raw API response for debugging
    """

    logger.debug(f"Checking status for job {job_id}")

    # Use the shared helper function
    return _make_results_request(job_id, api_url, auth_token)
