import requests
import time
from typing import Any

from .config import logger
from .exceptions import CyteTypeAPIError, CyteTypeTimeoutError, CyteTypeJobError


def submit_annotation_job(
    query: dict[str, Any],
    api_url: str,
    model_config: list[dict[str, Any]] | None = None,
) -> str:
    """Submits the annotation job to the API and returns the job ID."""

    submit_url = f"{api_url}/annotate"
    logger.debug(f"Submitting annotation job to {submit_url}")
    try:
        headers = {"Content-Type": "application/json"}

        payload = (
            query.copy()
        )  # Create a copy to avoid modifying the original query dict
        if model_config:
            payload["modelConfig"] = model_config

        response = requests.post(submit_url, json=payload, headers=headers, timeout=60)

        response.raise_for_status()

        job_id = response.json().get("job_id")
        if not job_id:
            raise ValueError("API response did not contain a 'job_id'.")
        logger.debug(f"Annotation job submitted successfully. Job ID: {job_id}")
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
    job_id: str, api_url: str, poll_interval: int, timeout: int
) -> dict[str, Any]:
    """Polls the API for results for a given job ID."""

    time.sleep(10)  # Initial delay before first poll

    retrieve_url = f"{api_url}/retrieve/{job_id}"
    logs_url = f"{api_url}/display_logs/{job_id}"  # Added log URL
    logger.debug(f"Polling for results for job {job_id} at {retrieve_url}")
    logger.debug(f"Fetching logs for job {job_id} from {logs_url}")  # Added log message
    start_time = time.time()
    last_logs = ""  # Initialize variable to store last fetched logs

    while True:
        elapsed_time = time.time() - start_time
        if elapsed_time > timeout:
            raise CyteTypeTimeoutError("Timeout while fetching results")

        logger.debug(
            f"Polling attempt for job {job_id}. Elapsed time: {elapsed_time:.1f}s"
        )
        try:
            response = requests.get(retrieve_url, timeout=30)
            response.raise_for_status()
            data = response.json()
            status = data.get("status")

            if status == "completed":
                logger.info(f"Annotation job {job_id} completed successfully.")
                result_data = data.get("result")
                if (
                    not isinstance(result_data, dict)
                    or "annotations" not in result_data
                ):
                    raise CyteTypeAPIError("Invalid response while parsing results")
                return result_data
            elif status == "error":
                error_message = data.get("message", "Unknown error")
                logger.debug(
                    f"Annotation job {job_id} failed on the server: {error_message}"
                )
                raise CyteTypeJobError("Server error: job failed")
            elif status in ["processing", "pending"]:
                logger.debug(
                    f"Job {job_id} status: {status}. Checking logs and waiting {poll_interval}s..."
                )
                try:
                    log_response = requests.get(
                        logs_url, timeout=10
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
            else:
                logger.warning(
                    f"Job {job_id} has unknown status: '{status}'. Continuing to poll."
                )
                time.sleep(poll_interval)

        except requests.exceptions.Timeout as e:
            logger.debug(
                f"Timeout during a polling request for {job_id}: {e}. Retrying..."
            )
            time.sleep(min(poll_interval, 5))

        except (ValueError, KeyError, requests.exceptions.JSONDecodeError) as e:
            # Catch JSON/Key/Value errors before broad RequestException
            raise CyteTypeAPIError("Invalid response while fetching results") from e

        except requests.exceptions.RequestException as e:
            error_details = "No response details available"
            if e.response is not None:
                try:
                    error_details = e.response.json()
                except requests.exceptions.JSONDecodeError:
                    error_details = e.response.text
            logger.debug(
                f"Network error during polling request for job {job_id}: {e}. Details: {error_details}"
            )
            raise CyteTypeAPIError("Network error while fetching results") from e
