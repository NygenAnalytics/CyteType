import requests
import time
from typing import Any

from .config import logger
from .exceptions import CyteTypeAPIError, CyteTypeTimeoutError, CyteTypeJobError
from .display import display_cluster_status
from .api import make_results_request


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

            # Check for authentication/authorization errors
            if e.response.status_code == 401:
                logger.error("❌ Authentication failed: Invalid or expired auth token")
                raise CyteTypeAPIError(
                    "Authentication failed during job submission"
                ) from e
            elif e.response.status_code == 403:
                logger.error("❌ Authorization failed: Access denied")
                logger.error("Your auth token doesn't have permission to submit jobs.")
                raise CyteTypeAPIError(
                    "Authorization failed during job submission"
                ) from e
            elif e.response.status_code == 422:
                if auth_token and "auth" in str(error_details).lower():
                    logger.error(
                        "❌ Authentication may have failed (server returned validation error)"
                    )
                    raise CyteTypeAPIError(
                        "Possible authentication failure during job submission"
                    ) from e
                else:
                    logger.error(f"❌ Validation error from server: {error_details}")
                    raise CyteTypeAPIError(
                        f"Validation error during job submission: {error_details}"
                    ) from e
            elif e.response.status_code == 429:
                logger.error("❌ Rate limit exceeded")
                raise CyteTypeAPIError(
                    "Rate limit exceeded. Rate limit is 5 annotation jobs every 24hrs."
                ) from e

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
    show_progress: bool = True,
) -> dict[str, Any]:
    """Polls the API for results for a given job ID.

    Args:
        job_id: The job ID to poll for results
        api_url: The API base URL
        poll_interval: How often to poll for results (in seconds)
        timeout: Maximum time to wait for results (in seconds)
        auth_token: Bearer token for API authentication
        show_progress: Whether to display progress updates (spinner and cluster status)

    Returns:
        The result data from the API when the job completes
    """

    logger.info(f"CyteType job (id: {job_id}) submitted. Polling for results...")

    time.sleep(5)

    # Create report URL with auth token if available
    report_url = f"{api_url}/report/{job_id}"
    if auth_token:
        logger.info(
            f"Token secured report (updates automatically) available at: {report_url}"
        )
    else:
        logger.info(f"Report (updates automatically) available at: {report_url}")

    logger.info(
        "If network disconnects, the results can still be fetched:\n`results = annotator.get_results()`"
    )

    start_time = time.time()
    last_cluster_status: dict[str, str] = {}
    spinner_frame = 0
    consecutive_not_found_count = 0  # Track consecutive 404 responses

    while True:
        elapsed_time = time.time() - start_time
        if elapsed_time > timeout:
            print()  # Add newline to clean up progress display
            raise CyteTypeTimeoutError("Timeout while fetching results")

        try:
            # Use the shared helper function for the results request
            status_response = make_results_request(job_id, api_url, auth_token)
            status = status_response["status"]

            # Reset consecutive not found count when we get a valid response
            if status != "not_found":
                consecutive_not_found_count = 0

            # Extract cluster status for all cases
            raw_response = status_response.get("raw_response", {})
            current_cluster_status = raw_response.get("clusterStatus", {})

            if status == "completed":
                # Show final cluster status if it wasn't already shown
                if show_progress:
                    if (
                        current_cluster_status
                        and current_cluster_status != last_cluster_status
                    ):
                        display_cluster_status(
                            current_cluster_status,
                            job_id,
                            is_final=True,
                            spinner_frame=spinner_frame,
                        )
                    elif current_cluster_status:
                        # If status didn't change, just add a newline to complete the progress line
                        print()

                logger.info(f"Job {job_id} completed successfully.")
                result = status_response["result"]
                # Ensure we return a proper dict[str, Any] instead of Any
                if not isinstance(result, dict):
                    raise CyteTypeAPIError(
                        f"Expected dict result from API, got {type(result)}"
                    )
                return result

            elif status == "failed":
                # Show final cluster status if it wasn't already shown
                if show_progress:
                    if (
                        current_cluster_status
                        and current_cluster_status != last_cluster_status
                    ):
                        display_cluster_status(
                            current_cluster_status,
                            job_id,
                            is_final=True,
                            spinner_frame=spinner_frame,
                        )
                    elif current_cluster_status:
                        # If status didn't change, just add a newline to complete the progress line
                        print()

                raise CyteTypeJobError(f"Server error: {status_response['message']}")

            elif status in ["processing", "pending"]:
                logger.debug(
                    f"Job {job_id} status: {status}. Checking cluster status and waiting {poll_interval}s..."
                )

                if show_progress:
                    # Only show updates if cluster status changed
                    if current_cluster_status != last_cluster_status:
                        if current_cluster_status:
                            display_cluster_status(
                                current_cluster_status,
                                job_id,
                                is_final=False,
                                spinner_frame=spinner_frame,
                            )
                        last_cluster_status = current_cluster_status.copy()
                    else:
                        # Even if status didn't change, update the spinner to show activity
                        if current_cluster_status:
                            display_cluster_status(
                                current_cluster_status,
                                job_id,
                                is_final=False,
                                spinner_frame=spinner_frame,
                            )
                else:
                    # Update last_cluster_status even when not showing progress
                    if current_cluster_status != last_cluster_status:
                        last_cluster_status = current_cluster_status.copy()

                # Sleep for poll_interval seconds, updating spinner every 0.5 seconds
                for _ in range(poll_interval * 2):
                    time.sleep(0.5)
                    spinner_frame += 1

                    # Update spinner during sleep if showing progress
                    if show_progress and current_cluster_status:
                        display_cluster_status(
                            current_cluster_status,
                            job_id,
                            is_final=False,
                            spinner_frame=spinner_frame,
                        )

            elif status == "not_found":
                consecutive_not_found_count += 1

                # If we have an auth token and get many consecutive 404s, it might be auth failure
                if auth_token and consecutive_not_found_count >= 3:
                    logger.warning(
                        "⚠️  Getting consecutive 404 responses with auth token. "
                        "This might indicate authentication issues."
                    )
                    logger.warning(
                        "Please verify your auth_token is valid and has proper permissions."
                    )
                    logger.warning(
                        "If you're using a shared server, contact your administrator."
                    )
                    # Reset counter to avoid spam
                    consecutive_not_found_count = 0

                logger.debug(
                    f"Results endpoint not ready yet for job {job_id} (404). Waiting {poll_interval}s..."
                )
                # Sleep for poll_interval seconds, updating spinner every 0.5 seconds
                for _ in range(poll_interval * 2):
                    time.sleep(0.5)
                    spinner_frame += 1

                    # Show spinner animation during wait if showing progress
                    if show_progress and current_cluster_status:
                        display_cluster_status(
                            current_cluster_status,
                            job_id,
                            is_final=False,
                            spinner_frame=spinner_frame,
                        )

            else:
                logger.warning(
                    f"Job {job_id} has unknown status: '{status}'. Continuing to poll."
                )
                # Sleep for poll_interval seconds, updating spinner every 0.5 seconds
                for _ in range(poll_interval * 2):
                    time.sleep(0.5)
                    spinner_frame += 1

                    # Show spinner animation during wait if showing progress
                    if show_progress and current_cluster_status:
                        display_cluster_status(
                            current_cluster_status,
                            job_id,
                            is_final=False,
                            spinner_frame=spinner_frame,
                        )

        except CyteTypeAPIError as e:
            error_msg = str(e)

            # Check for authentication/authorization errors and show them to users
            if (
                "Authentication failed" in error_msg
                or "Authorization failed" in error_msg
            ):
                print()  # Add newline to clean up progress display
                logger.error(f"❌ {error_msg}")
                logger.error(
                    "Please check your auth_token and ensure it's valid and has the necessary permissions."
                )
                logger.error(
                    "If you're using a shared server, contact your administrator for a valid token."
                )
                raise

            # Handle other network errors with retry logic
            elif "Network error" in error_msg:
                logger.debug(
                    f"Network error during polling request for {job_id}: {e}. Retrying..."
                )
                # Sleep for retry interval with spinner animation
                retry_interval = min(poll_interval, 5)
                for _ in range(retry_interval * 2):
                    time.sleep(0.5)
                    spinner_frame += 1

                    # Show spinner animation during retry if showing progress
                    if (
                        show_progress
                        and "current_cluster_status" in locals()
                        and current_cluster_status
                    ):
                        display_cluster_status(
                            current_cluster_status,
                            job_id,
                            is_final=False,
                            spinner_frame=spinner_frame,
                        )
            else:
                # Clean up progress display before re-raising other errors
                print()  # Add newline to clean up progress display
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
    return make_results_request(job_id, api_url, auth_token)
