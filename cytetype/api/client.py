import time
from typing import Any

from .transport import HTTPTransport
from .progress import ProgressDisplay
from .exceptions import JobFailedError, TimeoutError, APIError
from ..config import logger


def submit_annotation_job(
    base_url: str,
    auth_token: str | None,
    payload: dict[str, Any],
) -> str:
    """Submit annotation job and return job_id."""

    transport = HTTPTransport(base_url, auth_token)

    try:
        status_code, response = transport.post("annotate", payload, timeout=60)
        job_id = response.get("job_id")

        if not job_id:
            raise APIError("API response missing job_id", error_code="INVALID_RESPONSE")

        logger.debug(f"Job submitted successfully. Job ID: {job_id}")
        return str(job_id)

    except APIError:
        raise  # Re-raise API errors as-is
    except Exception as e:
        logger.error(f"Failed to submit job: {e}")
        raise


def get_job_status(
    base_url: str,
    auth_token: str | None,
    job_id: str,
) -> dict[str, Any]:
    """Get current status of a job."""

    transport = HTTPTransport(base_url, auth_token)
    status_code, data = transport.get(f"status/{job_id}")

    if status_code == 404:
        return {"jobStatus": "not_found", "clusterStatus": {}}

    return data


def fetch_job_results(
    base_url: str,
    auth_token: str | None,
    job_id: str,
) -> dict[str, Any]:
    """Fetch and transform results for completed job."""
    transport = HTTPTransport(base_url, auth_token)
    status_code, results_data = transport.get(f"results/{job_id}")

    if status_code == 404:
        raise APIError(
            f"Results not found for job {job_id}", error_code="JOB_NOT_FOUND"
        )

    # Validate response format
    if not isinstance(results_data, dict) or "annotations" not in results_data:
        raise APIError(
            "Invalid response format from API", error_code="INVALID_RESPONSE"
        )

    annotations_dict = results_data.get("annotations", {})
    if not isinstance(annotations_dict, dict):
        raise APIError(
            "Invalid annotations format from API", error_code="INVALID_RESPONSE"
        )

    # Convert dictionary format to list format for backward compatibility
    annotations_list = []
    for cluster_id, cluster_data in annotations_dict.items():
        if isinstance(cluster_data, dict) and "latest" in cluster_data:
            latest_data = cluster_data["latest"]
            if isinstance(latest_data, dict) and "annotation" in latest_data:
                annotation_data = latest_data["annotation"]
                if isinstance(annotation_data, dict):
                    # Transform to expected format
                    transformed_annotation = {
                        "clusterId": annotation_data.get("clusterId", cluster_id),
                        "annotation": annotation_data.get("annotation", "Unknown"),
                        "ontologyTerm": annotation_data.get(
                            "cellOntologyTermName", "Unknown"
                        ),
                        "ontologyTermID": annotation_data.get(
                            "cellOntologyTerm", "Unknown"
                        ),
                        "granularAnnotation": annotation_data.get(
                            "granularAnnotation", ""
                        ),
                        "cellState": annotation_data.get("cellState", ""),
                        "justification": annotation_data.get("justification", ""),
                        "supportingMarkers": annotation_data.get(
                            "supportingMarkers", []
                        ),
                        "conflictingMarkers": annotation_data.get(
                            "conflictingMarkers", []
                        ),
                        "missingExpression": annotation_data.get(
                            "missingExpression", ""
                        ),
                        "unexpectedExpression": annotation_data.get(
                            "unexpectedExpression", ""
                        ),
                    }
                    annotations_list.append(transformed_annotation)

    return {
        "annotations": annotations_list,
        "summary": results_data.get("summary", {}),
        "clusterCategories": results_data.get("clusterCategories", []),
        "studyContext": results_data.get("studyContext", ""),
        "raw_annotations": annotations_dict,
    }


def _sleep_with_spinner(
    seconds: int,
    progress: ProgressDisplay | None,
    cluster_status: dict[str, str],
) -> None:
    """Sleep for specified seconds while updating spinner animation.

    Args:
        seconds: Number of seconds to sleep
        progress: ProgressDisplay instance (if showing progress)
        cluster_status: Current cluster status for display
    """
    for _ in range(seconds * 2):
        if progress:
            progress.update(cluster_status)
        time.sleep(0.5)


def wait_for_completion(
    base_url: str,
    auth_token: str | None,
    job_id: str,
    poll_interval: int,
    timeout: int,
    show_progress: bool = True,
) -> dict[str, Any]:
    """Poll job until completion and return results."""
    progress = ProgressDisplay() if show_progress else None
    start_time = time.time()

    logger.info(f"CyteType job (id: {job_id}) submitted. Polling for results...")

    # Initial delay
    time.sleep(5)

    # Show report URL
    report_url = f"{base_url}/report/{job_id}"
    logger.info(f"Report (updates automatically) available at: {report_url}")
    logger.info(
        "If network disconnects, the results can still be fetched:\n"
        "`results = annotator.get_results()`"
    )

    consecutive_not_found = 0

    while (time.time() - start_time) < timeout:
        try:
            status_data = get_job_status(base_url, auth_token, job_id)
            job_status = status_data.get("jobStatus")
            cluster_status = status_data.get("clusterStatus", {})

            # Reset 404 counter on valid response
            if job_status != "not_found":
                consecutive_not_found = 0

            if job_status == "completed":
                if progress:
                    progress.finalize(cluster_status)
                logger.info(f"Job {job_id} completed successfully.")
                return fetch_job_results(base_url, auth_token, job_id)

            elif job_status == "failed":
                if progress:
                    progress.finalize(cluster_status)
                raise JobFailedError(f"Job {job_id} failed")

            elif job_status in ["processing", "pending"]:
                logger.debug(
                    f"Job {job_id} status: {job_status}. Waiting {poll_interval}s..."
                )
                _sleep_with_spinner(poll_interval, progress, cluster_status)

            elif job_status == "not_found":
                consecutive_not_found += 1

                # Warn about consecutive 404s with auth token
                if auth_token and consecutive_not_found >= 3:
                    logger.warning(
                        "⚠️  Getting consecutive 404 responses. "
                        "This might indicate authentication issues."
                    )
                    consecutive_not_found = 0  # Reset to avoid spam

                logger.debug(
                    f"Status endpoint not ready for job {job_id}. "
                    f"Waiting {poll_interval}s..."
                )
                _sleep_with_spinner(poll_interval, progress, cluster_status)

            else:
                logger.warning(f"Unknown job status: '{job_status}'. Continuing...")
                _sleep_with_spinner(poll_interval, progress, cluster_status)

        except APIError:
            # Let API errors (auth, etc.) bubble up immediately
            if progress:
                progress.finalize({})
            raise
        except Exception as e:
            # Network errors - log and retry
            logger.debug(f"Error during polling: {e}. Retrying...")
            retry_interval = min(poll_interval, 5)
            _sleep_with_spinner(retry_interval, progress, cluster_status)

    # Timeout reached
    if progress:
        progress.finalize({})
    raise TimeoutError(f"Job {job_id} did not complete within {timeout}s")
