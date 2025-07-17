"""
API request utilities for CyteType client.
"""

import requests
from typing import Any, Dict, Optional, Tuple

from .config import logger
from .exceptions import CyteTypeAPIError


def _check_job_status(
    job_id: str, api_url: str, auth_token: Optional[str] = None
) -> Tuple[int, Dict[str, Any]]:
    """Check the status of a job.

    Args:
        job_id: The job ID to check
        api_url: The API base URL
        auth_token: Bearer token for API authentication

    Returns:
        Tuple of (status_code, response_data)

    Raises:
        CyteTypeAPIError: For network or API response errors
    """
    status_url = f"{api_url}/status/{job_id}"

    headers = {}
    if auth_token:
        headers["Authorization"] = f"Bearer {auth_token}"

    try:
        response = requests.get(status_url, headers=headers, timeout=30)

        if response.status_code == 404:
            return 404, {}

        # Check for authentication errors before calling raise_for_status
        if response.status_code == 401:
            logger.debug(f"Authentication failed for job {job_id}: {response.text}")
            raise CyteTypeAPIError(
                "Authentication failed: Invalid or expired auth token"
            )
        elif response.status_code == 403:
            logger.debug(f"Authorization failed for job {job_id}: {response.text}")
            raise CyteTypeAPIError("Authorization failed: Access denied")

        response.raise_for_status()
        return response.status_code, response.json()

    except requests.exceptions.RequestException as e:
        # Don't wrap auth errors that we've already handled above
        if (
            hasattr(e, "response")
            and e.response is not None
            and e.response.status_code in [401, 403]
        ):
            raise  # Re-raise auth errors without wrapping
        logger.debug(f"Network error during status check for job {job_id}: {e}")
        raise CyteTypeAPIError(f"Network error while checking job status: {e}") from e
    except (ValueError, KeyError, requests.exceptions.JSONDecodeError) as e:
        logger.debug(f"Error processing status response for job {job_id}: {e}")
        raise CyteTypeAPIError(
            f"Invalid response while checking job status: {e}"
        ) from e


def _fetch_results(
    job_id: str, api_url: str, auth_token: Optional[str] = None
) -> Tuple[int, Dict[str, Any]]:
    """Fetch results for a completed job.

    Args:
        job_id: The job ID to fetch results for
        api_url: The API base URL
        auth_token: Bearer token for API authentication

    Returns:
        Tuple of (status_code, response_data)

    Raises:
        CyteTypeAPIError: For network or API response errors
    """
    results_url = f"{api_url}/results/{job_id}"

    headers = {}
    if auth_token:
        headers["Authorization"] = f"Bearer {auth_token}"

    try:
        response = requests.get(results_url, headers=headers, timeout=30)

        if response.status_code == 404:
            return 404, {}

        # Check for authentication errors before calling raise_for_status
        if response.status_code == 401:
            logger.debug(f"Authentication failed for job {job_id}: {response.text}")
            raise CyteTypeAPIError(
                "Authentication failed: Invalid or expired auth token"
            )
        elif response.status_code == 403:
            logger.debug(f"Authorization failed for job {job_id}: {response.text}")
            raise CyteTypeAPIError("Authorization failed: Access denied")

        response.raise_for_status()
        return response.status_code, response.json()

    except requests.exceptions.RequestException as e:
        # Don't wrap auth errors that we've already handled above
        if (
            hasattr(e, "response")
            and e.response is not None
            and e.response.status_code in [401, 403]
        ):
            raise  # Re-raise auth errors without wrapping
        logger.debug(f"Network error during results fetch for job {job_id}: {e}")
        raise CyteTypeAPIError(f"Network error while fetching results: {e}") from e
    except (ValueError, KeyError, requests.exceptions.JSONDecodeError) as e:
        logger.debug(f"Error processing results response for job {job_id}: {e}")
        raise CyteTypeAPIError(f"Invalid response while fetching results: {e}") from e


def _transform_results(results_data: Dict[str, Any]) -> Dict[str, Any]:
    """Transform API results data to expected format.

    Args:
        results_data: Raw results data from API

    Returns:
        Transformed results data

    Raises:
        CyteTypeAPIError: For invalid data format
    """
    # Validate results format
    if not isinstance(results_data, dict) or "annotations" not in results_data:
        raise CyteTypeAPIError("Invalid response format from API")

    # Transform new format to expected format
    annotations_dict = results_data.get("annotations", {})
    if not isinstance(annotations_dict, dict):
        raise CyteTypeAPIError("Invalid annotations format from API")

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
                            "cellOntologyTerm", "Unknown"
                        ),
                        # Include additional fields from new format
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

    # Build result in expected format
    return {
        "annotations": annotations_list,
        "summary": results_data.get("summary", {}),
        "semanticOrder": results_data.get("semanticOrder", []),
        "studyContext": results_data.get("studyContext", ""),
        # Include raw annotations for advanced usage
        "raw_annotations": annotations_dict,
    }


def make_results_request(
    job_id: str, api_url: str, auth_token: Optional[str] = None
) -> Dict[str, Any]:
    """Make a request to check job results status.

    This function checks job status, then fetches results if job is completed.

    Args:
        job_id: The job ID to check
        api_url: The API base URL
        auth_token: Bearer token for API authentication

    Returns:
        A dictionary containing:
        - 'status': The job status ('completed', 'processing', 'pending', 'failed', 'not_found')
        - 'result': The result data if status is 'completed'
        - 'message': Status message or error message
        - 'raw_response': The raw API response for debugging

    Raises:
        CyteTypeAPIError: For network or API response errors
    """
    # Check job status
    status_code, status_data = _check_job_status(job_id, api_url, auth_token)

    if status_code == 404:
        return {
            "status": "not_found",
            "result": None,
            "message": "Job not found",
            "raw_response": None,
        }

    job_status = status_data.get("jobStatus")

    if job_status == "completed":
        # Job is completed, fetch the actual results
        results_code, results_data = _fetch_results(job_id, api_url, auth_token)

        if results_code == 404:
            return {
                "status": "failed",
                "result": None,
                "message": "Job completed but results not available",
                "raw_response": status_data,
            }

        # Transform results data
        try:
            result_data = _transform_results(results_data)
        except CyteTypeAPIError as e:
            return {
                "status": "failed",
                "result": None,
                "message": str(e),
                "raw_response": status_data,
            }

        return {
            "status": "completed",
            "result": result_data,
            "message": "Job completed successfully",
            "raw_response": status_data,  # Return status data which contains cluster status
        }

    elif job_status == "failed":
        return {
            "status": "failed",
            "result": None,
            "message": "Job failed",
            "raw_response": status_data,
        }

    elif job_status in ["processing", "pending"]:
        return {
            "status": job_status,
            "result": None,
            "message": f"Job is {job_status}",
            "raw_response": status_data,
        }

    else:
        return {
            "status": "unknown",
            "result": None,
            "message": f"Unknown job status: {job_status}",
            "raw_response": status_data,
        }
