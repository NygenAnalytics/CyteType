import math
import time
import threading
from pathlib import Path
from typing import Any
from concurrent.futures import ThreadPoolExecutor

from .transport import HTTPTransport
from .progress import ProgressDisplay
from .exceptions import JobFailedError, TimeoutError, APIError, NetworkError
from .schemas import UploadResponse, UploadFileKind
from ..config import logger


MAX_UPLOAD_BYTES: dict[UploadFileKind, int] = {
    "obs_duckdb": 100 * 1024 * 1024,  # 100MB
    "vars_h5": 10 * 1024 * 1024 * 1024,  # 10GB
}

_CHUNK_RETRY_DELAYS = (1, 5, 20)
_RETRYABLE_API_ERROR_CODES = frozenset({"INTERNAL_ERROR", "HTTP_ERROR"})


def _upload_file(
    base_url: str,
    auth_token: str | None,
    file_kind: UploadFileKind,
    file_path: str,
    timeout: float | tuple[float, float] = (30.0, 3600.0),
    max_workers: int = 4,
) -> UploadResponse:
    path_obj = Path(file_path)
    if not path_obj.is_file():
        raise FileNotFoundError(f"Upload file not found: {path_obj}")

    size_bytes = path_obj.stat().st_size
    max_size = MAX_UPLOAD_BYTES[file_kind]
    if size_bytes > max_size:
        raise ValueError(
            f"{file_kind} exceeds upload limit: {size_bytes} bytes > {max_size} bytes"
        )

    transport = HTTPTransport(base_url, auth_token)

    # Step 1 – Initiate chunked upload
    _, init_data = transport.post_empty(f"upload/{file_kind}/initiate", timeout=timeout)
    upload_id: str = init_data["upload_id"]
    chunk_size: int = init_data["chunk_size_bytes"]

    server_max = init_data.get("max_size_bytes")
    if server_max is not None and size_bytes > server_max:
        raise ValueError(
            f"{file_kind} exceeds server upload limit: "
            f"{size_bytes} bytes > {server_max} bytes"
        )

    n_chunks = math.ceil(size_bytes / chunk_size) if size_bytes > 0 else 0

    # Step 2 – Upload chunks in parallel.
    # Each worker thread gets its own HTTPTransport (and thus its own
    # requests.Session / connection pool) for thread safety.
    # Memory is bounded to ~max_workers × chunk_size because each thread
    # reads its chunk on demand via seek+read.
    _tls = threading.local()
    _progress_lock = threading.Lock()
    _chunks_done = [0]

    def _upload_chunk(chunk_idx: int) -> None:
        if not hasattr(_tls, "transport"):
            _tls.transport = HTTPTransport(base_url, auth_token)
        offset = chunk_idx * chunk_size
        read_size = min(chunk_size, size_bytes - offset)
        with path_obj.open("rb") as f:
            f.seek(offset)
            chunk_data = f.read(read_size)

        last_exc: Exception | None = None
        for attempt in range(1 + len(_CHUNK_RETRY_DELAYS)):
            try:
                _tls.transport.put_binary(
                    f"upload/{upload_id}/chunk/{chunk_idx}",
                    data=chunk_data,
                    timeout=timeout,
                )
                with _progress_lock:
                    _chunks_done[0] += 1
                    done = _chunks_done[0]
                pct = 100 * done / n_chunks
                print(
                    f"\r  Uploading: {done}/{n_chunks} chunks ({pct:.0f}%)",
                    end="",
                    flush=True,
                )
                return
            except (NetworkError, TimeoutError) as exc:
                last_exc = exc
            except APIError as exc:
                if exc.error_code in _RETRYABLE_API_ERROR_CODES:
                    last_exc = exc
                else:
                    raise

            if attempt < len(_CHUNK_RETRY_DELAYS):
                delay = _CHUNK_RETRY_DELAYS[attempt]
                logger.warning(
                    "Chunk %d/%d upload failed (attempt %d/%d), retrying in %ds: %s",
                    chunk_idx + 1,
                    n_chunks,
                    attempt + 1,
                    1 + len(_CHUNK_RETRY_DELAYS),
                    delay,
                    last_exc,
                )
                time.sleep(delay)

        raise last_exc  # type: ignore[misc]

    if n_chunks > 0:
        effective_workers = min(max_workers, n_chunks)
        try:
            with ThreadPoolExecutor(max_workers=effective_workers) as pool:
                list(pool.map(_upload_chunk, range(n_chunks)))
            print(f"\r  \033[92m✓\033[0m Uploaded {n_chunks}/{n_chunks} chunks (100%)")
        except BaseException:
            print()  # ensure newline on failure
            raise

    # Step 3 – Complete upload (returns same UploadResponse shape as before)
    _, complete_data = transport.post_empty(
        f"upload/{upload_id}/complete", timeout=timeout
    )
    return UploadResponse(**complete_data)


def upload_obs_duckdb(
    base_url: str,
    auth_token: str | None,
    file_path: str,
    timeout: float | tuple[float, float] = (30.0, 3600.0),
    max_workers: int = 4,
) -> UploadResponse:
    return _upload_file(
        base_url,
        auth_token,
        "obs_duckdb",
        file_path,
        timeout=timeout,
        max_workers=max_workers,
    )


def upload_vars_h5(
    base_url: str,
    auth_token: str | None,
    file_path: str,
    timeout: float | tuple[float, float] = (30.0, 3600.0),
    max_workers: int = 4,
) -> UploadResponse:
    return _upload_file(
        base_url,
        auth_token,
        "vars_h5",
        file_path,
        timeout=timeout,
        max_workers=max_workers,
    )


def submit_annotation_job(
    base_url: str,
    auth_token: str | None,
    payload: dict[str, Any],
) -> str:
    """Submit annotation job and return job_id."""

    transport = HTTPTransport(base_url, auth_token)

    try:
        status_code, response = transport.post("annotate", payload, timeout=180)
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
                logger.success(f"Job {job_id} completed successfully.")
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
