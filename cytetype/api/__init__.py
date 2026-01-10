"""CyteType API client for server communication."""

from .client import (
    submit_annotation_job,
    get_job_status,
    fetch_job_results,
    wait_for_completion,
)
from .schemas import LLMModelConfig, InputData
from .exceptions import (
    CyteTypeError,
    APIError,
    AuthenticationError,
    AuthorizationError,
    RateLimitError,
    QuotaExceededError,
    JobNotFoundError,
    JobFailedError,
    TimeoutError,
    NetworkError,
)

__all__ = [
    # Client functions
    "submit_annotation_job",
    "get_job_status",
    "fetch_job_results",
    "wait_for_completion",
    # Schemas
    "LLMModelConfig",
    "InputData",
    # Exceptions
    "CyteTypeError",
    "APIError",
    "AuthenticationError",
    "AuthorizationError",
    "RateLimitError",
    "QuotaExceededError",
    "JobNotFoundError",
    "JobFailedError",
    "TimeoutError",
    "NetworkError",
]
