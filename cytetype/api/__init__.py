"""CyteType API client for server communication."""

from .client import (
    submit_annotation_job,
    get_job_status,
    fetch_job_results,
    wait_for_completion,
)
from .schemas import (
    LLMModelConfig,
    InputData,
    AnnotateRequest,
    UploadResponse,
    UploadedFiles,
    UploadFileKind,
)
from .exceptions import (
    CyteTypeError,
    APIError,
    AuthenticationError,
    AuthorizationError,
    RateLimitError,
    QuotaExceededError,
    JobNotFoundError,
    JobFailedError,
    LLMValidationError,
    ClientDisconnectedError,
    TimeoutError,
    NetworkError,
)

# The same class under the name that `cytetype` exports, which does not shadow the
# builtin TimeoutError.
from .exceptions import TimeoutError as CyteTypeTimeoutError

__all__ = [
    # Client functions
    "submit_annotation_job",
    "get_job_status",
    "fetch_job_results",
    "wait_for_completion",
    # Schemas
    "LLMModelConfig",
    "InputData",
    "AnnotateRequest",
    "UploadResponse",
    "UploadedFiles",
    "UploadFileKind",
    # Exceptions
    "CyteTypeError",
    "APIError",
    "AuthenticationError",
    "AuthorizationError",
    "RateLimitError",
    "QuotaExceededError",
    "JobNotFoundError",
    "JobFailedError",
    "LLMValidationError",
    "ClientDisconnectedError",
    # `from cytetype.api import *` has always exported TimeoutError, shadowing the
    # builtin, and code that catches it expects CyteType's class, so it stays.
    "TimeoutError",
    "CyteTypeTimeoutError",
    "NetworkError",
]
