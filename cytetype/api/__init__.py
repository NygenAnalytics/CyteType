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
    NetworkError,
)

# Importable by name, but left out of __all__ so `from cytetype.api import *`
# does not shadow the builtin TimeoutError; CyteTypeTimeoutError is the alias.
from .exceptions import TimeoutError as TimeoutError
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
    "CyteTypeTimeoutError",
    "NetworkError",
]
