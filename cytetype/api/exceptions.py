class CyteTypeError(Exception):
    """Base exception for all CyteType errors."""

    def __init__(self, message: str, error_code: str | None = None):
        super().__init__(message)
        self.error_code = error_code
        self.message = message


class APIError(CyteTypeError):
    """Base class for API-related errors."""

    pass


# API exceptions coupled to specific error codes
class AuthenticationError(APIError):
    """Authentication failed - INVALID_TOKEN, TOKEN_INACTIVE."""

    pass


class AuthorizationError(APIError):
    """Authorization failed - ACCESS_DENIED, LICENSE_EXPIRED, VPN_BLOCKED."""

    pass


class RateLimitError(APIError):
    """Rate limit exceeded - RATE_LIMIT_EXCEEDED."""

    pass


class QuotaExceededError(APIError):
    """Quota exceeded - QUOTA_EXCEEDED."""

    pass


class JobNotFoundError(APIError):
    """Job not found - JOB_NOT_FOUND."""

    pass


class JobFailedError(APIError):
    """Job failed - JOB_FAILED."""

    pass


class LLMValidationError(APIError):
    """LLM validation failed - LLM_VALIDATION_FAILED."""

    pass


# Client-side errors with default messages
class TimeoutError(CyteTypeError):
    """Client-side timeout waiting for results."""

    def __init__(
        self, message: str = "Request timed out", error_code: str | None = None
    ):
        super().__init__(message, error_code)


class NetworkError(CyteTypeError):
    """Network connectivity issue."""

    def __init__(
        self, message: str = "Network connection failed", error_code: str | None = None
    ):
        super().__init__(message, error_code)


# Map server error codes to exception classes
ERROR_CODE_TO_EXCEPTION = {
    "INVALID_TOKEN": AuthenticationError,
    "TOKEN_INACTIVE": AuthenticationError,
    "ACCESS_DENIED": AuthorizationError,
    "LICENSE_EXPIRED": AuthorizationError,
    "VPN_BLOCKED": AuthorizationError,
    "RATE_LIMIT_EXCEEDED": RateLimitError,
    "QUOTA_EXCEEDED": QuotaExceededError,
    "JOB_NOT_FOUND": JobNotFoundError,
    "JOB_FAILED": JobFailedError,
    "LLM_VALIDATION_FAILED": LLMValidationError,
    "JOB_PROCESSING": APIError,  # Generic - expected during polling
    "JOB_NOT_COMPLETED": APIError,  # Generic
    "HTTP_ERROR": APIError,  # Generic
    "INTERNAL_ERROR": APIError,  # Generic
}


def create_api_exception(error_code: str, message: str) -> APIError:
    """Create appropriate exception based on server's error_code.

    Args:
        error_code: Error code from server (e.g., "RATE_LIMIT_EXCEEDED")
        message: Human-readable error message from server

    Returns:
        Specific exception instance with server's message and error_code
    """
    exception_class = ERROR_CODE_TO_EXCEPTION.get(error_code, APIError)
    return exception_class(message, error_code=error_code)
