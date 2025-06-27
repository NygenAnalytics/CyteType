"""Tests for custom exceptions in the CyteType package."""

from cytetype.exceptions import (
    CyteTypeError,
    CyteTypeAPIError,
    CyteTypeTimeoutError,
    CyteTypeJobError,
)


def test_cytetype_error_base_class() -> None:
    """Test the base CyteTypeError exception."""
    msg = "Base error message"
    error = CyteTypeError(msg)
    assert str(error) == msg
    assert isinstance(error, Exception)


def test_cytetype_api_error() -> None:
    """Test CyteTypeAPIError exception."""
    msg = "API communication failed"
    error = CyteTypeAPIError(msg)
    assert str(error) == msg
    assert isinstance(error, CyteTypeError)
    assert isinstance(error, Exception)


def test_cytetype_timeout_error() -> None:
    """Test CyteTypeTimeoutError exception."""
    msg = "Request timed out"
    error = CyteTypeTimeoutError(msg)
    assert str(error) == msg
    assert isinstance(error, CyteTypeAPIError)
    assert isinstance(error, CyteTypeError)


def test_cytetype_job_error() -> None:
    """Test CyteTypeJobError exception."""
    msg = "Job processing failed"
    error = CyteTypeJobError(msg)
    assert str(error) == msg
    assert isinstance(error, CyteTypeAPIError)
    assert isinstance(error, CyteTypeError)


def test_exception_inheritance_chain() -> None:
    """Test that exception inheritance is properly set up."""
    # All custom exceptions should inherit from CyteTypeError
    assert issubclass(CyteTypeAPIError, CyteTypeError)
    assert issubclass(CyteTypeTimeoutError, CyteTypeAPIError)
    assert issubclass(CyteTypeJobError, CyteTypeAPIError)

    # All should ultimately inherit from Exception
    assert issubclass(CyteTypeError, Exception)
    assert issubclass(CyteTypeAPIError, Exception)
    assert issubclass(CyteTypeTimeoutError, Exception)
    assert issubclass(CyteTypeJobError, Exception)


def test_exception_with_chained_cause() -> None:
    """Test exception chaining works properly."""
    original_error = ValueError("Original error")

    try:
        raise CyteTypeAPIError("Wrapped error") from original_error
    except CyteTypeAPIError as e:
        assert str(e) == "Wrapped error"
        assert e.__cause__ is original_error


def test_exception_without_message() -> None:
    """Test exceptions work without explicit messages."""
    error = CyteTypeError()
    assert isinstance(error, CyteTypeError)

    api_error = CyteTypeAPIError()
    assert isinstance(api_error, CyteTypeAPIError)
