"""Custom exceptions for the CyteType package."""


class CyteTypeError(Exception):
    """Base class for exceptions in the CyteType package."""

    pass


class CyteTypeAPIError(CyteTypeError):
    """Raised for errors during communication with the CyteType API."""

    pass


class CyteTypeTimeoutError(CyteTypeAPIError):
    """Raised when a timeout occurs interacting with the CyteType API."""

    pass


class CyteTypeJobError(CyteTypeAPIError):
    """Raised when the CyteType API reports an error for a specific job."""

    pass
