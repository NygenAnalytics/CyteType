import requests
from typing import Any, BinaryIO

from .exceptions import create_api_exception, NetworkError, TimeoutError
from .schemas import ErrorResponse


class HTTPTransport:
    """Handles HTTP requests with automatic error parsing."""

    def __init__(self, base_url: str, auth_token: str | None = None):
        self.base_url = base_url.rstrip("/")
        self.auth_token = auth_token
        self.session = requests.Session()

    def _build_headers(self, content_type: str | None = None) -> dict[str, str]:
        """Build request headers with auth token."""
        headers = {}
        if self.auth_token:
            headers["Authorization"] = f"Bearer {self.auth_token}"
        if content_type:
            headers["Content-Type"] = content_type
        return headers

    def _parse_error(self, response: requests.Response) -> None:
        """Parse server error response and raise appropriate exception."""
        try:
            error_data = response.json()

            # Check if error is nested under "detail" key
            if "detail" in error_data and isinstance(error_data["detail"], dict):
                error_data = error_data["detail"]

            error = ErrorResponse(**error_data)
            raise create_api_exception(error.error_code, error.message)
        except ValueError:
            # Not JSON or missing fields - fallback to HTTP status
            raise NetworkError(
                f"HTTP {response.status_code}: {response.text}",
                error_code="HTTP_ERROR",
            )

    def _handle_request_error(self, e: requests.RequestException) -> None:
        """Handle request exceptions uniformly."""
        if isinstance(e, requests.Timeout):
            raise TimeoutError("Request timed out") from e
        elif isinstance(e, requests.ConnectionError):
            raise NetworkError(f"Connection failed: {e}") from e
        elif e.response is not None:
            self._parse_error(e.response)
        else:
            raise NetworkError(f"Request failed: {e}") from e

    def post(
        self, endpoint: str, data: dict[str, Any], timeout: int = 60
    ) -> tuple[int, dict[str, Any]]:
        """Make POST request with error handling."""
        url = f"{self.base_url}/{endpoint.lstrip('/')}"

        try:
            response = self.session.post(
                url,
                json=data,
                headers=self._build_headers(content_type="application/json"),
                timeout=timeout,
            )

            if not response.ok:
                self._parse_error(response)

            return response.status_code, response.json()

        except requests.RequestException as e:
            self._handle_request_error(e)
            raise  # For type checker

    def post_binary(
        self,
        endpoint: str,
        data: bytes | BinaryIO,
        timeout: float | tuple[float, float] = (30.0, 3600.0),
    ) -> tuple[int, dict[str, Any]]:
        """Make POST request with raw binary body (application/octet-stream)."""
        url = f"{self.base_url}/{endpoint.lstrip('/')}"

        try:
            response = self.session.post(
                url,
                data=data,
                headers=self._build_headers(content_type="application/octet-stream"),
                timeout=timeout,
            )

            if not response.ok:
                self._parse_error(response)

            return response.status_code, response.json()

        except requests.RequestException as e:
            self._handle_request_error(e)
            raise  # For type checker

    def get(self, endpoint: str, timeout: int = 30) -> tuple[int, dict[str, Any]]:
        """Make GET request and return (status_code, data)."""
        url = f"{self.base_url}/{endpoint.lstrip('/')}"

        try:
            response = self.session.get(
                url, headers=self._build_headers(), timeout=timeout
            )

            # 404 is expected for status/results endpoints (job not ready)
            if response.status_code == 404:
                return 404, {}

            if not response.ok:
                self._parse_error(response)

            return response.status_code, response.json()

        except requests.RequestException as e:
            self._handle_request_error(e)
            raise  # For type checker
