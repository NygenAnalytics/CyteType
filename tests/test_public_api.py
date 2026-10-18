"""Public names, provider list and server error-code mapping."""

from typing import Any

import pytest
import requests
from pydantic import ValidationError

import cytetype
import cytetype.api
from cytetype.api import exceptions
from cytetype.api.schemas import LLMModelConfig
from cytetype.api.transport import HTTPTransport

EXCEPTION_NAMES = [
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
    "NetworkError",
]


@pytest.mark.parametrize("name", EXCEPTION_NAMES)
def test_exceptions_are_exported_at_top_level(name: str) -> None:
    assert name in cytetype.__all__
    assert getattr(cytetype, name) is getattr(exceptions, name)


def test_timeout_error_is_exported_without_shadowing_the_builtin() -> None:
    assert cytetype.CyteTypeTimeoutError is exceptions.TimeoutError
    assert "CyteTypeTimeoutError" in cytetype.__all__
    assert "TimeoutError" not in cytetype.__all__

    namespace: dict[str, Any] = {}
    exec("from cytetype import *", namespace)
    assert "TimeoutError" not in namespace
    assert namespace["CyteTypeTimeoutError"] is exceptions.TimeoutError


def test_all_public_names_resolve() -> None:
    for name in cytetype.__all__:
        assert hasattr(cytetype, name), name
    for name in cytetype.api.__all__:
        assert hasattr(cytetype.api, name), name


def test_client_disconnected_error_is_exported_from_api() -> None:
    assert "ClientDisconnectedError" in cytetype.api.__all__
    assert cytetype.api.ClientDisconnectedError is exceptions.ClientDisconnectedError


def test_api_timeout_error_is_importable_but_not_star_exported() -> None:
    from cytetype.api import TimeoutError as ApiTimeoutError

    assert ApiTimeoutError is exceptions.TimeoutError
    assert cytetype.api.CyteTypeTimeoutError is exceptions.TimeoutError
    assert "TimeoutError" not in cytetype.api.__all__
    assert "CyteTypeTimeoutError" in cytetype.api.__all__

    namespace: dict[str, Any] = {}
    exec("from cytetype.api import *", namespace)
    assert "TimeoutError" not in namespace
    assert namespace["CyteTypeTimeoutError"] is exceptions.TimeoutError


def test_authentication_required_maps_to_authentication_error() -> None:
    error = exceptions.create_api_exception(
        "AUTHENTICATION_REQUIRED", "Authentication is required."
    )

    assert type(error) is exceptions.AuthenticationError
    assert error.error_code == "AUTHENTICATION_REQUIRED"
    assert error.message == "Authentication is required."


def test_transport_raises_authentication_error_for_authentication_required() -> None:
    response = requests.Response()
    response.status_code = 401
    response._content = (
        b'{"detail": {"error_code": "AUTHENTICATION_REQUIRED", '
        b'"message": "Authentication is required."}}'
    )

    with pytest.raises(exceptions.AuthenticationError) as error_info:
        HTTPTransport("https://api.test")._parse_error(response)

    assert error_info.value.error_code == "AUTHENTICATION_REQUIRED"


@pytest.mark.parametrize("provider", ["deepseek", "nebius", "together"])
def test_new_llm_providers_are_accepted(provider: str) -> None:
    config = LLMModelConfig(provider=provider, name="model", apiKey="key")  # type: ignore[arg-type]

    assert config.provider == provider


def test_unknown_llm_provider_is_rejected() -> None:
    with pytest.raises(ValidationError):
        LLMModelConfig(provider="unknown", name="model", apiKey="key")  # type: ignore[arg-type]
