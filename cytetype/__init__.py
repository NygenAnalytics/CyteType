__version__ = "0.17.0"

import requests

from .config import logger
from .main import CyteType

__all__ = ["CyteType"]

_PYPI_JSON_URL = "https://pypi.org/pypi/cytetype/json"


def _parse_simple_version(value: str) -> tuple[int, ...] | None:
    cleaned = value.lstrip("v")
    parts = cleaned.split(".")
    if not parts or any(not part.isdigit() for part in parts):
        return None
    return tuple(int(part) for part in parts)


def _warn_if_update_available() -> None:
    current_version = _parse_simple_version(__version__)
    if current_version is None:
        return

    try:
        response = requests.get(_PYPI_JSON_URL, timeout=2)
        response.raise_for_status()
        latest_raw = response.json().get("info", {}).get("version")
        if not isinstance(latest_raw, str):
            return
        latest_version = _parse_simple_version(latest_raw)
        if latest_version is None:
            return
        if latest_version > current_version:
            logger.warning(
                "A newer CyteType version is available on PyPI ({} -> {}). "
                "Run `pip install -U cytetype` to update.",
                __version__,
                latest_raw,
            )
    except Exception:
        # Best-effort check: stay silent on all failures (offline, timeout, parsing, etc).
        return


_warn_if_update_available()
