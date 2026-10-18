from typing import Any

import pytest

from cytetype.config import _log_format


def _record(level: str) -> Any:
    return {"level": type("Level", (), {"name": level})()}


def test_success_lines_are_colored_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("NO_COLOR", raising=False)
    assert _log_format(_record("SUCCESS")) == "\033[92m✓\033[0m {message}\n"


def test_no_color_removes_ansi_codes(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("NO_COLOR", "1")
    assert _log_format(_record("SUCCESS")) == "✓ {message}\n"
    assert _log_format(_record("INFO")) == "{message}\n"
