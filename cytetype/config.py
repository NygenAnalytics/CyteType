from __future__ import annotations

import sys
from typing import TYPE_CHECKING

from loguru import logger

if TYPE_CHECKING:
    from loguru import Record

logger.remove()


def _log_format(record: Record) -> str:
    if record["level"].name == "WARNING":
        return "⚠️  {message}\n"
    if record["level"].name == "SUCCESS":
        return "\033[92m✓\033[0m {message}\n"
    return "{message}\n"


logger.add(
    sys.stdout,
    level="INFO",
    format=_log_format,
)
