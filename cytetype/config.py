import json
import os
import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

from loguru import logger
from pydantic import BaseModel, Field, ValidationError

if TYPE_CHECKING:
    from loguru import Record

DEFAULT_API_URL = "https://cytetype.nygen.io"
DEFAULT_DASHBOARD_URL = "https://cytetype.nygen.io/dashboard"


class StoredCredentials(BaseModel):
    apiUrl: str
    dashboardUrl: str = DEFAULT_DASHBOARD_URL
    apiToken: str = Field(repr=False)
    tokenId: str
    userId: str
    email: str


def normalize_api_url(api_url: str) -> str:
    return api_url.strip().rstrip("/")


def get_default_api_url() -> str:
    configured_api_url = os.environ.get("CYTETYPE_API_URL")
    if configured_api_url is None:
        return DEFAULT_API_URL
    return normalize_api_url(configured_api_url)


def get_credentials_path() -> Path:
    config_home = os.environ.get("XDG_CONFIG_HOME")
    if config_home:
        return Path(config_home).expanduser() / "cytetype" / "credentials.json"
    if os.name == "nt" and os.environ.get("APPDATA"):
        return Path(os.environ["APPDATA"]) / "cytetype" / "credentials.json"
    return Path.home() / ".config" / "cytetype" / "credentials.json"


def save_credentials(credentials: StoredCredentials) -> Path:
    path = get_credentials_path()
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=".credentials-",
        text=True,
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            descriptor = -1
            json.dump(credentials.model_dump(), handle)
            handle.write("\n")
        temporary_path.replace(path)
        path.chmod(0o600)
    except Exception:
        if descriptor >= 0:
            os.close(descriptor)
        temporary_path.unlink(missing_ok=True)
        raise
    return path


def load_credentials(api_url: str | None = None) -> StoredCredentials | None:
    path = get_credentials_path()
    if not path.exists():
        return None
    try:
        credentials = StoredCredentials.model_validate_json(
            path.read_text(encoding="utf-8")
        )
    except (OSError, ValidationError) as error:
        raise ValueError(f"Invalid CyteType credentials file: {path}") from error

    if api_url and credentials.apiUrl != normalize_api_url(api_url):
        return None
    return credentials


def delete_credentials() -> bool:
    path = get_credentials_path()
    if not path.exists():
        return False
    path.unlink()
    return True


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

WRITE_MEM_BUDGET: int = 4 * 1024 * 1024 * 1024  # 4 GB
