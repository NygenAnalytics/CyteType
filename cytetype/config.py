from loguru import logger
import sys

logger.remove()

# Apply minimal format
logger.add(
    sys.stdout,
    level="INFO",
    format="{message}",
)


DEFAULT_API_URL = "https://prod.cytetype.nygen.io"
DEFAULT_POLL_INTERVAL = 10
DEFAULT_TIMEOUT = 7200
