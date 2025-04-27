from loguru import logger
import sys

logger.remove()

# Apply minimal format
logger.add(
    sys.stdout,
    level="INFO",
    format="<level>{level: <8}</level> <level>{message}</level>",
)

DEFAULT_API_URL = "https://nygen-labs-prod--cell-annotation-agent-fastapi-app.modal.run"
DEFAULT_POLL_INTERVAL = 10
DEFAULT_TIMEOUT = 1200
