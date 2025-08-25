# Development

## Setup
```bash
git clone https://github.com/NygenAnalytics/cytetype.git
cd cytetype
uv sync --all-extras
uv run pip install -e .
```

## Testing & Linting
```bash
uv run pytest
uv run ruff check .
uv run ruff format .
uv run mypy .
```

## Exceptions
- `CyteTypeError`: Base exception
- `CyteTypeAPIError`: API communication issues
- `CyteTypeTimeoutError`: API polling timed out
- `CyteTypeJobError`: Upstream job reported error 