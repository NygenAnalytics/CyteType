# Development

## Setup

```bash
git clone https://github.com/NygenAnalytics/CyteType.git
cd CyteType
uv sync
```

`uv sync` installs the package in editable mode with the `dev` dependency group. The version comes
from the git tags (setuptools-scm), so clone with tags to get a meaningful local version.

## Testing and linting

```bash
uv run pytest
uv run ruff check .
uv run ruff format .
uv run mypy .
```

CI runs the tests on Python 3.12 with `uv sync --locked --all-extras`. After changing dependencies,
run `uv lock` and commit `uv.lock`.

## Documentation

The documentation site at https://docs.nygen.io/CyteType/ is built from `docs/`. See
[docs/AGENTS.md](AGENTS.md) for how to build, preview and test it, and how to execute the
quickstart notebook.

## Exceptions

The exception classes are listed in the
[errors reference](https://docs.nygen.io/CyteType/reference/errors/). They are importable from
`cytetype` and `cytetype.api`.
