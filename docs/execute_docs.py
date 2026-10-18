"""Execute selected documentation notebooks locally and save their outputs.

Pages run one at a time in a fresh kernel, from a temporary working directory, so
downloads and CyteType's upload files never land in the repository. Executing a page
that calls ``CyteType.run()`` submits a real annotation job: it needs a saved CyteType
sign-in (``cytetype setup``), counts toward the account's limits, and waits for the job
to finish on CyteType's servers.

Outputs are published, so the kernel does not get environment variables whose names
look secret, and a page whose outputs contain a CyteType API key or another known
secret is not saved.

Run it from the repository root with the documentation dependencies, for example
``uv run --group docs python -m docs.execute_docs --pages "quickstart"``, or through
``make -C docs execute-docs PAGES="quickstart"``.
"""

import argparse
import os
import re
from collections.abc import Iterator, Mapping
from datetime import UTC, datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

DOCS_ROOT = Path(__file__).resolve().parent
SOURCE_DIR = DOCS_ROOT / "source"

# The annotation cell waits for the server job, so allow well beyond the SDK's
# default two-hour wait before giving up on a page.
TIMEOUT_SECONDS = 4 * 60 * 60

# Outputs are published, so keep download and upload progress bars and terminal color
# codes out of them, and render figures inline as PNGs.
KERNEL_ENVIRONMENT = {
    "TQDM_DISABLE": "1",
    "NO_COLOR": "1",
    "MPLBACKEND": "module://matplotlib_inline.backend_inline",
}

# The kernel does not get variables with these words in their names, so a page cannot
# print them. It reads the CyteType sign-in from the credentials file under HOME or
# XDG_CONFIG_HOME, which it keeps, as it keeps PATH and TMPDIR.
SECRET_NAME = re.compile("KEY|TOKEN|SECRET|PASSWORD|PASSWD|CREDENTIAL", re.IGNORECASE)
# Variables that the notebooks need and that are not secret, whatever their names.
KERNEL_ENVIRONMENT_ALLOWLIST = frozenset({"CYTETYPE_API_URL"})
# CyteType API keys: organization keys start with cyt_ and personal keys with cyt_p_.
API_KEY = re.compile(r"cyt_[A-Za-z0-9_-]{16,}")
# Shorter values, such as flags, are too common to look for in outputs.
MIN_SECRET_LENGTH = 8


def page_path(name: str) -> Path:
    """Return the source notebook for ``name``, such as ``quickstart``.

    A bare name may also match a notebook in a section folder, such as
    ``getting-started/quickstart.ipynb``; name the folder when several match.
    """
    relative = Path(name)
    if relative.suffix not in {"", ".ipynb"}:
        raise ValueError(f"Expected a documentation notebook, not {name}")
    relative = relative.with_suffix(".ipynb")
    source = SOURCE_DIR.resolve()

    def notebooks(candidates: list[Path]) -> list[Path]:
        paths = {candidate.resolve() for candidate in candidates}
        return sorted(
            path for path in paths if path.is_relative_to(source) and path.is_file()
        )

    matches = notebooks([SOURCE_DIR / relative])
    if not matches and len(relative.parts) == 1:
        matches = notebooks(list(SOURCE_DIR.glob(f"*/{relative.name}")))
    if not matches:
        raise ValueError(f"Unknown documentation page: {name}")
    if len(matches) > 1:
        found = ", ".join(path.relative_to(source).as_posix() for path in matches)
        raise ValueError(f"Ambiguous documentation page {name}: {found}")
    return matches[0]


def add_execution_timings(notebook) -> None:
    import nbformat

    # The build and its render plugin recognize timing footers by this metadata key,
    # which the documentation toolchain shares with CyteArc and CyteBase.
    for index, cell in enumerate(notebook.cells, start=1):
        if cell.cell_type != "code":
            continue
        cell.outputs = [
            output
            for output in cell.outputs
            if output.get("metadata", {}).get("cytearc_execution_timing") is not True
        ]
        if not cell.source.strip():
            continue
        if type(cell.execution_count) is not int or cell.execution_count < 1:
            raise ValueError(
                f"Cannot show execution timing for unexecuted cell {index}"
            )
        timing = cell.metadata.get("execution", {})
        try:
            started = datetime.fromisoformat(timing["iopub.status.busy"])
            finished = datetime.fromisoformat(timing["iopub.status.idle"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(
                f"Missing or invalid execution timing in cell {index}"
            ) from error
        if started.utcoffset() is None or finished.utcoffset() is None:
            raise ValueError(
                f"Execution timing in cell {index} must include a timezone"
            )
        elapsed = (finished - started).total_seconds()
        if elapsed < 0:
            raise ValueError(f"Execution timing in cell {index} has a negative runtime")
        stamp = started.astimezone(UTC).isoformat().replace("+00:00", "Z")
        text = f"Started: {stamp} | Wall time: {elapsed:.3f} s"
        cell.outputs.append(
            nbformat.v4.new_output(
                "display_data",
                data={"text/plain": text, "text/html": f"<small>{text}</small>"},
                metadata={"cytearc_execution_timing": True},
            )
        )


def is_secret_name(name: str) -> bool:
    """Return whether the environment variable ``name`` may hold a secret."""
    return name not in KERNEL_ENVIRONMENT_ALLOWLIST and bool(SECRET_NAME.search(name))


def kernel_environment(environment: Mapping[str, str]) -> dict[str, str]:
    """Return ``environment`` without its secret-looking variables, for the kernel."""
    kept = {
        name: value for name, value in environment.items() if not is_secret_name(name)
    }
    return {**kept, **KERNEL_ENVIRONMENT}


def saved_api_key() -> str | None:
    """Return the API key that ``cytetype setup`` saved, or None if none can be read."""
    try:
        from cytetype.config import load_credentials

        credentials = load_credentials()
    except Exception:
        return None
    return None if credentials is None else credentials.apiToken


def _strings(value: object) -> Iterator[str]:
    """Yield the strings in a notebook field: text, lists of lines and JSON data."""
    if isinstance(value, str):
        yield value
    elif isinstance(value, list):
        # Notebook files can store text as a list of lines.
        yield "".join(item for item in value if isinstance(item, str))
        for item in value:
            if not isinstance(item, str):
                yield from _strings(item)
    elif isinstance(value, dict):
        for key, item in value.items():
            yield str(key)
            yield from _strings(item)


def output_texts(cell: Mapping[str, Any]) -> Iterator[str]:
    """Yield the text of a cell's outputs: streams, errors and every data format."""
    streams: dict[str, list[str]] = {}
    for output in cell.get("outputs", []):
        if output.get("output_type") == "stream":
            chunks = streams.setdefault(output.get("name", ""), [])
            chunks.extend(_strings(output.get("text")))
        for field in ("text", "ename", "evalue", "traceback", "data"):
            yield from _strings(output.get(field))
    # A stream can arrive split over several outputs, so also check each one whole.
    for chunks in streams.values():
        yield "".join(chunks)


def check_outputs_for_secrets(
    notebook: Any, environment: Mapping[str, str] | None = None
) -> None:
    """Raise ValueError if outputs hold a secret, naming where but never the secret.

    Looks for CyteType API keys, the values of secret-looking variables in
    ``environment`` (by default, this process's environment) and the saved API key,
    which the kernel can read from the credentials file.
    """
    if environment is None:
        environment = os.environ
    # Map each secret to the description that the error shows instead of it.
    descriptions = {
        value: f"the value of {name}"
        for name, value in environment.items()
        if is_secret_name(name) and len(value) >= MIN_SECRET_LENGTH
    }
    api_key = saved_api_key()
    if api_key and len(api_key) >= MIN_SECRET_LENGTH:
        descriptions[api_key] = "the saved CyteType API key"
    places = [
        (f"cell {index}", output_texts(cell))
        for index, cell in enumerate(notebook.cells, start=1)
    ]
    # Widget state, such as the text of progress bars, is saved in the metadata.
    places.append(("the widget state", _strings(notebook.metadata.get("widgets"))))
    leaks = []
    for place, texts in places:
        found = set()
        for text in texts:
            if API_KEY.search(text):
                found.add("a CyteType API key")
            found.update(
                description
                for secret, description in descriptions.items()
                if secret in text
            )
        leaks.extend(f"{description} in {place}" for description in sorted(found))
    if leaks:
        raise ValueError(
            "Refusing to save outputs that contain secrets: " + "; ".join(leaks)
        )


def execute_notebook(source: str) -> str:
    import nbformat
    from jupyter_client.kernelspec import KernelSpecManager
    from jupyter_client.manager import AsyncKernelManager
    from nbclient import NotebookClient

    notebook = nbformat.reads(source, as_version=4)
    nbformat.validate(notebook)
    kernelspec = notebook.metadata.get("kernelspec")
    if not kernelspec or not any(cell.cell_type == "code" for cell in notebook.cells):
        raise ValueError("The selected page is not an executable notebook")
    notebook.metadata.pop("executed_at", None)
    notebook.metadata.pop("executed_code_sha256", None)
    for cell in notebook.cells:
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
            cell.metadata.pop("execution", None)
    # Start this interpreter's own Python kernel, which runs the checkout's cytetype with
    # the locked dependencies, even if another python3 kernelspec is installed.
    kernel_manager = AsyncKernelManager(
        kernel_name=kernelspec.get("name", ""),
        kernel_spec_manager=KernelSpecManager(kernel_dirs=[]),
    )
    with TemporaryDirectory(prefix="cytetype-docs-") as working_directory:
        NotebookClient(
            notebook,
            km=kernel_manager,
            timeout=TIMEOUT_SECONDS,
            allow_errors=False,
            record_timing=True,
            resources={"metadata": {"path": working_directory}},
        ).execute(cleanup_kc=True, env=kernel_environment(os.environ))
    notebook.metadata["executed_at"] = datetime.now(UTC).isoformat()
    return nbformat.writes(notebook)


def main(pages: list[str]) -> None:
    import nbformat

    from docs.build_docs import code_fingerprint, validate_execution

    sources = list(dict.fromkeys(page_path(name) for name in pages))
    if not sources:
        raise ValueError("Specify affected pages, for example PAGES='quickstart'")
    snapshots = []
    for path in sources:
        source = path.read_text()
        notebook = nbformat.reads(source, as_version=4)
        nbformat.validate(notebook)
        if not notebook.metadata.get("kernelspec") or not any(
            cell.cell_type == "code" for cell in notebook.cells
        ):
            raise ValueError(f"Not an executable notebook: {path}")
        snapshots.append((path, source, code_fingerprint(notebook)))

    failures = []
    for path, source, fingerprint in snapshots:
        relative = path.relative_to(SOURCE_DIR.resolve())
        print(f"Executing {relative}", flush=True)
        try:
            notebook = nbformat.reads(execute_notebook(source), as_version=4)
            nbformat.validate(notebook)
            if code_fingerprint(notebook) != fingerprint:
                raise ValueError("Executed code differs from the source notebook")
            notebook.metadata["executed_code_sha256"] = fingerprint
            validate_execution(notebook, relative)
            add_execution_timings(notebook)
            check_outputs_for_secrets(notebook)
            with TemporaryDirectory(prefix=f".{path.stem}-", dir=path.parent) as work:
                temporary = Path(work) / path.name
                nbformat.write(notebook, temporary)
                if path.read_text() != source:
                    raise ValueError(
                        "Notebook changed during execution; preserving current edits"
                    )
                temporary.replace(path)
            print(f"Saved {relative}: {notebook.metadata.executed_at}", flush=True)
        except Exception as error:
            failures.append(f"{relative}: {error}")
            print(f"FAILED {relative}: {error}", flush=True)
    if failures:
        raise SystemExit("Notebook execution failed:\n" + "\n".join(failures))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pages", required=True, help="Space-separated notebook names")
    main(parser.parse_args().pages.split())
