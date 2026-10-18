"""Execute selected documentation notebooks locally and save their outputs.

Pages run one at a time in a fresh kernel, from a temporary working directory, so
downloads and CyteType's upload files never land in the repository. Executing a page
that calls ``CyteType.run()`` submits a real annotation job: it needs a saved CyteType
sign-in (``cytetype setup``), counts toward the account's limits, and waits for the job
to finish on CyteType's servers.

Run it from the repository root with the documentation dependencies, for example
``uv run --group docs python -m docs.execute_docs --pages "quickstart"``, or through
``make -C docs execute-docs PAGES="quickstart"``.
"""

import argparse
import os
from datetime import UTC, datetime
from pathlib import Path
from tempfile import TemporaryDirectory

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
        ).execute(cleanup_kc=True, env={**os.environ, **KERNEL_ENVIRONMENT})
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
