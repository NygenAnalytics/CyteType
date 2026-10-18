"""Tests for the local notebook executor in docs/execute_docs.py.

They start real Jupyter kernels and need the docs dependency group; they skip without it.
"""

from collections.abc import Callable
from copy import deepcopy
from datetime import datetime
import os
from pathlib import Path
import sys
from typing import Any

import pytest

_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
_KERNELSPEC = {"name": "python3", "display_name": "Python 3", "language": "python"}
# Made up, in the format of a personal CyteType API key.
_FAKE_API_KEY = "cyt_p_" + "Fake0Key" * 5


def _require_docs_group() -> Any:
    # The executor validates its results with the build, which needs the docs group.
    pytest.importorskip("jupytext")
    pytest.importorskip("yaml")
    pytest.importorskip("nbclient")
    pytest.importorskip("ipykernel")
    return pytest.importorskip("nbformat")


@pytest.fixture(autouse=True)
def _no_local_secrets(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep this machine's saved sign-in and secret-looking variables out of the checks."""
    from docs import execute_docs

    for name in list(os.environ):
        if execute_docs.is_secret_name(name):
            monkeypatch.delenv(name)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))


@pytest.fixture(scope="module")
def executed_notebook() -> tuple[str, str]:
    nbformat = _require_docs_group()
    from docs import execute_docs

    notebook = nbformat.v4.new_notebook(
        cells=[
            nbformat.v4.new_markdown_cell("# Example"),
            nbformat.v4.new_code_cell('print("RESULT")'),
            nbformat.v4.new_code_cell(""),
            nbformat.v4.new_code_cell(
                "import os, sys\n"
                'print(os.environ.get("TQDM_DISABLE"))\n'
                'print(os.environ.get("MPLBACKEND"))\n'
                "print(os.path.realpath(sys.executable))\n"
                "print(os.getcwd())"
            ),
        ],
        metadata={
            "kernelspec": {
                "name": "python3",
                "display_name": "Python 3",
                "language": "python",
            },
            "executed_at": "old execution",
            "executed_code_sha256": "old fingerprint",
        },
    )
    source = nbformat.writes(notebook)
    environment = dict(os.environ)
    result = execute_docs.execute_notebook(source)
    assert dict(os.environ) == environment
    return source, result


@pytest.fixture
def build_docs() -> Any:
    _require_docs_group()
    from docs import build_docs

    return build_docs


@pytest.fixture
def execution_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, executed_notebook: tuple[str, str]
) -> Any:
    from docs import execute_docs

    source = tmp_path / "source"
    source.mkdir()
    (source / "getting-started").mkdir()
    (source / "example.ipynb").write_text(executed_notebook[0])
    monkeypatch.setattr(execute_docs, "SOURCE_DIR", source)
    return execute_docs


def _execute_results(
    monkeypatch: pytest.MonkeyPatch,
    runner: Any,
    results: list[str | Exception],
    *,
    before_result: Callable[[int, list[str]], None] | None = None,
) -> list[str]:
    """Replace kernel execution with saved results, recording the submitted sources."""
    submitted: list[str] = []

    def execute_notebook(source: str) -> str:
        index = len(submitted)
        submitted.append(source)
        if before_result is not None:
            before_result(index, submitted)
        result = results[index]
        if isinstance(result, Exception):
            raise result
        return result

    monkeypatch.setattr(runner, "execute_notebook", execute_notebook)
    return submitted


def test_execute_notebook_runs_json_notebook_and_records_success(
    executed_notebook: tuple[str, str],
) -> None:
    nbformat = pytest.importorskip("nbformat")

    notebook = nbformat.reads(executed_notebook[1], as_version=4)

    assert notebook.cells[1].execution_count == 1
    assert notebook.cells[1].outputs[0].text == "RESULT\n"
    assert notebook.cells[2].source == ""
    assert datetime.fromisoformat(notebook.metadata.executed_at).utcoffset() is not None
    assert "executed_code_sha256" not in notebook.metadata
    timing = notebook.cells[1].metadata.execution
    assert datetime.fromisoformat(
        timing["iopub.status.idle"]
    ) >= datetime.fromisoformat(timing["iopub.status.busy"])


def test_kernel_runs_this_interpreter_in_a_temporary_directory(
    executed_notebook: tuple[str, str],
) -> None:
    nbformat = pytest.importorskip("nbformat")

    notebook = nbformat.reads(executed_notebook[1], as_version=4)
    lines = notebook.cells[3].outputs[0].text.splitlines()

    assert lines[0] == "1"
    assert lines[1] == "module://matplotlib_inline.backend_inline"
    assert lines[2] == os.path.realpath(sys.executable)
    working_directory = Path(lines[3])
    assert working_directory.name.startswith("cytetype-docs-")
    assert not working_directory.is_relative_to(_REPOSITORY_ROOT)
    assert not working_directory.exists()


def test_execution_timings_preserve_results_and_are_idempotent(
    executed_notebook: tuple[str, str],
) -> None:
    nbformat = pytest.importorskip("nbformat")

    from docs.execute_docs import add_execution_timings

    notebook = nbformat.reads(executed_notebook[1], as_version=4)
    notebook.cells[1].metadata.execution.update(
        {
            "iopub.status.busy": "2026-10-08T12:34:56+02:00",
            "iopub.status.idle": "2026-10-08T12:34:57.250000+02:00",
        }
    )
    original = deepcopy(notebook)

    add_execution_timings(notebook)

    assert notebook.metadata == original.metadata
    assert notebook.cells[0] == original.cells[0]
    assert notebook.cells[2] == original.cells[2]
    cell = notebook.cells[1]
    assert cell.source == original.cells[1].source
    assert cell.execution_count == original.cells[1].execution_count
    assert cell.metadata == original.cells[1].metadata
    assert cell.outputs[:-1] == original.cells[1].outputs
    footer = cell.outputs[-1]
    assert footer.output_type == "display_data"
    assert footer.metadata.cytearc_execution_timing is True
    assert (
        footer.data["text/plain"]
        == "Started: 2026-10-08T10:34:56Z | Wall time: 1.250 s"
    )
    assert footer.data["text/plain"] in footer.data["text/html"]
    nbformat.validate(notebook)
    first = deepcopy(notebook)

    add_execution_timings(notebook)

    assert notebook == first


@pytest.mark.parametrize("failure", ["missing", "invalid", "naive", "negative"])
def test_execution_timings_reject_missing_or_invalid_records(
    executed_notebook: tuple[str, str], failure: str
) -> None:
    nbformat = pytest.importorskip("nbformat")

    from docs.execute_docs import add_execution_timings

    notebook = nbformat.reads(executed_notebook[1], as_version=4)
    timing = notebook.cells[1].metadata.execution
    if failure == "missing":
        timing.pop("iopub.status.idle")
    elif failure == "invalid":
        timing["iopub.status.busy"] = "invalid"
    elif failure == "naive":
        timing["iopub.status.busy"] = "2026-10-08T12:00:00"
    else:
        timing["iopub.status.busy"] = "2026-10-08T12:00:00Z"
        timing["iopub.status.idle"] = "2026-10-08T11:00:00Z"

    with pytest.raises(ValueError, match="timing"):
        add_execution_timings(notebook)

    assert notebook.cells[1].outputs[0].text == "RESULT\n"
    assert len(notebook.cells[1].outputs) == 1


def test_execute_notebook_does_not_return_failed_execution(
    executed_notebook: tuple[str, str],
) -> None:
    nbformat = pytest.importorskip("nbformat")
    from nbclient.exceptions import CellExecutionError

    from docs import execute_docs

    notebook = nbformat.reads(executed_notebook[0], as_version=4)
    notebook.cells[1].source = 'raise RuntimeError("execution failed")'

    with pytest.raises(CellExecutionError, match="execution failed"):
        execute_docs.execute_notebook(nbformat.writes(notebook))


def test_execute_notebook_rejects_unknown_kernels(
    executed_notebook: tuple[str, str],
) -> None:
    nbformat = pytest.importorskip("nbformat")
    from jupyter_client.kernelspec import NoSuchKernel

    from docs import execute_docs

    notebook = nbformat.reads(executed_notebook[0], as_version=4)
    notebook.metadata.kernelspec.name = "ir"

    with pytest.raises(NoSuchKernel):
        execute_docs.execute_notebook(nbformat.writes(notebook))


@pytest.mark.parametrize("change", ["empty", "skip_execution"])
def test_reexecution_clears_old_results_before_running_changed_cells(
    execution_source: Any,
    executed_notebook: tuple[str, str],
    monkeypatch: pytest.MonkeyPatch,
    build_docs: Any,
    change: str,
) -> None:
    nbformat = pytest.importorskip("nbformat")

    code_fingerprint = build_docs.code_fingerprint
    runner = execution_source
    path = runner.SOURCE_DIR / "example.ipynb"
    notebook = nbformat.reads(executed_notebook[1], as_version=4)
    notebook.metadata["executed_code_sha256"] = code_fingerprint(notebook)
    runner.add_execution_timings(notebook)
    cell = notebook.cells[1]
    cell.metadata["custom"] = {"preserve": True}
    cell.source = "" if change == "empty" else 'print("UPDATED")'
    if change == "skip_execution":
        cell.metadata["tags"] = ["skip-execution"]
    source = nbformat.writes(notebook)
    path.write_text(source)
    result = runner.execute_notebook(source)
    _execute_results(monkeypatch, runner, [result])

    if change == "empty":
        runner.main(["example"])
        saved = nbformat.read(path, as_version=4)
        assert saved.cells[1].outputs == []
        assert saved.cells[1].execution_count is None
        assert saved.metadata.executed_code_sha256 == code_fingerprint(saved)
    else:
        with pytest.raises(SystemExit, match="unexecuted cell 2"):
            runner.main(["example"])
        assert path.read_text() == source

    returned = nbformat.reads(result, as_version=4)
    assert returned.cells[1].outputs == []
    assert returned.cells[1].execution_count is None
    assert "execution" not in returned.cells[1].metadata
    assert returned.cells[1].metadata.custom == {"preserve": True}
    assert returned.cells[0] == notebook.cells[0]
    if change == "skip_execution":
        assert returned.cells[1].metadata.tags == ["skip-execution"]


@pytest.mark.parametrize(
    "name",
    [
        "example",
        "example.ipynb",
        "getting-started/example",
        "getting-started/example.ipynb",
    ],
)
def test_page_path_resolves_notebooks(execution_source: Any, name: str) -> None:
    runner = execution_source
    section = runner.SOURCE_DIR / "getting-started" / "example.ipynb"
    if name.startswith("getting-started/"):
        section.write_text((runner.SOURCE_DIR / "example.ipynb").read_text())

    expected = (
        section
        if name.startswith("getting-started/")
        else runner.SOURCE_DIR / "example.ipynb"
    )
    assert runner.page_path(name) == expected.resolve()


def test_page_path_resolves_a_notebook_in_a_section_by_name(
    execution_source: Any,
) -> None:
    runner = execution_source
    path = runner.SOURCE_DIR / "example.ipynb"
    destination = runner.SOURCE_DIR / "getting-started" / path.name
    path.rename(destination)

    assert runner.page_path("example") == destination.resolve()
    assert runner.page_path("example.ipynb") == destination.resolve()


def test_page_path_rejects_names_in_several_sections(execution_source: Any) -> None:
    runner = execution_source
    path = runner.SOURCE_DIR / "example.ipynb"
    (runner.SOURCE_DIR / "guides").mkdir()
    for section in ("getting-started", "guides"):
        (runner.SOURCE_DIR / section / path.name).write_text(path.read_text())

    # A notebook at the top level takes precedence over notebooks in sections.
    assert runner.page_path("example") == path.resolve()
    path.unlink()
    with pytest.raises(ValueError, match="Ambiguous") as error:
        runner.page_path("example")
    assert "getting-started/example.ipynb" in str(error.value)
    assert "guides/example.ipynb" in str(error.value)
    assert (
        runner.page_path("guides/example")
        == (runner.SOURCE_DIR / "guides" / "example.ipynb").resolve()
    )


@pytest.mark.parametrize("name", ["example.md", "getting-started/example.md"])
def test_page_path_rejects_markdown_aliases(execution_source: Any, name: str) -> None:
    runner = execution_source
    (runner.SOURCE_DIR / name).write_text("# Source-only page\n")

    with pytest.raises(ValueError, match="Expected a documentation notebook"):
        runner.page_path(name)


@pytest.mark.parametrize(
    "name", ["../outside.ipynb", "linked.ipynb", "missing", "getting-started"]
)
def test_page_path_rejects_notebooks_outside_source(
    execution_source: Any, name: str
) -> None:
    runner = execution_source
    outside = runner.SOURCE_DIR.parent / "outside.ipynb"
    outside.write_text((runner.SOURCE_DIR / "example.ipynb").read_text())
    (runner.SOURCE_DIR / "linked.ipynb").symlink_to(outside)

    with pytest.raises(ValueError, match="Unknown documentation page"):
        runner.page_path(name)


def test_main_validates_and_replaces_canonical_notebook(
    execution_source: Any,
    executed_notebook: tuple[str, str],
    monkeypatch: pytest.MonkeyPatch,
    build_docs: Any,
) -> None:
    nbformat = pytest.importorskip("nbformat")

    code_fingerprint = build_docs.code_fingerprint
    validate_execution = build_docs.validate_execution
    runner = execution_source
    path = runner.SOURCE_DIR / "example.ipynb"
    submitted = _execute_results(monkeypatch, runner, [executed_notebook[1]])

    runner.main(["example", "example.ipynb"])

    assert submitted == [executed_notebook[0]]
    notebook = nbformat.read(path, as_version=4)
    validate_execution(notebook, Path("example.ipynb"))
    assert notebook.metadata.executed_code_sha256 == code_fingerprint(notebook)
    assert notebook.cells[1].outputs[0].text == "RESULT\n"
    assert notebook.cells[1].outputs[-1].metadata.cytearc_execution_timing is True
    assert {entry.name for entry in runner.SOURCE_DIR.iterdir()} == {
        "example.ipynb",
        "getting-started",
    }


def test_main_continues_after_a_failed_page_and_preserves_its_source(
    execution_source: Any,
    executed_notebook: tuple[str, str],
    monkeypatch: pytest.MonkeyPatch,
    build_docs: Any,
) -> None:
    nbformat = pytest.importorskip("nbformat")

    runner = execution_source
    failed = runner.SOURCE_DIR / "example.ipynb"
    successful = runner.SOURCE_DIR / "getting-started" / "second.ipynb"
    successful.write_text(executed_notebook[0])
    submitted = _execute_results(
        monkeypatch,
        runner,
        [RuntimeError("kernel failed"), executed_notebook[1]],
    )

    with pytest.raises(SystemExit, match="example.ipynb: kernel failed"):
        runner.main(["example", "second"])

    assert len(submitted) == 2
    assert failed.read_text() == executed_notebook[0]
    assert (
        nbformat.read(successful, as_version=4).cells[1].outputs[0].text == "RESULT\n"
    )


def test_main_preserves_edits_made_during_execution(
    execution_source: Any,
    executed_notebook: tuple[str, str],
    monkeypatch: pytest.MonkeyPatch,
    build_docs: Any,
) -> None:
    nbformat = pytest.importorskip("nbformat")

    runner = execution_source
    path = runner.SOURCE_DIR / "example.ipynb"
    edited = nbformat.reads(executed_notebook[0], as_version=4)
    edited.cells[0].source = "# User edit while execution is running"
    text = nbformat.writes(edited)

    def edit_source(index: int, submitted: list[str]) -> None:
        path.write_text(text)

    _execute_results(
        monkeypatch, runner, [executed_notebook[1]], before_result=edit_source
    )

    with pytest.raises(SystemExit, match="Notebook changed during execution"):
        runner.main(["example"])

    assert path.read_text() == text
    assert not list(runner.SOURCE_DIR.glob(".example-*"))


@pytest.mark.parametrize("change", ["code", "empty_cell"])
def test_main_rejects_results_from_different_code(
    execution_source: Any,
    executed_notebook: tuple[str, str],
    monkeypatch: pytest.MonkeyPatch,
    build_docs: Any,
    change: str,
) -> None:
    nbformat = pytest.importorskip("nbformat")

    runner = execution_source
    result = nbformat.reads(executed_notebook[1], as_version=4)
    if change == "code":
        result.cells[1].source = 'print("OTHER")'
    else:
        result.cells.pop(2)
    _execute_results(monkeypatch, runner, [nbformat.writes(result)])

    with pytest.raises(SystemExit, match="Executed code differs"):
        runner.main(["example"])

    assert (runner.SOURCE_DIR / "example.ipynb").read_text() == executed_notebook[0]


@pytest.mark.parametrize("failure", ["timestamp", "unexecuted", "error", "timing"])
def test_main_rejects_incomplete_results_without_overwriting_source(
    execution_source: Any,
    executed_notebook: tuple[str, str],
    monkeypatch: pytest.MonkeyPatch,
    build_docs: Any,
    failure: str,
) -> None:
    nbformat = pytest.importorskip("nbformat")

    runner = execution_source
    notebook = nbformat.reads(executed_notebook[1], as_version=4)
    if failure == "timestamp":
        notebook.metadata.pop("executed_at")
    elif failure == "unexecuted":
        notebook.cells[1].execution_count = None
    elif failure == "timing":
        notebook.cells[1].metadata.pop("execution")
    else:
        notebook.cells[1].outputs = [
            nbformat.v4.new_output(
                "error", ename="RuntimeError", evalue="failed", traceback=["failed"]
            )
        ]
    _execute_results(monkeypatch, runner, [nbformat.writes(notebook)])

    with pytest.raises(SystemExit, match="Notebook execution failed"):
        runner.main(["example"])

    assert (runner.SOURCE_DIR / "example.ipynb").read_text() == executed_notebook[0]


def test_main_checks_all_sources_before_executing(
    execution_source: Any, monkeypatch: pytest.MonkeyPatch, build_docs: Any
) -> None:
    nbformat = pytest.importorskip("nbformat")

    runner = execution_source
    invalid = runner.SOURCE_DIR / "getting-started" / "invalid.ipynb"
    notebook = nbformat.read(runner.SOURCE_DIR / "example.ipynb", as_version=4)
    notebook.metadata.pop("kernelspec")
    nbformat.write(notebook, invalid)
    submitted = _execute_results(monkeypatch, runner, [])

    with pytest.raises(ValueError, match="Not an executable notebook"):
        runner.main(["example", "invalid"])

    assert submitted == []


def test_main_requires_pages(execution_source: Any, build_docs: Any) -> None:
    with pytest.raises(ValueError, match="Specify affected pages"):
        execution_source.main([])


def _notebook_with_outputs(nbformat: Any, outputs: list[Any]) -> Any:
    return nbformat.v4.new_notebook(
        cells=[
            nbformat.v4.new_markdown_cell("# Example"),
            nbformat.v4.new_code_cell(
                'print("RESULT")', execution_count=1, outputs=outputs
            ),
        ],
        metadata={"kernelspec": _KERNELSPEC},
    )


def test_kernel_environment_drops_secret_looking_variables() -> None:
    from docs import execute_docs

    kept = {
        "HOME": "/home/researcher",
        "XDG_CONFIG_HOME": "/home/researcher/.config",
        "PATH": "/usr/local/bin:/usr/bin",
        "TMPDIR": "/scratch/tmp",
        "CYTETYPE_API_URL": "https://cytetype.example",
    }
    dropped = {
        "OPENAI_API_KEY": "sk-example",
        "GITHUB_TOKEN": "ghp-example",
        "aws_secret_access_key": "example",
        "DB_PASSWORD": "example",
        "MYSQL_PASSWD": "example",
        "GOOGLE_APPLICATION_CREDENTIALS": "/home/researcher/service-account.json",
        "Hf_Token": "example",
    }
    environment = {**kept, **dropped, "TQDM_DISABLE": "0"}

    result = execute_docs.kernel_environment(environment)

    assert result == {**kept, **execute_docs.KERNEL_ENVIRONMENT}
    assert environment == {**kept, **dropped, "TQDM_DISABLE": "0"}


def test_kernel_starts_without_secret_looking_variables(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    nbformat = _require_docs_group()
    from docs import execute_docs

    monkeypatch.setenv("CYTETYPE_DOCS_TEST_TOKEN", "kept-out-of-the-kernel")
    monkeypatch.setenv("CYTETYPE_API_URL", "https://cytetype.example")
    names = ["CYTETYPE_DOCS_TEST_TOKEN", "CYTETYPE_API_URL", "XDG_CONFIG_HOME", "HOME"]
    notebook = nbformat.v4.new_notebook(
        cells=[
            nbformat.v4.new_code_cell(
                f"import os\nfor name in {names!r}:\n    print(os.environ.get(name))"
            )
        ],
        metadata={"kernelspec": _KERNELSPEC},
    )

    result = execute_docs.execute_notebook(nbformat.writes(notebook))

    # The kernel keeps what it needs to find the saved CyteType sign-in.
    lines = nbformat.reads(result, as_version=4).cells[0].outputs[0].text.splitlines()
    assert lines == [
        "None",
        "https://cytetype.example",
        os.environ["XDG_CONFIG_HOME"],
        str(os.environ.get("HOME")),
    ]


@pytest.mark.parametrize("place", ["stream", "split_stream", "html", "traceback"])
def test_outputs_with_a_cytetype_api_key_are_refused(place: str) -> None:
    nbformat = pytest.importorskip("nbformat")
    from docs import execute_docs

    new_output = nbformat.v4.new_output
    if place == "stream":
        outputs = [new_output("stream", name="stdout", text=f"{_FAKE_API_KEY}\n")]
    elif place == "split_stream":
        # Neither part matches on its own.
        outputs = [
            new_output("stream", name="stdout", text=f"Key: {_FAKE_API_KEY[:10]}"),
            new_output("stream", name="stdout", text=f"{_FAKE_API_KEY[10:]}\n"),
        ]
    elif place == "html":
        outputs = [
            new_output(
                "display_data",
                data={"text/plain": "Table", "text/html": f"<td>{_FAKE_API_KEY}</td>"},
            )
        ]
    else:
        outputs = [
            new_output(
                "error",
                ename="ValueError",
                evalue="Invalid API token",
                traceback=[f"auth_token = '{_FAKE_API_KEY}'"],
            )
        ]
    notebook = _notebook_with_outputs(nbformat, outputs)

    with pytest.raises(ValueError) as error:
        execute_docs.check_outputs_for_secrets(notebook, environment={})

    # Exact, so the message names the cell and never shows the key.
    assert str(error.value) == (
        "Refusing to save outputs that contain secrets: a CyteType API key in cell 2"
    )


def test_outputs_with_the_value_of_a_secret_variable_are_refused() -> None:
    nbformat = pytest.importorskip("nbformat")
    from docs import execute_docs

    secret = "correct-horse-battery"
    notebook = _notebook_with_outputs(
        nbformat,
        [
            nbformat.v4.new_output(
                "display_data", data={"text/plain": f"password: {secret}"}
            ),
            # Values of other variables, and values too short to check, may appear.
            nbformat.v4.new_output(
                "stream", name="stdout", text="/home/researcher 1234567\n"
            ),
        ],
    )
    environment = {
        "DATABASE_PASSWORD": secret,
        "HOME": "/home/researcher",
        "SHORT_TOKEN": "1234567",
    }

    with pytest.raises(ValueError) as error:
        execute_docs.check_outputs_for_secrets(notebook, environment)

    assert str(error.value) == (
        "Refusing to save outputs that contain secrets: "
        "the value of DATABASE_PASSWORD in cell 2"
    )


def test_outputs_with_the_saved_api_key_are_refused() -> None:
    nbformat = pytest.importorskip("nbformat")
    from cytetype.config import StoredCredentials, save_credentials
    from docs import execute_docs

    # Not in the cyt_ format, so only the saved key itself identifies it.
    key = "saved-key-0123456789abcdef"
    save_credentials(
        StoredCredentials(
            apiUrl="https://cytetype.nygen.io",
            apiToken=key,
            tokenId="token-id",
            userId="user-id",
            email="researcher@university.edu",
        )
    )
    notebook = _notebook_with_outputs(
        nbformat, [nbformat.v4.new_output("stream", name="stderr", text=f"{key}\n")]
    )

    with pytest.raises(ValueError) as error:
        execute_docs.check_outputs_for_secrets(notebook, environment={})

    assert str(error.value) == (
        "Refusing to save outputs that contain secrets: "
        "the saved CyteType API key in cell 2"
    )


def test_outputs_without_secrets_are_accepted(
    executed_notebook: tuple[str, str], monkeypatch: pytest.MonkeyPatch
) -> None:
    nbformat = pytest.importorskip("nbformat")
    from docs import execute_docs

    monkeypatch.setenv("CYTETYPE_DOCS_TEST_TOKEN", "not-in-any-output")
    notebook = nbformat.reads(executed_notebook[1], as_version=4)
    notebook.cells[1].outputs.append(
        nbformat.v4.new_output(
            "display_data",
            data={
                "text/plain": "cytetype_annotation_leiden\n"
                "https://cytetype.nygen.io/report/0123\n"
                "Personal API keys start with cyt_p_.",
                "text/html": "<th>cytetype_annotation_leiden</th>",
            },
        )
    )

    execute_docs.check_outputs_for_secrets(notebook)


def test_main_refuses_to_save_outputs_with_secrets(
    execution_source: Any,
    executed_notebook: tuple[str, str],
    monkeypatch: pytest.MonkeyPatch,
    build_docs: Any,
    capsys: pytest.CaptureFixture[str],
) -> None:
    nbformat = pytest.importorskip("nbformat")

    runner = execution_source
    notebook = nbformat.reads(executed_notebook[1], as_version=4)
    notebook.cells[1].outputs[0].text = f"RESULT {_FAKE_API_KEY}\n"
    _execute_results(monkeypatch, runner, [nbformat.writes(notebook)])

    with pytest.raises(SystemExit) as error:
        runner.main(["example"])

    assert str(error.value) == (
        "Notebook execution failed:\nexample.ipynb: Refusing to save outputs that "
        "contain secrets: a CyteType API key in cell 2"
    )
    assert _FAKE_API_KEY not in capsys.readouterr().out
    assert (runner.SOURCE_DIR / "example.ipynb").read_text() == executed_notebook[0]
