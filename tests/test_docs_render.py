"""Tests for the documentation build in docs/build_docs.py and docs/render_outputs.py.

The build tests need the docs dependency group and skip without it.
"""

from copy import deepcopy
from datetime import datetime
from html import unescape
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any
from xml.etree import ElementTree

import pytest

_DOCS = Path(__file__).resolve().parents[1] / "docs"


@pytest.fixture
def executable_docs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    pytest.importorskip("jupytext")
    nbformat = pytest.importorskip("nbformat")
    yaml = pytest.importorskip("yaml")
    from docs import build_docs

    source = tmp_path / "source"
    source.mkdir()
    (source / "myst.yml").write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "project": {
                    "github": "https://github.com/example/docs",
                    "toc": [{"file": "example.ipynb"}],
                },
            }
        )
    )
    saved = nbformat.v4.new_notebook(
        cells=[
            nbformat.v4.new_markdown_cell("# Example"),
            nbformat.v4.new_code_cell(
                'print("RESULT")',
                execution_count=1,
                outputs=[
                    nbformat.v4.new_output("stream", name="stdout", text="RESULT\n")
                ],
            ),
        ],
        metadata={"kernelspec": {"name": "python3", "display_name": "Python 3"}},
    )
    saved.metadata["executed_at"] = "2026-10-07T22:00:00+00:00"
    saved.metadata["executed_code_sha256"] = build_docs.code_fingerprint(saved)
    nbformat.write(saved, source / "example.ipynb")
    monkeypatch.setattr(build_docs, "SOURCE", source)
    monkeypatch.setattr(build_docs, "SITE", tmp_path / "site")
    return build_docs


def test_html_requires_saved_execution_for_every_executable_page(
    executable_docs: Any,
) -> None:
    nbformat = pytest.importorskip("nbformat")

    build_docs = executable_docs
    notebook = nbformat.read(build_docs.SOURCE / "example.ipynb", as_version=4)
    notebook.metadata.pop("executed_at")
    notebook.metadata.pop("executed_code_sha256")
    notebook.cells[1].execution_count = None
    notebook.cells[1].outputs = []
    for page in ("example.ipynb", "second.ipynb"):
        nbformat.write(notebook, build_docs.SOURCE / page)

    with pytest.raises(ValueError, match="Saved execution is missing") as error:
        build_docs.prepare_site()

    assert "example.ipynb" in str(error.value)
    assert "second.ipynb" in str(error.value)
    assert "--preview" in str(error.value)


@pytest.mark.parametrize("preview", [False, True])
def test_code_cells_require_kernel_metadata(
    executable_docs: Any, preview: bool
) -> None:
    nbformat = pytest.importorskip("nbformat")

    build_docs = executable_docs
    source = build_docs.SOURCE / "example.ipynb"
    notebook = nbformat.read(source, as_version=4)
    notebook.metadata.pop("kernelspec")
    nbformat.write(notebook, source)

    with pytest.raises(ValueError, match="example.ipynb has no kernelspec metadata"):
        build_docs.prepare_site(preview=preview)


def test_plain_python_examples_need_no_execution(executable_docs: Any) -> None:
    yaml = pytest.importorskip("yaml")

    build_docs = executable_docs
    (build_docs.SOURCE / "example.ipynb").unlink()
    (build_docs.SOURCE / "example.md").write_text(
        '# Example\n\n```python\nprint("RESULT")\n```\n'
    )
    config_path = build_docs.SOURCE / "myst.yml"
    config = yaml.safe_load(config_path.read_text())
    config["project"]["toc"] = [{"file": "example.md"}]
    config_path.write_text(yaml.safe_dump(config))

    site = build_docs.prepare_site()

    assert (site / "example.md").is_file()
    assert not (site / "example.ipynb").exists()


def test_source_and_edit_links_target_master(executable_docs: Any) -> None:
    nbformat = pytest.importorskip("nbformat")
    yaml = pytest.importorskip("yaml")

    build_docs = executable_docs
    (build_docs.SOURCE / "guide.md").write_text("# Guide\n")

    site = build_docs.prepare_site()
    notebook = nbformat.read(site / "example.ipynb", as_version=4)
    guide = yaml.safe_load((site / "guide.md").read_text().split("---", 2)[1])

    for page, metadata in (("example.ipynb", notebook.metadata), ("guide.md", guide)):
        assert metadata["source_url"] == (
            f"https://github.com/example/docs/blob/master/docs/source/{page}"
        )
        assert metadata["edit_url"] == (
            f"https://github.com/example/docs/edit/master/docs/source/{page}"
        )


def test_saved_notebooks_restore_complete_outputs_in_a_clean_checkout(
    executable_docs: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    nbformat = pytest.importorskip("nbformat")
    yaml = pytest.importorskip("yaml")

    build_docs = executable_docs
    held = build_docs.SOURCE / "held.ipynb"
    held.write_text("Excluded outputs")
    config_path = build_docs.SOURCE / "myst.yml"
    config = yaml.safe_load(config_path.read_text())
    config["project"]["exclude"] = ["held.ipynb"]
    config_path.write_text(yaml.safe_dump(config))
    checkout = tmp_path / "checkout" / "docs"
    shutil.copytree(build_docs.SOURCE, checkout / "source")
    monkeypatch.setattr(build_docs, "SOURCE", checkout / "source")
    monkeypatch.setattr(build_docs, "SITE", checkout / "build" / "site")
    assert not (checkout / "build").exists()
    source = build_docs.SOURCE / "example.ipynb"
    current = nbformat.read(source, as_version=4)
    current.cells[0].source = "# Updated"
    nbformat.write(current, source)
    source_text = source.read_text()

    site = build_docs.prepare_site()
    notebook = nbformat.read(site / "example.ipynb", as_version=4)

    assert notebook.cells[0].source == "# Updated"
    assert notebook.cells[1].outputs[0].text == "RESULT\n"
    assert notebook.cells[1].execution_count == 1
    assert notebook.metadata["executed_at"] == "2026-10-07T22:00:00+00:00"
    assert (
        notebook.metadata["executed_code_sha256"]
        == (current.metadata["executed_code_sha256"])
    )
    assert not (site / "held.ipynb").exists()
    assert source.read_text() == source_text

    current.cells[1].source = 'print("CHANGED")'
    nbformat.write(current, source)
    with pytest.raises(ValueError, match="Code changed"):
        build_docs.prepare_site()


def test_preview_labels_missing_execution_without_fabricating_outputs(
    executable_docs: Any,
) -> None:
    nbformat = pytest.importorskip("nbformat")

    build_docs = executable_docs
    source = build_docs.SOURCE / "example.ipynb"
    current = nbformat.read(source, as_version=4)
    current.metadata.pop("executed_at")
    current.metadata.pop("executed_code_sha256")
    current.cells[1].execution_count = None
    current.cells[1].outputs = []
    current.cells[1].metadata["figure_descriptions"] = [
        "An embedding that will be generated when this page is executed."
    ]
    nbformat.write(current, source)
    source_text = source.read_text()

    site = build_docs.prepare_site(preview=True)
    notebook = nbformat.read(site / "example.ipynb", as_version=4)

    assert notebook.cells[0].cell_type == "markdown"
    assert notebook.cells[0].source.startswith(":::{warning} Source preview\n")
    assert "Figures and results are not available" in notebook.cells[0].source
    assert "executed_at" not in notebook.metadata
    assert "executed_code_sha256" not in notebook.metadata
    assert all(
        cell.execution_count is None and cell.outputs == []
        for cell in notebook.cells
        if cell.cell_type == "code"
    )
    assert notebook.cells[2].metadata.figure_descriptions == (
        current.cells[1].metadata.figure_descriptions
    )
    assert source.read_text() == source_text


@pytest.mark.parametrize("preview", [False, True])
@pytest.mark.parametrize(
    ("failure", "message"),
    [
        ("timestamp", "no completed execution timestamp"),
        ("fingerprint", "fingerprint"),
        ("unexecuted", "unexecuted cell"),
        ("error", "has an error in cell"),
        ("changed_code", "Code changed"),
    ],
)
def test_build_rejects_invalid_saved_execution(
    executable_docs: Any, preview: bool, failure: str, message: str
) -> None:
    nbformat = pytest.importorskip("nbformat")

    build_docs = executable_docs
    saved_path = build_docs.SOURCE / "example.ipynb"
    saved = nbformat.read(saved_path, as_version=4)
    if failure == "timestamp":
        saved.metadata.pop("executed_at")
    elif failure == "fingerprint":
        saved.metadata.pop("executed_code_sha256")
    elif failure == "unexecuted":
        saved.cells[1].execution_count = None
    elif failure == "error":
        saved.cells[1].outputs = [
            nbformat.v4.new_output(
                "error", ename="RuntimeError", evalue="failed", traceback=["failed"]
            )
        ]
    else:
        saved.cells[1].source = 'print("OLD RESULT")'
    nbformat.write(saved, saved_path)

    with pytest.raises(ValueError, match=message):
        build_docs.prepare_site(preview=preview)


def test_completed_cells_need_not_produce_outputs(executable_docs: Any) -> None:
    nbformat = pytest.importorskip("nbformat")

    build_docs = executable_docs
    saved_path = build_docs.SOURCE / "example.ipynb"
    saved = nbformat.read(saved_path, as_version=4)
    saved.cells[1].source = "value = 1"
    saved.cells[1].outputs = []
    saved.metadata["executed_code_sha256"] = build_docs.code_fingerprint(saved)
    nbformat.write(saved, saved_path)

    site = build_docs.prepare_site()
    notebook = nbformat.read(site / "example.ipynb", as_version=4)

    assert notebook.cells[1].execution_count == 1
    assert notebook.cells[1].outputs == []
    assert "Source preview" not in notebook.cells[0].source


@pytest.mark.parametrize(
    ("exclude", "suffix"),
    [
        ("tutorials/held.ipynb", ".ipynb"),
        ("tutorials/**", ".ipynb"),
        ("**/held.ipynb", ".ipynb"),
        ("tutorials/held.md", ".md"),
    ],
)
def test_excluded_pages_need_no_execution_and_leave_no_generated_routes(
    executable_docs: Any, exclude: str, suffix: str
) -> None:
    nbformat = pytest.importorskip("nbformat")
    yaml = pytest.importorskip("yaml")

    build_docs = executable_docs
    held = build_docs.SOURCE / "tutorials" / f"held{suffix}"
    held.parent.mkdir()
    if suffix == ".md":
        held.write_text("# Held\n\n```{code-cell} python\nraise RuntimeError\n```\n")
    else:
        held.write_text("Invalid excluded notebook")
    config_path = build_docs.SOURCE / "myst.yml"
    config = yaml.safe_load(config_path.read_text())
    config["project"]["exclude"] = [exclude]
    config["project"]["toc"].append(
        {"title": "Tutorials", "children": [{"file": f"tutorials/held{suffix}"}]}
    )
    config_path.write_text(yaml.safe_dump(config))
    old_html = build_docs.SITE / "_build" / "html"
    old_html.mkdir(parents=True)
    (old_html / "held.html").write_text("Old held page")
    (old_html / "held.png").write_bytes(b"old figure")
    (build_docs.SITE / "held.ipynb").write_text("Old staged notebook")

    site = build_docs.prepare_site()
    notebook = nbformat.read(site / "example.ipynb", as_version=4)
    staged_config = yaml.safe_load((site / "myst.yml").read_text())

    assert notebook.cells[1].outputs[0].text == "RESULT\n"
    assert not (site / "tutorials" / "held.md").exists()
    assert not (site / "tutorials" / "held.ipynb").exists()
    assert not (site / "held.ipynb").exists()
    assert not old_html.exists()
    assert staged_config["project"]["toc"][1]["children"] == []
    assert held.is_file()
    assert yaml.safe_load(config_path.read_text()) == config


@pytest.mark.parametrize("change", ["edit", "reorder", "delete"])
@pytest.mark.parametrize("preview", [False, True])
def test_saved_outputs_follow_unchanged_code_and_current_prose(
    executable_docs: Any, change: str, preview: bool
) -> None:
    nbformat = pytest.importorskip("nbformat")

    build_docs = executable_docs
    source = build_docs.SOURCE / "example.ipynb"
    saved = nbformat.read(source, as_version=4)
    saved.cells.append(
        nbformat.v4.new_code_cell(
            'print("SECOND")',
            execution_count=2,
            outputs=[nbformat.v4.new_output("stream", name="stdout", text="SECOND\n")],
        )
    )
    saved.metadata["executed_code_sha256"] = build_docs.code_fingerprint(saved)
    current = deepcopy(saved)
    current.cells[0].source = "# Edited prose"
    current.cells.insert(2, nbformat.v4.new_markdown_cell("More explanation."))
    nbformat.write(current, source)

    site = build_docs.prepare_site(preview=preview)
    rendered = nbformat.read(site / "example.ipynb", as_version=4)

    assert rendered.cells[0].source == "# Edited prose"
    assert rendered.cells[2].source == "More explanation."
    assert [cell for cell in rendered.cells if cell.cell_type == "code"] == (
        [cell for cell in saved.cells if cell.cell_type == "code"]
    )
    assert rendered.metadata["executed_at"] == saved.metadata["executed_at"]
    assert (
        rendered.metadata["executed_code_sha256"]
        == (saved.metadata["executed_code_sha256"])
    )

    changed = deepcopy(current)
    if change == "edit":
        changed.cells[1].source = 'print("CHANGED")'
    elif change == "reorder":
        changed.cells[1], changed.cells[3] = changed.cells[3], changed.cells[1]
    else:
        del changed.cells[3]
    nbformat.write(changed, source)
    with pytest.raises(ValueError, match='execute-docs PAGES="example"'):
        build_docs.prepare_site(preview=preview)


def test_rerun_command_names_the_notebook_in_its_section(
    executable_docs: Any,
) -> None:
    nbformat = pytest.importorskip("nbformat")
    yaml = pytest.importorskip("yaml")

    build_docs = executable_docs
    section = build_docs.SOURCE / "getting-started"
    section.mkdir()
    notebook = nbformat.read(build_docs.SOURCE / "example.ipynb", as_version=4)
    notebook.cells[1].source = 'print("CHANGED")'
    nbformat.write(notebook, section / "quickstart.ipynb")
    (build_docs.SOURCE / "example.ipynb").unlink()
    config_path = build_docs.SOURCE / "myst.yml"
    config = yaml.safe_load(config_path.read_text())
    config["project"]["toc"] = [{"file": "getting-started/quickstart.ipynb"}]
    config_path.write_text(yaml.safe_dump(config))

    with pytest.raises(ValueError) as error:
        build_docs.prepare_site()

    assert 'make -C docs execute-docs PAGES="getting-started/quickstart"' in str(
        error.value
    )


@pytest.mark.parametrize("state", ["count", "output", "timestamp", "fingerprint"])
def test_preview_rejects_partial_execution(executable_docs: Any, state: str) -> None:
    nbformat = pytest.importorskip("nbformat")

    build_docs = executable_docs
    source = build_docs.SOURCE / "example.ipynb"
    current = nbformat.read(source, as_version=4)
    current.metadata.pop("executed_at")
    current.metadata.pop("executed_code_sha256")
    current.cells[1].execution_count = None
    current.cells[1].outputs = []
    if state == "count":
        current.cells[1].execution_count = 1
    elif state == "output":
        current.cells[1].outputs = [
            nbformat.v4.new_output("stream", name="stdout", text="Partial result\n")
        ]
    elif state == "timestamp":
        current.metadata["executed_at"] = "2026-10-07T22:00:00+00:00"
    else:
        current.metadata["executed_code_sha256"] = build_docs.code_fingerprint(current)
    nbformat.write(current, source)
    source_text = source.read_text()

    with pytest.raises(ValueError):
        build_docs.prepare_site(preview=True)

    assert source.read_text() == source_text


@pytest.mark.parametrize("preview", [False, True])
def test_executable_markdown_requires_notebook_source(
    executable_docs: Any, preview: bool
) -> None:
    yaml = pytest.importorskip("yaml")

    build_docs = executable_docs
    (build_docs.SOURCE / "example.ipynb").unlink()
    (build_docs.SOURCE / "example.md").write_text(
        "---\nkernelspec:\n  name: python3\n  display_name: Python 3\n---\n"
        '# Example\n\n```{code-cell} python\nprint("RESULT")\n```\n'
    )
    config_path = build_docs.SOURCE / "myst.yml"
    config = yaml.safe_load(config_path.read_text())
    config["project"]["toc"] = [{"file": "example.md"}]
    config_path.write_text(yaml.safe_dump(config))

    with pytest.raises(ValueError, match=r"\.ipynb") as error:
        build_docs.prepare_site(preview=preview)

    assert "example.md" in str(error.value)


@pytest.mark.parametrize(
    ("base_url", "prefix"),
    [
        (None, ""),
        ("", ""),
        ("/", ""),
        ("/CyteType", "/CyteType"),
        ("/CyteType/", "/CyteType"),
        ("/docs/CyteType/", "/docs/CyteType"),
    ],
)
def test_prepared_sidebar_urls_respect_site_base_url(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    base_url: str | None,
    prefix: str,
) -> None:
    pytest.importorskip("jupytext")
    pytest.importorskip("nbformat")
    yaml = pytest.importorskip("yaml")
    from docs import build_docs

    source = tmp_path / "source"
    source.mkdir()
    (source / "index.md").write_text("# Example\n")
    unchanged_urls = [
        "https://example.org/api/",
        "//example.org/api/",
        "mailto:docs@example.org",
        "relative.html",
        "#section",
    ]
    config = {
        "version": 1,
        "project": {
            "github": "https://github.com/example/docs",
            "toc": [
                {"file": "index.md"},
                {"url": "/quickstart#section"},
                {
                    "title": "Reference",
                    "children": [
                        {"url": "/api/index.html", "open_in_same_tab": True},
                        *({"url": url} for url in unchanged_urls),
                    ],
                },
            ],
        },
    }
    source_config = yaml.safe_dump(config)
    (source / "myst.yml").write_text(source_config)
    monkeypatch.setattr(build_docs, "SOURCE", source)
    monkeypatch.setattr(build_docs, "SITE", tmp_path / "site")
    if base_url is None:
        monkeypatch.delenv("BASE_URL", raising=False)
    else:
        monkeypatch.setenv("BASE_URL", base_url)

    site = build_docs.prepare_site()
    toc = yaml.safe_load((site / "myst.yml").read_text())["project"]["toc"]

    assert toc[0] == {"file": "index.md"}
    assert toc[1]["url"] == f"{prefix}/quickstart#section"
    assert toc[2]["children"][0] == {
        "url": f"{prefix}/api/index.html",
        "open_in_same_tab": True,
    }
    assert [entry["url"] for entry in toc[2]["children"][1:]] == unchanged_urls
    assert (source / "myst.yml").read_text() == source_config


@pytest.mark.parametrize(
    ("page", "link", "expected"),
    [
        (
            "guides/prepare-data.md",
            "../reference/api/index.html",
            "/api/index.html",
        ),
        (
            "guides/prepare-data.md",
            "../reference/api/cytetype.html#cytetype.CyteType.run",
            "/api/cytetype.html#cytetype.CyteType.run",
        ),
        (
            "results/report.md",
            "../reference/api/helpers.md#cytetype.marker_dotplot",
            "/api/helpers.html#cytetype.marker_dotplot",
        ),
        ("index.md", "reference/api/exceptions.html", "/api/exceptions.html"),
        ("index.md", "/api/index.html", "/api/index.html"),
        ("index.md", "reference/r-api.md", "reference/r-api.md"),
        (
            "guides/automation.md",
            "../reference/cli.md#setup",
            "../reference/cli.md#setup",
        ),
        ("index.md", "https://example.org/reference/api/x.html", None),
        ("index.md", "#reference", None),
    ],
)
def test_relative_api_links_point_to_the_sphinx_pages(
    page: str, link: str, expected: str | None
) -> None:
    pytest.importorskip("jupytext")
    pytest.importorskip("nbformat")
    pytest.importorskip("yaml")
    from docs import build_docs

    text = f"See [the API]({link}) and [more]({link}).\n"

    rewritten = build_docs.rewrite_api_links(text, Path(page))

    target = link if expected is None else expected
    assert rewritten == f"See [the API]({target}) and [more]({target}).\n"


def test_staged_pages_link_to_the_api_without_changing_sources(
    executable_docs: Any,
) -> None:
    nbformat = pytest.importorskip("nbformat")
    yaml = pytest.importorskip("yaml")

    build_docs = executable_docs
    guide = build_docs.SOURCE / "guides" / "guide.md"
    guide.parent.mkdir()
    guide_text = (
        "---\ndescription: A guide.\n---\n\n# Guide\n\n"
        "Read [run()](../reference/api/cytetype.html#cytetype.CyteType.run).\n"
    )
    guide.write_text(guide_text)
    source = build_docs.SOURCE / "example.ipynb"
    notebook = nbformat.read(source, as_version=4)
    notebook.cells[0].source = "# Example\n\nSee the [API](reference/api/index.html)."
    nbformat.write(notebook, source)
    notebook_text = source.read_text()

    site = build_docs.prepare_site()

    staged_guide = (site / "guides" / "guide.md").read_text()
    assert "[run()](/api/cytetype.html#cytetype.CyteType.run)" in staged_guide
    assert yaml.safe_load(staged_guide.split("---", 2)[1])["description"] == "A guide."
    staged = nbformat.read(site / "example.ipynb", as_version=4)
    assert staged.cells[0].source == "# Example\n\nSee the [API](/api/index.html)."
    assert guide.read_text() == guide_text
    assert source.read_text() == notebook_text


@pytest.mark.parametrize("figure_count", [1, 2])
def test_staged_figure_uses_authored_alt_without_changing_saved_output(
    executable_docs: Any,
    figure_count: int,
) -> None:
    nbformat = pytest.importorskip("nbformat")

    build_docs = executable_docs
    source = build_docs.SOURCE / "example.ipynb"
    notebook = nbformat.read(source, as_version=4)
    image = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+aD1sAAAAASUVORK5CYII="
    descriptions = [
        "UMAP of the PBMC clusters, colored by CyteType annotation.",
        "Dot plot of the marker genes that support each label.",
    ][:figure_count]
    notebook.cells[1].metadata["figure_descriptions"] = descriptions
    notebook.cells[1].outputs = [
        nbformat.v4.new_output(
            "display_data",
            data={
                "image/png": image,
                "text/plain": "<Figure size 400x300 with 1 Axes>",
            },
        )
        for _ in descriptions
    ]
    nbformat.write(notebook, source)
    original = source.read_bytes()

    site = build_docs.prepare_site()

    staged = nbformat.read(site / "example.ipynb", as_version=4)
    assert [
        output.data["text/plain"] for output in staged.cells[1].outputs
    ] == descriptions
    assert all(output.data["image/png"] == image for output in staged.cells[1].outputs)
    assert build_docs.code_fingerprint(staged) == notebook.metadata.executed_code_sha256
    assert source.read_bytes() == original


@pytest.mark.parametrize(
    "descriptions",
    [None, [], [""], ["   "], ["First", "Unexpected second"], "not a list"],
)
def test_staged_figure_rejects_missing_description(
    executable_docs: Any, descriptions: object
) -> None:
    nbformat = pytest.importorskip("nbformat")

    build_docs = executable_docs
    source = build_docs.SOURCE / "example.ipynb"
    notebook = nbformat.read(source, as_version=4)
    notebook.cells[1].outputs = [
        nbformat.v4.new_output(
            "display_data",
            data={"image/png": "aW1hZ2U=", "text/plain": "<Figure with 1 Axes>"},
        )
    ]
    if descriptions is not None:
        notebook.cells[1].metadata["figure_descriptions"] = descriptions
    nbformat.write(notebook, source)

    with pytest.raises(ValueError, match="example.ipynb.*metadata.figure_descriptions"):
        build_docs.prepare_site()


@pytest.fixture(params=["", "/CyteType"])
def rendered_discovery_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, request: pytest.FixtureRequest
) -> tuple[Path, Path, str]:
    pytest.importorskip("jupytext")
    pytest.importorskip("nbformat")
    pytest.importorskip("yaml")
    base_url = request.param
    monkeypatch.setenv("BASE_URL", base_url)
    html = tmp_path / "_build" / "html"
    api_text = tmp_path / "api-text"
    for directory in (html, api_text, html / "api" / "_sources"):
        directory.mkdir(parents=True)

    references = []
    for name, route in (("index", "/"), ("quickstart", "/quickstart")):
        references.append({"kind": "page", "url": route, "data": f"/{name}.json"})
        children = (
            [
                {
                    "type": "outputs",
                    "id": "saved-result",
                    "children": [
                        {
                            "type": "output",
                            "jupyter_data": {
                                "output_type": "execute_result",
                                "data": {
                                    "text/plain": {"content": "cell_count\n    1000"}
                                },
                            },
                        }
                    ],
                }
            ]
            if name == "quickstart"
            else []
        )
        (html / f"{name}.json").write_text(
            json.dumps(
                {
                    "frontmatter": {"title": name},
                    "mdast": {"type": "root", "children": children},
                }
            )
        )
        saved_output = (
            '<div data-mdast-node-id="saved-result">Loading...'
            '<img src="plot.png"/><p>LEAKED_OUTPUT</p></div>'
            if name == "quickstart"
            else ""
        )
        page = html / route.strip("/") / "index.html"
        page.parent.mkdir(parents=True, exist_ok=True)
        page.write_text(
            f'<html><head><title>{name}</title><meta name="description" '
            f'content="Read the {name} guide."></head>'
            f"<body><nav>Site navigation</nav><article><h1>{name}</h1>"
            "<p>Saved result: 42 cells.</p>"
            "<pre>def count():\n    return 42\n\nprint(count())\n42</pre>"
            "<table><tr><th>Group</th><th>Cells</th></tr>"
            "<tr><td>T cells</td><td>42</td></tr></table>"
            '<p>Read <a href="details#reading">the details</a>.</p>'
            '<img src="plot.png" alt="Cell embedding">'
            f"{saved_output}"
            "<p>Analysis continues after the saved output.</p>"
            "</article><footer>Footer navigation</footer></body></html>"
        )
    references.extend(
        [
            {
                "kind": "page",
                "identifier": "quickstart",
                "url": "/quickstart",
                "data": "/quickstart.json",
            },
            {
                "kind": "heading",
                "identifier": "results",
                "url": "/quickstart",
                "data": "/quickstart.json",
            },
        ]
    )
    (html / "myst.xref.json").write_text(json.dumps({"references": references}))
    (html / "api" / "_sources" / "cytetype.md").write_text(
        "# CyteType class\n\n```{autoclass} cytetype.CyteType\n:members:\n```\n"
    )
    (html / "api" / "cytetype.html").write_text(
        "<html><head><title>CyteType class</title>"
        '<meta name="description" content="Read the Python class reference.">'
        '<link rel="canonical" '
        'href="https://docs.nygen.io/CyteType/api/cytetype.html" /></head>'
        "<body>CyteType(adata, group_key, rank_key='rank_genes_groups')</body></html>"
    )
    (api_text / "cytetype.txt").write_text(
        "CyteType class\n**************\n\n"
        "class cytetype.CyteType(adata, group_key, rank_key='rank_genes_groups')\n"
        "   Prepare an AnnData object for annotation.\n"
    )
    for relative in (
        "build/api/cytetype.html",
        "api/search.html",
        "api/genindex.html",
    ):
        page = html / relative
        page.parent.mkdir(parents=True, exist_ok=True)
        page.write_text("<html><head></head><body>Extra generated page</body></html>")
    (html / "llms.txt").write_text(
        "# CyteType\n\n- [Home](./): overview\n"
        "\n## Getting started\n\n"
        "- [Quickstart (Python)](quickstart/): first analysis\n"
        "\n## Reference\n\n"
        "- [CyteType class](api/cytetype.html): public API\n"
    )
    (html / "robots.txt").write_text(
        "User-agent: *\nAllow: /\nSitemap: http://localhost:3000/sitemap.xml\n"
    )
    (html / "sitemap.xml").write_text(
        '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'
        "<url><loc>http://localhost:3000/quickstart</loc></url></urlset>"
    )
    return html, api_text, base_url


def _discovery_links(page: Path) -> list[dict[str, str | None]]:
    class LinkParser(HTMLParser):
        def __init__(self) -> None:
            super().__init__()
            self.links: list[dict[str, str | None]] = []

        def handle_starttag(
            self, tag: str, attrs: list[tuple[str, str | None]]
        ) -> None:
            if tag == "link":
                self.links.append(dict(attrs))

    parser = LinkParser()
    parser.feed(page.read_text())
    return parser.links


def test_discovery_metadata_uses_public_urls_and_local_export_links(
    rendered_discovery_files: tuple[Path, Path, str],
) -> None:
    from docs import build_docs

    html, api_text, base_url = rendered_discovery_files
    build_docs.write_discovery_files(html, api_text, revision="0123456789ab")

    public = "https://docs.nygen.io/CyteType"
    pages = [
        ("index.html", "/", "/index.txt", "text/plain"),
        (
            "quickstart/index.html",
            "/quickstart/",
            "/quickstart/index.txt",
            "text/plain",
        ),
        (
            "api/cytetype.html",
            "/api/cytetype.html",
            "/api/cytetype.txt",
            "text/plain",
        ),
    ]
    for page, route, export, media_type in pages:
        links = _discovery_links(html / page)
        canonicals = [link for link in links if link.get("rel") == "canonical"]
        assert len(canonicals) == 1
        assert canonicals[0]["href"] == f"{public}{route}"
        descriptions = [link for link in links if link.get("rel") == "describedby"]
        assert len(descriptions) == 1
        assert descriptions[0]["href"] == f"{base_url}/llms.txt"
        alternatives = [
            link
            for link in links
            if link.get("rel") == "alternate" and link.get("type") == media_type
        ]
        assert len(alternatives) == 1
        assert alternatives[0]["href"] == f"{base_url}{export}"
    assert "<h1>quickstart</h1>" in (html / "quickstart" / "index.html").read_text()

    sitemap = ElementTree.parse(html / "sitemap.xml")
    namespace = "{http://www.sitemaps.org/schemas/sitemap/0.9}"
    assert sitemap.getroot().tag == f"{namespace}urlset"
    locations = [
        element.text for element in sitemap.findall(f"{namespace}url/{namespace}loc")
    ]
    assert len(locations) == 3
    assert set(locations) == {f"{public}{route}" for _, route, _, _ in pages}
    robots = (html / "robots.txt").read_text()
    assert f"Sitemap: {public}/sitemap.xml" in robots
    assert "http://localhost" not in robots
    assert "User-agent: *" in robots
    assert "Allow: /" in robots


def test_discovery_exports_preserve_content_and_identify_revision(
    rendered_discovery_files: tuple[Path, Path, str],
) -> None:
    from docs import build_docs

    html, api_text, _ = rendered_discovery_files
    revision = "0123456789ab"
    build_docs.write_discovery_files(html, api_text, revision=revision)

    for relative in ("index.txt", "quickstart/index.txt", "api/cytetype.txt"):
        exported = (html / relative).read_text()
        assert revision in exported
        assert "{autoclass}" not in exported
    for relative in ("index.txt", "quickstart/index.txt"):
        exported = (html / relative).read_text()
        assert "Saved result: 42 cells." in exported
        assert "def count():\n    return 42\n\nprint(count())\n42" in exported
        assert all(value in exported for value in ("Group", "Cells", "T cells", "42"))
        assert "Site navigation" not in exported
        assert "Footer navigation" not in exported
    quickstart = (html / "quickstart" / "index.txt").read_text()
    assert "https://docs.nygen.io/CyteType/quickstart/details#reading" in quickstart
    assert "https://docs.nygen.io/CyteType/quickstart/plot.png" in quickstart
    assert "Cell embedding" in quickstart
    assert quickstart.count("cell_count\n    1000") == 1
    assert "Loading..." not in quickstart
    assert "LEAKED_OUTPUT" not in quickstart
    assert "Analysis continues after the saved output." in quickstart
    api = (html / "api" / "cytetype.txt").read_text()
    assert (
        "class cytetype.CyteType(adata, group_key, rank_key='rank_genes_groups')\n"
        "   Prepare an AnnData object for annotation." in api
    )
    llms = (html / "llms.txt").read_text()
    assert revision in llms
    links = dict(re.findall(r"\[([^\]]+)\]\(([^)]+)\)", llms))
    assert {
        label: links[label]
        for label in ("Home", "Quickstart (Python)", "CyteType class")
    } == {
        "Home": "https://docs.nygen.io/CyteType/index.txt",
        "Quickstart (Python)": "https://docs.nygen.io/CyteType/quickstart/index.txt",
        "CyteType class": "https://docs.nygen.io/CyteType/api/cytetype.txt",
    }


def test_discovery_rejects_missing_api_text(
    rendered_discovery_files: tuple[Path, Path, str],
) -> None:
    from docs import build_docs

    html, api_text, _ = rendered_discovery_files
    path = api_text / "cytetype.txt"
    path.unlink()

    with pytest.raises(FileNotFoundError, match=re.escape(path.name)):
        build_docs.write_discovery_files(html, api_text, revision="0123456789ab")


@pytest.mark.parametrize(
    "content", ["<p>Navigation without an article</p>", "<article> </article>"]
)
def test_discovery_rejects_missing_or_empty_article(
    rendered_discovery_files: tuple[Path, Path, str], content: str
) -> None:
    from docs import build_docs

    html, api_text, _ = rendered_discovery_files
    (html / "quickstart" / "index.html").write_text(
        f"<html><head></head><body>{content}</body></html>"
    )

    with pytest.raises(ValueError, match="quickstart"):
        build_docs.write_discovery_files(html, api_text, revision="0123456789ab")


def test_text_export_preserves_saved_markdown_instead_of_object_repr(
    rendered_discovery_files: tuple[Path, Path, str],
) -> None:
    from docs import build_docs

    html, api_text, _ = rendered_discovery_files
    document = json.loads((html / "quickstart.json").read_text())
    output = document["mdast"]["children"][0]["children"][0]["jupyter_data"]
    markdown = "| Dataset | Cells |\n| --- | ---: |\n| PBMC | 1000 |"
    output["data"] = {
        "text/plain": {"content": "<IPython.core.display.Markdown object>"},
        "text/markdown": {"content": markdown},
    }
    (html / "quickstart.json").write_text(json.dumps(document))

    build_docs.write_discovery_files(html, api_text, revision="0123456789ab")

    exported = (html / "quickstart" / "index.txt").read_text()
    assert markdown in exported
    assert "<IPython.core.display.Markdown object>" not in exported


def test_text_export_uses_job_progress_text_instead_of_object_repr(
    rendered_discovery_files: tuple[Path, Path, str],
) -> None:
    from docs import build_docs

    html, api_text, _ = rendered_discovery_files
    document = json.loads((html / "quickstart.json").read_text())
    outputs = document["mdast"]["children"][0]["children"]
    outputs.append(
        {
            "type": "block",
            "jupyter_data": {
                "output_type": "display_data",
                "data": {
                    "text/html": {
                        "content": "<pre style='margin: 0;'>[DONE] CyteType job "
                        "completed &amp; results &lt;stored&gt;.</pre>"
                    },
                    "text/plain": {"content": "<IPython.core.display.HTML object>"},
                },
            },
        }
    )
    (html / "quickstart.json").write_text(json.dumps(document))

    build_docs.write_discovery_files(html, api_text, revision="0123456789ab")

    exported = (html / "quickstart" / "index.txt").read_text()
    assert "[DONE] CyteType job completed & results <stored>." in exported
    assert "<IPython.core.display.HTML object>" not in exported
    assert "<pre" not in exported


@pytest.mark.parametrize("as_lines", [False, True])
def test_saved_html_tables_become_native_nodes_without_changing_output_data(
    as_lines: bool,
) -> None:
    from docs.render_outputs import render_outputs

    content = (
        '<div><style>.dataframe { text-align: right; }</style><table class="dataframe">'
        '<thead><tr><th></th><th colspan="2">Counts</th></tr></thead><tbody>'
        '<tr><th rowspan="2" valign="top">Donor A</th><td>1000</td><td>&lt;5</td></tr>'
        "<tr><td>1200</td><td>6</td></tr></tbody></table>"
        "<p>2 rows × 2 columns</p></div>"
    )
    output = {
        "output_type": "execute_result",
        "execution_count": 1,
        "metadata": {},
        "data": {
            "text/html": [content] if as_lines else content,
            "text/plain": "Counts\nDonor A 1000 <5\n        1200 6",
        },
    }
    node: dict[str, Any] = {
        "type": "output",
        "jupyter_data": deepcopy(output),
        "children": [],
    }

    render_outputs(node)

    assert node["type"] == "block"
    assert node["jupyter_data"] == output
    table, summary = node["children"]
    assert table["type"] == "table"
    assert len(table["children"]) == 3
    heading = table["children"][0]["children"][1]
    assert heading == {
        "type": "tableCell",
        "header": True,
        "colspan": 2,
        "children": [{"type": "text", "value": "Counts"}],
    }
    row = table["children"][1]["children"]
    assert row[0]["header"] is True
    assert row[0]["rowspan"] == 2
    assert row[1]["header"] is False
    assert row[1]["children"][0]["value"] == "1000"
    assert row[2]["children"][0]["value"] == "<5"
    assert summary == {
        "type": "paragraph",
        "children": [{"type": "text", "value": "2 rows × 2 columns"}],
    }


def test_saved_output_transform_preserves_streams_and_renders_timing() -> None:
    from docs.render_outputs import render_outputs

    timing = {
        "type": "output",
        "jupyter_data": {
            "output_type": "display_data",
            "metadata": {"cytearc_execution_timing": True},
            "data": {
                "text/html": "<small>Wall time: 1.234 s</small>",
                "text/plain": ["Wall time: ", "1.234 s"],
            },
        },
        "children": [],
    }
    stream = {
        "type": "output",
        "jupyter_data": {"output_type": "stream", "name": "stdout", "text": "Done"},
        "children": [],
    }
    tree: dict[str, Any] = {
        "type": "outputs",
        "children": [deepcopy(stream), deepcopy(timing)],
    }

    render_outputs(tree)

    assert tree["children"][0] == stream
    rendered = tree["children"][1]
    assert rendered["jupyter_data"] == timing["jupyter_data"]
    assert rendered["children"] == [
        {
            "type": "paragraph",
            "children": [{"type": "text", "value": "Wall time: 1.234 s"}],
        }
    ]


@pytest.mark.parametrize(
    ("content", "text"),
    [
        (
            "<pre style='margin: 0; white-space: pre-wrap; font-family: monospace;'>"
            "⠋ CyteType job running... 01:05 elapsed</pre>",
            "⠋ CyteType job running... 01:05 elapsed",
        ),
        ("<pre>[DONE] CyteType job completed.</pre>", "[DONE] CyteType job completed."),
        (
            "<pre>Cluster &lt;1&gt; &amp; &quot;2&quot; &#x27;3&#x27;</pre>",
            "Cluster <1> & \"2\" '3'",
        ),
        ("  \n<pre>\nFirst line\n  Indented\n</pre>\n", "First line\n  Indented\n"),
        (["<pre>Split ", "across lines</pre>"], "Split across lines"),
        ("<PRE>Upper case</PRE>", "Upper case"),
        ("<pre></pre>", ""),
    ],
)
def test_single_pre_output_becomes_a_text_code_block(
    content: str | list[str], text: str
) -> None:
    from docs.render_outputs import render_outputs

    output = {
        "output_type": "display_data",
        "metadata": {},
        "data": {
            "text/html": content,
            "text/plain": "<IPython.core.display.HTML object>",
        },
    }
    node: dict[str, Any] = {
        "type": "output",
        "jupyter_data": deepcopy(output),
        "children": [],
    }

    render_outputs(node)

    assert node["type"] == "block"
    assert node["jupyter_data"] == output
    assert node["children"] == [{"type": "code", "lang": "text", "value": text}]


def test_cytetype_progress_display_renders_as_its_message() -> None:
    pytest.importorskip("IPython")
    from cytetype.api.progress import _render_notebook_message
    from docs.render_outputs import render_outputs

    message = '[FAILED] CyteType job failed. <see "report"> & retry'
    html = _render_notebook_message(message)
    node: dict[str, Any] = {
        "type": "output",
        "jupyter_data": {
            "output_type": "display_data",
            "data": {"text/html": html.data, "text/plain": repr(html)},
        },
        "children": [],
    }

    render_outputs(node)

    assert node["children"] == [{"type": "code", "lang": "text", "value": message}]


@pytest.mark.parametrize(
    "content",
    [
        "<table><tr><td><script>unknown()</script></td></tr></table>",
        "<table><tr><td>Unclosed",
        '<table><tr><td rowspan="0">Invalid span</td></tr></table>',
        "<div>Unrecognized output</div>",
        "<pre>One</pre><pre>Two</pre>",
        "<pre>Nested <b>markup</b></pre>",
        "<pre>Text</pre> after",
        "Before <pre>Text</pre>",
        "<pre>Unclosed",
        "<pre>Text</pre><!-- comment -->",
        "<pre>Text</pre><br/>",
        "<pre>Table follows</pre><table><tr><td>1</td></tr></table>",
    ],
)
def test_saved_output_transform_rejects_unsupported_or_invalid_html(
    content: str,
) -> None:
    from docs.render_outputs import render_outputs

    node = {
        "type": "output",
        "jupyter_data": {"data": {"text/html": content}},
        "children": [],
    }

    with pytest.raises(ValueError):
        render_outputs(node)


def test_render_plugin_follows_the_myst_executable_protocol() -> None:
    plugin = _DOCS / "render_outputs.py"
    specification = subprocess.run(
        [sys.executable, str(plugin)], capture_output=True, text=True, check=True
    )
    assert json.loads(specification.stdout) == {
        "name": "Saved notebook outputs",
        "transforms": [{"name": "saved-outputs", "stage": "document"}],
    }
    tree = {
        "type": "root",
        "children": [
            {
                "type": "output",
                "jupyter_data": {
                    "output_type": "display_data",
                    "data": {"text/html": "<pre>[DONE] CyteType job completed.</pre>"},
                },
                "children": [],
            }
        ],
    }
    transformed = subprocess.run(
        [sys.executable, str(plugin), "--transform", "saved-outputs"],
        input=json.dumps(tree),
        capture_output=True,
        text=True,
        check=True,
    )
    assert json.loads(transformed.stdout)["children"][0]["children"] == [
        {"type": "code", "lang": "text", "value": "[DONE] CyteType job completed."}
    ]
    assert os.access(plugin, os.X_OK), "MyST runs the plugin as an executable"


def _git(repository: Path, *arguments: str, date: str | None = None) -> None:
    """Run git in a test repository, isolated from the user's configuration."""
    environment = {
        **os.environ,
        "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_NOSYSTEM": "1",
        "GIT_AUTHOR_NAME": "Documentation tests",
        "GIT_AUTHOR_EMAIL": "docs@example.org",
        "GIT_COMMITTER_NAME": "Documentation tests",
        "GIT_COMMITTER_EMAIL": "docs@example.org",
    }
    if date is not None:
        environment["GIT_AUTHOR_DATE"] = environment["GIT_COMMITTER_DATE"] = date
    subprocess.run(
        ["git", *arguments],
        cwd=repository,
        env=environment,
        capture_output=True,
        check=True,
    )


@pytest.fixture
def documentation_git_checkout(tmp_path: Path) -> Path:
    pytest.importorskip("jupytext")
    pytest.importorskip("nbformat")
    pytest.importorskip("yaml")
    repository = tmp_path / "checkout"
    repository.mkdir()
    _git(repository, "init", "--quiet")
    files = {
        "README.md": "# Example\n",
        "LICENSE": "Example license\n",
        "docs/source/reference/license.md": "# License\n",
        "docs/source/references.bib": "@article{paper, title={Paper}}\n",
    }
    # Each file has its own commit and date, so the latest one is unambiguous.
    for day, (relative, contents) in enumerate(files.items(), start=1):
        path = repository / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents)
        _git(repository, "add", relative)
        _git(
            repository,
            "commit",
            "--quiet",
            "--message",
            f"Add {relative}",
            date=f"2026-10-0{day}T12:00:00+02:00",
        )
    return repository


def test_last_modified_uses_source_and_include_history(
    documentation_git_checkout: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from docs import build_docs

    repository = documentation_git_checkout
    monkeypatch.setattr(build_docs, "SOURCE", repository / "docs/source")
    monkeypatch.setattr(build_docs, "SITE", repository / "docs/build/site")
    bibliography = repository / "docs/source/references.bib"
    source = repository / "docs/source/reference/license.md"
    dependencies = build_docs.myst_source_files(
        {
            "location": "/reference/license.md",
            "dependencies": [],
            "mdast": {
                "type": "root",
                "children": [
                    {"type": "include", "file": "../references.bib"},
                    {
                        "type": "include",
                        "file": "../../../../README.md",
                        "children": [{"type": "include", "file": "LICENSE"}],
                    },
                    {"type": "cite", "label": "paper"},
                ],
            },
        },
        {bibliography},
    )
    assert dependencies == {
        source,
        bibliography,
        repository / "README.md",
        repository / "LICENSE",
    }
    dates = [
        subprocess.run(
            ["git", "log", "-1", "--format=%cI", "--", str(path)],
            cwd=repository,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        for path in dependencies
    ]
    expected = max(dates, key=datetime.fromisoformat)
    assert build_docs.latest_source_commit(repository, dependencies) == expected
    for path in dependencies:
        path.touch()
    assert build_docs.latest_source_commit(repository, dependencies) == expected


@pytest.mark.parametrize(
    "change", ["unstaged", "staged", "untracked", "deleted", "outside"]
)
def test_last_modified_omits_uncommitted_or_untracked_dependencies(
    documentation_git_checkout: Path, change: str
) -> None:
    from docs import build_docs

    repository = documentation_git_checkout
    source = repository / "LICENSE"
    dependencies = {source}
    if change in {"unstaged", "staged"}:
        source.write_text(source.read_text() + "\nUncommitted change\n")
        if change == "staged":
            _git(repository, "add", "LICENSE")
    elif change == "untracked":
        include = repository / "untracked.txt"
        include.write_text("Untracked included content\n")
        dependencies.add(include)
    elif change == "deleted":
        source.unlink()
    else:
        dependencies.add(repository.parent / "outside.txt")
    assert build_docs.latest_source_commit(repository, dependencies) is None


def test_last_modified_ignores_unrelated_working_tree_changes(
    documentation_git_checkout: Path,
) -> None:
    from docs import build_docs

    repository = documentation_git_checkout
    dependencies = {repository / "LICENSE"}
    original = build_docs.latest_source_commit(repository, dependencies)
    assert original is not None
    (repository / "README.md").write_text("Unrelated local change\n")
    assert build_docs.latest_source_commit(repository, dependencies) == original


def test_last_modified_tracks_lazy_api_implementations() -> None:
    addnodes = pytest.importorskip("sphinx.addnodes")
    pytest.importorskip("jupytext")
    pytest.importorskip("nbformat")
    pytest.importorskip("yaml")
    from docs import build_docs

    doctree = addnodes.desc()
    for module, fullname, kind in (
        ("cytetype", "marker_dotplot", "function"),
        ("cytetype", "CyteType.run", "method"),
        ("cytetype", "CyteTypeError.message", "attribute"),
        ("cytetype", "CyteTypeTimeoutError", "exception"),
        ("cytetype.api", "LLMModelConfig.provider", "attribute"),
        ("cytetype", "subsample_by_group", "function"),
    ):
        description = addnodes.desc(objtype=kind)
        description += addnodes.desc_signature(module=module, fullname=fullname)
        doctree += description
    paths = build_docs.api_source_files(doctree)
    package = Path(__file__).resolve().parents[1] / "cytetype"
    assert paths == {
        package / "plotting" / "dotplot.py",
        package / "main.py",
        package / "api" / "exceptions.py",
        package / "api" / "schemas.py",
        package / "preprocessing" / "subsampling.py",
    }


def test_sitemap_includes_only_known_last_modified_dates(
    rendered_discovery_files: tuple[Path, Path, str],
) -> None:
    from docs import build_docs

    html, api_text, _ = rendered_discovery_files
    date = "2026-10-08T20:02:08+02:00"
    build_docs.write_discovery_files(
        html, api_text, revision="0123456789ab", last_modified={"/quickstart/": date}
    )
    namespace = "{http://www.sitemaps.org/schemas/sitemap/0.9}"
    entries = ElementTree.parse(html / "sitemap.xml").findall(f"{namespace}url")
    dated = {
        entry.findtext(f"{namespace}loc"): entry.findtext(f"{namespace}lastmod")
        for entry in entries
        if entry.find(f"{namespace}lastmod") is not None
    }
    assert dated == {"https://docs.nygen.io/CyteType/quickstart/": date}


def test_staged_config_paths_are_relative_to_the_site(executable_docs: Any) -> None:
    yaml = pytest.importorskip("yaml")

    # The theme embeds the project configuration in every page, so absolute paths
    # would publish the paths of the machine that built the site.
    build_docs = executable_docs
    (build_docs.SOURCE / "llms.txt").write_text("# Example\n")
    config_path = build_docs.SOURCE / "myst.yml"
    config = yaml.safe_load(config_path.read_text())
    config["project"]["static_files"] = ["llms.txt", "../build/api"]
    config["project"]["plugins"] = [
        {"type": "executable", "path": "../render_outputs.py"}
    ]
    config_path.write_text(yaml.safe_dump(config))

    site = build_docs.prepare_site()

    staged = yaml.safe_load((site / "myst.yml").read_text())["project"]
    paths = [*staged["static_files"], staged["plugins"][0]["path"]]
    originals = ["llms.txt", "../build/api", "../render_outputs.py"]
    for path, original in zip(paths, originals, strict=True):
        assert not Path(path).is_absolute(), path
        assert (site / path).resolve() == (build_docs.SOURCE / original).resolve()


@pytest.mark.parametrize("base_url", ["", "/CyteType"])
def test_api_pages_join_the_staged_guide_navigation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, base_url: str
) -> None:
    pytest.importorskip("jupytext")
    pytest.importorskip("nbformat")
    yaml = pytest.importorskip("yaml")
    from docs import build_docs

    api = tmp_path / "docs" / "build" / "api"
    api.mkdir(parents=True)
    (api / "index.html").write_text(
        '<nav class="docs-api-tree"><ul><li>'
        '<a href="cytetype.html"><span>CyteType class</span></a>'
        '</li><li><a href="helpers.html">Helpers &amp; utilities</a>'
        "</li></ul></nav>"
    )
    site = tmp_path / "site"
    site.mkdir()
    config = {
        "project": {
            "toc": [
                {"file": "index.md"},
                {
                    "title": "Reference",
                    "children": [
                        {"url": f"{base_url}/api/index.html", "title": "API reference"}
                    ],
                },
            ]
        }
    }
    (site / "myst.yml").write_text(yaml.safe_dump(config))
    monkeypatch.setattr(build_docs, "DOCS", api.parents[1])
    monkeypatch.setenv("BASE_URL", base_url)

    build_docs.include_api_navigation(site)

    toc = yaml.safe_load((site / "myst.yml").read_text())["project"]["toc"]
    assert toc[0] == {"file": "index.md"}
    assert toc[1]["children"][0]["children"] == [
        {
            "title": "CyteType class",
            "url": f"{base_url}/api/cytetype.html",
            "open_in_same_tab": True,
        },
        {
            "title": "Helpers & utilities",
            "url": f"{base_url}/api/helpers.html",
            "open_in_same_tab": True,
        },
    ]


@pytest.mark.parametrize("base_url", ["", "/CyteType"])
def test_api_pages_keep_the_guide_header_navigation_and_existing_content(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, base_url: str
) -> None:
    pytest.importorskip("jupytext")
    pytest.importorskip("nbformat")
    pytest.importorskip("yaml")
    from docs import build_docs

    monkeypatch.setenv("BASE_URL", base_url)
    html = tmp_path / "html"
    api = html / "api"
    api.mkdir(parents=True)
    header = (
        '<nav class="myst-top-nav-bar"><div>'
        '<button class="myst-top-nav-menu-button flex">Menu</button>'
        f'<a class="myst-home-link" href="{base_url}/"><img src="logo.svg"/>CyteType</a>'
        '</div><button aria-expanded="false" data-state="closed" class="myst-search-bar flex myst-search-bar-disabled">Search'
        "<script>window.headerScript = true;</script></button>"
        '<button class="myst-theme-button">Theme</button></nav>'
    )
    navigation = (
        '<nav class="myst-primary-sidebar-toc"><div class="myst-toc">'
        f'<a class="myst-toc-item myst-toc-item-exact font-bold p-2 my-1" aria-current="page" href="{base_url}/">Overview</a>'
        f'<a href="{base_url}/tutorials/remote-stores/">Remote stores</a>'
        f'<a href="{base_url}/api/index.html">API reference</a>'
        f'<a href="{base_url}/api/cytetype.html">CyteType</a>'
        "</div></nav>"
    )
    (html / "index.html").write_text(
        f"<html><head></head><body>{header}{navigation}</body></html>"
    )
    article = '<article><h1 id="cytetype">CyteType</h1><a href="cytetype.html#method">Method</a></article>'
    for name in ("index.html", "cytetype.html", "search.html", "genindex.html"):
        (api / name).write_text(
            '<html><head></head><body class="docs-api"><header>'
            '<nav class="docs-header-placeholder"></nav></header>'
            '<div class="myst-primary-sidebar-nav docs-api-navigation">'
            '<nav class="docs-api-tree"><ul><li>Original API menu</li></ul></nav>'
            f"</div>{article}</body></html>"
        )
    (api / "objects.inv").write_bytes(b"existing Sphinx inventory")
    static = api / "_static"
    static.mkdir()
    (static / "webpack-macros.html").write_text("{% macro assets() %}{% endmacro %}")

    build_docs.integrate_api_navigation(html)

    for page in api.glob("*.html"):
        contents = page.read_text()
        copied_header = build_docs.html_fragment(contents, "myst-top-nav-bar")
        copied_navigation = build_docs.html_fragment(
            contents, "myst-primary-sidebar-toc"
        )
        assert 'class="myst-search-bar search-button__button"' in copied_header
        assert (
            'class="myst-top-nav-menu-button sidebar-toggle primary-toggle"'
            in copied_header
        )
        assert "myst-search-bar-disabled" not in copied_header
        assert 'aria-label="Search documentation"' in copied_header
        assert 'aria-controls="pst-search-dialog"' in copied_header
        assert "aria-expanded" not in copied_header
        assert "data-state" not in copied_header
        assert "<script" not in copied_header
        assert "Original API menu" not in contents
        assert "font-bold" not in copied_navigation
        assert "myst-toc-item-exact" not in copied_navigation
        assert "p-2" not in copied_navigation
        assert "my-1" not in copied_navigation
        assert re.findall(r'href="([^"]+)"', copied_navigation) == re.findall(
            r'href="([^"]+)"', navigation
        )
        selected = re.findall(r'<a\b[^>]*aria-current="page"[^>]*>', copied_navigation)
        if page.name in {"index.html", "cytetype.html"}:
            assert len(selected) == 1
            assert f'href="{base_url}/api/{page.name}"' in selected[0]
        else:
            assert not selected
        assert article in contents
        assert contents.count('src="https://docs.nygen.io/assets/brand/docs.js"') == 1
        assert "nygen-docs-sidebar:" in contents
        assert "fonts/inter/latin.woff2" in contents
    assert (api / "objects.inv").read_bytes() == b"existing Sphinx inventory"
    assert (
        static / "webpack-macros.html"
    ).read_text() == "{% macro assets() %}{% endmacro %}"


def test_html_fragment_rejects_missing_or_duplicate_navigation() -> None:
    pytest.importorskip("jupytext")
    pytest.importorskip("nbformat")
    pytest.importorskip("yaml")
    from docs import build_docs

    for contents in (
        "<html></html>",
        '<nav class="menu"></nav><nav class="menu"></nav>',
    ):
        with pytest.raises(ValueError, match="Expected one complete .menu element"):
            build_docs.html_fragment(contents, "menu")


@pytest.mark.parametrize(
    ("route", "title", "expected"),
    [
        ("/", "CyteType - CyteType", "CyteType documentation"),
        (
            "/api/index.html",
            "CyteType Python API reference &#8212; CyteType API reference",
            "CyteType Python API reference",
        ),
        (
            "/guides/example/",
            "Data &amp; labels - CyteType",
            "Data &amp; labels - CyteType",
        ),
    ],
)
def test_page_metadata_preserves_content_and_completes_share_tags(
    route: str, title: str, expected: str
) -> None:
    from docs import build_docs

    class Metadata(HTMLParser):
        def __init__(self) -> None:
            super().__init__()
            self.values: dict[str, list[str]] = {}

        def handle_starttag(
            self, tag: str, attrs: list[tuple[str, str | None]]
        ) -> None:
            attributes = dict(attrs)
            if tag == "meta":
                name = attributes.get("name") or attributes.get("property")
                self.values.setdefault(name, []).append(attributes.get("content"))

    canonical = build_docs.PUBLIC_URL + route
    body = "<body><h1>Keep this title</h1><p>Keep this content.</p></body></html>"
    source = (
        f"<html><head><title>{title}</title>"
        '<meta content="Compare &quot;labels&quot; &amp; cell types." name="description">'
        '<meta property="og:title" content="Old title">'
        '<meta property="og:url" content="https://example.invalid/">'
        '<meta name="twitter:card" content="summary_large_image">'
        "</head>" + body
    )

    result = build_docs.update_page_metadata(source, canonical)

    assert result.endswith(body)
    assert f"<title>{expected}</title>" in result
    assert build_docs.update_page_metadata(result, canonical) == result
    metadata = Metadata()
    metadata.feed(result)
    values = metadata.values
    assert values["og:title"] == [unescape(expected)]
    assert all(len(value) == 1 for value in values.values())
    assert (
        values["description"]
        == values["og:description"]
        == ['Compare "labels" & cell types.']
    )
    assert values["og:url"] == [canonical]
    assert values["og:type"] == ["website"]
    image = "https://docs.nygen.io/assets/brand/logo/cytetype-dark.png"
    assert values["og:image"] == values["twitter:image"] == [image]
    assert values["og:image:width"] == values["og:image:height"] == ["512"]
    assert values["og:image:type"] == ["image/png"]
    assert values["og:image:alt"] == values["twitter:image:alt"] == ["CyteType logo"]
    assert values["twitter:card"] == ["summary"]


@pytest.mark.parametrize(
    "head",
    [
        "<title>Missing description</title>",
        '<meta name="description" content="Missing title">',
        '<title>Empty description</title><meta name="description" content="">',
    ],
)
def test_page_metadata_rejects_missing_authored_metadata(head: str) -> None:
    from docs import build_docs

    with pytest.raises(ValueError, match="Expected a title and description"):
        build_docs.update_page_metadata(
            f"<html><head>{head}</head><body>Content</body></html>",
            build_docs.PUBLIC_URL + "/",
        )
