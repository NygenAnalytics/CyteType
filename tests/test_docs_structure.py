"""Structure checks for the documentation sources in docs/source."""

import io
import json
import re
import shutil
import subprocess
import tokenize
from pathlib import Path
from typing import Any, get_args

import pytest

yaml = pytest.importorskip("yaml")

_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
_DOCS_SOURCE = _REPOSITORY_ROOT / "docs" / "source"
_API_SOURCE = _DOCS_SOURCE / "reference" / "api"
_MARKDOWN_LINK = re.compile(r"\[[^\]]*\]\(([^)\s]+)\)")
_AUTODOC_TARGET = re.compile(
    r"(?m)^(?:\.\. (?:auto\w+|py:\w+)::|[`:]{3,}\{(?:auto\w+|py:\w+)\})[ \t]*([^\s(]+)"
)
_TAB_SET = re.compile(
    r"(?ms)^(?P<fence>:{4,}|`{4,})\{tab-set\}[^\n]*\n(?P<body>.*?)^(?P=fence)[ \t]*$"
)
_TAB_ITEM = re.compile(
    r"(?m)^(?::{3,}|`{3,})\{tab-item\}[ \t]*(?P<title>[^\n]*)\n"
    r"(?P<options>(?:[ \t]*:[\w-]+:[^\n]*\n)*)"
)
_CLIENT_TABS = {"Python": "python", "R": "r"}
_SECTIONS = [
    "Getting started",
    "Prepare and run",
    "Results",
    "Platform",
    "Reference",
    "Support",
]


def _project() -> dict[str, Any]:
    project: dict[str, Any] = yaml.safe_load((_DOCS_SOURCE / "myst.yml").read_text())[
        "project"
    ]
    return project


def _notebook_cells(path: Path) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = json.loads(path.read_text())["cells"]
    for cell in cells:
        if not isinstance(cell["source"], str):
            cell["source"] = "".join(cell["source"])
    return cells


def _source_text(path: Path) -> str:
    if path.suffix != ".ipynb":
        return path.read_text()
    return "\n".join(cell["source"] for cell in _notebook_cells(path))


def _pages() -> list[Path]:
    return sorted(
        path
        for path in _DOCS_SOURCE.rglob("*")
        if path.suffix in {".md", ".ipynb"}
        and not {"_build", ".ipynb_checkpoints"} & set(path.parts)
    )


def _api_pages() -> list[Path]:
    return sorted(_API_SOURCE.glob("*.md"))


def _guide_pages() -> list[Path]:
    """Return the MyST pages; Sphinx builds the pages in reference/api."""
    return [page for page in _pages() if not page.is_relative_to(_API_SOURCE)]


def _toc_entries(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    found = []
    for entry in entries:
        found.append(entry)
        found.extend(_toc_entries(entry.get("children", [])))
    return found


def _toc_files() -> set[str]:
    return {
        entry["file"] for entry in _toc_entries(_project()["toc"]) if "file" in entry
    }


def _slug(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")


def _route(page: Path) -> str:
    """Return the MyST route of a page, such as ``getting-started/quickstart``."""
    relative = page.relative_to(_DOCS_SOURCE).with_suffix("")
    if relative.as_posix() == "index":
        return ""
    return "/".join(_slug(part) for part in relative.parts)


def _route_source(route: str) -> Path:
    """Return the source page published at ``route``, such as ``api/cytetype.html``."""
    route = route.split("#", maxsplit=1)[0]
    if route.startswith("api/"):
        name = route.removeprefix("api/")
        assert name.endswith(".html"), route
        return _API_SOURCE / name.replace(".html", ".md")
    routes = {_route(page): page for page in _guide_pages()}
    route = "" if route in {"", "./"} else route.strip("/")
    assert route in routes, f"No documentation page publishes {route!r}"
    return routes[route]


def _link_source(page: Path, target: str) -> Path:
    """Return the documentation source that a relative link from ``page`` opens."""
    path = target.split("#", maxsplit=1)[0]
    if path.startswith("/"):
        return _route_source(path.lstrip("/"))
    if page.parent == _API_SOURCE:
        # Sphinx pages are served at /api/<page>.html, beside the MyST routes.
        if path.startswith("../"):
            return _route_source(path.removeprefix("../"))
        return _API_SOURCE / path.replace(".html", ".md")
    resolved = (page.parent / path).resolve()
    if resolved.parent == _API_SOURCE.resolve():
        # The build points links to the Sphinx sources at their /api/ pages.
        assert resolved.suffix in {".html", ".md"}, (page, target)
        return resolved.with_suffix(".md")
    assert path.endswith((".md", ".ipynb")), (page, target)
    return resolved


def _python_sources(path: Path) -> list[str]:
    if path.suffix != ".ipynb":
        return []
    return [
        cell["source"]
        for cell in _notebook_cells(path)
        if cell["cell_type"] == "code" and cell["source"].strip()
    ]


def test_llms_index_links_to_every_published_page() -> None:
    contents = (_DOCS_SOURCE / "llms.txt").read_text()
    targets = _MARKDOWN_LINK.findall(contents)

    assert targets
    assert contents.startswith("# CyteType\n\n> ")
    for target in targets:
        assert "://" not in target, target
        # The build maps these exact routes to the plain-text exports.
        assert target == "./" or target.endswith(("/", ".html")), target
        source = _route_source(target)
        assert source.resolve().is_relative_to(_DOCS_SOURCE.resolve())
        assert source.is_file(), target
    listed = {_route_source(target).resolve() for target in targets}
    published = {(_DOCS_SOURCE / page).resolve() for page in _toc_files()}
    published |= {page.resolve() for page in _api_pages()}
    assert published - listed == set(), "Pages missing from llms.txt"


def test_navigation_publishes_every_page() -> None:
    project = _project()
    toctree = project["toc"]
    sections = {
        entry["title"]: entry["children"] for entry in toctree if "children" in entry
    }

    assert toctree[0] == {"file": "index.md"}
    assert list(sections) == _SECTIONS
    assert [entry["file"] for entry in sections["Getting started"][:3]] == [
        "getting-started/installation.md",
        "getting-started/sign-in.md",
        "getting-started/quickstart.ipynb",
    ]
    assert sections["Reference"][0] == {
        "title": "Python API",
        "url": "/api/index.html",
        "open_in_same_tab": True,
    }
    issue_trackers = {entry["url"] for entry in sections["Support"] if "url" in entry}
    assert issue_trackers == {
        "https://github.com/NygenAnalytics/CyteType/issues",
        "https://github.com/NygenAnalytics/CyteTypeR/issues",
    }
    assert "llms.txt" in project["static_files"]
    assert "../build/api" not in project["static_files"]
    assert "reference/api/**" in project["exclude"]

    published = _toc_files()
    for file in published:
        assert (_DOCS_SOURCE / file).is_file(), f"{file} is in the TOC but missing"
    for page in _guide_pages():
        relative = page.relative_to(_DOCS_SOURCE).as_posix()
        assert relative in published, f"{relative} is not in the TOC"


def test_relative_links_resolve_to_documentation_pages() -> None:
    for page in _pages():
        for target in _MARKDOWN_LINK.findall(_source_text(page)):
            if "://" in target or target.startswith(("#", "mailto:")):
                continue
            source = _link_source(page, target)
            assert source.is_file(), (page, target)
            assert source.resolve().is_relative_to(_DOCS_SOURCE.resolve())


def test_pages_have_a_description_and_one_title() -> None:
    for page in _guide_pages():
        if page.suffix == ".ipynb":
            metadata = json.loads(page.read_text())["metadata"]
            description = metadata.get("description")
            text = "\n".join(
                cell["source"]
                for cell in _notebook_cells(page)
                if cell["cell_type"] == "markdown"
            )
        else:
            text = page.read_text()
            frontmatter = re.match(r"\A---\n(.*?)\n---\n", text, re.DOTALL)
            assert frontmatter, f"{page} has no front matter"
            description = (yaml.safe_load(frontmatter[1]) or {}).get("description")
            text = text[frontmatter.end() :]
        assert isinstance(description, str) and description.strip(), page
        text = re.sub(r"(?ms)^(`{3,}|~{3,}).*?^\1", "", text)
        assert len(re.findall(r"(?m)^# ", text)) == 1, f"{page} needs one H1 title"


def test_autodoc_targets_use_the_public_cytetype_names() -> None:
    import cytetype
    import cytetype.api

    # Document objects under the names users import them by, not their modules.
    public = {"cytetype": cytetype.__all__, "cytetype.api": cytetype.api.__all__}
    targets = []
    for page in _api_pages():
        for target in _AUTODOC_TARGET.findall(page.read_text()):
            module, _, name = target.rpartition(".")
            assert name in public.get(module, []), (page, target)
            targets.append(target)
    assert "cytetype.CyteType" in targets


def test_client_tabs_are_synced_with_python_first() -> None:
    for page in _pages():
        text = _source_text(page)
        tab_sets = list(_TAB_SET.finditer(text))
        assert len(tab_sets) == text.count("{tab-set}"), (page, "unparsed tab-set")
        items = 0
        for tab_set in tab_sets:
            titles = []
            for item in _TAB_ITEM.finditer(tab_set["body"]):
                title = item["title"].strip()
                options = dict(
                    re.findall(
                        r"(?m)^[ \t]*:([\w-]+):[ \t]*(.*?)[ \t]*$", item["options"]
                    )
                )
                sync = options.get("sync")
                if title in _CLIENT_TABS or sync in _CLIENT_TABS.values():
                    assert sync == _CLIENT_TABS.get(title), (page, title, sync)
                titles.append(title)
                items += 1
            if {"Python", "R"} & set(titles):
                assert "Python" in titles and "R" in titles, (page, titles)
                assert titles.index("Python") < titles.index("R"), (page, titles)
        assert items == text.count("{tab-item}"), (page, "unparsed tab-item")


def test_notebook_comments_explain_each_step_on_their_own_lines() -> None:
    notebooks = [page for page in _pages() if page.suffix == ".ipynb"]
    assert notebooks
    for notebook in notebooks:
        for index, source in enumerate(_python_sources(notebook), start=1):
            lines = source.splitlines()
            assert lines[0].startswith("# "), (notebook.name, index)
            tokens = tokenize.generate_tokens(io.StringIO(source).readline)
            for token in tokens:
                if token.type != tokenize.COMMENT:
                    continue
                row = token.start[0]
                assert token.line.strip().startswith("#"), (notebook.name, token.line)
                previous = lines[row - 2].strip() if row > 1 else ""
                assert row == 1 or not previous or previous.startswith("#"), (
                    notebook.name,
                    index,
                    token.line,
                )


def test_pages_contain_no_html_comments() -> None:
    # The site publishes page sources, page data and the search index, so comments leak.
    for page in _pages():
        assert "<!--" not in _source_text(page), page


def test_notebook_outputs_contain_no_local_paths() -> None:
    # Warnings can print paths of the machine that executed a page; outputs are published.
    local_path = re.compile(r"/home/|/Users/|/tmp/|[A-Za-z]:\\\\Users\\\\")
    for page in _pages():
        if page.suffix != ".ipynb":
            continue
        for index, cell in enumerate(_notebook_cells(page), start=1):
            outputs = json.dumps(cell.get("outputs", []))
            assert not local_path.search(outputs), (page, index)


def test_notebooks_contain_no_api_keys_or_bearer_tokens() -> None:
    # Notebooks are published with their outputs. Report only where a match is, so that
    # a failure does not print the key.
    secret = re.compile(r"cyt_[A-Za-z0-9_-]{16,}|Bearer ")
    found = []
    for page in _pages():
        if page.suffix != ".ipynb":
            continue
        notebook = json.loads(page.read_text())
        if secret.search(json.dumps(notebook["metadata"])):
            found.append((page.name, "metadata"))
        for index, cell in enumerate(notebook["cells"], start=1):
            if secret.search(json.dumps(cell)):
                found.append((page.name, f"cell {index}"))
    assert found == []


def test_llm_configuration_reference_matches_the_schema() -> None:
    from cytetype.api.schemas import AgentType, LLMModelConfig, LLMProvider

    text = (_API_SOURCE / "llm-config.md").read_text()
    documented = re.findall(
        r"(?m)^   \.\. py:attribute:: (\w+)\n((?:      :\w+:.*\n)*)", text
    )
    fields = LLMModelConfig.model_fields

    assert [name for name, _ in documented] == list(fields)
    for name, options in documented:
        field = fields[name]
        values = dict(re.findall(r":(\w+): (.*)", options))
        if field.is_required():
            assert "value" not in values, name
        else:
            default = field.get_default(call_default_factory=True)
            assert values.get("value") == repr(default), name
    for provider in get_args(LLMProvider):
        assert f"``{provider}``" in text, provider
    for agent in get_args(AgentType):
        assert f"``{agent}``" in text, agent


def test_documentation_notebooks_are_tracked() -> None:
    git = shutil.which("git")
    if git is None or not (_REPOSITORY_ROOT / ".git").exists():
        pytest.skip("Not a git checkout")
    notebooks = [
        page.relative_to(_REPOSITORY_ROOT).as_posix()
        for page in _pages()
        if page.suffix == ".ipynb"
    ]
    ignored = subprocess.run(
        [git, "check-ignore", "--no-index", *notebooks],
        cwd=_REPOSITORY_ROOT,
        capture_output=True,
        text=True,
    )
    assert ignored.stdout == "", f"Ignored notebooks: {ignored.stdout}"
