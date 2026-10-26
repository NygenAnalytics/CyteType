"""Render current documentation with saved notebook outputs, without execution."""

import argparse
import hashlib
import inspect
import json
import os
import pickle
import re
import shutil
import subprocess
import tomllib
from datetime import datetime
from fnmatch import fnmatchcase
from html import escape
from html.parser import HTMLParser
from importlib import import_module
from pathlib import Path
from urllib.parse import urljoin
from xml.etree import ElementTree

import jupytext
import nbformat
import yaml

from docs.render_outputs import preformatted_text

DOCS = Path(__file__).resolve().parent
SOURCE = DOCS / "source"
SITE = DOCS / "build" / "site"
PUBLIC_URL = tomllib.loads((DOCS.parent / "pyproject.toml").read_text())["project"][
    "urls"
]["Documentation"].rstrip("/")


class ArticleText(HTMLParser):
    void_tags = {
        "area",
        "base",
        "br",
        "col",
        "embed",
        "hr",
        "img",
        "input",
        "link",
        "meta",
        "param",
        "source",
        "track",
        "wbr",
    }

    def __init__(self, canonical: str, base_url: str, outputs: dict[str, str]):
        super().__init__()
        self.canonical = canonical
        self.base_url = base_url
        self.outputs = outputs
        self.article = False
        self.started = False
        self.ignored = 0
        self.pre = False
        self.link = ""
        self.parts: list[str] = []
        self.lines: list[str] = []

    def absolute_url(self, url: str) -> str:
        if url.startswith("/") and not url.startswith("//"):
            if self.base_url and url.startswith(f"{self.base_url}/"):
                url = url[len(self.base_url) :]
            return f"{PUBLIC_URL}{url}"
        return urljoin(self.canonical, url)

    def flush(self) -> None:
        text = "".join(self.parts)
        if self.pre:
            text = text.strip("\n")
        else:
            text = re.sub(r"[^\S\t]+", " ", text).strip()
        if text.strip():
            self.lines.append(text)
        self.parts.clear()

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = dict(attrs)
        if tag == "article":
            self.article = True
        if not self.article:
            return
        if self.ignored:
            if tag not in self.void_tags:
                self.ignored += 1
            return
        output_id = attributes.get("data-mdast-node-id")
        if output_id in self.outputs:
            self.flush()
            self.lines.append(self.outputs.pop(output_id))
            self.ignored = 1
            return
        if (
            tag in {"script", "style", "nav", "button", "svg"}
            or "myst-footer-links" in (attributes.get("class") or "").split()
            or attributes.get("aria-label") == "Link to this Section"
        ):
            self.ignored = 1
            return
        if tag == "h1":
            self.started = True
        if not self.started:
            return
        if tag in {
            "p",
            "div",
            "li",
            "dt",
            "dd",
            "tr",
            "h1",
            "h2",
            "h3",
            "h4",
            "h5",
            "h6",
            "pre",
            "figure",
        }:
            self.flush()
        if tag == "pre":
            self.pre = True
        if tag == "a":
            self.link = attributes.get("href") or ""
        if tag == "img":
            self.parts.append(
                f" [Image: {attributes.get('alt') or 'Figure'} "
                f"<{self.absolute_url(attributes['src'])}>] "
            )
        if tag == "br":
            self.parts.append("\n" if self.pre else " ")

    def handle_endtag(self, tag: str) -> None:
        if not self.article or tag in self.void_tags:
            return
        if self.ignored:
            self.ignored -= 1
            return
        if tag == "a" and self.link:
            self.parts.append(f" <{self.absolute_url(self.link)}>")
            self.link = ""
        if tag in {"td", "th"}:
            self.parts.append("\t")
        if tag in {
            "article",
            "p",
            "div",
            "li",
            "dt",
            "dd",
            "tr",
            "h1",
            "h2",
            "h3",
            "h4",
            "h5",
            "h6",
            "pre",
            "figure",
        }:
            self.flush()
        if tag == "pre":
            self.pre = False
        if tag == "article":
            self.article = False

    def handle_data(self, data: str) -> None:
        if self.article and self.started and not self.ignored:
            self.parts.append(data)


def notebook_output_text(node: dict, parser: ArticleText) -> dict[str, str]:
    outputs = {}
    if node["type"] == "outputs":
        parts = []
        for child in node["children"]:
            output = child["jupyter_data"]
            if output.get("metadata", {}).get("cytearc_execution_timing"):
                continue
            if output["output_type"] == "stream":
                parts.append(output["text"])
            else:
                data = output["data"]
                # The render plugin shows a single <pre> element, such as CyteType's
                # job progress, as its text; text/plain only names the HTML object.
                html = data.get("text/html", {}).get("content")
                preformatted = None if html is None else preformatted_text(html)
                if preformatted is not None:
                    parts.append(preformatted)
                else:
                    text = data.get("text/markdown", data["text/plain"])
                    parts.append(text["content"])
                parts.extend(
                    f"Image: {parser.absolute_url(value['path'])}"
                    for mime, value in data.items()
                    if mime.startswith("image/")
                )
        outputs[node["id"]] = "\n\n".join(parts)
    else:
        for child in node.get("children", []):
            outputs.update(notebook_output_text(child, parser))
    return outputs


def latest_source_commit(repository: Path, dependencies: set[Path]) -> str | None:
    repository = repository.resolve()
    paths = {path.resolve() for path in dependencies}
    if not paths or any(not path.is_relative_to(repository) for path in paths):
        return None
    relative = sorted(path.relative_to(repository).as_posix() for path in paths)

    def git(*arguments: str, paths: list[str]) -> str:
        return subprocess.run(
            ["git", *arguments, "--", *paths],
            cwd=repository,
            capture_output=True,
            text=True,
            check=True,
        ).stdout

    tracked = set(git("ls-files", "-z", paths=relative).rstrip("\0").split("\0"))
    if tracked != set(relative) or git("status", "--porcelain", paths=relative):
        return None
    # Combined path histories can report merges that changed no individual file.
    dates = [
        git("log", "-1", "--format=%cI", paths=[path]).strip() for path in relative
    ]
    return max(dates, key=datetime.fromisoformat)


def api_source_files(doctree) -> set[Path]:
    paths = set()
    for node in doctree.findall():
        if node.tagname != "desc_signature" or not node.get("module"):
            continue
        obj = import_module(node["module"])
        names = node["fullname"].split(".")
        if node.parent["objtype"] in {"attribute", "data"}:
            names = names[:-1]
        for name in names:
            obj = getattr(obj, name)
        if isinstance(obj, property):
            obj = obj.fget
        source = (
            obj.__file__
            if inspect.ismodule(obj)
            else inspect.getsourcefile(inspect.unwrap(obj))
        )
        if source is None:
            raise ValueError(f"No source file for documented symbol {node['fullname']}")
        paths.add(Path(source))
        if inspect.isclass(obj):
            for member in vars(obj).values():
                if isinstance(member, classmethod | staticmethod):
                    member = member.__func__
                elif isinstance(member, property):
                    member = member.fget
                if inspect.isfunction(member):
                    source = inspect.unwrap(member).__globals__.get("__file__")
                    if source and Path(source).is_relative_to(DOCS.parent):
                        paths.add(Path(source))
    return paths


def myst_source_files(document: dict, bibliography: set[Path]) -> set[Path]:
    source = Path(document["location"].lstrip("/"))
    paths = {SOURCE / source}

    def include_files(node: dict, containing_file: Path) -> None:
        if node["type"] == "include":
            containing_file = (containing_file.parent / node["file"]).resolve()
            paths.add(
                SOURCE / containing_file.relative_to(SITE)
                if containing_file.is_relative_to(SITE)
                else containing_file
            )
        if node["type"] in {"cite", "citeGroup"}:
            paths.update(bibliography)
        for child in node.get("children", []):
            include_files(child, containing_file)

    include_files(document["mdast"], SITE / source)
    return paths


def sitemap_last_modified(html: Path, doctrees: Path) -> dict[str, str]:
    references = json.loads((html / "myst.xref.json").read_text())["references"]
    pages = {
        reference["url"].rstrip("/") + "/": reference["data"]
        for reference in references
        if reference["kind"] == "page"
    }
    project = yaml.safe_load((SOURCE / "myst.yml").read_text())["project"]
    bibliography = {SOURCE / path for path in project.get("bibliography", [])}
    dependencies = {}
    for route, data_url in pages.items():
        document = json.loads((html / data_url.lstrip("/")).read_text())
        dependencies[route] = myst_source_files(document, bibliography)

    environment = pickle.loads((doctrees / "environment.pickle").read_bytes())
    for source in sorted((html / "api" / "_sources").glob("*.md")):
        name = source.stem
        paths = {Path(environment.doc2path(name))}
        paths.update(
            Path(environment.srcdir) / path for path in environment.dependencies[name]
        )
        paths.update(
            Path(environment.doc2path(child))
            for child in environment.toctree_includes.get(name, [])
        )
        paths.update(
            api_source_files(pickle.loads((doctrees / f"{name}.doctree").read_bytes()))
        )
        dependencies[f"/api/{name}.html"] = paths

    return {
        route: date
        for route, paths in dependencies.items()
        if (date := latest_source_commit(DOCS.parent, paths)) is not None
    }


def write_discovery_files(
    html: Path,
    api_text: Path,
    *,
    revision: str,
    last_modified: dict[str, str] | None = None,
) -> None:
    base_url = os.environ.get("BASE_URL", "").rstrip("/")
    references = json.loads((html / "myst.xref.json").read_text())["references"]
    pages = {}
    for reference in references:
        if reference["kind"] == "page":
            route = reference["url"].rstrip("/") + "/"
            if route in pages and pages[route] != reference["data"]:
                raise ValueError(f"Conflicting documentation routes: {route}")
            pages[route] = reference["data"]
    if not pages:
        raise ValueError(f"No documentation pages in {html / 'myst.xref.json'}")

    exports = {}
    for route, data_url in pages.items():
        page = html / route.lstrip("/") / "index.html"
        canonical = f"{PUBLIC_URL}{route}"
        document = json.loads((html / data_url.lstrip("/")).read_text())
        parser = ArticleText(canonical, base_url, {})
        parser.outputs = notebook_output_text(document["mdast"], parser)
        parser.feed(page.read_text())
        if parser.outputs:
            raise ValueError(
                f"Notebook outputs missing from {page}: {list(parser.outputs)}"
            )
        text = "\n\n".join(parser.lines).strip()
        if not text:
            raise ValueError(f"No article content in {page}")
        export = page.with_name("index.txt")
        export.write_text(
            f"Source: {canonical}\nRevision: {revision}\n"
            f"Structured page: {PUBLIC_URL}{data_url}\n\n{text}\n"
        )
        exports[route] = export.relative_to(html).as_posix()

    api_sources = sorted((html / "api" / "_sources").glob("*.md"))
    if not api_sources:
        raise ValueError(f"No API source pages in {html / 'api' / '_sources'}")
    for source in api_sources:
        route = f"/api/{source.stem}.html"
        export = html / "api" / f"{source.stem}.txt"
        text = (api_text / export.name).read_text()
        if not text.strip():
            raise ValueError(f"Empty API text: {api_text / export.name}")
        export.write_text(
            f"Source: {PUBLIC_URL}{route}\nRevision: {revision}\n\n{text}"
        )
        exports[route] = export.relative_to(html).as_posix()

    sitemap = ElementTree.Element(
        "urlset", xmlns="http://www.sitemaps.org/schemas/sitemap/0.9"
    )
    for route, export in sorted(exports.items()):
        canonical = f"{PUBLIC_URL}{route}"
        page = html / route.lstrip("/")
        if route.endswith("/"):
            page /= "index.html"
        contents = page.read_text()
        canonical_tag = (
            f'<link rel="canonical" href="{escape(canonical, quote=True)}"/>'
        )
        canonical_pattern = r"""<link\b(?=[^>]*\brel=["']canonical["'])[^>]*>"""
        contents, count = re.subn(canonical_pattern, canonical_tag, contents)
        links = (
            f'<link rel="describedby" href="{base_url}/llms.txt"/>'
            f'<link rel="alternate" type="text/plain" href="{base_url}/{export}"/>'
        )
        if route.endswith("/"):
            script_url = "https://docs.nygen.io/assets/brand/docs.js"
            links += f'<script defer src="{script_url}"></script>'
        if not count:
            links = canonical_tag + links
        if contents.count("</head>") != 1:
            raise ValueError(f"Expected one HTML head in {page}")
        page.write_text(contents.replace("</head>", links + "</head>", 1))
        entry = ElementTree.SubElement(sitemap, "url")
        ElementTree.SubElement(entry, "loc").text = canonical
        if last_modified and route in last_modified:
            ElementTree.SubElement(entry, "lastmod").text = last_modified[route]
    ElementTree.indent(sitemap)
    ElementTree.ElementTree(sitemap).write(
        html / "sitemap.xml", encoding="utf-8", xml_declaration=True
    )
    (html / "robots.txt").write_text(
        f"User-agent: *\nAllow: /\n\nSitemap: {PUBLIC_URL}/sitemap.xml\n"
    )

    def export_link(match: re.Match[str]) -> str:
        route = "/" if match[1] == "./" else f"/{match[1]}"
        return f"]({PUBLIC_URL}/{exports[route]})"

    index = re.sub(r"\]\(([^)]+)\)", export_link, (html / "llms.txt").read_text())
    index = re.sub(
        r"(?m)^(?=## |\- \[)",
        f"Documentation revision: {revision}. Links below open plain-text documentation.\n\n",
        index,
        count=1,
    )
    (html / "llms.txt").write_text(index)


def code_fingerprint(notebook: nbformat.NotebookNode) -> str:
    sources = [cell.source for cell in notebook.cells if cell.cell_type == "code"]
    encoded = json.dumps(sources, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def validate_execution(notebook: nbformat.NotebookNode, page: Path) -> None:
    nbformat.validate(notebook)
    if not notebook.metadata.get("executed_at"):
        raise ValueError(f"Notebook {page} has no completed execution timestamp")
    fingerprint = notebook.metadata.get("executed_code_sha256")
    if not fingerprint:
        raise ValueError(f"Notebook {page} has no executed-code fingerprint")
    if fingerprint != code_fingerprint(notebook):
        raise ValueError(
            f"Code changed in {page} since execution. Rerun:\n"
            f'make -C docs execute-docs PAGES="{page.with_suffix("").as_posix()}"'
        )
    for index, cell in enumerate(notebook.cells, start=1):
        if cell.cell_type != "code":
            continue
        if any(output.output_type == "error" for output in cell.get("outputs", [])):
            raise ValueError(f"Notebook {page} has an error in cell {index}")
        count = cell.get("execution_count")
        if cell.source.strip() and (type(count) is not int or count < 1):
            raise ValueError(f"Notebook {page} has an unexecuted cell {index}")


def rewrite_includes(text: str, page: Path) -> str:
    def include_path(match: re.Match[str]) -> str:
        original = (SOURCE / page.parent / match[2]).resolve()
        return f"{match[1]}{original.relative_to(SITE / page.parent, walk_up=True)}"

    return re.sub(
        r"(?m)^((?:`{3,}|:{3,})\{include\}[ \t]+)(\S+)",
        include_path,
        text,
    )


def rewrite_api_links(text: str, page: Path) -> str:
    """Point relative links to the Sphinx API sources at their pages under /api/.

    Sphinx builds ``reference/api/<name>.md`` to ``/api/<name>.html``, beside the MyST
    routes, so ``../reference/api/cytetype.html#cytetype.CyteType.run`` (or ``.md``)
    becomes ``/api/cytetype.html#cytetype.CyteType.run``. MyST adds the base URL.
    """
    api = (SOURCE / "reference" / "api").resolve()

    def api_link(match: re.Match[str]) -> str:
        target, hash_sign, fragment = match[1].partition("#")
        if not target or "://" in target or target.startswith(("/", "mailto:")):
            return match[0]
        path = (SOURCE / page.parent / target).resolve()
        if path.parent != api or path.suffix not in {".html", ".md"}:
            return match[0]
        return f"](/api/{path.stem}.html{hash_sign}{fragment})"

    return re.sub(r"\]\(([^)\s]+)\)", api_link, text)


def prepare_site(*, preview: bool = False) -> Path:
    config = yaml.safe_load((SOURCE / "myst.yml").read_text())
    excludes = config["project"].get("exclude", [])

    def excluded(path: Path) -> bool:
        relative = path.as_posix()
        return any(
            fnmatchcase(relative, pattern)
            or (pattern.startswith("**/") and fnmatchcase(relative, pattern[3:]))
            for pattern in excludes
        )

    def ignore(directory: str, names: list[str]) -> set[str]:
        relative = Path(directory).relative_to(SOURCE)
        return shutil.ignore_patterns("_build", "__pycache__", ".ipynb_checkpoints")(
            directory, names
        ) | {name for name in names if excluded(relative / name)}

    if SITE.exists():
        shutil.rmtree(SITE)
    shutil.copytree(
        SOURCE,
        SITE,
        ignore=ignore,
    )
    pages = {
        path.relative_to(SITE): (
            nbformat.read(path, as_version=4)
            if path.suffix == ".ipynb"
            else jupytext.read(path, fmt="md:myst")
        )
        for path in sorted(SITE.rglob("*"))
        if path.suffix in {".md", ".ipynb"}
        and "_build" not in path.relative_to(SITE).parts
    }
    executable = {
        page
        for page, notebook in pages.items()
        if any(cell.cell_type == "code" for cell in notebook.cells)
    }
    for page in executable:
        if page.suffix != ".ipynb":
            raise ValueError(
                f"Executable Markdown page {page} must be an .ipynb notebook. "
                "Keep its code and outputs in one source file."
            )
        if not pages[page].metadata.get("kernelspec"):
            raise ValueError(
                f"Executable page {page} has no kernelspec metadata. "
                "Add the kernel configuration for its execution environment "
                "before building."
            )
    missing = []
    github = config["project"]["github"].rstrip("/")
    for page, notebook in pages.items():
        links = {
            key: f"{github}/{action}/master/docs/source/{page.as_posix()}"
            for key, action in (("source_url", "blob"), ("edit_url", "edit"))
        }
        staged = SITE / page
        if page.suffix == ".ipynb":
            nbformat.validate(notebook)
            code = [cell for cell in notebook.cells if cell.cell_type == "code"]
            pristine = (
                "executed_at" not in notebook.metadata
                and "executed_code_sha256" not in notebook.metadata
                and all(
                    cell.execution_count is None and not cell.outputs for cell in code
                )
            )
            if code and pristine:
                missing.append(page.as_posix())
                if preview:
                    notebook.cells.insert(
                        0,
                        nbformat.v4.new_markdown_cell(
                            ":::{warning} Source preview\n"
                            "This page has not been executed. Figures and results "
                            "are not available.\n:::"
                        ),
                    )
            elif code or "executed_code_sha256" in notebook.metadata:
                validate_execution(notebook, page)
            for cell in notebook.cells:
                if cell.cell_type == "markdown":
                    cell.source = rewrite_api_links(
                        rewrite_includes(cell.source, page), page
                    )
                if preview and pristine:
                    continue
                figures = [
                    output["data"]
                    for output in cell.get("outputs", [])
                    if any(mime.startswith("image/") for mime in output.get("data", {}))
                ]
                descriptions = cell.metadata.get("figure_descriptions", [])
                if (
                    not isinstance(descriptions, list)
                    or len(descriptions) != len(figures)
                    or any(
                        not isinstance(text, str) or not text.strip()
                        for text in descriptions
                    )
                ):
                    raise ValueError(
                        f"Figures in {page}, cell {cell.id}, need "
                        f"{len(figures)} nonempty metadata.figure_descriptions entries"
                    )
                for data, description in zip(figures, descriptions, strict=True):
                    data["text/plain"] = description
            notebook.metadata.pop("jupytext", None)
            notebook.metadata.update(links)
            nbformat.write(notebook, staged)
        else:
            text = rewrite_api_links(rewrite_includes(staged.read_text(), page), page)
            frontmatter = re.match(r"\A---\n(.*?)\n---\n", text, re.DOTALL)
            metadata = yaml.safe_load(frontmatter[1]) or {} if frontmatter else {}
            metadata.update(links)
            body = text[frontmatter.end() :] if frontmatter else text
            staged.write_text("---\n" + yaml.safe_dump(metadata) + "---\n" + body)

    base_url = os.environ.get("BASE_URL", "").rstrip("/")

    def update_toc(entries: list[dict]) -> None:
        entries[:] = [
            entry
            for entry in entries
            if "file" not in entry or not excluded(Path(entry["file"]))
        ]
        for entry in entries:
            url = entry.get("url", "")
            if url.startswith("/") and not url.startswith("//"):
                entry["url"] = f"{base_url}{url}"
            update_toc(entry.get("children", []))

    update_toc(config["project"]["toc"])

    # Point paths at the originals relative to the staged site: the theme embeds the
    # project configuration in every page, so absolute paths would publish build paths.
    def staged_path(path: str) -> str:
        return os.path.relpath((SOURCE / path).resolve(), SITE.resolve())

    config["project"]["static_files"] = [
        staged_path(path) for path in config["project"].get("static_files", [])
    ]
    for plugin in config["project"].get("plugins", []):
        plugin["path"] = staged_path(plugin["path"])
    (SITE / "myst.yml").write_text(yaml.safe_dump(config, sort_keys=False))
    if missing:
        if not preview:
            raise ValueError(
                "Saved execution is missing for:\n  "
                + "\n  ".join(missing)
                + "\nExecute these pages before building HTML, or use --preview "
                "for a labelled source preview."
            )
        print(
            "Unexecuted notebooks:\n  " + "\n  ".join(missing),
            flush=True,
        )
    return SITE


def link_shared_brand_assets(html: Path) -> None:
    source_options = yaml.safe_load((SOURCE / "myst.yml").read_text())["site"]["options"]
    config = json.loads((html / "config.json").read_text())
    urls = {
        config["options"][key]: source_options[key]
        for key in ("logo", "logo_dark", "favicon", "style")
    }
    nodes = [config["parts"]["navbar_end"]["mdast"]]
    symbols = []
    while nodes:
        node = nodes.pop()
        if node.get("type") == "image" and node.get("class") == "brand-symbol":
            symbols.append(node["url"])
        nodes.extend(node.get("children", []))
    if len(symbols) != 1:
        raise ValueError(f"Expected one Nygen symbol in {html / 'config.json'}")
    urls[symbols[0]] = "https://docs.nygen.io/assets/brand/logo/nygen-symbol.svg"
    base_url = os.environ.get("BASE_URL", "").rstrip("/")
    for page in html.rglob("*"):
        if page.suffix not in {".html", ".json"}:
            continue
        contents = page.read_text()
        for local, public in urls.items():
            contents = contents.replace(f"{base_url}/build{local}", public)
            contents = contents.replace(local, public)
        for name, public in (
            ("favicon.ico", source_options["favicon"]),
            ("myst-theme.css", source_options["style"]),
        ):
            contents = contents.replace(f"{base_url}/{name}", public)
        page.write_text(contents)
    roots = list((html / "build").glob("root-*.js"))
    if len(roots) != 1:
        raise ValueError(f"Expected one generated theme root in {html / 'build'}")
    contents = roots[0].read_text()
    for name, public in (
        ("favicon.ico", source_options["favicon"]),
        ("myst-theme.css", source_options["style"]),
    ):
        pattern = r'href:`\$\{[\w$]+\|\|""\}/' + re.escape(name) + r'`'
        contents, count = re.subn(pattern, f"href:{json.dumps(public)}", contents)
        if count != 1:
            raise ValueError(f"Expected one {name} link in {roots[0]}")
    roots[0].write_text(contents)
    renames = {
        roots[0].name: f"root-{hashlib.sha256(contents.encode()).hexdigest()[:12]}.js"
    }
    manifests = list((html / "build").glob("manifest-*.js"))
    if len(manifests) != 1:
        raise ValueError(f"Expected one generated manifest in {html / 'build'}")
    contents = manifests[0].read_text().replace(roots[0].name, renames[roots[0].name])
    renames[manifests[0].name] = (
        f"manifest-{hashlib.sha256(contents.encode()).hexdigest()[:12]}.js"
    )
    for page in html.rglob("*"):
        if page.suffix not in {".html", ".json", ".js"}:
            continue
        contents = page.read_text()
        updated = contents
        for old, new in renames.items():
            updated = updated.replace(old, new)
        if updated != contents:
            page.write_text(updated)
    for old, new in renames.items():
        (html / "build" / old).rename(html / "build" / new)
    for local in urls:
        (html / "build" / local.lstrip("/")).unlink()
    (html / "favicon.ico").unlink()
    (html / "myst-theme.css").unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--preview",
        action="store_true",
        help="Allow missing executions and label those pages as source previews",
    )
    args = parser.parse_args()
    subprocess.run(
        ["myst", "build", "--html", "--strict"],
        cwd=prepare_site(preview=args.preview),
        check=True,
    )
    revision = subprocess.run(
        ["git", "describe", "--always", "--dirty", "--abbrev=40"],
        cwd=DOCS.parent,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    write_discovery_files(
        SITE / "_build" / "html",
        DOCS / "build" / "api-text",
        revision=revision,
        last_modified=sitemap_last_modified(
            SITE / "_build" / "html", DOCS / "build" / "api" / ".doctrees"
        ),
    )
    link_shared_brand_assets(SITE / "_build" / "html")


if __name__ == "__main__":
    main()
