#!/usr/bin/env -S uv run --no-sync python
"""Render saved notebook HTML outputs as native MyST nodes."""

from html.parser import HTMLParser
import json
import sys


class OutputTable(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.root = {"type": "root", "children": []}
        self.stack = [("root", self.root)]
        self.style = False

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "style":
            self.style = True
            return
        if tag in {"div", "thead", "tbody"}:
            return
        kinds = {
            "table": "table",
            "tr": "tableRow",
            "th": "tableCell",
            "td": "tableCell",
            "p": "paragraph",
        }
        if tag not in kinds:
            raise ValueError(f"Unsupported saved table element: {tag}")
        node = {"type": kinds[tag], "children": []}
        if tag in {"th", "td"}:
            node["header"] = tag == "th"
            attributes = dict(attrs)
            for name in ("rowspan", "colspan"):
                if name in attributes:
                    value = int(attributes[name])
                    if value < 1:
                        raise ValueError(f"Invalid table {name}: {value}")
                    node[name] = value
        self.stack[-1][1]["children"].append(node)
        self.stack.append((tag, node))

    def handle_endtag(self, tag: str) -> None:
        if tag == "style":
            self.style = False
            return
        if tag in {"div", "thead", "tbody"}:
            return
        if len(self.stack) == 1 or self.stack[-1][0] != tag:
            raise ValueError(f"Unexpected saved table closing element: {tag}")
        self.stack.pop()

    def handle_data(self, data: str) -> None:
        if self.style or not data.strip():
            return
        if self.stack[-1][0] not in {"th", "td", "p"}:
            raise ValueError("Saved table text is outside a cell or summary")
        self.stack[-1][1]["children"].append({"type": "text", "value": data.strip()})

    def close(self) -> None:
        super().close()
        if len(self.stack) != 1 or self.style:
            raise ValueError("Unclosed saved table element")
        if not any(node["type"] == "table" for node in self.root["children"]):
            raise ValueError("Saved HTML output has no table")


class Preformatted(HTMLParser):
    """Collect the text of HTML that is exactly one ``<pre>`` element with only text.

    CyteType shows job progress in notebooks as ``<pre>`` HTML that it updates in place.
    """

    def __init__(self) -> None:
        super().__init__()
        self.state = "before"
        self.parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.state = "inside" if self.state == "before" and tag == "pre" else "other"

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.state = "other"

    def handle_endtag(self, tag: str) -> None:
        self.state = "after" if self.state == "inside" and tag == "pre" else "other"

    def handle_data(self, data: str) -> None:
        if self.state == "inside":
            self.parts.append(data)
        elif data.strip():
            self.state = "other"

    def handle_comment(self, data: str) -> None:
        self.state = "other"

    def handle_decl(self, decl: str) -> None:
        self.state = "other"

    def handle_pi(self, data: str) -> None:
        self.state = "other"

    def unknown_decl(self, data: str) -> None:
        self.state = "other"


def preformatted_text(html: str | list[str]) -> str | None:
    """Return the unescaped text of a single ``<pre>`` element, or None otherwise."""
    parser = Preformatted()
    parser.feed("".join(html))
    parser.close()
    if parser.state != "after":
        return None
    text = "".join(parser.parts)
    # Browsers ignore a newline that directly follows the <pre> start tag.
    return text.removeprefix("\n")


def render_outputs(node: dict) -> None:
    for child in node.get("children", []):
        render_outputs(child)
    if node.get("type") != "output":
        return
    output = node["jupyter_data"]
    data = output.get("data", {})
    if "text/html" not in data:
        return
    if output.get("metadata", {}).get("cytearc_execution_timing"):
        text = data["text/plain"]
        children = [
            {
                "type": "paragraph",
                "children": [{"type": "text", "value": "".join(text)}],
            }
        ]
    elif (text := preformatted_text(data["text/html"])) is not None:
        children = [{"type": "code", "lang": "text", "value": text}]
    else:
        parser = OutputTable()
        parser.feed("".join(data["text/html"]))
        parser.close()
        children = parser.root["children"]
    node["type"] = "block"
    node["children"] = children


def main() -> None:
    if len(sys.argv) == 1:
        result = {
            "name": "Saved notebook outputs",
            "transforms": [{"name": "saved-outputs", "stage": "document"}],
        }
    elif sys.argv[1:] == ["--transform", "saved-outputs"]:
        result = json.load(sys.stdin)
        render_outputs(result)
    else:
        raise ValueError(f"Unknown plugin arguments: {sys.argv[1:]}")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
