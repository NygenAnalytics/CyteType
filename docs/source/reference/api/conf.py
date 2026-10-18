import sys
from pathlib import Path

SOURCE_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SOURCE_ROOT.parents[1]))

import cytetype  # noqa: E402

# cytetype exports its TimeoutError as CyteTypeTimeoutError, so that star imports do not
# hide the builtin. Autodoc documents a class under another name as a bare alias, so give
# the class its exported name for this build. Its canonical name stays
# cytetype.api.exceptions.TimeoutError, which docstrings reference.
cytetype.CyteTypeTimeoutError.__name__ = "CyteTypeTimeoutError"

project = "CyteType API"
copyright = "2025-2026, Nygen Analytics AB"
author = "Nygen Analytics AB"
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
]
master_doc = "index"
exclude_patterns = ["_build", "**.ipynb_checkpoints"]
myst_enable_extensions = ["colon_fence"]
myst_all_links_external = True
pygments_style = "sphinx"
language = "en"
# Docstrings mark code with single backticks, as in Markdown.
default_role = "literal"

# Show each class's __init__ arguments with its description.
autoclass_content = "both"
autodoc_member_order = "bysource"
autodoc_inherit_docstrings = False
# Render "str, optional" in Google-style argument types as text, not as references.
napoleon_preprocess_types = True
napoleon_use_rtype = False
typehints_defaults = "comma"
# Write unions as "str | None", as the code does, instead of Optional[str].
always_use_bars_union = True

html_theme = "sphinx_book_theme"
html_title = "CyteType API reference"
html_baseurl = "https://docs.nygen.io/CyteType/api/"
html_theme_options = {
    "repository_url": "https://github.com/NygenAnalytics/CyteType",
    "path_to_docs": "docs/source/reference/api",
    "use_repository_button": True,
    "show_navbar_depth": 2,
    "navigation_with_keys": False,
}

# The reference names these types without linking to their documentation.
nitpick_ignore = [
    ("py:class", "anndata.AnnData"),
    ("py:class", "anndata._core.anndata.AnnData"),
    ("py:class", "pandas.DataFrame"),
    ("py:data", "typing.Any"),
    ("py:data", "typing.Literal"),
    ("py:exc", "pydantic.ValidationError"),
]
