# Documentation instructions

These instructions cover documentation sources, notebook execution, and site builds. The
site documents the Python SDK (`cytetype`) and the R package (`CyteTypeR`) together and is
published at `https://docs.nygen.io/CyteType/`.

- Author executable pages as tracked `.ipynb` files in `docs/source/`. Each notebook is the
  page's sole source and contains its prose, code, outputs, and execution metadata. The
  executed page today is `getting-started/quickstart.ipynb`.
- Keep prose guides and API references in Markdown. Code in Markdown pages is not executed,
  so check every function, argument, default and column name against the source.
- Explain each logical step in notebook code with a comment. Put comments on their own lines,
  with a blank line before each comment unless it starts the code cell.
- Give every code cell that draws figures a `figure_descriptions` list in its cell metadata,
  with one nonempty description per figure, in output order. The build uses them as alt text
  and fails without them. Cell metadata is not part of the code fingerprint, so adding or
  editing descriptions needs no re-execution.
- MyST builds guides and notebooks. Sphinx builds only `docs/source/reference/api/`, which
  the site serves at `/api/`.
- Import the SDK as `import cytetype` or `from cytetype import ...`, and use `cytetype.*`
  autodoc targets, the names under which users import each object.
- Start every page with front matter holding one `description:` line, then a label and one
  H1 title. Add every page to the table of contents in `docs/source/myst.yml` and to
  `docs/source/llms.txt`, with a one-line summary.

## Python and R

Show client code for both clients in synced tabs, Python first, so that a reader's choice
carries across pages:

````markdown
::::{tab-set}
:::{tab-item} Python
:sync: python
```python
adata = annotator.run(study_context="...")
```
:::
:::{tab-item} R
:sync: r
```r
obj <- CyteTypeR(obj, prepped, study_context = "...")
```
:::
::::
````

The pages document CyteTypeR 0.6.0 and the Python client 0.21.0. When a client release is
planned that changes documented behaviour, keep documenting the released version and list each
affected passage under "Update for new CyteTypeR releases" in `docs/MAINTENANCE.md`; update the
listed passages when the release ships. Regenerate `reference/r-api.md` from CyteTypeR's help
pages with `Rscript docs/render_r_reference.R ../CyteTypeR` after each CyteTypeR release.

Never leave HTML comments in the pages. The site publishes each page's source file, page data
and search index, so comments become public. Track pending work, such as screenshots to add,
in `docs/MAINTENANCE.md` instead.

## Links

- Link other pages with relative file links, such as `[Sign in](../getting-started/sign-in.md)`
  or `[Quickstart](../getting-started/quickstart.ipynb)`. Relative links must resolve to a
  page in `docs/source/`.
- Link the Python API through its Sphinx source with an `.html` (or `.md`) suffix, such as
  `[run()](../reference/api/cytetype.html#cytetype.CyteType.run)`. The build points these
  links at `/api/cytetype.html#cytetype.CyteType.run`. Anchors are the documented names,
  such as `#cytetype.CyteType.run` or `#cytetype.QuotaExceededError`. The absolute form
  `/api/<page>.html` also works.
- In the Sphinx pages, link guides relative to `/api/`, such as
  `[Python quickstart](../getting-started/quickstart)`.

## Setup and rendering

Use Node.js 24 and install the Python tools:

```bash
uv sync --group docs
make -C docs html
```

The site is written to ignored `docs/build/site/_build/html/`. Rendering does not execute code.
Every executable notebook needs a completed execution timestamp, a matching code fingerprint,
executed code cells, and no errors. Prose edits can reuse its outputs. Changed code requires
re-execution and the updated notebook in the same change. For a lightweight source preview,
use `make -C docs preview`; unexecuted pages display a warning that figures and results are
unavailable. Preview builds still reject changed code and incomplete or failed executions.
The render plugin, `docs/render_outputs.py`, turns saved HTML outputs into native page
content: pandas tables, the timing footers, and CyteType's job progress display, which is a
single `<pre>` element and renders as a text block. Any other HTML output fails the build;
show such objects as text, Markdown or a figure instead. To view a build, serve the output
directory, for example with
`python -m http.server -d docs/build/site/_build/html 8000`, and open
`http://localhost:8000/`.

Pages excluded by `project.exclude` in `docs/source/myst.yml` remain source-only and need no
saved execution. Builds clear the generated site and API output so excluded pages leave no old
routes.

Pull requests and the Pages workflow validate and render the notebooks from the checkout
with `make -C docs html`. CI does not execute them. Missing, failed, or stale
executions fail the build. GitHub Pages publishes from `master` at
`https://docs.nygen.io/CyteType/`. The workflow builds with `BASE_URL=/CyteType` and uploads
the generated site as the Pages artifact; the organization's Pages site serves it under
`/CyteType/` and owns the root `robots.txt`. Keep rendered site files outside Git.

## Execute notebooks locally

Executing a notebook runs every cell. A page that calls `CyteType.run()`, such as the
quickstart, uploads data and submits a real annotation job: it needs a saved CyteType sign-in
on this machine, counts toward the account's job limits, and waits for the job to finish.
Ask before executing unless the current request already authorizes it.

```bash
uv run cytetype setup
make -C docs execute-docs PAGES="quickstart"
```

`cytetype setup` is needed once per machine; `cytetype login` saves an existing API key
instead. `PAGES` takes notebook names, with the section folder when a name is ambiguous, such
as `getting-started/quickstart`. Pages run one at a time in a new kernel of the `uv`
environment, which uses the checkout's `cytetype`, from a temporary working directory, so
downloads, `vars.h5`, `obs.duckdb` and `query.json` never land in the repository. The kernel
runs with `TQDM_DISABLE=1`, so progress bars stay out of the published outputs.

Outputs are published, so the kernel does not get environment variables whose names contain
`KEY`, `TOKEN`, `SECRET`, `PASSWORD`, `PASSWD` or `CREDENTIAL`; it reads the CyteType sign-in
from the credentials file under `HOME` or `XDG_CONFIG_HOME`. The runner also refuses to save a
page whose outputs contain a CyteType API key (`cyt_...`), the saved API key or the value of
such a variable in your environment. The error names the cell, never the secret.

Each successful execution atomically replaces the same `.ipynb` file in `docs/source/` and
records `metadata.executed_at` and `metadata.executed_code_sha256`. The latter identifies the
executed code; do not edit this fingerprint by hand. A failed page leaves its previous notebook
intact and the command exits with an error after the remaining pages; rerun the failed pages
directly. Do not edit or render a notebook while it is being executed: the runner keeps your
edits and discards its result.

The runner records each code cell's UTC start time and elapsed wall time and displays them
below its output. CPU time is not measured. Keep `%%time` out of notebook cells.

The published outputs show the job's report link. The quickstart says that its report is
public, so after executing it, set the job's visibility to public in the dashboard. The job
also records the client version of the checkout, so execute from the released commit when the
version shown in the report matters.

## Checks

```bash
uv run --group docs pytest tests/test_docs_structure.py tests/test_docs_render.py \
    tests/test_docs_execution.py
make -C docs check-reference
```

The API check is a nitpicky, warnings-as-errors Sphinx build. Its inventory must document every
`cytetype.__all__` export under its public `cytetype.` name, and every public `CyteType`
member exactly once. The structure tests check the table of contents, `llms.txt`, relative
links, tabs, that pages hold no HTML comments, notebook comments, that notebooks hold no API
keys, and that the LLM configuration reference matches `cytetype.api.LLMModelConfig`. Inspect the rendered pages after changing execution or
presentation.
