---
description: Install the CyteType Python package or the CyteTypeR R package, check the installed version and keep it up to date.
---

(installation)=
# Installation

CyteType has two clients:

- `cytetype`, the Python package, for AnnData and Scanpy workflows;
- `CyteTypeR`, the R package, for Seurat workflows.

Both prepare your data on your own computer and send it to the CyteType service at `https://cytetype.nygen.io`, where the annotation runs. Install the client for the environment you analyze your data in. You can install both; they share the same sign-in.

## Requirements

| | Python | R |
|---|---|---|
| Language version | Python 3.12 or newer | R 4.2 or newer (see [R notes](#installation-r-notes)) |
| Package source | PyPI | GitHub, with dependencies from CRAN and Bioconductor |
| Single-cell toolkit | AnnData, installed with `cytetype` | Seurat, installed with `CyteTypeR` |

You also need internet access to `https://cytetype.nygen.io` and a CyteType account, which needs an institutional email address (see [Sign in and API keys](sign-in.md)).

## Install the client

::::{tab-set}
:::{tab-item} Python
:sync: python
With pip, in the virtual environment or conda environment you use for your analysis:

```bash
pip install cytetype
```

With uv:

```bash
# Add CyteType to a uv project
uv add cytetype

# Or install it into the active virtual environment
uv pip install cytetype
```

The package also installs the `cytetype` command, which you use to sign in.
:::
:::{tab-item} R
:sync: r
With devtools, as in the CyteTypeR README:

```r
install.packages("devtools")
devtools::install_github("NygenAnalytics/CyteTypeR")
```

Or with pak:

```r
install.packages("pak")
pak::pak("NygenAnalytics/CyteTypeR")
```

Both installers also install the packages that CyteTypeR requires, including Seurat, duckdb and the Bioconductor package rhdf5 (see [R notes](#installation-r-notes)).
:::
::::

(installation-python-notes)=
### Python notes

- **Check your Python version first** with `python --version`. pip only offers releases that support the Python you run it with. On Python 3.11, `pip install cytetype` installs release 0.12.0, the last one for 3.11, without any error. On older versions it finds no release at all.
- **Prefer a dedicated environment.** `cytetype` requires narrow version ranges of some common packages (anndata 0.12.x, pydantic 2.12.x and requests 2.32.x), so installing it into an existing environment can upgrade or downgrade them.

(installation-r-notes)=
### R notes

- **R 4.2 or newer.** CyteTypeR requires it, and so does the duckdb package it uses. Check your version with `R.version.string`, and update R first if it is older.
- **duckdb** writes the cell metadata file that is uploaded with each job, and both installers add it for you. On Linux, CRAN installs duckdb from source, which can take a while.
- **rhdf5 from Bioconductor.** CyteTypeR writes the expression file with rhdf5, a Bioconductor package, and both installers add it for you. pak installs it from the Bioconductor repositories, which provide binary packages for Windows and macOS. devtools looks it up on Bioconductor's Git server instead, which needs Git (or the git2r package) on your computer, and when it installs rhdf5 it builds it from source, which needs a compiler (Rtools on Windows, the Xcode Command Line Tools on macOS). If you do not have these tools, use pak. If rhdf5 fails to build with devtools, install it from Bioconductor first and tell devtools to keep it:

  ```r
  install.packages("BiocManager")
  BiocManager::install("rhdf5")

  # Keep the installed rhdf5 instead of rebuilding it
  devtools::install_github("NygenAnalytics/CyteTypeR", upgrade = "never")
  ```

- **Seurat and the other dependencies,** such as SeuratObject and withr, are installed from CRAN when they are missing, so the first installation can take some time.

Optional packages that CyteTypeR uses when they are installed (the installers do not add them):

| Package | What it is for |
|---|---|
| `furrr` and `progressr` | With furrr, CyteTypeR uploads large files in parallel chunks (up to `upload_max_workers`, 6 by default); with the default `show_progress = TRUE` it also needs progressr. Without them, chunks are uploaded one at a time. |

## Install the packages for the quickstart

The quickstarts use the standard Scanpy and Seurat workflows to cluster the example data before annotating it.

::::{tab-set}
:::{tab-item} Python
:sync: python
Scanpy is not installed with `cytetype`. Install it with the `leiden` extra:

```bash
pip install "scanpy[leiden]"

# or, with uv
uv add "scanpy[leiden]"
```

The extra installs igraph and leidenalg, the two Leiden implementations Scanpy can use; the quickstart's `sc.tl.leiden(flavor="igraph")` needs igraph. CyteType's plotting helper, `cytetype.marker_dotplot()`, also needs Scanpy.
:::
:::{tab-item} R
:sync: r
The R quickstart uses Seurat and dplyr, which are installed with CyteTypeR. No other packages are needed.
:::
::::

## Check the installation

::::{tab-set}
:::{tab-item} Python
:sync: python
```bash
cytetype --version
python -c "import cytetype; print(cytetype.__version__)"
```

Both print the installed version, for example `cytetype 0.21.0` and `0.21.0`.

If your shell cannot find the `cytetype` command, the environment's scripts folder is not on your `PATH`. Activate the environment, or run the command through Python instead, for example `python -m cytetype.cli --version`.
:::
:::{tab-item} R
:sync: r
```r
packageVersion("CyteTypeR")
#> [1] '0.6.0'
```
:::
::::

## Update

::::{tab-set}
:::{tab-item} Python
:sync: python
```bash
pip install -U cytetype

# or, with uv in a virtual environment
uv pip install -U cytetype

# or, with uv in a project
uv lock --upgrade-package cytetype
uv sync
```

When Python imports `cytetype`, and when you run the `cytetype` command, the package asks PyPI for the latest version, waiting at most 2 seconds. If a newer release exists, it prints a warning such as:

```text
⚠️  A newer CyteType version is available on PyPI (<installed> -> <latest>). Run `pip install -U cytetype` to update.
```

The check stays silent when you are offline or PyPI cannot be reached.
:::
:::{tab-item} R
:sync: r
```r
devtools::install_github("NygenAnalytics/CyteTypeR", force = TRUE)

# or, with pak
pak::pak("NygenAnalytics/CyteTypeR")
```

`force = TRUE` reinstalls the package even if devtools thinks it is up to date.

When you load the package with `library(CyteTypeR)`, it asks GitHub for the latest CyteTypeR release, waiting at most 2 seconds. If a newer release exists, it prints a startup message such as:

```text
A newer version of CyteTypeR is available (<installed> -> <latest>). Run devtools::install_github('NygenAnalytics/CyteTypeR') to update.
```

The check stays silent when you are offline or GitHub cannot be reached.
:::
::::

After updating, restart your Python kernel or R session so that the new version is loaded.

(installation-cytetyper)=
## Optional: the `cytetyper` command for R

Python users get the `cytetype` terminal command with the package. R users can install an equivalent command, `cytetyper`, which runs the CyteTypeR sign-in functions from a terminal:

```r
library(CyteTypeR)
InstallCyteTypeRCli()
```

The function prints and returns the path of the installed command. `cytetyper` has the same commands as the Python `cytetype` command:

```bash
cytetyper setup           # sign in with your browser (also: cytetyper get-key)
cytetyper login           # save an existing API key
cytetyper logout          # remove the saved key from this computer
cytetyper dashboard       # open the dashboard
cytetyper view <JOB_ID>   # open the report of a job
cytetyper --version
```

`cytetyper` only signs you in and opens pages. Annotation jobs always run from an R session with `CyteTypeR()`. You do not need the command at all: the R functions `SetupCyteTypeR()`, `LoginCyteTypeR()` and `LogoutCyteTypeR()` do the same things.

Before you rely on it, check these points:

- **It must stay next to `Rscript`.** The command runs the `Rscript` in its own folder, so `InstallCyteTypeRCli()` copies it into the folder of the first `Rscript` on your R session's `PATH`, after following symbolic links (with Homebrew or the CRAN installer for macOS, that is R's own program folder). Keep the default `destination`; a copy in a folder without `Rscript` does not work. If you have several R versions, check that this `Rscript` belongs to the R where CyteTypeR is installed. If the function stops with `Could not locate Rscript on PATH`, add R's `bin` folder to `PATH` and try again.
- **That folder must be on your `PATH`.** Otherwise a new terminal cannot find `cytetyper`. With the R installers from CRAN for macOS and Windows, R's own program folder is usually not on `PATH`: add it, or run the command with the full path that `InstallCyteTypeRCli()` printed. On Windows the installed file is `cytetyper.cmd`.
- **A system-wide R may need administrator rights.** If you cannot write to R's program folder, the function stops with `Could not install cytetyper in <folder>. Check that the directory is writable.` Ask your administrator, or use the R functions instead.
- **It ignores your R startup files.** The command starts R with `--vanilla`, which skips `.Rprofile` and `.Renviron`. If CyteTypeR is installed in a project library (for example with renv) or in a library that only those files point to, the command fails with `Error: Unable to load CyteTypeR`.
- **Reinstall it for a new R version.** A new R installation has its own program folder; run `InstallCyteTypeRCli()` again from the new version. Use `overwrite = TRUE` to replace an existing copy.

## Next step

Sign in once on each computer you run CyteType from: [Sign in and API keys](sign-in.md). If the installation fails, see [Troubleshooting](../support/troubleshooting.md).
