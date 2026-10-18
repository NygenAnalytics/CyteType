---
description: Reference for the cytetype (Python) and cytetyper (R) command-line tools, with every command, option, default, environment variable, output and exit code, and the matching R functions.
---

(cli)=
# Command-line reference

CyteType has two small command-line tools. Both sign you in and open CyteType pages in your browser:

- `cytetype` is installed with the Python package.
- `cytetyper` is an optional launcher for the R package. Every `cytetyper` command is also an R function, so you can do the same things from an R session without installing it.

Neither tool runs annotations. You submit jobs from Python with `CyteType.run()` and from R with `CyteTypeR()` (see the [Python quickstart](../getting-started/quickstart.ipynb) and the [R quickstart](../getting-started/quickstart-r.md)). For a step-by-step guide to signing in, see [Sign in and API keys](../getting-started/sign-in.md).

(cli-commands)=
## Commands at a glance

| Task | Python CLI | R function | R CLI |
|---|---|---|---|
| Sign in with the browser and save a personal API key | `cytetype setup [--api-url URL] [--force]` | `SetupCyteTypeR(api_url = NULL, force = FALSE)` | `cytetyper setup [--api-url URL] [--force]` |
| Same as `setup` | `cytetype get-key [--api-url URL] [--force]` | Not available | `cytetyper get-key [--api-url URL] [--force]` |
| Save an API key you already have | `cytetype login [--api-url URL]` | `LoginCyteTypeR(api_token = NULL, api_url = NULL)` | `cytetyper login [--api-url URL]` |
| Delete the saved key from this computer | `cytetype logout` | `LogoutCyteTypeR()` | `cytetyper logout` |
| Open the dashboard | `cytetype dashboard` | `OpenCyteTypeDashboard()` | `cytetyper dashboard` |
| Open the report of a job | `cytetype view JOB_ID` | `ViewCyteTypeJob(job_id, api_url = NULL)` | `cytetyper view JOB_ID [--api-url URL]` |
| Print the installed version | `cytetype --version` | `packageVersion("CyteTypeR")` | `cytetyper --version` |
| Show help | `cytetype --help`, `cytetype COMMAND --help` | `?SetupCyteTypeR` and the other help pages | `cytetyper --help`, `cytetyper COMMAND --help` |
| Install the R CLI | Not applicable | `InstallCyteTypeRCli(destination = NULL, overwrite = FALSE)` | Not applicable |

All of these share one saved sign-in (see [Saved credentials](#cli-credentials)). Signing in with either tool, or from either package, signs in both clients.

The same commands in the two clients:

::::{tab-set}
:::{tab-item} Python
:sync: python
```bash
cytetype setup            # sign in with the browser
cytetype view <JOB_ID>    # open a report
cytetype logout           # delete the saved key on this computer
```
:::
:::{tab-item} R
:sync: r
```r
library(CyteTypeR)
SetupCyteTypeR()             # sign in with the browser
ViewCyteTypeJob("<JOB_ID>")  # open a report
LogoutCyteTypeR()            # delete the saved key on this computer
```

With the optional launcher, `cytetyper setup`, `cytetyper view <JOB_ID>` and `cytetyper logout` do the same from a terminal.
:::
::::

(cli-cytetype)=
## `cytetype` (Python)

The `cytetype` command is installed with the Python package. If your shell cannot find it, activate the environment where `cytetype` is installed, or run the same commands through Python:

```bash
python -m cytetype.cli setup
```

Usage:

```text
cytetype [-h] [--version] {setup,get-key,login,dashboard,view,logout} ...
```

Running `cytetype` without a command prints the help. Options go after the command name: `cytetype setup --force`, not `cytetype --force setup`.

(cli-setup)=
### `cytetype setup`

```text
cytetype setup [--api-url URL] [--force]
```

Signs you in with your browser and saves a personal API key on this computer.

| Option | Default | Description |
|---|---|---|
| `--api-url URL` | `CYTETYPE_API_URL` if set, else `https://cytetype.nygen.io` | The CyteType server to sign in to. Must be a server origin (see [Server URL](#cli-server-url)). |
| `--force` | Off | Skip the check of the saved key and start a new browser sign-in. Use it when the server no longer accepts the saved key, for example after you disabled it in the dashboard. |

What the command does:

1. If a key for this server is already saved and you did not pass `--force`, it checks the key with the server. If the key is valid, it prints `CyteType is already configured for <email>.` and your dashboard link, and opens no browser. If the server rejects the key as invalid or disabled, the command stops with the server's message followed by "Run `cytetype setup --force` to re-authenticate with a new key".
2. Otherwise it starts a temporary listener on `127.0.0.1`, on a free port, and opens the CyteType sign-in page in your default browser. It also prints the same link, so you can copy it if no browser opens.
3. In the browser, sign in with the one-time code sent to your email (if that browser is not signed in yet), then select **Authorize CLI**.
4. The browser returns to the listener. The command exchanges a short-lived code for your API key, saves it and prints `API key saved for <email> in <path>.` and your dashboard link. The browser shows "CyteType setup is complete." and moves on to the dashboard after a few seconds.

The command waits at most 5 minutes for you to finish in the browser. The API key never appears in the browser's address bar or in the terminal.

The key you get is your personal key named "CyteType CLI". If that key exists and is active, setup returns the same key, so every computer you sign in to this way shares it. If you disabled it, the next `setup --force` creates a new one.

Example output (the banner is left out):

```text
Opening CyteType sign-in in your browser.
If it does not open automatically, use this URL:

<sign-in link>

API key saved for you@university.edu in /home/you/.config/cytetype/credentials.json.
Dashboard:  https://cytetype.nygen.io/dashboard
```

The browser is opened with `open` on macOS, `xdg-open` on Linux and the Windows default browser under WSL. If none is available, the command prints `Was not able to launch web browser` and keeps waiting: copy the printed link into a browser on the same computer.

### `cytetype get-key`

Another name for `cytetype setup`, with the same options and behavior.

(cli-login)=
### `cytetype login`

```text
cytetype login [--api-url URL]
```

Saves an API key that you already have, for example one created on the **API keys** page of the [dashboard](../platform/dashboard.md). The command asks for the key at a hidden `API key:` prompt, checks it with the server, and saves it only if the server accepts it. If the check fails, nothing is saved and any key saved earlier stays in place. Personal keys (starting with `cyt_p_`) and organization keys (starting with `cyt_`) both work.

| Option | Default | Description |
|---|---|---|
| `--api-url URL` | `CYTETYPE_API_URL` if set, else `https://cytetype.nygen.io` | The server that issued the key. |

Output on success:

```text
Signed in as you@university.edu.
API key saved in /home/you/.config/cytetype/credentials.json.
Dashboard:  https://cytetype.nygen.io/dashboard
```

The prompt reads the key from your terminal. A notebook has no terminal to type into, so open one (in JupyterLab: **File**, **New**, **Terminal**) or pass the key to `run()` as `auth_token` instead (see [Sign in and API keys](../getting-started/sign-in.md)).

### `cytetype logout`

```text
cytetype logout
```

Deletes the saved credentials file on this computer and prints `Local API key removed. Revoke it from the dashboard if needed.`, or `No local CyteType API key was found.` if there was none. The key itself keeps working wherever else it is saved or used. To stop it from working everywhere, disable it in the [dashboard](../platform/dashboard.md).

### `cytetype dashboard`

```text
cytetype dashboard
```

Prints the dashboard address of the server you are signed in to and opens it in your browser. Without a saved sign-in, it opens `https://cytetype.nygen.io/dashboard`. The command takes no options and does not use `CYTETYPE_API_URL`.

### `cytetype view`

```text
cytetype view JOB_ID
```

| Argument | Description |
|---|---|
| `JOB_ID` | The job ID, as printed when the job starts and listed in the dashboard. |

Prints and opens the report of a job. The link goes through the sign-in page of the server you are signed in to (without a saved sign-in, `https://cytetype.nygen.io`), so you may be asked to sign in first. If the browser is already signed in, you go straight to the report. The command does not check that the job exists or that you may open it; the report page does that. It takes no `--api-url` option and does not use `CYTETYPE_API_URL`.

### `cytetype --version` and `--help`

`cytetype --version` prints `cytetype <version>`, for example `cytetype 0.21.0`. `cytetype --help` (or `-h`) lists the commands, and `cytetype COMMAND --help` describes one command.

(cli-cytetyper)=
## `cytetyper` (R)

### Install the launcher

```r
library(CyteTypeR)
InstallCyteTypeRCli()
```

`InstallCyteTypeRCli()` copies a launcher named `cytetyper` (`cytetyper.cmd` on Windows) into a folder, prints `Installed <path>.` and returns the path invisibly.

| Argument | Default | Description |
|---|---|---|
| `destination` | `NULL`: the folder that contains the `Rscript` found on your `PATH` | Where to install the launcher. The launcher runs the `Rscript` in its own folder, so keep the default unless the folder you choose also contains `Rscript`. |
| `overwrite` | `FALSE` | Replace an existing launcher. Without it, the function stops with `The cytetyper command already exists at <path>. Set overwrite = TRUE to replace it.` |

Other errors you may see:

| Message | Meaning |
|---|---|
| `Could not locate Rscript on PATH` | R's `bin` folder is not on the `PATH` of your R session. Pass the folder of the running R's `Rscript`, `InstallCyteTypeRCli(destination = R.home("bin"))`, or add the folder to `PATH`. |
| `CLI destination directory does not exist` | The `destination` folder does not exist. |
| `Could not install cytetyper in <folder>. Check that the directory is writable.` | You cannot write to that folder, which is common with a system-wide R. Ask your administrator, or use the R functions instead. |
| `Could not make the cytetyper command executable` | The copied file could not be made executable (Linux and macOS). |

The launcher starts R with `Rscript --vanilla`, which skips your `.Rprofile` and `.Renviron` files. CyteTypeR must therefore be installed in a library that R finds without them (not only in a project library such as renv), and `CYTETYPE_API_URL` must be set in the shell, not in `.Renviron`. See [Installation](../getting-started/installation.md#installation-cytetyper) for the other requirements.

### Usage

```text
cytetyper setup [--api-url URL] [--force]
cytetyper get-key [--api-url URL] [--force]
cytetyper login [--api-url URL]
cytetyper logout
cytetyper dashboard
cytetyper view JOB_ID [--api-url URL]
cytetyper --version
```

`cytetyper` without arguments, `cytetyper -h` and `cytetyper --help` print this help; `cytetyper COMMAND --help` describes one command. Give options as `--api-url URL` or `--api-url=URL`, each at most once.

Each command calls the R function shown in [Commands at a glance](#cli-commands), so it behaves like the Python command of the same name described above. The differences:

- `cytetyper login` reads the key from standard input. Type or paste it at the `API key:` prompt (on Linux and macOS, what you type is not shown), or pipe it in, for example `printf '%s\n' "$CYTETYPE_API_TOKEN" | cytetyper login`. An empty key stops the command with `Error: API key is required`.
- `cytetyper view` accepts `--api-url`. Without it, it uses the server you are signed in to, then `CYTETYPE_API_URL`, then `https://cytetype.nygen.io`.
- `cytetyper --version` prints `cytetyper <version>`, for example `cytetyper 0.6.0`.
- Usage errors (an unknown command or option, a missing `JOB_ID`) end with exit status 1, not 2.
- If `setup` cannot open a browser, R warns `Could not launch a web browser.` instead of printing `Was not able to launch web browser`. Copy the printed link in either case.

### The R functions

| Function | Arguments | Returns (invisibly) |
|---|---|---|
| `SetupCyteTypeR()` | `api_url = NULL`: the server, defaulting as for `--api-url`. `force = FALSE`: skip the check of the saved key. | The saved account details: server, dashboard address, email, user and key IDs. Never the key itself. |
| `LoginCyteTypeR()` | `api_token = NULL`: the key; when `NULL`, a hidden prompt asks for it in an interactive session, and otherwise the key is read from the first line of standard input. `api_url = NULL`: the server. | The saved account details, as above. |
| `LogoutCyteTypeR()` | None | `TRUE` if a saved key was deleted, otherwise `FALSE`. |
| `OpenCyteTypeDashboard()` | None | The dashboard address it opened. |
| `ViewCyteTypeJob()` | `job_id`: the job ID. `api_url = NULL`: the server. | The address it opened. |
| `InstallCyteTypeRCli()` | `destination = NULL`, `overwrite = FALSE` | The path of the installed launcher. |

Because `LoginCyteTypeR()` accepts the key as an argument, it is useful on machines where a prompt is inconvenient:

```r
LoginCyteTypeR(api_token = Sys.getenv("CYTETYPE_API_TOKEN"))
```

The variable name is your choice: CyteTypeR does not read it on its own.

(cli-server-url)=
## Server URL and `CYTETYPE_API_URL`

The default server is `https://cytetype.nygen.io`. You only need another address if you were given one, for example for a development server.

`setup`, `get-key` and `login` (and the R functions `SetupCyteTypeR()` and `LoginCyteTypeR()`) choose the server in this order:

1. `--api-url` (R: the `api_url` argument);
2. the `CYTETYPE_API_URL` environment variable;
3. `https://cytetype.nygen.io`.

The same variable sets the default server of `CyteType()` in Python and of the CyteTypeR functions. In an R session, when the variable is not set, the R option `cytetype.default.api.url` is used if you set it.

```bash
export CYTETYPE_API_URL=https://cytetype.example.org
cytetype setup
```

Both clients accept only a server origin:

- a scheme, a host name and an optional port, with no path, query string, fragment or user name (a trailing `/` is removed);
- HTTPS, except for `http://localhost` and `http://127.0.0.1`.

Anything else fails with `API URL must be an HTTP or HTTPS server origin` or `Non-local API URLs must use HTTPS`. In Python, a `CYTETYPE_API_URL` that is set but empty is also rejected; unset it instead (`unset CYTETYPE_API_URL`).

The saved sign-in belongs to one server. Signing in to another server replaces it, and the clients do not send a key saved for one server to another.

(cli-credentials)=
## Saved credentials

`setup`, `get-key` and `login` write a file named `credentials.json`: in `$XDG_CONFIG_HOME/cytetype/` when that variable is set, otherwise in `%APPDATA%\cytetype\` on Windows and in `~/.config/cytetype/` on Linux and macOS. [Where the key is saved](../getting-started/sign-in.md#sign-in-credentials) describes how the file is protected.

- There is one file per user account on the computer. The Python and R clients and both command-line tools share it.
- It holds the sign-in for one server: the API key as plain text, the server and dashboard addresses, your email, and your user and key IDs. Do not share it, commit it or copy it into notebooks.
- No environment variable holds the API key itself. In scripts and CI, pass the key as `auth_token` to `run()` or `CyteTypeR()` (see [Sign in and API keys](../getting-started/sign-in.md)).

## Environment variables

| Variable | Read by | Effect |
|---|---|---|
| `CYTETYPE_API_URL` | `setup`, `get-key` and `login` in both tools; `cytetyper view` and `ViewCyteTypeJob()` without a saved sign-in; `CyteType()` and the CyteTypeR functions | The default server origin. |
| `XDG_CONFIG_HOME` | Both clients | The folder that holds `cytetype/credentials.json`. |
| `APPDATA` | Both clients, on Windows | The folder that holds `cytetype\credentials.json` when `XDG_CONFIG_HOME` is not set. |

## Output and exit codes

`cytetype` prints results to standard output and errors to standard error, as `Error: <message>`. `cytetyper` also prints errors as `Error: <message>`, and prints its other messages, such as links and the signed-in email, to standard error too.

| Exit status | `cytetype` | `cytetyper` |
|---|---|---|
| 0 | Success, including `--help` and `--version` | Success, including help and `--version` |
| 1 | The command failed: for example a rejected key, a timeout, a network error, an invalid server URL or a damaged credentials file | The command failed, including usage errors and `Error: Unable to load CyteTypeR: ...` |
| 2 | Usage error: an unknown command or option, or a missing `JOB_ID` | Not used |
| 130 | Interrupted with Ctrl+C; prints `Cancelled.` | Interrupted; prints `Cancelled.` |

Messages you may see, and what to do:

| Message | Cause | What to do |
|---|---|---|
| `CyteType setup timed out` | The browser step was not finished within 5 minutes, or the browser could not reach the listener on `127.0.0.1`. | Run the command again. See [Troubleshooting](../support/troubleshooting.md#troubleshooting-sign-in). |
| `CyteType authorization was denied` | **Cancel** was selected on the authorization page. | Run the command again and select **Authorize CLI**. |
| `Invalid API token` or `API token has expired or been deactivated`, followed by "Run `cytetype setup --force` ..." | The server rejected the saved key, usually because it was disabled. | `cytetype setup --force` (R: `SetupCyteTypeR(force = TRUE)`). |
| `Invalid CyteType credentials file: <path>` | The saved file is damaged or incomplete. | `cytetype logout`, then `cytetype setup`. |
| `Could not validate the saved API key` or `Could not validate the API key` | The server could not be reached. | Check your network or proxy, then try again. |
| `API key is required` | The `login` prompt received an empty key. | Run `cytetype login` again and paste the key. |
| `API URL must be an HTTP or HTTPS server origin` or `Non-local API URLs must use HTTPS` | `--api-url` or `CYTETYPE_API_URL` is not a valid server origin. | See [Server URL](#cli-server-url). |
| `CyteType credentials directory is not owned by the current user: <path>` | The credentials folder belongs to another user, for example after running the command with `sudo`. | Fix the folder's owner, or set `XDG_CONFIG_HOME` to a folder of your own. |

For errors that appear while you run jobs, see [Errors and error codes](errors.md).
