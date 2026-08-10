# CLI and Authentication

CyteType requires authentication before submitting jobs or fetching remote results. The recommended setup is browser-based sign-in through the CyteType CLI.

## Browser-Based Setup

Run:

```bash
cytetype setup
```

The command:

1. Starts a temporary callback server on `127.0.0.1` using an available port.
2. Opens the CyteType authorization page in your browser and prints the same URL in the terminal.
3. Verifies the callback state and exchanges the one-time authorization code using PKCE.
4. Saves the returned API credentials locally.

The API key is not included in the browser URL or printed in the terminal. If the browser does not open automatically, copy the printed URL into a browser. The command times out after five minutes if authorization is not completed.

Running `cytetype setup` again for the same server validates the saved key before reporting the configured account. A valid key does not open another browser. If the key is invalid or inactive, the command fails without changing the saved credentials. Run `cytetype setup --force` to skip validation and authenticate with a new key.

`cytetype get-key` is an alias for `cytetype setup`.

## Use Saved Credentials from Python

No authentication argument is needed after setup:

```python
from cytetype import CyteType

annotator = CyteType(
    adata,
    group_key="leiden",
)
adata = annotator.run(study_context="Human PBMC from a healthy donor")
```

When `run()` starts, CyteType loads the saved API key that matches the selected API server. `get_results()` uses the server saved with the job and resolves credentials for that server when a remote fetch is needed.

If no matching credentials are available, CyteType raises an authentication error and asks you to run `cytetype setup`.

## Use an Existing API Key

If you already have a personal API key, save and validate it with:

```bash
cytetype login
```

The key is entered through a hidden prompt. CyteType validates it with the selected server before replacing any saved credentials. A failed login leaves existing credentials unchanged.

## Commands

| Command | Purpose |
| --- | --- |
| `cytetype setup [--force]` | Validate saved credentials or sign in through a browser |
| `cytetype get-key` | Alias for `cytetype setup` |
| `cytetype login` | Validate and save an existing API key |
| `cytetype dashboard` | Open the dashboard for the saved server |
| `cytetype view <job_id>` | Open a job report through the saved server's sign-in flow |
| `cytetype logout` | Delete the locally saved credentials |
| `cytetype --version` | Print the installed CyteType version |
| `cytetype --help` | Show all available commands |

`cytetype logout` only removes the local credentials file. Revoke the key from the dashboard if it should no longer be accepted by the server.

## Custom Servers

Pass a server origin directly:

```bash
cytetype setup --api-url https://cytetype.example.org
```

Or set the default server for CLI and Python usage:

```bash
export CYTETYPE_API_URL=https://cytetype.example.org
cytetype setup
```

An explicit `--api-url` takes precedence over `CYTETYPE_API_URL`. An explicit `api_url` passed to `CyteType` or `run()` takes precedence in Python.

The API URL must be a server origin containing only the scheme and host, with an optional port. Paths, credentials, query strings, and fragments are rejected. Non-local servers must use HTTPS. `http://localhost` and `http://127.0.0.1` are allowed for local development.

Saved credentials are tied to the selected API origin. Run setup or login against the same origin used by Python:

```python
annotator = CyteType(
    adata,
    group_key="leiden",
    api_url="https://cytetype.example.org",
)
```

CyteType stores one credential set at a time. Completing setup or login for another server replaces the previously saved set.

## Credential Storage

Credentials are stored in `credentials.json` at:

| Platform | Default location |
| --- | --- |
| Linux and macOS | `~/.config/cytetype/credentials.json` |
| Linux and macOS with `XDG_CONFIG_HOME` | `$XDG_CONFIG_HOME/cytetype/credentials.json` |
| Windows | `%APPDATA%\cytetype\credentials.json` |

On POSIX systems, CyteType sets the directory to mode `0700` and the credentials file to mode `0600`. It also refuses to write into a credentials directory owned by another user.

The file contains the API key in plain JSON so the client can use it. Do not share it, commit it, or copy it into notebooks.

## Direct Tokens for CI, Remote Notebooks, and Managed Environments

Browser-based setup requires the authorization callback to reach `127.0.0.1` in the environment where the CLI is running. It is preferred for local use, but it may not work from a remote notebook, an SSH session without port forwarding, or CI.

If the remote environment has an interactive terminal and you already have an API key, use `cytetype login`. For non-interactive environments, read a token from the platform's secret store and pass it explicitly:

```python
import os

from cytetype import CyteType

annotator = CyteType(
    adata,
    group_key="leiden",
    api_url="https://cytetype.example.org",
    auth_token=os.environ["CYTETYPE_API_TOKEN"],
)
adata = annotator.run(study_context="Human PBMC from a healthy donor")
```

`CYTETYPE_API_TOKEN` in this example is a user-managed secret. CyteType does not read it automatically.

Authentication is resolved in this order:

1. An `auth_token` passed directly to `run()`.
2. An `auth_token` previously supplied to the `CyteType` instance for the same API origin.
3. Saved CLI credentials matching the API origin.

Tokens are not reused when the API origin changes. Pass a token for the new origin or run CLI setup against that origin.

For common setup failures, see [Troubleshooting](./troubleshooting.md).
