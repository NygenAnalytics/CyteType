---
description: Create a CyteType account, sign in from Python or R, use API keys on remote machines and in CI, and manage the saved key.
---

(sign-in)=
# Sign in and API keys

CyteType needs an API key to upload your data, submit jobs and fetch results. On your own computer you sign in once with your browser, and the client saves a personal API key that it then uses automatically, from both Python and R. On a machine where you cannot open a browser, such as a remote server or a CI runner, you create a key in the dashboard and give it to the client yourself.

## Your CyteType account

You sign in at `https://cytetype.nygen.io/login` with your email address. There is no password: CyteType emails you an 8-digit verification code, which you enter on the sign-in page. Your account is created the first time you sign in.

- **Use an institutional email address.** CyteType accepts university, institute and company addresses. Addresses from personal email providers (such as Gmail or Outlook) are rejected with the message "Please use your university/institute/company email".
- **Organization accounts.** If your email domain is registered to an organization with a CyteType license, your new account joins that organization automatically. The dashboard then shows an organization workspace next to your **Personal** workspace.
- **No code?** Check your spam or junk folder, or select **Resend code** on the sign-in page.

CyteType is free for academic and non-commercial research within the [usage limits](../platform/limits.md). Commercial use needs a license; contact [contact@nygen.io](mailto:contact@nygen.io).

(sign-in-browser)=
## Sign in with your browser

On a computer where you can open a web browser, run:

::::{tab-set}
:::{tab-item} Python
:sync: python
In a terminal, in the environment where `cytetype` is installed:

```bash
cytetype setup
```
:::
:::{tab-item} R
:sync: r
In an R session:

```r
library(CyteTypeR)
SetupCyteTypeR()
```

If you installed the [`cytetyper` command](installation.md#installation-cytetyper), `cytetyper setup` in a terminal does the same.
:::
::::

What happens next:

1. The command prints a sign-in link and opens it in your default browser. If no browser opens, copy the link into a browser on the same computer. The command keeps waiting.
2. If you are not signed in to CyteType in that browser, enter your email address and the 8-digit code from your inbox.
3. A page titled **Authorize CyteType CLI?** shows the email of the account that will receive the key. Select **Authorize CLI**. This page appears every time, even if you are already signed in; only approve requests that you started yourself.
4. The browser shows "CyteType setup is complete." and moves on to your dashboard. The terminal prints your email address and the location of the saved key.

Finish these steps within 5 minutes. After that the command stops with "CyteType setup timed out", and you can run it again. If you select **Cancel**, it stops with "CyteType authorization was denied".

Running the command again on a computer that is already signed in does not open the browser. It checks the saved key with the server and prints "CyteType is already configured for" followed by your email. If the server rejects the saved key, for example because the key was disabled, the command stops and asks you to sign in again with `--force` (see [Switch accounts or servers](#sign-in-switch)).

**The key it creates.** Browser sign-in gives you a personal API key named "CyteType CLI", listed under **API keys** in the **Personal** workspace of the [dashboard](https://cytetype.nygen.io/dashboard). Later browser sign-ins, on this or another computer, reuse that key while it is active.

(sign-in-remote)=
## Remote servers, SSH and headless machines

Browser sign-in needs a browser on the same computer as the command, because the sign-in page hands the key back through a temporary local address on that computer (`http://127.0.0.1:<port>/callback`). On a remote server, an HPC login node, an SSH session or a hosted notebook, create a key in the dashboard and save it on that machine instead:

1. On your own computer, open the [dashboard](https://cytetype.nygen.io/dashboard), select **API keys**, then **Create API key**.
2. Name the key (for example after the server) and select **Create API key**. Copy the key: it is shown only once.
3. On the remote machine, run the login command and paste the key at the prompt. Nothing is shown while you paste.

::::{tab-set}
:::{tab-item} Python
:sync: python
```bash
cytetype login
```
:::
:::{tab-item} R
:sync: r
```r
library(CyteTypeR)
LoginCyteTypeR()
```

`LoginCyteTypeR()` also accepts the key as an argument, for example `LoginCyteTypeR(api_token = Sys.getenv("CYTETYPE_API_TOKEN"))`. Read it from an environment variable or a secret store rather than typing it into a script.

If you installed the [`cytetyper` command](installation.md#installation-cytetyper), `cytetyper login` in a terminal does the same.
:::
::::

The client checks the key with the server before saving it, then prints the account email and where the key was saved. If the server rejects the key, nothing is saved and any key saved earlier stays in place.

In an organization workspace, only organization admins can create keys. Switch to your **Personal** workspace to create a personal key, or ask an admin for an organization key. Personal keys start with `cyt_p_` and organization keys with `cyt_`. Keys do not expire, but they can be disabled.

## Pass a key in code

Instead of saving a key, you can pass it to the client in code. This is the usual choice for CI and scheduled jobs, and for shared machines where you do not want a key stored in your home directory.

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
import os

from cytetype import CyteType

annotator = CyteType(adata, group_key="leiden", rank_key="rank_genes_leiden")
adata = annotator.run(
    study_context="Human PBMC from a healthy donor",
    auth_token=os.environ["CYTETYPE_API_TOKEN"],
)
```

You can also pass `auth_token` to `CyteType()` instead of `run()`. If you pass it only to `run()` on a machine without a saved key, `CyteType()` prints a note that no saved sign-in was found; you can ignore it.
:::
:::{tab-item} R
:sync: r
```r
pbmc <- CyteTypeR(
  obj = pbmc,
  prepped_data = prepped_data,
  study_context = "Human PBMC from a healthy donor",
  auth_token = Sys.getenv("CYTETYPE_API_TOKEN")
)
```

`GetResults()` also accepts `auth_token`, for example `GetResults(job_id = "<job_id>", auth_token = Sys.getenv("CYTETYPE_API_TOKEN"))`. Results already stored in a Seurat object are returned without contacting the server.
:::
::::

`CYTETYPE_API_TOKEN` is only an example name: the clients never read it on their own. Keep keys out of scripts, notebooks and version control.

## CI and scheduled jobs

Store the key as a secret in your CI system or scheduler, expose it to the job as an environment variable, and pass it as `auth_token` as shown above. Browser sign-in is not possible there. [Automation, HPC and CI](../guides/automation.md) covers batch jobs in more detail.

(sign-in-key-order)=
## How the clients choose a key

**Python** uses the first key it finds, in this order:

1. the `auth_token` passed to `run()`;
2. an `auth_token` given earlier to `CyteType()`, or to an earlier `run()` of the same object, for the same server;
3. the saved key for that server. It is read again at every `run()`, so a key you replace with `cytetype setup --force` is used without restarting Python.

`get_results()` uses the server recorded with the job.

**R** uses:

1. the `auth_token` passed to `CyteTypeR()` or `GetResults()`;
2. the saved key for that server.

If there is no key, Python raises `AuthenticationError` with the message "CyteType sign-in is required. Run `cytetype setup` first." R stops with "CyteType sign-in is required. Run SetupCyteTypeR() first." Python also checks for a saved key as soon as you create `CyteType`, before the slower data preparation, and prints a note if it finds none.

(sign-in-credentials)=
## Where the key is saved

| System | Location |
|---|---|
| Any, when `XDG_CONFIG_HOME` is set | `$XDG_CONFIG_HOME/cytetype/credentials.json` |
| Windows | `%APPDATA%\cytetype\credentials.json` |
| Linux and macOS | `~/.config/cytetype/credentials.json` |

- Python and R use the same file, so signing in with one client signs in the other.
- The file holds the key for one server at a time.
- On Linux and macOS, the folder (mode `0700`) and the file (mode `0600`) are readable only by you, and the clients refuse to write into a folder that belongs to another user.
- The key is stored as plain text so that the clients can use it. Do not share, commit or copy the file.

(sign-in-switch)=
## Switch accounts or servers

**Saved key rejected.** If your saved key was disabled, sign in again with `--force`. It skips the check of the saved key and starts a new browser sign-in:

::::{tab-set}
:::{tab-item} Python
:sync: python
```bash
cytetype setup --force
```
:::
:::{tab-item} R
:sync: r
```r
SetupCyteTypeR(force = TRUE)
```
:::
::::

**Another account.** Browser sign-in uses the account you are signed in to in the browser. Sign out there first (dashboard, **Account**, **Sign out**), then run the `--force` command above and sign in with the other email address. Alternatively, save a key from the other account with `cytetype login` (R: `LoginCyteTypeR()`). Either way, the new key replaces the saved one.

**Another server.** If you were given the address of a different CyteType server, pass it when you sign in (`cytetype setup --api-url <address>`, R: `SetupCyteTypeR(api_url = "<address>")`) and to the client (`api_url`), or set the `CYTETYPE_API_URL` environment variable for both. The address must use HTTPS; plain HTTP is accepted only for `localhost` and `127.0.0.1`. Because the file holds one server at a time, signing in to another server replaces the saved key.

## Sign out or revoke a key

Removing a key from your computer and revoking it are different things:

- **Sign out on this computer.** `cytetype logout` (R: `LogoutCyteTypeR()`) deletes the saved key file. The key itself keeps working anywhere else it is used.
- **Revoke the key.** Open the [dashboard](https://cytetype.nygen.io/dashboard), select **API keys** and **Disable** the key. The server then rejects it everywhere. Disabling cannot be undone.

Because browser sign-in reuses the "CyteType CLI" key, disabling that key signs out every computer that saved it; the next browser sign-in creates a new one. If you think a key has leaked, disable it and create a new one.

## Commands at a glance

| Task | Python (terminal) | R | R (terminal, optional) |
|---|---|---|---|
| Sign in with the browser | `cytetype setup` | `SetupCyteTypeR()` | `cytetyper setup` |
| Replace a rejected key | `cytetype setup --force` | `SetupCyteTypeR(force = TRUE)` | `cytetyper setup --force` |
| Save an existing key | `cytetype login` | `LoginCyteTypeR()` | `cytetyper login` |
| Remove the saved key | `cytetype logout` | `LogoutCyteTypeR()` | `cytetyper logout` |
| Open the dashboard | `cytetype dashboard` | `OpenCyteTypeDashboard()` | `cytetyper dashboard` |

`cytetyper login` reads the key from standard input, so a script can also pipe it in, for example `printf '%s\n' "$CYTETYPE_API_TOKEN" | cytetyper login`. All options are listed in the [command-line reference](../reference/cli.md).

## Next step

Run your first annotation: the [Python quickstart](quickstart.ipynb) or the [R quickstart](quickstart-r.md).
