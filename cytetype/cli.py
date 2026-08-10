import argparse
import base64
import getpass
import hashlib
import os
import secrets
import shutil
import subprocess
import sys
import time
import webbrowser
from html import escape
from http.server import BaseHTTPRequestHandler, HTTPServer
from importlib.resources import files
from string import Template
from urllib.parse import parse_qs, quote, urlencode, urlparse

import requests
from pydantic import ValidationError  # pyright: ignore[reportMissingImports]

from . import __version__
from .config import (
    DEFAULT_API_URL,
    DEFAULT_DASHBOARD_URL,
    StoredCredentials,
    delete_credentials,
    get_default_api_url,
    load_credentials,
    resolve_dashboard_url,
    save_credentials,
    validate_api_url,
)

_SETUP_TIMEOUT_SECONDS = 300
_ANSI_RED = "\033[91m"
_ANSI_BLUE = "\033[94m"
_ANSI_RESET = "\033[0m"
_NYGEN_GLYPHS = {
    "n": ("      ", "# ### ", "##   #", "#    #", "#    #", "#    #", "#    #"),
    "y": ("#   #", "#   #", "#   #", " ####", "    #", "#   #", " ### "),
    "g": (" ### ", "#   #", "#   #", "#   #", " ####", "    #", " ### "),
    "e": ("#### ", "#   #", "#####", "#    ", "#    ", "#   #", " ### "),
}
_NYGEN_WORDMARK_PIXELS = tuple(
    " ".join(_NYGEN_GLYPHS[letter][row] for letter in "nygen") for row in range(7)
)
_NYGEN_BANNER = tuple(
    "".join("██" if pixel == "#" else "  " for pixel in row).rstrip()
    for row in _NYGEN_WORDMARK_PIXELS
)
_CALLBACK_PAGE_TEMPLATE = Template(
    files("cytetype")
    .joinpath("templates", "cli_callback.html")
    .read_text(encoding="utf-8")
)


def _render_callback_page(
    message: str,
    credentials: StoredCredentials | None = None,
) -> bytes:
    context = {
        "message": escape(message),
        "meta_refresh": "",
        "state_class": "error",
        "icon": "!",
        "success_hidden": "hidden",
        "error_hidden": "",
        "email": "",
        "dashboard_url": "",
    }
    if credentials is not None:
        dashboard_url = escape(
            resolve_dashboard_url(credentials),
            quote=True,
        )
        context.update(
            {
                "meta_refresh": (
                    f'<meta http-equiv="refresh" content="5;url={dashboard_url}">'
                ),
                "state_class": "success",
                "icon": "✓",
                "success_hidden": "",
                "error_hidden": "hidden",
                "email": escape(credentials.email),
                "dashboard_url": dashboard_url,
            }
        )

    return _CALLBACK_PAGE_TEMPLATE.substitute(context).encode("utf-8")


def _create_pkce_pair() -> tuple[str, str]:
    verifier = secrets.token_urlsafe(64)
    digest = hashlib.sha256(verifier.encode("ascii")).digest()
    challenge = base64.urlsafe_b64encode(digest).rstrip(b"=").decode("ascii")
    return verifier, challenge


def _open_browser_silently(url: str) -> bool:
    """
    Handle the case of WSL. Browser opening is set to default windows browser.
    """
    command: list[str] | None = None
    is_wsl = sys.platform == "linux" and bool(
        os.environ.get("WSL_DISTRO_NAME") or os.environ.get("WSL_INTEROP")
    )

    if is_wsl:
        launcher = shutil.which("rundll32.exe")
        if launcher is None:
            return False
        command = [launcher, "url.dll,FileProtocolHandler", url]
    elif sys.platform == "linux":
        launcher = shutil.which("xdg-open")
        if launcher is None:
            return False
        command = [launcher, url]
    elif sys.platform == "darwin":
        launcher = shutil.which("open")
        if launcher is None:
            return False
        command = [launcher, url]

    if command is not None:
        try:
            subprocess.Popen(
                command,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
        except OSError:
            return False
        return True

    stderr_fd = 2
    saved_stderr_fd: int | None = None
    sink_fd: int | None = None

    try:
        saved_stderr_fd = os.dup(stderr_fd)
        sink_fd = os.open(os.devnull, os.O_WRONLY)
        os.dup2(sink_fd, stderr_fd)
    except OSError:
        if saved_stderr_fd is not None:
            os.close(saved_stderr_fd)
        if sink_fd is not None:
            os.close(sink_fd)
        saved_stderr_fd = None
        sink_fd = None

    try:
        return bool(webbrowser.open(url))
    except (OSError, webbrowser.Error):
        return False
    finally:
        if saved_stderr_fd is not None:
            os.dup2(saved_stderr_fd, stderr_fd)
            os.close(saved_stderr_fd)
        if sink_fd is not None:
            os.close(sink_fd)


def _print_setup_banner() -> None:
    print()
    print("\n".join(_NYGEN_BANNER))
    print()
    print(f"Nygen Analytics:      {_ANSI_BLUE}https://nygen.io{_ANSI_RESET}")
    print()


def _run_setup(api_url: str, force: bool = False) -> StoredCredentials:
    _print_setup_banner()
    existing = None if force else load_credentials(api_url)
    if existing is not None:
        try:
            response = requests.get(
                f"{api_url}/auth/cli/credentials",
                headers={"Authorization": f"Bearer {existing.apiToken}"},
                timeout=30,
            )
        except requests.RequestException as error:
            raise RuntimeError("Could not validate the saved API key") from error

        if not response.ok:
            detail: object | None = None
            error_code: str | None = None
            try:
                data = response.json()
                if isinstance(data, dict):
                    detail = data.get("detail")
                    if isinstance(detail, dict):
                        response_error_code = detail.get("error_code")
                        if isinstance(response_error_code, str):
                            error_code = response_error_code
                        detail = detail.get("message")
            except ValueError:
                pass
            message = str(detail) if detail else "Saved API key validation failed"
            if error_code in {"INVALID_TOKEN", "TOKEN_INACTIVE"}:
                message += (
                    ". Run `cytetype setup --force` to re-authenticate with a new key"
                )
            raise RuntimeError(message)

        try:
            data = response.json()
            if not isinstance(data, dict):
                raise TypeError
            existing = StoredCredentials(
                apiUrl=api_url,
                dashboardUrl=data["dashboardUrl"],
                apiToken=existing.apiToken,
                tokenId=data["tokenId"],
                userId=data["userId"],
                email=data["email"],
            )
        except (KeyError, TypeError, ValidationError, ValueError) as error:
            raise RuntimeError("Server returned invalid CLI credentials") from error

        save_credentials(existing)
        print(f"CyteType is already configured for {existing.email}.")
        dashboard_url = resolve_dashboard_url(existing)
        print(f"Dashboard:  {_ANSI_BLUE}{dashboard_url}{_ANSI_RESET}")
        return existing

    state = secrets.token_urlsafe(32)
    verifier, challenge = _create_pkce_pair()
    result: StoredCredentials | Exception | None = None
    redirect_uri = ""

    class CallbackHandler(BaseHTTPRequestHandler):
        def log_message(self, format: str, *args: object) -> None:
            return

        def respond(
            self,
            status: int,
            message: str,
            credentials: StoredCredentials | None = None,
        ) -> None:
            body = _render_callback_page(message, credentials)
            self.send_response(status)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Referrer-Policy", "no-referrer")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("X-Frame-Options", "DENY")
            self.send_header(
                "Content-Security-Policy",
                "default-src 'none'; style-src 'unsafe-inline'; "
                "base-uri 'none'; form-action 'none'; frame-ancestors 'none'",
            )
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self) -> None:
            nonlocal result
            parsed = urlparse(self.path)
            if parsed.path != "/callback":
                self.respond(404, "Unknown callback path.")
                return

            query = parse_qs(parsed.query)
            returned_state = query.get("state", [""])[0]
            if not secrets.compare_digest(returned_state, state):
                self.respond(400, "Authorization state did not match.")
                return

            callback_error = query.get("error", [""])[0]
            if callback_error:
                result = RuntimeError("CyteType authorization was denied")
                self.respond(400, "CyteType authorization failed.")
                return

            code = query.get("code", [""])[0]
            if not code:
                self.respond(400, "Authorization code was missing.")
                return

            try:
                response = requests.post(
                    f"{api_url}/auth/cli/token",
                    json={
                        "code": code,
                        "codeVerifier": verifier,
                        "redirectUri": redirect_uri,
                    },
                    timeout=30,
                )
                if not response.ok:
                    try:
                        detail = response.json().get("detail")
                    except ValueError:
                        detail = None
                    raise RuntimeError(
                        str(detail) if detail else "Token exchange failed"
                    )
                data = response.json()
                if not isinstance(data, dict):
                    raise TypeError
                result = StoredCredentials(
                    apiUrl=api_url,
                    dashboardUrl=data.get(
                        "dashboardUrl",
                        DEFAULT_DASHBOARD_URL,
                    ),
                    apiToken=data["apiToken"],
                    tokenId=data["tokenId"],
                    userId=data["userId"],
                    email=data["email"],
                )
            except (KeyError, TypeError, ValidationError, ValueError):
                result = RuntimeError("Server returned invalid CLI credentials")
            except requests.RequestException:
                result = RuntimeError("Could not exchange the authorization code")
            except RuntimeError as exchange_error:
                result = exchange_error

            if isinstance(result, StoredCredentials):
                self.respond(
                    200,
                    "CyteType setup is complete.",
                    credentials=result,
                )
            else:
                self.respond(400, "CyteType setup failed.")

    server = HTTPServer(("127.0.0.1", 0), CallbackHandler)
    try:
        redirect_uri = f"http://127.0.0.1:{server.server_port}/callback"
        query = urlencode(
            {
                "redirectUri": redirect_uri,
                "state": state,
                "codeChallenge": challenge,
            }
        )
        authorize_url = f"{api_url}/auth/cli/authorize?{query}"
        print("Opening CyteType sign-in in your browser.")
        print("If it does not open automatically, use this URL:")
        print()
        print(f"{_ANSI_BLUE}{authorize_url}{_ANSI_RESET}")
        print(flush=True)
        if not _open_browser_silently(authorize_url):
            print(f"{_ANSI_RED}Was not able to launch web browser{_ANSI_RESET}")

        deadline = time.monotonic() + _SETUP_TIMEOUT_SECONDS
        while result is None and time.monotonic() < deadline:
            server.timeout = min(1.0, max(0.0, deadline - time.monotonic()))
            server.handle_request()
    finally:
        server.server_close()

    if result is None:
        raise TimeoutError("CyteType setup timed out")
    if isinstance(result, Exception):
        raise result

    path = save_credentials(result)
    print(f"API key saved for {result.email} in {path}.")
    dashboard_url = resolve_dashboard_url(result)
    print(f"Dashboard:  {_ANSI_BLUE}{dashboard_url}{_ANSI_RESET}")
    return result


def _run_login(api_url: str) -> StoredCredentials:
    api_token = getpass.getpass("API key: ").strip()
    if not api_token:
        raise ValueError("API key is required")

    try:
        response = requests.get(
            f"{api_url}/auth/cli/credentials",
            headers={"Authorization": f"Bearer {api_token}"},
            timeout=30,
        )
    except requests.RequestException as error:
        raise RuntimeError("Could not validate the API key") from error

    if not response.ok:
        detail: object | None = None
        try:
            data = response.json()
            if isinstance(data, dict):
                detail = data.get("detail")
                if isinstance(detail, dict):
                    detail = detail.get("message")
        except ValueError:
            pass
        raise RuntimeError(str(detail) if detail else "API key validation failed")

    try:
        data = response.json()
        if not isinstance(data, dict):
            raise TypeError
        credentials = StoredCredentials(
            apiUrl=api_url,
            dashboardUrl=data["dashboardUrl"],
            apiToken=api_token,
            tokenId=data["tokenId"],
            userId=data["userId"],
            email=data["email"],
        )
    except (KeyError, TypeError, ValidationError, ValueError) as error:
        raise RuntimeError("Server returned invalid CLI credentials") from error

    path = save_credentials(credentials)
    print(f"Signed in as {credentials.email}.")
    print(f"API key saved in {path}.")
    dashboard_url = resolve_dashboard_url(credentials)
    print(f"Dashboard:  {_ANSI_BLUE}{dashboard_url}{_ANSI_RESET}")
    return credentials


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cytetype",
        description="Authenticate with CyteType and open your jobs.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    commands = parser.add_subparsers(dest="command")

    setup = commands.add_parser(
        "setup",
        aliases=["get-key"],
        help="Sign in and save a personal API key.",
    )
    setup.add_argument(
        "--api-url",
        default=get_default_api_url(),
        help=(
            "CyteType server origin. Defaults to CYTETYPE_API_URL or "
            f"{DEFAULT_API_URL}."
        ),
    )
    setup.add_argument(
        "--force",
        action="store_true",
        help="Re-authenticate without validating saved credentials.",
    )
    login = commands.add_parser(
        "login",
        help="Save and validate an existing API key.",
    )
    login.add_argument(
        "--api-url",
        default=get_default_api_url(),
        help=(
            "CyteType server origin. Defaults to CYTETYPE_API_URL or "
            f"{DEFAULT_API_URL}."
        ),
    )
    commands.add_parser("dashboard", help="Open the CyteType dashboard.")
    view = commands.add_parser("view", help="Open a CyteType job report.")
    view.add_argument("job_id", help="Job identifier to open.")
    commands.add_parser("logout", help="Remove the locally saved API key.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        if args.command in {"setup", "get-key"}:
            _run_setup(validate_api_url(args.api_url), force=args.force)
            return 0

        if args.command == "login":
            _run_login(validate_api_url(args.api_url))
            return 0

        if args.command == "logout":
            if delete_credentials():
                print("Local API key removed. Revoke it from the dashboard if needed.")
            else:
                print("No local CyteType API key was found.")
            return 0

        if args.command == "dashboard":
            credentials = load_credentials()
            target = (
                resolve_dashboard_url(credentials)
                if credentials
                else DEFAULT_DASHBOARD_URL
            )
            print(f"{_ANSI_BLUE}{target}{_ANSI_RESET}")
            _open_browser_silently(target)
            return 0

        if args.command == "view":
            credentials = load_credentials()
            api_url = credentials.apiUrl if credentials else DEFAULT_API_URL
            redirect = f"/report/{quote(args.job_id, safe='')}"
            target = f"{api_url}/login?{urlencode({'redirect': redirect})}"
            print(target)
            _open_browser_silently(target)
            return 0

        parser.print_help()
        return 0
    except (EOFError, OSError, RuntimeError, TimeoutError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nCancelled.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
