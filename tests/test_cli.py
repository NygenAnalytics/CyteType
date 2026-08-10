import os
import stat
import threading
from pathlib import Path
from types import SimpleNamespace
from urllib.error import HTTPError
from urllib.parse import parse_qs, urlparse
from urllib.request import urlopen

import pytest

import cytetype.cli as cli
from cytetype.config import (
    DEFAULT_DASHBOARD_URL,
    StoredCredentials,
    delete_credentials,
    get_credentials_path,
    load_credentials,
    resolve_dashboard_url,
    save_credentials,
    validate_api_url,
)


@pytest.fixture
def credentials() -> StoredCredentials:
    return StoredCredentials(
        apiUrl="https://dev.cytetype.example",
        dashboardUrl="https://dashboard.cytetype.example/dashboard",
        apiToken="cyt_p_secret",
        tokenId="token-id",
        userId="user-id",
        email="researcher@university.edu",
    )


def test_credentials_round_trip_is_private_and_server_specific(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    credentials: StoredCredentials,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))

    path = save_credentials(credentials)

    assert path == tmp_path / "cytetype" / "credentials.json"
    assert stat.S_IMODE(path.parent.stat().st_mode) == 0o700
    assert stat.S_IMODE(path.stat().st_mode) == 0o600
    assert load_credentials("https://dev.cytetype.example/") == credentials
    assert load_credentials("https://other.example") is None
    assert delete_credentials() is True
    assert delete_credentials() is False


@pytest.mark.skipif(os.name != "posix", reason="POSIX permissions required")
def test_save_credentials_repairs_existing_directory_permissions(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    credentials: StoredCredentials,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    directory = get_credentials_path().parent
    directory.mkdir(parents=True)
    directory.chmod(0o777)

    path = save_credentials(credentials)

    assert stat.S_IMODE(directory.stat().st_mode) == 0o700
    assert path.exists()


@pytest.mark.skipif(os.name != "posix", reason="POSIX ownership required")
def test_save_credentials_rejects_directory_owned_by_another_user(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    credentials: StoredCredentials,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    path = get_credentials_path()
    path.parent.mkdir(parents=True)
    owner_uid = path.parent.stat().st_uid
    monkeypatch.setattr("cytetype.config.os.getuid", lambda: owner_uid + 1)

    with pytest.raises(PermissionError, match="not owned by the current user"):
        save_credentials(credentials)

    assert not path.exists()


def test_invalid_credentials_file_has_clear_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    path = get_credentials_path()
    path.parent.mkdir(parents=True)
    path.write_text('{"apiUrl": "https://dev.example"}', encoding="utf-8")

    with pytest.raises(ValueError, match="Invalid CyteType credentials file"):
        load_credentials()


def test_existing_credentials_default_to_production_dashboard(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    path = get_credentials_path()
    path.parent.mkdir(parents=True)
    path.write_text(
        (
            '{"apiUrl":"https://dev.cytetype.example",'
            '"apiToken":"cyt_p_secret","tokenId":"token-id",'
            '"userId":"user-id","email":"researcher@university.edu"}'
        ),
        encoding="utf-8",
    )

    credentials = load_credentials()

    assert credentials is not None
    assert credentials.dashboardUrl == DEFAULT_DASHBOARD_URL


def test_help_and_get_key_alias(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["--help"])

    assert exit_info.value.code == 0
    output = capsys.readouterr().out
    assert "setup" in output
    assert "get-key" in output
    assert "login" in output
    assert "dashboard" in output
    assert "view" in output
    assert cli._build_parser().parse_args(["get-key"]).command == "get-key"
    assert cli._build_parser().parse_args(["setup"]).force is False
    assert cli._build_parser().parse_args(["get-key", "--force"]).force is True


def test_setup_api_url_argument_overrides_environment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv(
        "CYTETYPE_API_URL",
        "https://dev.cytetype.example/",
    )
    parser = cli._build_parser()

    assert parser.parse_args(["setup"]).api_url == "https://dev.cytetype.example"
    assert (
        parser.parse_args(["setup", "--api-url", "https://explicit.example"]).api_url
        == "https://explicit.example"
    )
    assert parser.parse_args(["login"]).api_url == "https://dev.cytetype.example"


def test_get_key_force_reaches_setup(
    monkeypatch: pytest.MonkeyPatch,
    credentials: StoredCredentials,
) -> None:
    calls: list[tuple[str, bool]] = []

    def fake_setup(api_url: str, force: bool = False) -> StoredCredentials:
        calls.append((api_url, force))
        return credentials

    monkeypatch.setattr(cli, "_run_setup", fake_setup)

    assert (
        cli.main(
            [
                "get-key",
                "--api-url",
                "https://dev.cytetype.example",
                "--force",
            ]
        )
        == 0
    )
    assert calls == [("https://dev.cytetype.example", True)]


def test_dashboard_url_uses_api_origin_only_when_origins_differ(
    credentials: StoredCredentials,
) -> None:
    assert (
        resolve_dashboard_url(credentials)
        == "https://dev.cytetype.example/dashboard"
    )

    same_origin = credentials.model_copy(
        update={
            "apiUrl": "https://api.cytetype.example",
            "dashboardUrl": "https://api.cytetype.example/custom-dashboard",
        }
    )
    assert (
        resolve_dashboard_url(same_origin)
        == "https://api.cytetype.example/custom-dashboard"
    )


def test_setup_always_shows_nygen_banner_and_links(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    credentials: StoredCredentials,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    save_credentials(credentials)

    response = SimpleNamespace(
        ok=True,
        json=lambda: {
            "tokenId": "refreshed-token-id",
            "userId": "refreshed-user-id",
            "email": "current@university.edu",
            "dashboardUrl": "https://dev.cytetype.example/dashboard",
        },
    )

    def fake_get(
        url: str,
        headers: dict[str, str],
        timeout: int,
    ) -> SimpleNamespace:
        assert url == "https://dev.cytetype.example/auth/cli/credentials"
        assert headers == {"Authorization": f"Bearer {credentials.apiToken}"}
        assert timeout == 30
        return response

    def fail_browser_open(url: str) -> bool:
        pytest.fail(f"Existing setup should not open a browser: {url}")

    monkeypatch.setattr(cli.requests, "get", fake_get)
    monkeypatch.setattr(cli, "_open_browser_silently", fail_browser_open)

    result = cli._run_setup(credentials.apiUrl)

    output = capsys.readouterr().out
    assert result.apiToken == credentials.apiToken
    assert result.tokenId == "refreshed-token-id"
    assert result.userId == "refreshed-user-id"
    assert result.email == "current@university.edu"
    assert load_credentials() == result
    assert cli._NYGEN_GLYPHS["n"] == (
        "      ",
        "# ### ",
        "##   #",
        "#    #",
        "#    #",
        "#    #",
        "#    #",
    )
    assert len(cli._NYGEN_BANNER) == 7
    assert all(len(line) <= 72 for line in cli._NYGEN_BANNER)
    assert "\n".join(cli._NYGEN_BANNER) in output
    assert (
        f"Nygen Analytics:      {cli._ANSI_BLUE}https://nygen.io{cli._ANSI_RESET}"
        in output
    )
    assert (
        "Dashboard:  "
        f"{cli._ANSI_BLUE}https://dev.cytetype.example/dashboard"
        f"{cli._ANSI_RESET}" in output
    )
    assert credentials.apiToken not in output


@pytest.mark.parametrize("error_code", ["INVALID_TOKEN", "TOKEN_INACTIVE"])
def test_setup_rejects_invalid_saved_credentials_without_changing_them(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    credentials: StoredCredentials,
    error_code: str,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    save_credentials(credentials)

    response = SimpleNamespace(
        ok=False,
        json=lambda: {
            "detail": {
                "error_code": error_code,
                "message": "Saved API key is no longer valid",
            }
        },
    )

    monkeypatch.setattr(cli.requests, "get", lambda *args, **kwargs: response)
    monkeypatch.setattr(
        cli,
        "_open_browser_silently",
        lambda url: pytest.fail(f"Invalid setup should not open a browser: {url}"),
    )

    with pytest.raises(RuntimeError, match=r"cytetype setup --force"):
        cli._run_setup(credentials.apiUrl)

    assert load_credentials() == credentials


def test_setup_reports_validation_network_errors_without_changing_credentials(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    credentials: StoredCredentials,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    save_credentials(credentials)

    def fail_validation(*args: object, **kwargs: object) -> None:
        raise cli.requests.ConnectionError("offline")

    monkeypatch.setattr(cli.requests, "get", fail_validation)
    monkeypatch.setattr(
        cli,
        "_open_browser_silently",
        lambda url: pytest.fail(f"Failed validation should not open a browser: {url}"),
    )

    with pytest.raises(RuntimeError, match="Could not validate the saved API key"):
        cli._run_setup(credentials.apiUrl)

    assert load_credentials() == credentials


def test_setup_completes_callback_exchange_without_exposing_key(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    credentials: StoredCredentials,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    save_credentials(credentials)
    opened_urls: list[str] = []
    callback_thread: threading.Thread | None = None
    exchange_body: dict[str, str] = {}
    callback_page: dict[str, str] = {}

    class FakeResponse:
        ok = True

        @staticmethod
        def json() -> dict[str, str]:
            return {
                "apiToken": "cyt_p_returned_secret",
                "tokenId": "server-token-id",
                "userId": "server-user-id",
                "email": "researcher@university.edu",
                "dashboardUrl": "https://dashboard.cytetype.example/dashboard",
            }

    def fake_post(
        url: str,
        json: dict[str, str],
        timeout: int,
    ) -> FakeResponse:
        assert url == "https://dev.cytetype.example/auth/cli/token"
        assert timeout == 30
        exchange_body.update(json)
        return FakeResponse()

    def fake_browser_open(url: str) -> bool:
        nonlocal callback_thread
        opened_urls.append(url)
        query = parse_qs(urlparse(url).query)
        callback_url = (
            f"{query['redirectUri'][0]}?code=signed-code&state={query['state'][0]}"
        )

        def call_back() -> None:
            with urlopen(callback_url, timeout=5) as response:
                assert response.status == 200
                callback_page["body"] = response.read().decode("utf-8")
                callback_page["cacheControl"] = response.headers["Cache-Control"]
                callback_page["referrerPolicy"] = response.headers["Referrer-Policy"]
                callback_page["contentSecurityPolicy"] = response.headers[
                    "Content-Security-Policy"
                ]

        callback_thread = threading.Thread(target=call_back)
        callback_thread.start()
        return True

    monkeypatch.setattr(cli.requests, "post", fake_post)
    monkeypatch.setattr(
        cli.requests,
        "get",
        lambda *args, **kwargs: pytest.fail(
            "Forced setup must not validate saved credentials"
        ),
    )
    monkeypatch.setattr(cli, "_open_browser_silently", fake_browser_open)

    result = cli._run_setup("https://dev.cytetype.example", force=True)

    assert callback_thread is not None
    callback_thread.join(timeout=5)
    assert not callback_thread.is_alive()
    assert result.apiToken == "cyt_p_returned_secret"
    assert exchange_body["code"] == "signed-code"
    assert exchange_body["redirectUri"].startswith("http://127.0.0.1:")
    assert exchange_body["codeVerifier"] not in opened_urls[0]
    assert "CyteType setup is complete." in callback_page["body"]
    assert "researcher@university.edu" in callback_page["body"]
    assert (
        'content="5;url=https://dev.cytetype.example/dashboard"'
        in callback_page["body"]
    )
    assert 'href="https://dev.cytetype.example/dashboard"' in callback_page["body"]
    assert "cyt_p_returned_secret" not in callback_page["body"]
    assert "signed-code" not in callback_page["body"]
    assert callback_page["cacheControl"] == "no-store"
    assert callback_page["referrerPolicy"] == "no-referrer"
    assert "default-src 'none'" in callback_page["contentSecurityPolicy"]
    assert load_credentials() == result
    output = capsys.readouterr().out
    assert "cyt_p_returned_secret" not in output
    assert result.dashboardUrl == "https://dashboard.cytetype.example/dashboard"
    assert (
        "Dashboard:  "
        f"{cli._ANSI_BLUE}https://dev.cytetype.example/dashboard"
        f"{cli._ANSI_RESET}" in output
    )
    assert "If it does not open automatically, use this URL:" in output
    assert f"{cli._ANSI_BLUE}https://dev.cytetype.example/auth/cli/authorize?" in output


def test_setup_times_out_without_callback(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    monkeypatch.setattr(cli, "_SETUP_TIMEOUT_SECONDS", 0.01)
    monkeypatch.setattr(cli, "_open_browser_silently", lambda url: False)

    with pytest.raises(TimeoutError, match="timed out"):
        cli._run_setup("https://dev.cytetype.example")

    output = capsys.readouterr().out
    assert (
        f"{cli._ANSI_RED}Was not able to launch web browser{cli._ANSI_RESET}" in output
    )
    assert "If it does not open automatically, use this URL:" in output
    assert f"{cli._ANSI_BLUE}https://dev.cytetype.example/auth/cli/authorize?" in output


def test_login_validates_hidden_api_key_before_saving(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    api_token = "cyt_p_existing_secret"
    request_data: dict[str, object] = {}

    class FakeResponse:
        ok = True

        @staticmethod
        def json() -> dict[str, str]:
            return {
                "tokenId": "existing-token-id",
                "userId": "existing-user-id",
                "email": "researcher@university.edu",
                "dashboardUrl": "https://dashboard.cytetype.example/dashboard",
            }

    def fake_get(
        url: str,
        headers: dict[str, str],
        timeout: int,
    ) -> FakeResponse:
        request_data.update(url=url, headers=headers, timeout=timeout)
        return FakeResponse()

    monkeypatch.setattr(cli.getpass, "getpass", lambda prompt: api_token)
    monkeypatch.setattr(cli.requests, "get", fake_get)

    credentials = cli._run_login("https://dev.cytetype.example")

    assert request_data == {
        "url": "https://dev.cytetype.example/auth/cli/credentials",
        "headers": {"Authorization": f"Bearer {api_token}"},
        "timeout": 30,
    }
    assert load_credentials() == credentials
    assert credentials.email == "researcher@university.edu"
    output = capsys.readouterr().out
    assert api_token not in output
    assert "Signed in as researcher@university.edu." in output
    assert (
        "Dashboard:  "
        f"{cli._ANSI_BLUE}https://dev.cytetype.example/dashboard"
        f"{cli._ANSI_RESET}" in output
    )


def test_login_failure_preserves_existing_credentials(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    credentials: StoredCredentials,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    save_credentials(credentials)

    class FakeResponse:
        ok = False

        @staticmethod
        def json() -> dict[str, dict[str, str]]:
            return {"detail": {"message": "Invalid API token"}}

    def fake_get(*args: object, **kwargs: object) -> FakeResponse:
        return FakeResponse()

    monkeypatch.setattr(cli.getpass, "getpass", lambda prompt: "invalid-token")
    monkeypatch.setattr(cli.requests, "get", fake_get)

    with pytest.raises(RuntimeError, match="Invalid API token"):
        cli._run_login("https://dev.cytetype.example")

    assert load_credentials() == credentials


def test_browser_launcher_suppresses_diagnostics(
    monkeypatch: pytest.MonkeyPatch,
    capfd: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.delenv("WSL_DISTRO_NAME", raising=False)
    monkeypatch.delenv("WSL_INTEROP", raising=False)
    monkeypatch.setattr(cli.sys, "platform", "win32")

    def noisy_browser_open(url: str) -> bool:
        os.write(2, f"launcher failed for {url}\n".encode())
        return False

    monkeypatch.setattr(cli.webbrowser, "open", noisy_browser_open)

    assert cli._open_browser_silently("https://dev.cytetype.example") is False
    assert capfd.readouterr().err == ""


def test_wsl_browser_launcher_is_detached_from_terminal(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    launched: list[tuple[list[str], dict[str, object]]] = []

    def fake_which(command: str) -> str | None:
        assert command == "rundll32.exe"
        return "/mnt/c/WINDOWS/system32/rundll32.exe"

    def fake_popen(command: list[str], **options: object) -> None:
        launched.append((command, options))

    monkeypatch.setenv("WSL_DISTRO_NAME", "Ubuntu-24.04")
    monkeypatch.setattr(cli.shutil, "which", fake_which)
    monkeypatch.setattr(cli.subprocess, "Popen", fake_popen)

    url = "https://dev.cytetype.example/auth/cli/authorize?state=test"
    assert cli._open_browser_silently(url) is True
    assert launched == [
        (
            [
                "/mnt/c/WINDOWS/system32/rundll32.exe",
                "url.dll,FileProtocolHandler",
                url,
            ],
            {
                "stdin": cli.subprocess.DEVNULL,
                "stdout": cli.subprocess.DEVNULL,
                "stderr": cli.subprocess.DEVNULL,
                "start_new_session": True,
            },
        )
    ]


def test_wsl_browser_does_not_fall_back_to_terminal_browser(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_browser_open(url: str) -> bool:
        pytest.fail(f"WSL must not fall back to a terminal browser: {url}")

    monkeypatch.setenv("WSL_DISTRO_NAME", "Ubuntu-24.04")
    monkeypatch.setattr(cli.shutil, "which", lambda command: None)
    monkeypatch.setattr(cli.webbrowser, "open", fail_browser_open)

    assert cli._open_browser_silently("https://dev.cytetype.example") is False


def test_setup_keyboard_interrupt_exits_without_traceback(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def interrupt_setup(api_url: str, force: bool = False) -> None:
        raise KeyboardInterrupt

    monkeypatch.setattr(cli, "_run_setup", interrupt_setup)

    assert cli.main(["setup"]) == 130
    assert capsys.readouterr().err == "\nCancelled.\n"


def test_setup_rejects_callback_with_wrong_state(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    monkeypatch.setattr(cli, "_SETUP_TIMEOUT_SECONDS", 0.1)
    callback_thread: threading.Thread | None = None
    callback_page: dict[str, str] = {}

    def fail_post(*args: object, **kwargs: object) -> None:
        pytest.fail("A callback with the wrong state must not exchange a token")

    def fake_browser_open(url: str) -> bool:
        nonlocal callback_thread
        query = parse_qs(urlparse(url).query)
        callback_url = f"{query['redirectUri'][0]}?code=signed-code&state=wrong"

        def call_back() -> None:
            with pytest.raises(HTTPError) as error_info:
                urlopen(callback_url, timeout=5)
            response = error_info.value
            assert response.code == 400
            callback_page["body"] = response.read().decode("utf-8")
            response.close()

        callback_thread = threading.Thread(target=call_back)
        callback_thread.start()
        return True

    monkeypatch.setattr(cli.requests, "post", fail_post)
    monkeypatch.setattr(cli, "_open_browser_silently", fake_browser_open)

    with pytest.raises(TimeoutError, match="timed out"):
        cli._run_setup("https://dev.cytetype.example")

    assert callback_thread is not None
    callback_thread.join(timeout=5)
    assert not callback_thread.is_alive()
    assert "Authorization state did not match." in callback_page["body"]
    assert 'class="card error"' in callback_page["body"]
    assert 'http-equiv="refresh"' not in callback_page["body"]


def test_dashboard_view_and_logout_use_saved_server(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    credentials: StoredCredentials,
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    save_credentials(credentials)
    opened_urls: list[str] = []

    def open_url(url: str) -> bool:
        opened_urls.append(url)
        return True

    monkeypatch.setattr(cli, "_open_browser_silently", open_url)

    assert cli.main(["dashboard"]) == 0
    assert cli.main(["view", "job/with space"]) == 0

    report_redirect = parse_qs(urlparse(opened_urls[1]).query)["redirect"][0]
    assert opened_urls[0] == "https://dev.cytetype.example/dashboard"
    assert report_redirect == "/report/job%2Fwith%20space"
    assert "cyt_p_secret" not in "".join(opened_urls)
    assert (
        f"{cli._ANSI_BLUE}https://dev.cytetype.example/dashboard"
        f"{cli._ANSI_RESET}" in capsys.readouterr().out
    )

    assert cli.main(["logout"]) == 0
    assert load_credentials() is None


@pytest.mark.parametrize(
    "api_url",
    [
        "http://example.com",
        "https://user@example.com",
        "https://example.com/path",
        "ftp://example.com",
    ],
)
def test_validate_api_url_rejects_unsafe_api_url(api_url: str) -> None:
    with pytest.raises(ValueError):
        validate_api_url(api_url)
