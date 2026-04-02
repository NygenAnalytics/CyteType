import io

import pytest

from cytetype.api import client
from cytetype.api import progress as progress_module
from cytetype.api.progress import ProgressDisplay


class _FakeStream(io.StringIO):
    def isatty(self) -> bool:
        return False


class _CapturedLogger:
    def __init__(self) -> None:
        self.messages: list[tuple[str, str]] = []

    def info(self, message: str) -> None:
        self.messages.append(("info", message))

    def success(self, message: str) -> None:
        self.messages.append(("success", message))

    def warning(self, message: str) -> None:
        self.messages.append(("warning", message))

    def debug(self, message: str) -> None:
        self.messages.append(("debug", message))


class _FakeDisplayHandle:
    def __init__(self, message: str) -> None:
        self.messages = [message]

    def update(self, message: str) -> None:
        self.messages.append(message)


def test_progress_display_plain_output_avoids_repeated_lines(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progress_module, "_in_notebook", lambda: False)

    stream = _FakeStream()
    progress = ProgressDisplay(stream=stream)

    progress.update("processing")
    progress.update("processing")
    progress.finalize("completed")
    progress.finalize("completed")

    assert stream.getvalue().splitlines() == [
        "CyteType job running...",
        "[DONE] CyteType job completed.",
    ]


def test_progress_display_updates_single_notebook_output(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progress_module, "_in_notebook", lambda: True)
    monkeypatch.setattr(
        progress_module, "_render_notebook_message", lambda message: message
    )

    display_handle: _FakeDisplayHandle | None = None

    def _create_handle(message: str) -> _FakeDisplayHandle:
        nonlocal display_handle
        display_handle = _FakeDisplayHandle(message)
        return display_handle

    monkeypatch.setattr(
        progress_module,
        "_create_notebook_display_handle",
        _create_handle,
    )

    progress = ProgressDisplay(stream=_FakeStream())
    progress.update("processing")
    progress.update("processing")
    progress.finalize("completed")

    assert display_handle is not None
    assert len(display_handle.messages) == 3
    assert display_handle.messages[0].startswith("⠋ CyteType job running...")
    assert display_handle.messages[1].startswith("⠙ CyteType job running...")
    assert display_handle.messages[2] == "[DONE] CyteType job completed."


def test_wait_for_completion_logs_report_cta_before_polling(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured_logger = _CapturedLogger()
    sleep_calls: list[int | float] = []
    result_payload: dict[str, list[object]] = {"annotations": []}

    monkeypatch.setattr(client, "logger", captured_logger)
    monkeypatch.setattr(
        client.time, "sleep", lambda seconds: sleep_calls.append(seconds)
    )
    monkeypatch.setattr(
        client,
        "get_job_status",
        lambda *args, **kwargs: {
            "jobStatus": "completed",
            "clusterStatus": {"0": "completed"},
        },
    )
    monkeypatch.setattr(
        client,
        "fetch_job_results",
        lambda *args, **kwargs: result_payload,
    )

    result = client.wait_for_completion(
        "https://example.com/",
        None,
        "job-123",
        poll_interval=5,
        timeout=30,
        show_progress=False,
    )

    assert result == result_payload
    assert sleep_calls == [5]
    assert captured_logger.messages[:3] == [
        ("info", "CyteType job submitted."),
        (
            "info",
            "If your session disconnects, results can still be fetched later with:\n"
            "`results = annotator.get_results()`",
        ),
        ("info", "\n[TRACK PROGRESS]\nhttps://example.com/report/job-123"),
    ]
    assert (
        "success",
        "Job job-123 completed successfully.",
    ) in captured_logger.messages
