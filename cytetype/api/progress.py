import sys
import time
from html import escape
from typing import Any, Callable, TextIO, cast


def _in_notebook() -> bool:
    try:
        from IPython import get_ipython
    except ImportError:
        return False

    shell_getter = cast(Callable[[], Any | None], get_ipython)
    shell = shell_getter()
    return bool(shell and shell.__class__.__name__ == "ZMQInteractiveShell")


def _create_notebook_display_handle(message: str) -> Any | None:
    try:
        from IPython.display import display
    except ImportError:
        return None

    _display: Callable[..., Any] = display
    return _display(_render_notebook_message(message), display_id=True)


def _render_notebook_message(message: str) -> Any:
    try:
        from IPython.display import HTML
    except ImportError:
        return message

    html_cls = cast(Callable[[str], Any], HTML)
    return html_cls(
        "<pre style='margin: 0; white-space: pre-wrap; font-family: monospace;'>"
        f"{escape(message)}"
        "</pre>"
    )


class ProgressDisplay:
    """Manages terminal and notebook progress display during job polling."""

    COLORS = {"failed": "\033[91m", "reset": "\033[0m"}
    SPINNER_CHARS = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]

    def __init__(self, stream: TextIO | None = None) -> None:
        self.stream = stream or sys.stdout
        self._interactive = bool(
            hasattr(self.stream, "isatty") and self.stream.isatty()
        )
        self._use_notebook_display = not self._interactive and _in_notebook()
        self._display_handle: Any | None = None
        self._finalized = False
        self._last_plain_status: str | None = None
        self._start_time = time.monotonic()
        self.spinner_frame = 0

    def update(self, job_status: str) -> None:
        """Update progress display with the overall job status."""
        if self._finalized:
            return

        if self._interactive:
            message = self._build_running_line(job_status)
            print(f"\r{message}", end="", file=self.stream, flush=True)
        elif self._use_notebook_display:
            self._update_notebook_display(self._build_running_line(job_status))
        else:
            message = self._build_plain_line(job_status)
            if message != self._last_plain_status:
                print(message, file=self.stream, flush=True)
                self._last_plain_status = message

        self.spinner_frame += 1

    def finalize(
        self,
        final_status: str | None = None,
        cluster_status: dict[str, str] | None = None,
    ) -> None:
        """Show final status and cleanup."""
        if self._finalized:
            return
        self._finalized = True

        if final_status is None:
            if self._interactive:
                print(file=self.stream, flush=True)
            return

        message = self._build_final_line(final_status)
        if self._interactive:
            print(f"\r{message}", file=self.stream, flush=True)
        elif self._use_notebook_display:
            self._update_notebook_display(message)
        else:
            print(message, file=self.stream, flush=True)

        if final_status == "failed" and cluster_status:
            self._show_failed_clusters(cluster_status)

    def _update_notebook_display(self, message: str) -> None:
        """Update a single notebook output cell instead of printing many lines."""
        if self._display_handle is None:
            self._display_handle = _create_notebook_display_handle(message)
            if self._display_handle is None:
                self._use_notebook_display = False
                print(message, file=self.stream, flush=True)
            return

        self._display_handle.update(_render_notebook_message(message))

    def _build_running_line(self, job_status: str) -> str:
        spinner = self.SPINNER_CHARS[self.spinner_frame % len(self.SPINNER_CHARS)]
        elapsed = self._format_elapsed()
        return f"{spinner} {self._status_message(job_status)} {elapsed} elapsed"

    def _build_plain_line(self, job_status: str) -> str:
        return self._status_message(job_status)

    @staticmethod
    def _build_final_line(final_status: str) -> str:
        if final_status == "completed":
            return "[DONE] CyteType job completed."
        if final_status == "failed":
            return "[FAILED] CyteType job failed."
        if final_status == "timed_out":
            return "[TIMEOUT] CyteType job timed out."
        return "[STOPPED] CyteType job stopped."

    @staticmethod
    def _status_message(job_status: str) -> str:
        if job_status == "pending":
            return "CyteType job queued..."
        if job_status == "processing":
            return "CyteType job running..."
        if job_status == "not_found":
            return "Waiting for CyteType job to start..."
        return "Waiting for CyteType results..."

    def _format_elapsed(self) -> str:
        elapsed = int(time.monotonic() - self._start_time)
        minutes, seconds = divmod(elapsed, 60)
        return f"{minutes:02d}:{seconds:02d}"

    def _show_failed_clusters(self, cluster_status: dict[str, str]) -> None:
        """Show details of failed clusters."""
        failed_details = []
        for cluster_id in self._sorted_cluster_ids(cluster_status):
            if cluster_status[cluster_id] != "failed":
                continue

            if self._interactive:
                failed_details.append(
                    f"{self.COLORS['failed']}✗ Cluster {cluster_id}{self.COLORS['reset']}"
                )
            else:
                failed_details.append(f"✗ Cluster {cluster_id}")

        for i in range(0, len(failed_details), 4):
            print(f"  {' | '.join(failed_details[i : i + 4])}", file=self.stream)

    @staticmethod
    def _sorted_cluster_ids(cluster_status: dict[str, str]) -> list[str]:
        """Sort cluster IDs naturally (numeric if possible)."""
        return sorted(cluster_status.keys(), key=lambda x: int(x) if x.isdigit() else x)
