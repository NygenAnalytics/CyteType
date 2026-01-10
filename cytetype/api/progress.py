import sys


class ProgressDisplay:
    """Manages terminal progress display during job polling."""

    # Class constants
    COLORS = {
        "completed": "\033[92m",
        "processing": "\033[93m",
        "pending": "\033[94m",
        "failed": "\033[91m",
        "reset": "\033[0m",
    }
    SYMBOLS = {"completed": "✓", "processing": "⟳", "pending": "○", "failed": "✗"}
    SPINNER_CHARS = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]

    def __init__(self) -> None:
        self.spinner_frame = 0
        self.last_status: dict[str, str] = {}

    def update(self, cluster_status: dict[str, str]) -> None:
        """Update progress display with current cluster status."""
        if not cluster_status:
            return

        # Always render to keep spinner animating
        self._render(cluster_status, is_final=False)

        # Track last status for potential future use
        if cluster_status != self.last_status:
            self.last_status = cluster_status.copy()

        # Always increment spinner to show activity
        self.spinner_frame += 1

    def finalize(self, cluster_status: dict[str, str]) -> None:
        """Show final status and cleanup."""
        if cluster_status:
            self._render(cluster_status, is_final=True)
        print()  # Ensure newline

    def _render(self, cluster_status: dict[str, str], is_final: bool) -> None:
        """Render status to terminal."""
        status_counts = self._count_statuses(cluster_status)
        progress_bar = self._build_progress_bar(cluster_status)
        status_line = self._build_status_line(
            progress_bar, status_counts, is_final=is_final
        )

        # Print status line
        if is_final:
            print(f"\r{status_line}{self.COLORS['reset']}")
            sys.stdout.flush()
            self._show_failed_clusters(cluster_status, status_counts["failed"])
        else:
            print(f"\r{status_line}{self.COLORS['reset']}", end="", flush=True)

    def _count_statuses(self, cluster_status: dict[str, str]) -> dict[str, int]:
        """Count occurrences of each status."""
        counts = {"completed": 0, "failed": 0}
        for status in cluster_status.values():
            counts[status] = counts.get(status, 0) + 1
        return counts

    def _build_progress_bar(self, cluster_status: dict[str, str]) -> str:
        """Build colored progress bar from cluster statuses."""
        progress_units = []
        for cluster_id in self._sorted_cluster_ids(cluster_status):
            status = cluster_status[cluster_id]
            color = self.COLORS.get(status, self.COLORS["reset"])
            symbol = self.SYMBOLS.get(status, "?")
            progress_units.append(f"{color}{symbol}{self.COLORS['reset']}")
        return "".join(progress_units)

    def _build_status_line(
        self, progress_bar: str, counts: dict[str, int], is_final: bool
    ) -> str:
        """Build status line with progress bar and counts."""
        total = sum(counts.values())
        completed = counts["completed"]

        if is_final:
            status_line = f"[DONE] [{progress_bar}] {completed}/{total}"
            if counts["failed"] > 0:
                status_line += f" ({counts['failed']} failed)"
            elif completed == total:
                status_line += " completed"
        else:
            spinner = self.SPINNER_CHARS[self.spinner_frame % len(self.SPINNER_CHARS)]
            status_line = f"{spinner} [{progress_bar}] {completed}/{total} completed"

        return status_line

    def _show_failed_clusters(
        self, cluster_status: dict[str, str], failed_count: int
    ) -> None:
        """Show details of failed clusters."""
        if failed_count == 0:
            return

        failed_details = []
        for cluster_id in self._sorted_cluster_ids(cluster_status):
            if cluster_status[cluster_id] == "failed":
                color = self.COLORS["failed"]
                symbol = self.SYMBOLS["failed"]
                failed_details.append(
                    f"{color}{symbol} Cluster {cluster_id}{self.COLORS['reset']}"
                )

        # Group into lines of 4
        for i in range(0, len(failed_details), 4):
            print(f"  {' | '.join(failed_details[i : i + 4])}")

    @staticmethod
    def _sorted_cluster_ids(cluster_status: dict[str, str]) -> list[str]:
        """Sort cluster IDs naturally (numeric if possible)."""
        return sorted(cluster_status.keys(), key=lambda x: int(x) if x.isdigit() else x)
