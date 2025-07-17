"""
Display utilities for CyteType client progress and status visualization.
"""

import sys
from typing import Dict


def display_cluster_status(
    cluster_status: Dict[str, str],
    job_id: str,
    is_final: bool = False,
    spinner_frame: int = 0,
) -> None:
    """Display cluster status with colors and progress indicators.

    Args:
        cluster_status: Dictionary mapping cluster IDs to their status
        job_id: The job ID for context
        is_final: Whether this is the final status update (shows cluster details)
        spinner_frame: Current frame of the spinner animation (0-3)
    """
    # Color codes for terminal output
    colors = {
        "completed": "\033[92m",  # Green
        "processing": "\033[93m",  # Yellow
        "pending": "\033[94m",  # Blue
        "failed": "\033[91m",  # Red
        "reset": "\033[0m",  # Reset
    }

    # Status symbols
    symbols = {"completed": "✓", "processing": "⟳", "pending": "○", "failed": "✗"}

    if not cluster_status:
        return

    # Count statuses
    status_counts: dict[str, int] = {}
    for status in cluster_status.values():
        status_counts[status] = status_counts.get(status, 0) + 1

    total_clusters = len(cluster_status)
    completed_count = status_counts.get("completed", 0)
    failed_count = status_counts.get("failed", 0)

    # Create progress bar with one unit per cluster in specific positions
    progress_units = []
    for cluster_id in sorted(
        cluster_status.keys(), key=lambda x: int(x) if x.isdigit() else x
    ):
        status = cluster_status[cluster_id]
        color = colors.get(status, colors["reset"])
        symbol = symbols.get(status, "?")
        progress_units.append(f"{color}{symbol}{colors['reset']}")

    progress_bar = "".join(progress_units)

    # Create spinner animation
    spinner_chars = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]
    spinner = spinner_chars[spinner_frame % len(spinner_chars)]

    # Create status line
    if is_final:
        # No spinner for final status
        status_line = f"[DONE] [{progress_bar}] {completed_count}/{total_clusters}"
        if total_clusters > completed_count and failed_count > 0:
            status_line += f" ({failed_count} failed)"
        elif total_clusters == completed_count:
            status_line += " completed"
    else:
        # Include spinner for progress updates
        status_line = (
            f"{spinner} [{progress_bar}] {completed_count}/{total_clusters} completed"
        )

    # Use carriage return to overwrite the same line unless it's final
    if is_final:
        # Final status gets a new line
        print(f"\r{status_line}{colors['reset']}")  # Ensure color reset
        sys.stdout.flush()

        # Show individual cluster details only if there are failures
        if failed_count > 0:
            cluster_details = []
            for cluster_id in sorted(
                cluster_status.keys(), key=lambda x: int(x) if x.isdigit() else x
            ):
                status = cluster_status[cluster_id]
                if status == "failed":  # Only show failed clusters
                    color = colors.get(status, colors["reset"])
                    symbol = symbols.get(status, "?")
                    cluster_details.append(
                        f"{color}{symbol} Cluster {cluster_id}{colors['reset']}"
                    )

            if cluster_details:
                # Group details into lines of 4 for better readability
                for i in range(0, len(cluster_details), 4):
                    line_details = cluster_details[i : i + 4]
                    print(f"  {' | '.join(line_details)}")
    else:
        # Progress updates overwrite the same line
        print(f"\r{status_line}{colors['reset']}", end="", flush=True)
