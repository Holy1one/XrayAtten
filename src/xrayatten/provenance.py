"""Provenance helpers for reproducible outputs."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import subprocess

from . import __version__
from .local_nist import PROJECT_ROOT


def generated_at_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def git_commit(root: Path = PROJECT_ROOT) -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return "unavailable"
    return result.stdout.strip() or "unavailable"


def runtime_metadata() -> dict[str, str]:
    return {
        "generated_at_utc": generated_at_utc(),
        "package_version": __version__,
        "git_commit": git_commit(),
    }

