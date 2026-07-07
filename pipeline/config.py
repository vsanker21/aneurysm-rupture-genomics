"""Pipeline configuration and path helpers."""
from __future__ import annotations

import os
import sys
from pathlib import Path

PIPELINE_VERSION = "1.0.0"

# Ensure repo root is importable when running `python -m pipeline` or `run_pipeline.py`.
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from project_paths import get_project_root  # noqa: E402


def get_python() -> str:
    return os.environ.get("PIPELINE_PYTHON", sys.executable)


def results_dir() -> Path:
    return get_project_root() / "results"


def data_dir() -> Path:
    return get_project_root() / "data"
