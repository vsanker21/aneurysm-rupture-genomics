"""
Resolve the project root for portable paths (local clone, CI, or custom location).

Set ANEURYSM_GENOMICS_ROOT to override. Otherwise uses the directory containing
this file (repository root).
"""
from __future__ import annotations

import os
from pathlib import Path


def get_project_root() -> Path:
    env = os.environ.get("ANEURYSM_GENOMICS_ROOT", "").strip()
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parent


BASE = get_project_root()
