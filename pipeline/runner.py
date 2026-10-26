"""Execute pipeline stages as subprocesses with prerequisite checks."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from pipeline.config import get_project_root, get_python
from pipeline.stages import STAGES, Stage, check_requirements, should_skip_stage


class PipelineError(RuntimeError):
    pass


def _run_script(root: Path, script: str, extra_args: list[str] | None = None) -> None:
    path = root / script
    if not path.is_file():
        raise PipelineError(f"Script not found: {path}")
    cmd = [get_python(), str(path)]
    if extra_args:
        cmd.extend(extra_args)
    print(f"\n>>> {' '.join(cmd)}", flush=True)
    env = os.environ.copy()
    env["ANEURYSM_GENOMICS_ROOT"] = str(root)
    subprocess.run(cmd, cwd=str(root), env=env, check=True)


def run_stage(stage: Stage, *, root: Path | None = None, dry_run: bool = False) -> None:
    root = root or get_project_root()
    skip_reason = should_skip_stage(root, stage)
    if skip_reason:
        print(f"[skip] {stage.name}: {skip_reason}")
        return

    missing = check_requirements(root, stage)
    if missing and not stage.optional:
        raise PipelineError(
            f"Stage '{stage.name}' missing required inputs:\n  "
            + "\n  ".join(missing)
        )
    if missing and stage.optional:
        print(f"[skip] {stage.name}: missing {', '.join(missing)}")
        return

    print(f"\n=== Stage: {stage.name} — {stage.description} ===", flush=True)
    if dry_run:
        for script in stage.scripts:
            args = stage.script_args.get(script, [])
            print(f"  would run: {script} {' '.join(args)}".strip())
        return

    for script in stage.scripts:
        args = list(stage.script_args.get(script, []))
        _run_script(root, script, args)


def run_pipeline(
    stage_names: list[str],
    *,
    root: Path | None = None,
    dry_run: bool = False,
    stop_on_error: bool = True,
) -> None:
    root = root or get_project_root()
    os.environ["ANEURYSM_GENOMICS_ROOT"] = str(root)
    print(f"Project root: {root}")
    print(f"Stages: {', '.join(stage_names)}")

    for name in stage_names:
        stage = STAGES[name]
        try:
            run_stage(stage, root=root, dry_run=dry_run)
        except subprocess.CalledProcessError as exc:
            msg = f"Stage '{name}' failed (exit {exc.returncode})"
            if stop_on_error:
                raise PipelineError(msg) from exc
            print(f"[error] {msg}", file=sys.stderr)
        except PipelineError:
            if stop_on_error:
                raise
            print(f"[error] stage '{name}' failed", file=sys.stderr)

    print("\nPipeline finished.", flush=True)
