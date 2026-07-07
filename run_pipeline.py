#!/usr/bin/env python3
"""
End-to-end rupture genomics analysis pipeline (no manuscript generation).

Usage:
  python run_pipeline.py --list-stages
  python run_pipeline.py --stage all
  python run_pipeline.py --stage deseq2,phases_g,blood_bridge,supplementary
  python run_pipeline.py --stage all --dry-run
  python run_pipeline.py --stage supplementary --skip-docx   # alias; docx always skipped here

Outputs land under results/ (tables, figures, shortlist, panel performance).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from pipeline.config import PIPELINE_VERSION, get_project_root
from pipeline.runner import PipelineError, run_pipeline
from pipeline.stages import STAGE_ORDER, STAGES, resolve_stage_list


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Aneurysm rupture genomics analysis pipeline (analysis outputs only)."
    )
    ap.add_argument(
        "--stage",
        default="all",
        help=f"Comma-separated stage names or 'all'. Available: {', '.join(STAGE_ORDER)}",
    )
    ap.add_argument("--list-stages", action="store_true", help="List stages and exit")
    ap.add_argument("--dry-run", action="store_true", help="Print commands without running")
    ap.add_argument(
        "--project-root",
        type=Path,
        default=None,
        help="Override project root (default: directory containing run_pipeline.py)",
    )
    ap.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue after a stage fails",
    )
    args = ap.parse_args()

    if args.list_stages:
        print(f"Pipeline v{PIPELINE_VERSION}\n")
        for name in STAGE_ORDER:
            s = STAGES[name]
            opt = " [optional]" if s.optional else ""
            print(f"  {name:20} {s.description}{opt}")
        print("\nManuscript Word generation is NOT part of this pipeline.")
        print("Outputs: results/phases_g/, results/blood_bridge_panel/, results/supplementary/")
        return 0

    root = args.project_root.resolve() if args.project_root else get_project_root()
    try:
        stages = resolve_stage_list(args.stage)
        run_pipeline(
            stages,
            root=root,
            dry_run=args.dry_run,
            stop_on_error=not args.continue_on_error,
        )
    except (PipelineError, ValueError) as exc:
        print(f"Pipeline error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
