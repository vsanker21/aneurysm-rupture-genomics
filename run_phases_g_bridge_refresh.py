"""
One command: export GSE122897 FDR genes → symbol TSV, link to PHASES-G shortlist, write bridge + summary JSON.

Does not regenerate the Word document (run create_document_with_inferences.py after).

Usage:
  python run_phases_g_bridge_refresh.py
  python run_phases_g_bridge_refresh.py --max-padj 0.05 --skip-export
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
PY = sys.executable
COHORT_DE = BASE / "results/phases_g/cohort_de_gse122897_for_phases_g_bridge.tsv"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-padj", type=float, default=0.05)
    ap.add_argument("--skip-export", action="store_true", help="Only run link_tissue_de_to_phases_g.py")
    ap.add_argument("--skip-mygene", action="store_true", help="Pass to export (local symbol maps only).")
    args = ap.parse_args()

    if not args.skip_export:
        cmd = [
            PY,
            str(BASE / "export_gse122897_cohort_de_for_bridge.py"),
            "--max-padj",
            str(args.max_padj),
        ]
        if args.skip_mygene:
            cmd.append("--skip-mygene")
        r = subprocess.run(cmd, cwd=str(BASE))
        if r.returncode != 0:
            return r.returncode

    if not COHORT_DE.is_file():
        print(f"Missing {COHORT_DE}; run without --skip-export.", file=sys.stderr)
        return 1

    r = subprocess.run(
        [PY, str(BASE / "link_tissue_de_to_phases_g.py"), "--de-tsv", str(COHORT_DE)],
        cwd=str(BASE),
    )
    if r.returncode != 0:
        return r.returncode

    print("\nNext: python run_pipeline.py --stage blood_bridge,supplementary")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
