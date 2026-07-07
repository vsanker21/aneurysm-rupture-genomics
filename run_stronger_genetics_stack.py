"""
Execute the stronger-genetics layer (executable without R):

  1) prepare_stronger_genetics_followup.py --from-shortlist --no-ensembl
     (regional Bakker GWAS + manifest; fast sentinel-centered windows)
  2) stronger_genetics_eqtl_fetch.py
     (per-locus blood eQTL rows in the same windows)
  3) stronger_genetics_zscreen.py
     (Spearman correlation of GWAS z vs eQTL z at overlapping rsIDs — screen only)

Formal COLOC / SMR: install R + coloc, then see RUNBOOK_stronger_genetics.txt

Usage:
  python run_stronger_genetics_stack.py
  python run_stronger_genetics_stack.py --max-loci-prep 300 --max-loci-eqtl 50
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

BASE = Path(r"e:\Aneurysms Genomics Project")
PY = sys.executable


def run(script: str, args: list[str] | None = None) -> None:
    cmd = [PY, str(BASE / script)]
    if args:
        cmd.extend(args)
    print("\n>>>", " ".join(cmd))
    subprocess.run(cmd, cwd=BASE, check=True)


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--max-loci-prep",
        type=int,
        default=300,
        help="Rows from shortlist with rsID to consider (Bakker filter inside prep).",
    )
    ap.add_argument("--window-kb", type=int, default=500)
    ap.add_argument(
        "--max-loci-eqtl",
        type=int,
        default=40,
        help="Cap Catalogue fetches (each locus hits 3 blood datasets).",
    )
    ap.add_argument(
        "--max-pages-per-dataset",
        type=int,
        default=4,
        help="Pagination cap per blood dataset (lower = faster, fewer variants).",
    )
    ap.add_argument("--skip-eqtl-fetch", action="store_true")
    ap.add_argument("--skip-zscreen", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument(
        "--tier1-tsv",
        type=str,
        default="results/phases_g/tier1_science_supplement.tsv",
        help="When set, eQTL fetch targets tier-1 loci from this TSV (rank order).",
    )
    ap.add_argument(
        "--tier1-top-n",
        type=int,
        default=10,
        help="When --tier1-tsv is used, number of tier-1 loci to fetch eQTL regions for.",
    )
    args = ap.parse_args()

    run(
        "prepare_stronger_genetics_followup.py",
        [
            "--from-shortlist",
            "--no-ensembl",
            "--max-loci",
            str(args.max_loci_prep),
            "--window-kb",
            str(args.window_kb),
        ],
    )

    if not args.skip_eqtl_fetch:
        run(
            "stronger_genetics_eqtl_fetch.py",
            [
                "--max-loci",
                str(args.max_loci_eqtl),
                "--tier1-tsv",
                str(args.tier1_tsv),
                "--tier1-top-n",
                str(args.tier1_top_n),
                "--max-pages-per-dataset",
                str(args.max_pages_per_dataset),
            ],
        )

    if not args.skip_zscreen:
        run("stronger_genetics_zscreen.py", [])

    if args.verify:
        min_gw = min(10, max(3, args.max_loci_prep // 10))
        # Tier-1-targeted fetch cannot exceed --max-loci-eqtl; verify against the smaller cap.
        min_eq = min(max(1, int(args.tier1_top_n)), max(1, int(args.max_loci_eqtl)))
        run(
            "verify_stronger_genetics_completion.py",
            [
                "--min-loci-gwas",
                str(min_gw),
                "--min-eqtl-files",
                str(min_eq),
            ],
        )

    print("\nDone. See results/phases_g/stronger_genetics/")


if __name__ == "__main__":
    main()
