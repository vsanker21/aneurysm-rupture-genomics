"""
Run PHASES-G evidence downloads + immune-adjusted DE + integration (single command).

Steps:
  1. download_phases_g_public_data.py  (Bakker GWAS, eQTL blood index, GSE122897 TPM+raw;
     verifies gzip readability; exits non-zero on failure)
  2. prepare_cibersortx_mixture_gse122897.py  (if mixture matrix missing)
  3. immune_deconvolution_lm22_gse122897.py  (if lm22_fractions missing)
  4. deg_immune_adjusted_gse122897.py  (requires results/gse122897_deseq2_results.tsv)
  5. run_phases_g_evidence_pipeline.py
  6. append_mouse_scrna_stopgap_to_shortlist.py  (optional stopgap; adds mouse compartment evidence)

Does not run CIBERSORTx web jobs.

Usage:
  python run_phases_g_full_stack.py
  python run_phases_g_full_stack.py --max-genes 200
  python run_phases_g_full_stack.py --verify   # run verify_phases_g_objectives.py at end

After a successful full run (with eQTL), refresh tier-1 table and check objectives:
  python phases_g_tier1_supplement.py --top-n 40
  python verify_phases_g_objectives.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

try:
    from project_paths import get_project_root
    BASE = get_project_root()
except ImportError:
    BASE = Path(__file__).resolve().parent
PY = sys.executable


def run(script: str, args: list[str] | None = None) -> None:
    cmd = [PY, str(BASE / script)]
    if args:
        cmd.extend(args)
    print("\n>>>", " ".join(cmd))
    subprocess.run(cmd, cwd=BASE, check=True)


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--max-genes", type=int, default=300)
    ap.add_argument("--skip-download", action="store_true")
    ap.add_argument(
        "--include-lineage-ge-contexts",
        action="store_true",
        help="More blood-lineage GE datasets in eQTL index (download step).",
    )
    ap.add_argument(
        "--eqtl-pages-per-dataset",
        type=int,
        default=4,
        help="Pass through to evidence pipeline.",
    )
    ap.add_argument(
        "--eqtl-timeout",
        type=int,
        default=240,
        help="Seconds per eQTL Catalogue request (matches evidence pipeline default).",
    )
    ap.add_argument(
        "--allow-missing-gwas",
        action="store_true",
        help="If Bakker sumstats not present, skip IA overlap.",
    )
    ap.add_argument(
        "--skip-eqtl",
        action="store_true",
        help="Evidence step: skip eQTL Catalogue (fast; replication + markers only).",
    )
    ap.add_argument(
        "--verify",
        action="store_true",
        help="Run verify_phases_g_objectives.py after the evidence step (exit non-zero if incomplete).",
    )
    args = ap.parse_args()

    if not args.skip_download:
        dl_args = (
            ["--include-lineage-ge-contexts"]
            if args.include_lineage_ge_contexts
            else None
        )
        run("download_phases_g_public_data.py", dl_args)

    deseq = BASE / "results/gse122897_deseq2_results.tsv"
    if not deseq.is_file():
        raise SystemExit(
            f"Missing {deseq}. Generate it with run_deseq2_gse122897.R (or your DESeq2 pipeline) "
            "after raw counts and sample metadata are in place."
        )

    mixture = BASE / "results/cibersortx/GSE122897_CIBERSORTx_mixture_TPM.txt"
    if not mixture.is_file():
        run("prepare_cibersortx_mixture_gse122897.py")

    lm22 = BASE / "results/cibersortx/lm22_fractions_gse122897.tsv"
    if not lm22.is_file():
        run("immune_deconvolution_lm22_gse122897.py")

    run("deg_immune_adjusted_gse122897.py", ["--max-genes", "15000"])
    ev_args = [
        "--max-genes",
        str(args.max_genes),
        "--eqtl-pages-per-dataset",
        str(args.eqtl_pages_per_dataset),
        "--eqtl-timeout",
        str(args.eqtl_timeout),
        "--append-mouse-stopgap",
    ]
    if args.allow_missing_gwas:
        ev_args.append("--allow-missing-gwas")
    if args.skip_eqtl:
        ev_args.append("--skip-eqtl")
    run("run_phases_g_evidence_pipeline.py", ev_args)
    # Optional stopgap: append mouse scRNA compartment evidence (GSE193533) to the shortlist TSV.
    # This is supportive only (species/model/coarse annotation).
    try:
        run("append_mouse_scrna_stopgap_to_shortlist.py")
        print(
            "\nDone. See results/phases_g/shortlist_evidence.tsv and "
            "results/phases_g/shortlist_evidence_with_mouse_gse193533_stopgap.tsv"
        )
    except subprocess.CalledProcessError:
        print(
            "\nDone. See results/phases_g/shortlist_evidence.tsv "
            "(mouse stopgap layer was not appended)."
        )
    if args.verify:
        run("verify_phases_g_objectives.py")


if __name__ == "__main__":
    main()
