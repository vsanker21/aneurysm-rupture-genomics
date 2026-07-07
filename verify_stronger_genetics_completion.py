"""
Verify stronger-genetics stack outputs after run_stronger_genetics_stack.py.

Requires:
  - loci_manifest.json with multiple loci carrying gwas_region_file
  - eqtl_fetch_summary.json and eqtl_region_*.tsv files
  - zscreen_correlation.tsv with at least one row

Exit 1 if checks fail. Writes results/phases_g/stronger_genetics/COMPLETION_VERIFICATION.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
SG = BASE / "results/phases_g/stronger_genetics"


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--min-loci-gwas", type=int, default=10)
    ap.add_argument("--min-eqtl-files", type=int, default=5)
    args = ap.parse_args()

    report: dict = {"checks": [], "errors": []}

    def ok(msg: str) -> None:
        report["checks"].append({"ok": True, "detail": msg})

    def err(msg: str) -> None:
        report["errors"].append(msg)

    man = SG / "loci_manifest.json"
    if not man.is_file():
        err("Missing loci_manifest.json — run prepare_stronger_genetics_followup.py")
    else:
        loci = json.loads(man.read_text(encoding="utf-8"))
        n_gw = sum(1 for L in loci if L.get("gwas_region_file"))
        if n_gw < args.min_loci_gwas:
            err(f"Too few loci with gwas_region_file: {n_gw} < {args.min_loci_gwas}")
        else:
            ok(f"loci_manifest: {n_gw} loci with regional GWAS TSV")

    eqtl_sum = SG / "eqtl_fetch_summary.json"
    if not eqtl_sum.is_file():
        err("Missing eqtl_fetch_summary.json — run stronger_genetics_eqtl_fetch.py")
    else:
        sm = json.loads(eqtl_sum.read_text(encoding="utf-8"))
        n_eq = sum(1 for x in sm if int(x.get("n_eqtl_rows", 0)) > 0)
        if len(sm) < args.min_eqtl_files:
            err(f"eqtl_fetch_summary has too few loci: {len(sm)} < {args.min_eqtl_files}")
        elif n_eq < 1:
            err("No locus has n_eqtl_rows > 0 (Catalogue/API or window issue)")
        else:
            ok(f"eqtl_fetch: {len(sm)} loci fetched, {n_eq} with >=1 eQTL row")

    zs = SG / "zscreen_correlation.tsv"
    if not zs.is_file():
        err("Missing zscreen_correlation.tsv — run stronger_genetics_zscreen.py")
    else:
        nlines = sum(1 for _ in open(zs, encoding="utf-8")) - 1
        if nlines < 1:
            err("zscreen_correlation.tsv empty")
        else:
            ok(f"zscreen: {nlines} loci screened (Spearman on z-scores; not formal COLOC)")

    rb = SG / "RUNBOOK_stronger_genetics.txt"
    if rb.is_file():
        ok("RUNBOOK_stronger_genetics.txt present (R coloc / SMR / eQTLGen)")

    report["all_ok"] = len(report["errors"]) == 0
    out = SG / "COMPLETION_VERIFICATION.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["all_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
