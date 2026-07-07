"""
Re-extract one or more regional GWAS TSVs from Bakker sumstats using manifest coordinates
plus an optional symmetric expansion. Use after widening the eQTL fetch window so COLOC
has GWAS rows in the same genomic span.

Example:
  python reslice_gwas_region.py --locus-ids FN3K_rs1046896 --expand-window-kb 1000
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent
DEFAULT_MANIFEST = BASE / "results/phases_g/stronger_genetics/loci_manifest.json"
DEFAULT_GWAS = BASE / "data/gwas/bakker_2020/downloaded/IA.GWAS.BakkerMK.2020.sumstats.Stage_2.txt.gz"
OUT_DIR = BASE / "results/phases_g/stronger_genetics"


def norm_chr_label(v) -> str:
    s = str(v).strip().upper().replace("CHR", "")
    if s in ("23", "X"):
        return "X"
    try:
        return str(int(float(s)))
    except ValueError:
        return s


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--gwas", type=Path, default=DEFAULT_GWAS)
    ap.add_argument(
        "--locus-ids",
        type=str,
        required=True,
        help="Comma-separated locus_id values from the manifest.",
    )
    ap.add_argument(
        "--expand-window-kb",
        type=int,
        default=0,
        help="Extend region_start/end by this many kb on each side (must match eQTL fetch).",
    )
    ap.add_argument("--max-snps-per-locus", type=int, default=50_000)
    args = ap.parse_args()

    if not args.manifest.is_file():
        print(f"Missing manifest {args.manifest}", file=__import__("sys").stderr)
        return 1
    if not args.gwas.is_file():
        print(f"Missing GWAS {args.gwas}", file=__import__("sys").stderr)
        return 1

    want = {x.strip() for x in args.locus_ids.split(",") if x.strip()}
    loci = json.loads(args.manifest.read_text(encoding="utf-8"))
    intervals: dict[str, tuple[str, int, int]] = {}
    for L in loci:
        lid = L.get("locus_id")
        if lid not in want:
            continue
        chrom = norm_chr_label(L["sentinel_chr"])
        lo = max(1, int(L["region_start"]) - args.expand_window_kb * 1000)
        hi = int(L["region_end"]) + args.expand_window_kb * 1000
        intervals[lid] = (chrom, lo, hi)
        print(f"[reslice] {lid} chr{chrom}:{lo}-{hi}", flush=True)

    missing = want - set(intervals)
    if missing:
        print(f"Unknown locus_id(s) in manifest: {sorted(missing)}", file=__import__("sys").stderr)
        return 1

    hits: dict[str, list] = defaultdict(list)
    for chunk in pd.read_csv(
        args.gwas, sep="\t", compression="gzip", chunksize=400_000, low_memory=False
    ):
        chunk["_chrn"] = chunk["CHR"].map(norm_chr_label)
        for lid, (cwant, lo, hi) in intervals.items():
            sub = chunk[chunk["_chrn"] == cwant]
            if sub.empty:
                continue
            bp = sub["BP"].astype(int)
            m = (bp >= lo) & (bp <= hi)
            if m.any():
                hits[lid].append(sub.loc[m].drop(columns=["_chrn"], errors="ignore"))

    for lid in sorted(intervals):
        parts = hits.get(lid) or []
        if not parts:
            print(f"[reslice] WARNING: no SNPs for {lid}", flush=True)
            continue
        reg = pd.concat(parts, ignore_index=True).drop_duplicates(subset=["SNP"])
        if len(reg) > args.max_snps_per_locus:
            reg = reg.nsmallest(args.max_snps_per_locus, "P")
        fn = OUT_DIR / f"gwas_region_{lid}.tsv"
        reg.to_csv(fn, sep="\t", index=False)
        print(f"[reslice] Wrote {fn} ({len(reg)} SNPs)", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
