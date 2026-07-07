"""
Z-score correlation screen between regional IA GWAS (Bakker) and regional blood eQTL
(Catalogue export) at overlapping rsIDs.

This is a **fast, non-colocalization** screen (not coloc.abf PP4). Use R `coloc`
for formal shared-variant posterior inference.

Output: results/phases_g/stronger_genetics/zscreen_correlation.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(r"e:\Aneurysms Genomics Project")
DEFAULT_MANIFEST = BASE / "results/phases_g/stronger_genetics/loci_manifest.json"
OUT_DIR = BASE / "results/phases_g/stronger_genetics"
OUT_TSV = OUT_DIR / "zscreen_correlation.tsv"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--min-snps", type=int, default=15)
    args = ap.parse_args()

    if not args.manifest.is_file():
        print("Missing loci_manifest.json", file=sys.stderr)
        return 1

    import json

    loci = json.loads(args.manifest.read_text(encoding="utf-8"))
    rows_out = []

    for L in loci:
        lid = L.get("locus_id")
        gwas_p = L.get("gwas_region_file")
        if not lid or not gwas_p:
            continue
        gwas_p = Path(gwas_p)
        eqtl_p = OUT_DIR / f"eqtl_region_{lid}.tsv"
        if not gwas_p.is_file() or not eqtl_p.is_file():
            continue
        if gwas_p.stat().st_size == 0 or eqtl_p.stat().st_size == 0:
            rows_out.append(
                {
                    "locus_id": lid,
                    "symbol": L.get("symbol"),
                    "n_overlap_snps": 0,
                    "spearman_r": np.nan,
                    "spearman_p": np.nan,
                    "note": "empty_gwas_or_eqtl_tsv",
                }
            )
            continue
        try:
            gwas = pd.read_csv(gwas_p, sep="\t")
            eqtl = pd.read_csv(eqtl_p, sep="\t")
        except pd.errors.EmptyDataError:
            rows_out.append(
                {
                    "locus_id": lid,
                    "symbol": L.get("symbol"),
                    "n_overlap_snps": 0,
                    "spearman_r": np.nan,
                    "spearman_p": np.nan,
                    "note": "empty_parse",
                }
            )
            continue
        if "SNP" not in gwas.columns or "BETA" not in gwas.columns or "SE" not in gwas.columns:
            continue
        if "rsid" not in eqtl.columns or "beta" not in eqtl.columns or "se" not in eqtl.columns:
            continue
        gwas = gwas.copy()
        gwas["SNP"] = gwas["SNP"].astype(str)
        eqtl = eqtl.copy()
        eqtl["rsid"] = eqtl["rsid"].astype(str)
        gwas["z_ia"] = pd.to_numeric(gwas["BETA"], errors="coerce") / pd.to_numeric(
            gwas["SE"], errors="coerce"
        )
        eqtl["z_eq"] = pd.to_numeric(eqtl["beta"], errors="coerce") / pd.to_numeric(
            eqtl["se"], errors="coerce"
        )
        m = gwas.merge(
            eqtl,
            left_on="SNP",
            right_on="rsid",
            how="inner",
            suffixes=("_drop", "_y"),
        )
        m = m[np.isfinite(m["z_ia"]) & np.isfinite(m["z_eq"])]
        n = len(m)
        if n < args.min_snps:
            rows_out.append(
                {
                    "locus_id": lid,
                    "symbol": L.get("symbol"),
                    "n_overlap_snps": n,
                    "spearman_r": np.nan,
                    "spearman_p": np.nan,
                    "note": f"n<{args.min_snps}",
                }
            )
            continue
        r_s, p_s = stats.spearmanr(m["z_ia"].values, m["z_eq"].values)
        rows_out.append(
            {
                "locus_id": lid,
                "symbol": L.get("symbol"),
                "n_overlap_snps": n,
                "spearman_r": float(r_s),
                "spearman_p": float(p_s),
                "note": "zscore_screen_not_coloc",
            }
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_df = pd.DataFrame(rows_out)
    out_df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"Wrote {OUT_TSV} ({len(out_df)} loci screened)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
