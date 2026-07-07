#!/usr/bin/env python3
"""Export a regional Bakker GWAS TSV to GCTA-SMR .ma format (outcome / GWAS side).

SMR expects (tab-separated): SNP A1 A2 freq b se p N
We map: A1=A_EFF, A2=A_NONEFF, freq=Freq_EFF, b=BETA, se=SE, p=P, N=Neff.

Example:
  python export_smr_gwas_ma.py --gwas-tsv results/phases_g/stronger_genetics/gwas_region_OXER1_rs2043912.tsv \\
    --out results/phases_g/stronger_genetics/smr/OXER1_rs2043912_outcome.ma
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser(description="Regional GWAS TSV -> SMR .ma")
    ap.add_argument("--gwas-tsv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    req = ["SNP", "A_EFF", "A_NONEFF", "Freq_EFF", "BETA", "SE", "P", "Neff"]
    df = pd.read_csv(args.gwas_tsv, sep="\t")
    miss = [c for c in req if c not in df.columns]
    if miss:
        raise SystemExit(f"Missing columns in GWAS TSV: {miss}")

    out = pd.DataFrame(
        {
            "SNP": df["SNP"].astype(str),
            "A1": df["A_EFF"].astype(str).str.upper(),
            "A2": df["A_NONEFF"].astype(str).str.upper(),
            "freq": pd.to_numeric(df["Freq_EFF"], errors="coerce"),
            "b": pd.to_numeric(df["BETA"], errors="coerce"),
            "se": pd.to_numeric(df["SE"], errors="coerce"),
            "p": pd.to_numeric(df["P"], errors="coerce"),
            "N": pd.to_numeric(df["Neff"], errors="coerce"),
        }
    )
    out = out.dropna(subset=["SNP", "A1", "A2", "b", "se", "p", "N"], how="any")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False, na_rep="NA")
    print(f"Wrote {args.out} ({len(out)} SNPs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
