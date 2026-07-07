"""
Append stopgap mouse scRNA compartment evidence (GSE193533) to a PHASES-G shortlist TSV.

This is a SUPPORTIVE layer only:
  - species difference (mouse vs human)
  - disease model (elastase-induced mouse aneurysm)
  - coarse EC/SMC/fibroblast/immune annotation

Inputs:
  - results/phases_g/shortlist_evidence*.tsv
  - results/phases_g/gse193533_mouse_celltype_gene_expression.tsv

Outputs:
  - results/phases_g/shortlist_evidence_with_mouse_gse193533_stopgap.tsv (default)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(r"e:\Aneurysms Genomics Project")
DEFAULT_IN = BASE / r"results\phases_g\shortlist_evidence.tsv"
MOUSE_EXPR = BASE / r"results\phases_g\gse193533_mouse_celltype_gene_expression.tsv"
DEFAULT_OUT = BASE / r"results\phases_g\shortlist_evidence_with_mouse_gse193533_stopgap.tsv"


def build_mouse_gene_summary(mouse_expr_tsv: Path) -> pd.DataFrame:
    mex = pd.read_csv(mouse_expr_tsv, sep="\t")
    mex["mean_log1p"] = pd.to_numeric(mex["mean_log1p"], errors="coerce")
    mex["gene_human"] = mex["gene_human"].astype(str)

    # Best compartment (max mean_log1p across coarse cell_type)
    best = (
        mex.dropna(subset=["mean_log1p"])
        .sort_values(["gene_human", "mean_log1p"], ascending=[True, False])
        .drop_duplicates("gene_human")
        .rename(
            columns={
                "gene_human": "symbol",
                "cell_type": "mouse_gse193533_best_cell_type",
                "mean_log1p": "mouse_gse193533_best_mean_log1p",
            }
        )[["symbol", "mouse_gse193533_best_cell_type", "mouse_gse193533_best_mean_log1p"]]
    )

    # Wall minus immune (max wall - immune)
    wall = {"EC", "SMC", "FIBROBLAST"}
    rows = []
    for gene, sub in mex.groupby("gene_human"):
        sub = sub.dropna(subset=["mean_log1p"])
        if sub.empty:
            rows.append((gene, np.nan))
            continue
        immune = (
            float(sub.loc[sub["cell_type"] == "IMMUNE", "mean_log1p"].iloc[0])
            if (sub["cell_type"] == "IMMUNE").any()
            else np.nan
        )
        wall_max = (
            float(sub.loc[sub["cell_type"].isin(wall), "mean_log1p"].max())
            if (sub["cell_type"].isin(wall)).any()
            else np.nan
        )
        rows.append((gene, wall_max - immune))
    wmi = pd.DataFrame(rows, columns=["symbol", "mouse_gse193533_wall_minus_immune"])

    return best.merge(wmi, on="symbol", how="outer")


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--in-tsv", type=Path, default=DEFAULT_IN)
    ap.add_argument("--mouse-expr", type=Path, default=MOUSE_EXPR)
    ap.add_argument("--out-tsv", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    if not args.in_tsv.is_file():
        raise FileNotFoundError(f"Missing input shortlist TSV: {args.in_tsv}")
    if not args.mouse_expr.is_file():
        raise FileNotFoundError(
            f"Missing mouse expression summary TSV: {args.mouse_expr}\n"
            "Run: python analyze_gse193533_celltype_expression.py"
        )

    df = pd.read_csv(args.in_tsv, sep="\t")
    if "symbol" not in df.columns:
        raise ValueError(f"Input TSV missing 'symbol' column: {args.in_tsv}")

    m = build_mouse_gene_summary(args.mouse_expr)
    out = df.merge(m, on="symbol", how="left")

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"Wrote {args.out_tsv}")


if __name__ == "__main__":
    main()

