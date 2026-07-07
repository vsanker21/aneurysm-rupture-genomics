"""
Summarize stopgap cell-compartment evidence for genes lacking human scRNA specificity assignment.

Inputs:
  - results/phases_g/tier1_supplement_with_scrna_specificity_cvae186_v2.tsv
  - results/phases_g/gse193533_mouse_celltype_gene_expression.tsv

Outputs:
  - results/phases_g/nan_scrna_stopgap_mouse_summary.tsv
  - results/phases_g/nan_scrna_stopgap_mouse_report.md

Notes:
  This uses mouse aneurysm-model scRNA-seq (GSE193533) and a very coarse annotation.
  Treat this as supportive evidence only; not definitive for human IA wall cell-of-origin.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(r"e:\Aneurysms Genomics Project")
TIER1 = BASE / r"results\phases_g\tier1_supplement_with_scrna_specificity_cvae186_v2.tsv"
MOUSE_EXPR = BASE / r"results\phases_g\gse193533_mouse_celltype_gene_expression.tsv"
OUT_TSV = BASE / r"results\phases_g\nan_scrna_stopgap_mouse_summary.tsv"
OUT_MD = BASE / r"results\phases_g\nan_scrna_stopgap_mouse_report.md"


def main() -> None:
    df = pd.read_csv(TIER1, sep="\t")
    mex = pd.read_csv(MOUSE_EXPR, sep="\t")
    mex["mean_log1p"] = pd.to_numeric(mex["mean_log1p"], errors="coerce")

    # Candidate genes: missing human scRNA assignment in tier1, but present in tier1 table.
    missing = df[df.get("scrna_best_cell_type").isna() | (df.get("scrna_best_cell_type").astype(str).str.strip() == "")]
    genes = missing["symbol"].astype(str).tolist()

    if not genes:
        OUT_MD.write_text("No missing `scrna_best_cell_type` genes found in Tier-1 table.\n", encoding="utf-8")
        return

    wall = {"EC", "SMC", "FIBROBLAST"}
    rows = []
    for g in genes:
        sub = mex[mex["gene_human"].astype(str) == g].dropna(subset=["mean_log1p"])
        if sub.empty:
            rows.append(
                {
                    "symbol": g,
                    "mouse_best_cell_type": np.nan,
                    "mouse_best_mean_log1p": np.nan,
                    "mouse_wall_minus_immune": np.nan,
                }
            )
            continue
        s = sub.sort_values("mean_log1p", ascending=False).iloc[0]
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
        rows.append(
            {
                "symbol": g,
                "mouse_best_cell_type": str(s["cell_type"]),
                "mouse_best_mean_log1p": float(s["mean_log1p"]),
                "mouse_wall_minus_immune": float(wall_max - immune),
            }
        )

    out = pd.DataFrame(rows).sort_values("mouse_wall_minus_immune", ascending=False)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    md = []
    md.append("## Stopgap cell-compartment evidence (mouse scRNA; GSE193533)\n")
    md.append(
        "These genes lacked **human IA wall scRNA marker-table specificity** in the Tier-1 table. "
        "As a stopgap, we summarize which **broad compartment** expresses the gene in a **mouse aneurysm model** scRNA dataset (GSE193533). "
        "Interpret cautiously (species + model + coarse annotation).\n"
    )
    md.append("\n### Summary table\n")
    md.append(out.to_markdown(index=False))
    md.append("\n\n### Interpretation guidance\n")
    md.append(
        "- **mouse_wall_minus_immune**: (max(mean_log1p in EC/SMC/fibroblast) - mean_log1p in immune). "
        "Positive suggests wall-biased expression in this mouse dataset.\n"
    )
    md.append("- **mouse_best_cell_type**: coarse call from canonical markers; not a definitive cell-of-origin assignment.\n")
    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")

    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()

