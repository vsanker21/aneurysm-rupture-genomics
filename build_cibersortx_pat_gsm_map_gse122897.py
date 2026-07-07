"""
Build Pat## -> GSM mapping when CIBERSORTx mixture used patient labels.

If you upload a mixture whose column headers are Pat02, Pat126, ... instead of
GSM3487853, ..., the job results use those labels. This script matches each Pat
row to the GSM row with the most similar LM22 fraction profile from
`lm22_fractions_gse122897.tsv` (same biology, different algorithm than CIBERSORTx).

**Validate** the table against your own upload order / clinical key before publication.
If the mixture column order matched GEO GSM order, you can instead build a simple
two-column map by hand.

Usage:
  python build_cibersortx_pat_gsm_map_gse122897.py
  python build_cibersortx_pat_gsm_map_gse122897.py --results results/cibersortx/CIBERSORTx_Job5_Results.txt

Default output: results/cibersortx/cibersortx_pat_to_gsm_map.tsv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

BASE = Path(r"e:\Aneurysms Genomics Project")
DEFAULT_CX = BASE / "results/cibersortx/CIBERSORTx_Job3_Results.txt"
LM = BASE / "results/cibersortx/lm22_fractions_gse122897.tsv"
OUT = BASE / "results/cibersortx/cibersortx_pat_to_gsm_map.tsv"


def normc(x: str) -> str:
    return str(x).strip().lower()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--results",
        type=Path,
        default=DEFAULT_CX,
        help="CIBERSORTx fractions output (tab-separated)",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=OUT,
        help="Output mapping TSV path",
    )
    args = ap.parse_args()
    cx_path = Path(args.results)
    if not cx_path.is_file():
        raise FileNotFoundError(cx_path)

    cx = pd.read_csv(cx_path, sep="\t")
    lm = pd.read_csv(LM, sep="\t")
    qc_like = {"p-value", "p value", "p_value", "pvalue", "correlation", "rmse"}
    ct = [
        c
        for c in cx.columns
        if c != "Mixture"
        and normc(c) not in qc_like
        and c in lm.columns
        and pd.api.types.is_numeric_dtype(cx[c])
    ]

    pat = cx["Mixture"].astype(str).values
    gsm = lm["Sample"].astype(str).values
    Xp = cx[ct].values.astype(float)
    Xg = lm[ct].values.astype(float)

    def row_corr(a, b):
        if np.std(a) < 1e-12 or np.std(b) < 1e-12:
            return 0.0
        return float(np.corrcoef(a, b)[0, 1])

    C = np.zeros((len(pat), len(gsm)))
    for i in range(len(pat)):
        for j in range(len(gsm)):
            C[i, j] = row_corr(Xp[i], Xg[j])
    ri, cj = linear_sum_assignment(1.0 - C)
    out = pd.DataFrame(
        {
            "Mixture": pat[ri],
            "Sample": gsm[cj],
            "profile_correlation": C[ri, cj],
        }
    ).sort_values("Mixture")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"Source CIBERSORTx results: {cx_path.resolve()}")
    print(out.to_string(index=False))
    print(f"\nmean profile_correlation={out['profile_correlation'].mean():.3f}")
    print(f"Wrote {args.out.resolve()}")


if __name__ == "__main__":
    main()
