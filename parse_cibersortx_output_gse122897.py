"""
Parse CIBERSORTx job output and test ruptured vs unruptured fractions (GSE122897).

Usage:
  python parse_cibersortx_output_gse122897.py --input path/to/CIBERSORTx_Results.txt

Expects a table where rows are samples (mixtures) and columns include immune cell
fractions, typically ending with columns like 'P-value', 'Correlation', 'RMSE' —
these are auto-excluded from testing.

Outputs:
  results/cibersortx/cibersortx_merged_with_group.tsv
  results/cibersortx/cibersortx_ruptured_vs_unruptured.tsv
  results/cibersortx/cibersortx_qc_per_sample.tsv
  results/cibersortx/cibersortx_parse_summary.json
  results/cibersortx/figures_ultra_hd/cibersortx_fractions_by_group_ultra_hd.png
"""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

BASE = Path(r"e:\Aneurysms Genomics Project")
DEFAULT_OUT = BASE / "results/cibersortx"
META_PATH = BASE / "results/gse122897_sample_metadata_supplementary.tsv"

SKIP_COLS = {
    "p-value",
    "p value",
    "p_value",
    "pvalue",
    "correlation",
    "pearson correlation",
    "rmse",
    "rms",
    "absolute score (abs)",
    "absolute score",
    "mixture",
}
QC_COL_ALIASES = {
    "p-value": "CIBERSORTx_P_value",
    "p value": "CIBERSORTx_P_value",
    "p_value": "CIBERSORTx_P_value",
    "pvalue": "CIBERSORTx_P_value",
    "correlation": "CIBERSORTx_Correlation",
    "pearson correlation": "CIBERSORTx_Correlation",
    "rmse": "CIBERSORTx_RMSE",
    "rms": "CIBERSORTx_RMSE",
}


def normalize_col(c: str) -> str:
    return str(c).strip().lower()


def is_fraction_column(series: pd.Series, col_name: str) -> bool:
    """True if column looks like numeric cell-type fractions (not sample IDs)."""
    if normalize_col(col_name) in SKIP_COLS:
        return False
    num = pd.to_numeric(series, errors="coerce")
    ok = num.notna().sum()
    if ok < max(2, int(0.85 * len(series))):
        return False
    # Fractions are usually in [0, 1]; allow slight numerical overshoot
    v = num.dropna()
    return bool(v.min() >= -0.05 and v.max() <= 1.05)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="CIBERSORTx output .txt or .csv")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: results/cibersortx)",
    )
    ap.add_argument(
        "--sample-map",
        type=Path,
        default=None,
        help=(
            "Optional TSV: column1 = CIBERSORTx mixture label (e.g. Pat02), "
            "column2 = GEO sample ID (e.g. GSM3487853). "
            "Use when mixture columns were renamed for upload. "
            "See results/cibersortx/cibersortx_pat_to_gsm_map.tsv (validate before citing)."
        ),
    )
    args = ap.parse_args()
    inp = Path(args.input)
    out_dir = Path(args.out_dir) if args.out_dir is not None else DEFAULT_OUT
    fig_dir = out_dir / "figures_ultra_hd"
    if not inp.exists():
        raise FileNotFoundError(inp)

    sep = "\t" if inp.suffix.lower() in (".txt", ".tsv") else ","
    raw = pd.read_csv(inp, sep=sep, low_memory=False)
    # First column often 'Mixture' with sample IDs
    first = raw.columns[0]
    meta = pd.read_csv(META_PATH, sep="\t")
    gmap = dict(zip(meta["Sample"].astype(str), meta["Group"]))

    if normalize_col(first) in ("mixture", "sample", "sample name", "id"):
        raw = raw.rename(columns={first: "Sample"})
    else:
        raw = raw.rename(columns={first: "Sample"})

    raw["Sample"] = raw["Sample"].astype(str).str.strip()
    if args.sample_map:
        smp = pd.read_csv(args.sample_map, sep="\t")
        c0, c1 = smp.columns[0], smp.columns[1]
        mmap = dict(
            zip(
                smp[c0].astype(str).str.strip(),
                smp[c1].astype(str).str.strip(),
            )
        )
        raw["CIBERSORTx_mixture"] = raw["Sample"]
        raw["Sample"] = raw["Sample"].map(mmap)
        n_miss = int(raw["Sample"].isna().sum())
        if n_miss:
            print(
                f"[warn] {n_miss} mixture rows have no entry in --sample-map "
                "(will drop when merging groups)."
            )
    raw["Group"] = raw["Sample"].map(gmap)
    raw = raw.dropna(subset=["Group"])
    raw = raw[raw["Group"].isin(["Ruptured", "Unruptured"])]

    cell_cols = []
    for c in raw.columns:
        if c in ("Sample", "Group"):
            continue
        if normalize_col(c) in SKIP_COLS:
            continue
        if raw[c].dtype in (float, int) or pd.api.types.is_numeric_dtype(raw[c]):
            cell_cols.append(c)
        elif is_fraction_column(raw[c], c):
            cell_cols.append(c)

    r = raw[raw["Group"] == "Ruptured"]
    u = raw[raw["Group"] == "Unruptured"]
    rows = []
    for ct in cell_cols:
        a = pd.to_numeric(r[ct], errors="coerce").dropna().values
        b = pd.to_numeric(u[ct], errors="coerce").dropna().values
        if len(a) < 2 or len(b) < 2:
            continue
        tstat, pval = ttest_ind(a, b, equal_var=False)
        if not (np.isfinite(tstat) and np.isfinite(pval)):
            tstat, pval = np.nan, np.nan
        rows.append(
            {
                "Cell_type": ct,
                "Mean_Ruptured": float(np.nanmean(a)),
                "Mean_Unruptured": float(np.nanmean(b)),
                "Diff_R_minus_U": float(np.nanmean(a) - np.nanmean(b)),
                "T_statistic": float(tstat) if np.isfinite(tstat) else np.nan,
                "P_value": float(pval) if np.isfinite(pval) else np.nan,
            }
        )

    stat_df = pd.DataFrame(rows)
    if len(stat_df):
        stat_df["FDR_BH"] = np.nan
        pv_ok = stat_df["P_value"].notna() & np.isfinite(stat_df["P_value"])
        if pv_ok.any():
            _, fdr, _, _ = multipletests(stat_df.loc[pv_ok, "P_value"], method="fdr_bh")
            stat_df.loc[pv_ok, "FDR_BH"] = fdr
        stat_df = stat_df.sort_values("P_value", na_position="last")

    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    qc_parts = {"Sample": raw["Sample"], "Group": raw["Group"]}
    for c in raw.columns:
        nc = normalize_col(c)
        if nc in QC_COL_ALIASES:
            qc_parts[QC_COL_ALIASES[nc]] = pd.to_numeric(raw[c], errors="coerce")
    qc_df = pd.DataFrame(qc_parts)
    qc_df.to_csv(out_dir / "cibersortx_qc_per_sample.tsv", sep="\t", index=False)

    raw.to_csv(out_dir / "cibersortx_merged_with_group.tsv", sep="\t", index=False)
    stat_df.to_csv(out_dir / "cibersortx_ruptured_vs_unruptured.tsv", sep="\t", index=False)

    if len(cell_cols) and len(raw):
        plot_df = raw.melt(
            id_vars=["Sample", "Group"],
            value_vars=cell_cols,
            var_name="Cell_type",
            value_name="Fraction",
        )
        plot_df["Fraction"] = pd.to_numeric(plot_df["Fraction"], errors="coerce")
        plot_df = plot_df.dropna(subset=["Fraction"])
        plt.figure(figsize=(14, 8), dpi=200)
        sns.barplot(
            data=plot_df,
            x="Cell_type",
            y="Fraction",
            hue="Group",
            palette={"Ruptured": "#c0392b", "Unruptured": "#2980b9"},
            errorbar="se",
        )
        plt.xticks(rotation=45, ha="right")
        plt.title("CIBERSORTx fractions: ruptured vs unruptured (GSE122897)")
        plt.tight_layout()
        pfig = fig_dir / "cibersortx_fractions_by_group_ultra_hd.png"
        plt.savefig(pfig, dpi=600, bbox_inches="tight")
        plt.close()
        print(f"Wrote {pfig}")

    qc_summary = {}
    for std_name in (
        "CIBERSORTx_P_value",
        "CIBERSORTx_Correlation",
        "CIBERSORTx_RMSE",
    ):
        if std_name not in qc_df.columns:
            continue
        s = qc_df[std_name].dropna()
        if len(s) == 0:
            continue
        qc_summary[std_name] = {
            "n": int(len(s)),
            "median": float(s.median()),
            "min": float(s.min()),
            "max": float(s.max()),
        }
        by_g = qc_df.groupby("Group")[std_name].median()
        qc_summary[std_name]["median_by_group"] = {
            str(k): float(v) for k, v in by_g.items() if pd.notna(v)
        }

    summary_payload = {
        "input_file": str(inp.resolve()),
        "n_samples_parsed": int(len(raw)),
        "cell_type_columns": cell_cols,
        "qc_summary": qc_summary,
    }
    if args.sample_map:
        summary_payload["sample_map_file"] = str(Path(args.sample_map).resolve())

    (out_dir / "cibersortx_parse_summary.json").write_text(
        json.dumps(summary_payload, indent=2, allow_nan=False),
        encoding="utf-8",
    )
    print(f"Merged: {out_dir / 'cibersortx_merged_with_group.tsv'}")
    print(f"Stats:  {out_dir / 'cibersortx_ruptured_vs_unruptured.tsv'}")
    print(f"QC:     {out_dir / 'cibersortx_qc_per_sample.tsv'}")


if __name__ == "__main__":
    main()
