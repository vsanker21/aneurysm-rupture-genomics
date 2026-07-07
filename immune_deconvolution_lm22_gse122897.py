"""
Local immune deconvolution (LM22 + non-negative least squares) for GSE122897 bulk RNA-seq.

Purpose:
  - Address reviewer questions on bulk vs cell-type-specific signal with a standard,
    citable approach (Newman et al., Nat Methods 2015; LM22 signature).
  - Complements (does not replace) CIBERSORTx: same biological intent; use CIBERSORTx
    output when available for manuscript text that names "CIBERSORTx" specifically.

Method:
  - Intersect LM22 genes with mixture genes (HGNC symbols).
  - Transform mixture TPM with log2(TPM+1); LM22 with log2(x+1).
  - Per sample: scipy.optimize.nnls(signature_matrix, mixture_vector); normalize to sum=1.
  - Compare Ruptured vs Unruptured (Welch t-test on proportions; FDR Benjamini-Hochberg).

LM22 source: public mirror used for automation; for publication, verify against the
matrix obtained through https://cibersort.stanford.edu/ per their terms.

Outputs:
  results/cibersortx/lm22_fractions_gse122897.tsv
  results/cibersortx/lm22_ruptured_vs_unruptured.tsv
  results/cibersortx/figures_ultra_hd/lm22_fractions_by_group_ultra_hd.png
  results/cibersortx/figures_ultra_hd/lm22_heatmap_samples_ultra_hd.png
"""
import io
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from scipy.optimize import nnls
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

try:
    from project_paths import get_project_root
    BASE = get_project_root()
except ImportError:
    BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "results/cibersortx"
FIG_DIR = OUT_DIR / "figures_ultra_hd"
FIG_DIR.mkdir(parents=True, exist_ok=True)

LM22_URL = (
    "https://raw.githubusercontent.com/mdozmorov/Immuno_notes/master/data/Cibersoft/LM22.txt"
)
MIXTURE_PATH = OUT_DIR / "GSE122897_CIBERSORTx_mixture_TPM.txt"
META_PATH = BASE / "results/gse122897_sample_metadata_supplementary.tsv"


def load_lm22(path_or_url=LM22_URL):
    if Path(path_or_url).exists():
        txt = Path(path_or_url).read_text(encoding="utf-8", errors="ignore")
    else:
        r = requests.get(path_or_url, timeout=120)
        r.raise_for_status()
        txt = r.text
    df = pd.read_csv(io.StringIO(txt), sep="\t")
    df = df.rename(columns={df.columns[0]: "GeneSymbol"})
    df["GeneSymbol"] = df["GeneSymbol"].astype(str).str.strip()
    df = df.set_index("GeneSymbol")
    return np.log2(df.astype(float) + 1.0), list(df.columns)


def deconvolve_sample(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    coef, _ = nnls(X, y)
    s = coef.sum()
    if s <= 0:
        return np.zeros(X.shape[1])
    return coef / s


def main():
    if not MIXTURE_PATH.exists():
        raise FileNotFoundError(
            f"Missing {MIXTURE_PATH}. Run: python prepare_cibersortx_mixture_gse122897.py"
        )

    meta = pd.read_csv(META_PATH, sep="\t")
    samples = meta["Sample"].astype(str).tolist()
    group_map = dict(zip(meta["Sample"].astype(str), meta["Group"]))

    mix = pd.read_csv(MIXTURE_PATH, sep="\t")
    if mix.columns[0] != "NAME":
        mix = mix.rename(columns={mix.columns[0]: "NAME"})
    mix = mix.set_index("NAME")
    mix = mix[[c for c in samples if c in mix.columns]]

    sig, cell_types = load_lm22(LM22_URL)
    common = sorted(set(sig.index) & set(mix.index))
    print(f"Intersecting genes (LM22 and mixture): {len(common)}")
    if len(common) < 100:
        raise RuntimeError("Too few overlapping genes; check symbol mapping.")

    X = sig.loc[common, cell_types].values.astype(float)
    Y = np.log2(mix.loc[common].values.astype(float) + 1.0)

    fracs = []
    for j, sid in enumerate(mix.columns):
        w = deconvolve_sample(Y[:, j], X)
        row = {"Sample": sid, "Group": group_map.get(sid, "Unknown")}
        for i, ct in enumerate(cell_types):
            row[ct] = w[i]
        fracs.append(row)

    frac_df = pd.DataFrame(fracs)
    frac_df.to_csv(OUT_DIR / "lm22_fractions_gse122897.tsv", sep="\t", index=False)

    # Group tests
    r = frac_df[frac_df["Group"] == "Ruptured"]
    u = frac_df[frac_df["Group"] == "Unruptured"]
    rows = []
    for ct in cell_types:
        a, b = r[ct].values.astype(float), u[ct].values.astype(float)
        tstat, pval = ttest_ind(a, b, equal_var=False, nan_policy="omit")
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
    stat_df["FDR_BH"] = np.nan
    pv_ok = stat_df["P_value"].notna() & np.isfinite(stat_df["P_value"])
    if pv_ok.any():
        _, fdr, _, _ = multipletests(stat_df.loc[pv_ok, "P_value"], method="fdr_bh")
        stat_df.loc[pv_ok, "FDR_BH"] = fdr
    stat_df = stat_df.sort_values("P_value", na_position="last")
    stat_df.to_csv(OUT_DIR / "lm22_ruptured_vs_unruptured.tsv", sep="\t", index=False)
    print(stat_df.head(10).to_string(index=False))

    # --- Figures ---
    plot_df = frac_df.melt(
        id_vars=["Sample", "Group"],
        value_vars=cell_types,
        var_name="Cell_type",
        value_name="Fraction",
    )

    plt.figure(figsize=(14, 8), dpi=200)
    order = (
        stat_df.sort_values("Mean_Ruptured", ascending=False)["Cell_type"].tolist()
    )
    sns.barplot(
        data=plot_df,
        x="Cell_type",
        y="Fraction",
        hue="Group",
        order=order,
        palette={"Ruptured": "#c0392b", "Unruptured": "#2980b9"},
        errorbar="se",
    )
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Estimated immune fraction (LM22 + NNLS)")
    plt.title(
        f"GSE122897 immune deconvolution (n={len(r)} ruptured, n={len(u)} unruptured)\n"
        f"{len(common)} genes intersected with LM22"
    )
    plt.tight_layout()
    p1 = FIG_DIR / "lm22_fractions_by_group_ultra_hd.png"
    plt.savefig(p1, dpi=600, bbox_inches="tight")
    plt.close()
    print(f"Wrote {p1}")

    cell_mat = frac_df.set_index("Sample")[cell_types].T
    # column order: group then sample
    meta_ord = meta.set_index("Sample").loc[[s for s in meta["Sample"] if s in cell_mat.columns]]
    cols = meta_ord.index.tolist()
    cell_mat = cell_mat[cols]
    grp_colors = meta_ord["Group"].map({"Ruptured": "#c0392b", "Unruptured": "#2980b9"})

    g = sns.clustermap(
        cell_mat,
        col_colors=grp_colors,
        cmap="viridis",
        figsize=(16, 8),
        dendrogram_ratio=(0.12, 0.05),
        cbar_kws={"label": "Fraction"},
        xticklabels=True,
        yticklabels=True,
    )
    g.fig.suptitle("LM22 cell-type fractions (samples x cell types)", y=1.02, fontsize=14)
    p2 = FIG_DIR / "lm22_heatmap_samples_ultra_hd.png"
    g.savefig(p2, dpi=600, bbox_inches="tight")
    plt.close()
    print(f"Wrote {p2}")

    # Summary JSON for manuscript snippets (no NaN: invalid JSON)
    sig_n = int((stat_df["FDR_BH"].notna() & (stat_df["FDR_BH"] < 0.05)).sum())
    top_up = (
        stat_df.sort_values("Diff_R_minus_U", ascending=False)
        .head(3)[["Cell_type", "Diff_R_minus_U", "P_value", "FDR_BH"]]
        .replace({np.nan: None})
        .to_dict(orient="records")
    )
    top_dn = (
        stat_df.sort_values("Diff_R_minus_U", ascending=True)
        .head(3)[["Cell_type", "Diff_R_minus_U", "P_value", "FDR_BH"]]
        .replace({np.nan: None})
        .to_dict(orient="records")
    )
    summary = {
        "n_genes_intersect_lm22": len(common),
        "n_ruptured": len(r),
        "n_unruptured": len(u),
        "cell_types": len(cell_types),
        "fdr_significant_cell_types_005": sig_n,
        "top_up_in_ruptured": top_up,
        "top_down_in_ruptured": top_dn,
    }
    import json

    (OUT_DIR / "lm22_deconvolution_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8"
    )
    print("Wrote results/cibersortx/lm22_deconvolution_summary.json")


if __name__ == "__main__":
    main()
