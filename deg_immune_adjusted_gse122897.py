"""
GSE122897: differential expression signal adjusted for inferred immune cell fractions.

CIBERSORTx (or local LM22+NNLS) gives per-sample immune proportions. This script fits,
for each gene, a linear model on log2(TPM+1):

    expression ~ rupture_status + selected immune fraction covariates

so you can report which genes stay associated with rupture after accounting for
composition (reviewer-facing bulk confound check).

Important framing for the manuscript:
  - Adjusted p-values support "the expression–rupture association is not fully explained
    by these LM22-inferred immune proportions under a linear additive model."
  - This is NOT proof of cell-intrinsic programs: bulk tissue, blood-derived reference,
    and model assumptions still apply. Pair with pathway/endothelial biomarker logic
    and acknowledge limitations.

Inputs (defaults):
  - TPM: data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
  - DESeq2 results: results/gse122897_deseq2_results.tsv (defines which genes to test)
  - Fractions: results/cibersortx/lm22_fractions_gse122897.tsv
    Use CIBERSORTx merged output when available:
      results/cibersortx/cibersortx_merged_with_group.tsv
    and pass: --fractions path --fraction-cols auto

Outputs:
  - results/cibersortx/deg_immune_covariate_adjusted.tsv
  - results/cibersortx/deg_immune_adjustment_summary.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

try:
    from project_paths import get_project_root
    BASE = get_project_root()
except ImportError:
    BASE = Path(__file__).resolve().parent
TPM_DEFAULT = BASE / "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
META_DEFAULT = BASE / "results/gse122897_sample_metadata_supplementary.tsv"
DESEQ_DEFAULT = BASE / "results/gse122897_deseq2_results.tsv"
FRAC_DEFAULT = BASE / "results/cibersortx/lm22_fractions_gse122897.tsv"
OUT_TSV = BASE / "results/cibersortx/deg_immune_covariate_adjusted.tsv"
OUT_JSON = BASE / "results/cibersortx/deg_immune_adjustment_summary.json"

SKIP_FRAC_COLS = {
    "sample",
    "group",
    "cibersortx_p_value",
    "cibersortx_correlation",
    "cibersortx_rmse",
}


def load_fraction_table(path: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_csv(path, sep="\t")
    df.columns = [str(c).strip() for c in df.columns]
    if "Sample" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Sample"})
    df["Sample"] = df["Sample"].astype(str).str.strip()

    frac_cols = []
    for c in df.columns:
        if c == "Sample" or c == "Group":
            continue
        if c.strip().lower() in SKIP_FRAC_COLS:
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            frac_cols.append(c)
        else:
            test = pd.to_numeric(df[c], errors="coerce")
            if test.notna().sum() >= 0.85 * len(df):
                df[c] = test
                frac_cols.append(c)

    return df, frac_cols


def design_matrix(frac_wide: pd.DataFrame, frac_cols: list[str], meta: pd.DataFrame):
    """Align samples to meta (ruptured/unruptured only); return X, sample list, rupture 0/1."""
    m = meta[meta["Group"].isin(["Ruptured", "Unruptured"])].copy()
    samples = m["Sample"].astype(str).tolist()
    fr = frac_wide[frac_wide["Sample"].isin(samples)].set_index("Sample").loc[samples]

    active = []
    for c in frac_cols:
        if c not in fr.columns:
            continue
        v = pd.to_numeric(fr[c], errors="coerce").astype(float).values
        if np.nanstd(v) > 1e-10:
            active.append(c)

    if len(active) < 1:
        raise RuntimeError("No variable fraction columns after filtering.")

    # Sum-to-one across LM22-like sets: drop last active column to reduce collinearity with intercept.
    F = fr[active].astype(float).values
    F = np.nan_to_num(F, nan=0.0)
    use_cols = active[:-1]
    F_red = F[:, : len(use_cols)]

    ruptured = (m["Group"].values == "Ruptured").astype(float)
    X = np.column_stack([np.ones(len(samples)), ruptured, F_red])
    colnames = ["const", "Ruptured"] + use_cols
    return X, samples, colnames, ruptured


def main():
    ap = argparse.ArgumentParser(description="DE vs rupture adjusted for immune fractions (GSE122897).")
    ap.add_argument("--tpm", type=Path, default=TPM_DEFAULT)
    ap.add_argument("--meta", type=Path, default=META_DEFAULT)
    ap.add_argument("--deseq2", type=Path, default=DESEQ_DEFAULT)
    ap.add_argument("--fractions", type=Path, default=FRAC_DEFAULT)
    ap.add_argument("--max-genes", type=int, default=15000, help="Cap genes tested (sorted by DESeq2 padj).")
    ap.add_argument("--min-basemean", type=float, default=20.0)
    args = ap.parse_args()

    meta = pd.read_csv(args.meta, sep="\t")
    frac_df, frac_cols = load_fraction_table(args.fractions)
    X, samples, xnames, _ = design_matrix(frac_df, frac_cols, meta)

    # Full column rank check
    rank = np.linalg.matrix_rank(X)
    if rank < X.shape[1]:
        raise RuntimeError(
            f"Design matrix rank-deficient (rank={rank}, cols={X.shape[1]}). "
            "Drop redundant fraction columns or samples."
        )

    de = pd.read_csv(args.deseq2, sep="\t")
    de.columns = [c.strip().replace('"', "") for c in de.columns]
    gid = de.columns[0] if "Gene_ID" not in de.columns else "Gene_ID"
    de["Gene_ID"] = de[gid].astype(str).str.replace(".0", "", regex=False)
    de["padj"] = pd.to_numeric(de["padj"], errors="coerce")
    de["baseMean"] = pd.to_numeric(de["baseMean"], errors="coerce")
    de = de[np.isfinite(de["padj"])].sort_values("padj")
    de = de[de["baseMean"] >= args.min_basemean].head(args.max_genes)

    print(f"Loading TPM ({args.tpm.name})...")
    tpm = pd.read_csv(args.tpm, sep="\t", compression="gzip", low_memory=False)
    tpm = tpm.rename(columns={tpm.columns[0]: "GeneID"})
    tpm["GeneID"] = tpm["GeneID"].astype(str).str.replace(".0", "", regex=False)
    tpm = tpm.set_index("GeneID")
    for s in samples:
        if s not in tpm.columns:
            raise ValueError(f"Sample {s} missing from TPM matrix")
    tpm_sub = tpm[samples].apply(pd.to_numeric, errors="coerce").astype(float)

    rows = []
    coef_idx = 1  # Ruptured
    for _, dr in de.iterrows():
        gid = dr["Gene_ID"]
        if gid not in tpm_sub.index:
            continue
        y = np.log2(tpm_sub.loc[gid].values + 1.0)
        if not np.all(np.isfinite(y)):
            continue
        try:
            fit = sm.OLS(y, X).fit()
        except Exception:
            continue
        if fit.df_resid < 1:
            continue
        p_rup = float(fit.pvalues[coef_idx])
        b_rup = float(fit.params[coef_idx])
        rows.append(
            {
                "Gene_ID": gid,
                "DESeq2_log2FC": float(dr.get("log2FoldChange", np.nan)),
                "DESeq2_padj": float(dr["padj"]),
                "baseMean": float(dr["baseMean"]),
                "LM_immune_adj_log2_coef_Ruptured": b_rup,
                "LM_immune_adj_P_Ruptured": p_rup,
            }
        )

    out = pd.DataFrame(rows)
    if len(out) == 0:
        raise RuntimeError("No genes produced results; check ID alignment and inputs.")

    _, out["LM_immune_adj_Padj_BH"], _, _ = multipletests(
        out["LM_immune_adj_P_Ruptured"], method="fdr_bh"
    )

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    sig_de = out["DESeq2_padj"] < 0.05
    sig_adj = out["LM_immune_adj_Padj_BH"] < 0.05
    nom_rup = out["LM_immune_adj_P_Ruptured"] < 0.05
    summary = {
        "n_samples": len(samples),
        "fraction_file": str(args.fractions.resolve()),
        "immune_covariates_in_model": [c for c in xnames if c not in ("const", "Ruptured")],
        "n_genes_tested": int(len(out)),
        "n_DESeq2_FDR_lt_0.05_among_tested": int(sig_de.sum()),
        "n_rupture_coef_nominal_P_lt_0.05": int(nom_rup.sum()),
        "n_rupture_coef_FDR_BH_lt_0.05_on_tested_set": int(sig_adj.sum()),
        "n_DESeq2_sig_and_rupture_coef_nominal_P_lt_0.05": int((sig_de & nom_rup).sum()),
        "n_DESeq2_sig_and_rupture_coef_FDR_BH_lt_0.05": int((sig_de & sig_adj).sum()),
        "n_DESeq2_sig_but_rupture_coef_FDR_ge_0.05": int((sig_de & ~sig_adj).sum()),
        "note": (
            "Linear adjustment on log2(TPM+1); immune fractions from deconvolution reference. "
            "BH-FDR is computed across the genes in this run only (not full transcriptome). "
            "Does not establish cell-intrinsic biology without single-cell or orthogonal validation."
        ),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(out.head(15).to_string(index=False))
    print(f"\nWrote {OUT_TSV}")
    print(f"Wrote {OUT_JSON}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
