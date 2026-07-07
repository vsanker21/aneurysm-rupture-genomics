"""
Compare two analysis methodologies on the SAME local expression matrices:

A) Pipeline method (blood_bridge_panel_gene_attribution.py):
   - log2(CPM+1) or preprocessed log2fpkm
   - Welch t-test from saved DE tables
   - Single-gene AUROC: cohort z-score x discovery (GSE122897) sign, no flip
   - LOO: discovery-signed row z-scores, mean panel

B) User GEOparse-script method (applied locally):
   - Same expression matrices
   - Welch t-test on sample means (recomputed)
   - Single-gene AUROC: max(AUC, 1-AUC) direction-agnostic
   - LOO: cohort sign flip if log2FC<0, sklearn StandardScaler, mean
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from blood_bridge_panel import (
    BASE,
    _auc_from_scores,
    _load_gse159610_expression,
    _load_simple_expression,
)
from blood_bridge_panel_gene_attribution import (
    PANEL_GENES,
    _gene_score,
    _load_panel_signs,
    _panel_score,
    _prepare_cohort,
)


OUT_DIR = BASE / "results" / "blood_bridge_panel"
COHORTS = {
    "GSE159610": ("Aneurysm", "Control"),
    "GSE106520": ("Aneurysm", "Control"),
}


def _load_cohort_expr(base: Path, cohort: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if cohort == "GSE159610":
        return _load_gse159610_expression(base)
    return _load_simple_expression(
        base / "results" / "external_validation_gse106520" / "gse106520_gene_level_log2fpkm.tsv",
        base / "results" / "external_validation_gse106520" / "gse106520_sample_metadata.tsv",
        sample_col="Sample",
        group_col="Group",
    )


def _auroc_flip(y: np.ndarray, scores: np.ndarray) -> float:
    if len(np.unique(y)) < 2:
        return float("nan")
    auc = roc_auc_score(y, scores)
    return float(max(auc, 1 - auc))


def _method_user(expr: pd.DataFrame, labels: pd.Series, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    de_rows = []
    auc_rows = []
    y = labels.to_numpy(dtype=bool)

    gene_vals = {}
    for gene in genes:
        if gene not in expr.index:
            continue
        v = pd.to_numeric(expr.loc[gene], errors="coerce")
        gene_vals[gene] = v
        pos = v[labels].values
        neg = v[~labels].values
        log2fc = float(np.nanmean(pos) - np.nanmean(neg))
        t_stat, p_val = stats.ttest_ind(pos, neg, equal_var=False, nan_policy="omit")
        de_rows.append({"Gene": gene, "Log2FC": log2fc, "T_statistic": t_stat, "P_value": p_val})
        auc_rows.append({"Gene": gene, "Single_gene_AUROC": _auroc_flip(y.astype(int), v.values)})

    de_df = pd.DataFrame(de_rows)
    auc_df = pd.DataFrame(auc_rows)

    panel = pd.DataFrame(gene_vals).dropna(axis=1)
    available = list(panel.columns)
    for gene in available:
        if de_df.loc[de_df["Gene"] == gene, "Log2FC"].iloc[0] < 0:
            panel[gene] = -panel[gene]
    scaler = StandardScaler()
    scaled = pd.DataFrame(scaler.fit_transform(panel), columns=available, index=panel.index)
    full_auc = _auroc_flip(y.astype(int), scaled.mean(axis=1).values)

    loo_rows = []
    for gene in available:
        rem = [g for g in available if g != gene]
        loo_auc = _auroc_flip(y.astype(int), scaled[rem].mean(axis=1).values)
        loo_rows.append({"Gene": gene, "LOO_AUROC": loo_auc, "AUROC_drop": full_auc - loo_auc})
    loo_df = pd.DataFrame(loo_rows)
    perf = auc_df.merge(loo_df, on="Gene", how="outer")
    return de_df, perf, full_auc


def _method_pipeline(expr: pd.DataFrame, labels: pd.Series, genes: list[str], signs: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    auc_rows = []
    y = labels

    for gene in genes:
        if gene not in expr.index:
            continue
        scores = _gene_score(expr, gene, float(signs[gene]))
        auc = _auc_from_scores(y.to_numpy(dtype=bool), scores.to_numpy(dtype=float))
        auc_rows.append({"Gene": gene, "Single_gene_AUROC": auc})

    full_scores = _panel_score(expr, genes, signs)
    full_auc = _auc_from_scores(y.to_numpy(dtype=bool), full_scores.to_numpy(dtype=float))

    loo_rows = []
    for gene in genes:
        if gene not in expr.index:
            continue
        rem = [g for g in genes if g != gene]
        loo_auc = _auc_from_scores(
            y.to_numpy(dtype=bool),
            _panel_score(expr, rem, signs).to_numpy(dtype=float),
        )
        loo_rows.append({"Gene": gene, "LOO_AUROC": loo_auc, "AUROC_drop": full_auc - loo_auc})

    return pd.DataFrame(), pd.DataFrame(auc_rows).merge(pd.DataFrame(loo_rows), on="Gene"), full_auc


def main() -> None:
    signs = _load_panel_signs(BASE)
    genes = [g for g in PANEL_GENES if g in signs.index]

    rows = []
    for cohort, (pos, neg) in COHORTS.items():
        expr, meta = _load_cohort_expr(BASE, cohort)
        expr, labels = _prepare_cohort(expr, meta, pos, neg)

        de_u, perf_u, full_u = _method_user(expr, labels, genes)
        _, perf_p, full_p = _method_pipeline(expr, labels, genes, signs)

        de_saved = pd.read_csv(OUT_DIR / "blood_bridge_panel_within_cohort_de.tsv", sep="\t")
        de_p = de_saved[de_saved["Cohort"] == cohort][
            ["Gene", "Log2FC_aneurysm_vs_control", "T_statistic", "P_value"]
        ].rename(columns={"Log2FC_aneurysm_vs_control": "Log2FC"})

        de_u = de_u.assign(Cohort=cohort)
        de_p = de_p.assign(Cohort=cohort)
        comp = de_u.merge(de_p, on=["Cohort", "Gene"], suffixes=("_user", "_pipeline"))
        rows.append(comp)

        perf_u = perf_u.assign(Cohort=cohort, Method="User_script_logic", Full_panel_AUROC=full_u)
        perf_p = perf_p.assign(Cohort=cohort, Method="Pipeline", Full_panel_AUROC=full_p)

        print(f"\n=== {cohort} ===")
        print(f"Full panel AUROC  User-script logic: {full_u:.3f}  Pipeline: {full_p:.3f}")
        print("\nMPP1 / MARCO — User script logic:")
        print(perf_u[perf_u["Gene"].isin(["MPP1", "MARCO"])][["Gene", "Single_gene_AUROC", "AUROC_drop"]].to_string(index=False))
        print("\nMPP1 / MARCO — Pipeline:")
        print(perf_p[perf_p["Gene"].isin(["MPP1", "MARCO"])][["Gene", "Single_gene_AUROC", "AUROC_drop"]].to_string(index=False))

    comparison = pd.concat(rows, ignore_index=True)
    out = OUT_DIR / "blood_bridge_panel_method_comparison.tsv"
    comparison.to_csv(out, sep="\t", index=False)
    print(f"\nWrote: {out}")


if __name__ == "__main__":
    main()
