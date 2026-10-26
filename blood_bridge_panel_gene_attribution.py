"""
Gene-level attribution for the 8-gene blood bridge panel.

Outputs:
  - within-cohort DE stats (log2FC, t, P) for GSE159610 and GSE106520
  - single-gene AUROC per bridge cohort
  - leave-one-out AUROC drop from the full 8-gene panel
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm

from blood_bridge_panel import (
    BASE,
    _auc_from_scores,
    _load_gse159610_expression,
    _load_simple_expression,
    _zscore_rows,
)


OUT_DIR = BASE / "results" / "blood_bridge_panel"
PANEL_GENES = [
    "DYNLT1",
    "OASL",
    "MPP1",
    "TMEM130",
    "MARCO",
    "MCOLN1",
    "RNASET2",
    "ANK3",
]


def _load_panel_signs(base: Path) -> pd.Series:
    panel = pd.read_csv(OUT_DIR / "blood_bridge_panel_gene_summary.tsv", sep="\t")
    signs = np.sign(pd.to_numeric(panel.set_index("Gene")["DESeq2_log2FC"], errors="coerce"))
    return signs.reindex(PANEL_GENES).astype(float)


def _cohort_frames(base: Path) -> dict[str, tuple[pd.DataFrame, pd.DataFrame, str, str]]:
    g159610_expr, g159610_meta = _load_gse159610_expression(base)
    g106520_expr, g106520_meta = _load_simple_expression(
        base / "results" / "external_validation_gse106520" / "gse106520_gene_level_log2fpkm.tsv",
        base / "results" / "external_validation_gse106520" / "gse106520_sample_metadata.tsv",
        sample_col="Sample",
        group_col="Group",
    )
    return {
        "GSE159610": (g159610_expr, g159610_meta, "Aneurysm", "Control"),
        "GSE106520": (g106520_expr, g106520_meta, "Aneurysm", "Control"),
    }


def _prepare_cohort(expr: pd.DataFrame, meta: pd.DataFrame, positive: str, negative: str) -> tuple[pd.DataFrame, pd.Series]:
    available = [c for c in expr.columns if c in set(meta["Sample"])]
    expr = expr.loc[:, available].copy()
    meta = meta[meta["Sample"].isin(available)].copy()
    meta = meta[meta["Group"].isin([positive, negative])].copy()
    meta = meta.drop_duplicates(subset=["Sample"]).set_index("Sample")
    expr = expr.loc[:, expr.columns.isin(meta.index)].copy()
    meta = meta.loc[expr.columns]
    labels = (meta["Group"] == positive).astype(bool)
    return expr, labels


def _gene_score(expr: pd.DataFrame, gene: str, sign: float) -> pd.Series:
    if gene not in expr.index or sign == 0 or np.isnan(sign):
        return pd.Series(dtype=float)
    row = expr.loc[[gene]].apply(pd.to_numeric, errors="coerce")
    z = _zscore_rows(row).iloc[0]
    return z * sign


def _panel_score(expr: pd.DataFrame, genes: list[str], signs: pd.Series) -> pd.Series:
    usable = [g for g in genes if g in expr.index and signs.get(g, 0) not in (0, np.nan)]
    if not usable:
        return pd.Series(dtype=float)
    panel_expr = expr.loc[usable].apply(pd.to_numeric, errors="coerce")
    panel_z = _zscore_rows(panel_expr).dropna(how="all")
    usable = [g for g in usable if g in panel_z.index]
    signed = panel_z.loc[usable].mul(signs.loc[usable], axis=0)
    return signed.mean(axis=0, skipna=True)


def _auc_for_labels(labels: pd.Series, scores: pd.Series) -> float:
    aligned = scores.reindex(labels.index).dropna()
    y = labels.loc[aligned.index].to_numpy(dtype=bool)
    return _auc_from_scores(y, aligned.to_numpy(dtype=float))


def _within_cohort_de(base: Path) -> pd.DataFrame:
    g159 = pd.read_csv(base / "results" / "gse159610_full_de_results_annotated.tsv", sep="\t")
    g159 = g159[g159["Gene_Symbol"].isin(PANEL_GENES)].rename(
        columns={
            "Gene_Symbol": "Gene",
            "Log2FC": "Log2FC_aneurysm_vs_control",
            "T_statistic": "T_statistic",
            "P_value": "P_value",
        }
    )[["Gene", "Log2FC_aneurysm_vs_control", "T_statistic", "P_value"]]
    g159["Cohort"] = "GSE159610"

    g106 = pd.read_csv(
        base / "results" / "external_validation_gse106520" / "gse106520_gene_level_context.tsv",
        sep="\t",
    )
    g106 = g106[g106["Gene"].isin(PANEL_GENES)].rename(
        columns={
            "Log2FC_AvC": "Log2FC_aneurysm_vs_control",
            "P_value": "P_value",
        }
    )[["Gene", "Log2FC_aneurysm_vs_control", "T_statistic", "P_value"]]
    g106["Cohort"] = "GSE106520"

    bridge = pd.read_csv(OUT_DIR / "blood_bridge_panel_gene_summary.tsv", sep="\t")
    bridge = bridge[["Gene", "Bridge_Stouffer_Z", "Bridge_meta_P"]]

    out = pd.concat([g159, g106], ignore_index=True)
    out = out.merge(bridge, on="Gene", how="left")
    return out.sort_values(["Gene", "Cohort"]).reset_index(drop=True)


def _single_gene_auroc(base: Path, signs: pd.Series) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cohort, (expr, meta, pos, neg) in _cohort_frames(base).items():
        expr, labels = _prepare_cohort(expr, meta, pos, neg)
        for gene in PANEL_GENES:
            scores = _gene_score(expr, gene, float(signs[gene]))
            auc = _auc_for_labels(labels, scores)
            pos_scores = scores[labels]
            neg_scores = scores[~labels]
            mw_p = (
                float(mannwhitneyu(pos_scores, neg_scores, alternative="two-sided").pvalue)
                if len(pos_scores) and len(neg_scores)
                else float("nan")
            )
            rows.append(
                {
                    "Gene": gene,
                    "Cohort": cohort,
                    "Single_gene_AUROC": auc,
                    "Mann_Whitney_P": mw_p,
                    "N_positive": int(labels.sum()),
                    "N_negative": int((~labels).sum()),
                }
            )
    return pd.DataFrame(rows)


def _leave_one_out(base: Path, signs: pd.Series) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cohort, (expr, meta, pos, neg) in _cohort_frames(base).items():
        expr, labels = _prepare_cohort(expr, meta, pos, neg)
        full_auc = _auc_for_labels(labels, _panel_score(expr, PANEL_GENES, signs))
        rows.append(
            {
                "Gene": "FULL_PANEL",
                "Cohort": cohort,
                "Panel_AUROC": full_auc,
                "Leave_one_out_AUROC": full_auc,
                "AUROC_drop": 0.0,
                "Pct_drop": 0.0,
            }
        )
        for gene in PANEL_GENES:
            reduced = [g for g in PANEL_GENES if g != gene]
            loo_auc = _auc_for_labels(labels, _panel_score(expr, reduced, signs))
            drop = full_auc - loo_auc
            rows.append(
                {
                    "Gene": gene,
                    "Cohort": cohort,
                    "Panel_AUROC": full_auc,
                    "Leave_one_out_AUROC": loo_auc,
                    "AUROC_drop": drop,
                    "Pct_drop": (drop / full_auc * 100.0) if full_auc else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def _write_summary_md(de: pd.DataFrame, single: pd.DataFrame, loo: pd.DataFrame) -> None:
    mpp1_marco_de = de[de["Gene"].isin(["MPP1", "MARCO"])]
    mpp1_marco_auc = single[single["Gene"].isin(["MPP1", "MARCO"])]
    loo_focus = loo[loo["Gene"].isin(["MPP1", "MARCO", "FULL_PANEL"])]

    lines = [
        "# Blood bridge panel gene attribution",
        "",
        "## MPP1 and MARCO — within-cohort aneurysm vs control",
        "",
        mpp1_marco_de.to_markdown(index=False),
        "",
        "## MPP1 and MARCO — single-gene AUROC",
        "",
        mpp1_marco_auc.to_markdown(index=False),
        "",
        "## MPP1 and MARCO — leave-one-out panel contribution",
        "",
        loo_focus.to_markdown(index=False),
        "",
        "## All eight genes — single-gene AUROC",
        "",
        single.pivot(index="Gene", columns="Cohort", values="Single_gene_AUROC").reindex(PANEL_GENES).to_markdown(),
        "",
        "## All eight genes — leave-one-out AUROC drop",
        "",
        loo[loo["Gene"] != "FULL_PANEL"].pivot(index="Gene", columns="Cohort", values="AUROC_drop").reindex(PANEL_GENES).to_markdown(),
        "",
        "## Interpretation threshold",
        "",
        "- Individual gene AUROC > 0.60 is treated here as a minimal meaningful single-gene discriminator.",
        "- Larger leave-one-out AUROC drop indicates greater contribution to the composite panel score.",
    ]
    (OUT_DIR / "blood_bridge_panel_gene_attribution_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    signs = _load_panel_signs(BASE)

    de = _within_cohort_de(BASE)
    single = _single_gene_auroc(BASE, signs)
    loo = _leave_one_out(BASE, signs)

    de_path = OUT_DIR / "blood_bridge_panel_within_cohort_de.tsv"
    single_path = OUT_DIR / "blood_bridge_panel_single_gene_auroc.tsv"
    loo_path = OUT_DIR / "blood_bridge_panel_leave_one_out_auroc.tsv"

    de.to_csv(de_path, sep="\t", index=False)
    single.to_csv(single_path, sep="\t", index=False)
    loo.to_csv(loo_path, sep="\t", index=False)
    _write_summary_md(de, single, loo)

    summary = {
        "mpp1_marco_within_cohort_de": de[de["Gene"].isin(["MPP1", "MARCO"])].to_dict(orient="records"),
        "mpp1_marco_single_gene_auroc": single[single["Gene"].isin(["MPP1", "MARCO"])].to_dict(orient="records"),
        "mpp1_marco_leave_one_out": loo[loo["Gene"].isin(["MPP1", "MARCO", "FULL_PANEL"])].to_dict(orient="records"),
        "full_panel_auroc": loo[loo["Gene"] == "FULL_PANEL"][["Cohort", "Panel_AUROC"]].to_dict(orient="records"),
    }
    (OUT_DIR / "blood_bridge_panel_gene_attribution_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("Within-cohort DE:\n", de.to_string(index=False))
    print("\nSingle-gene AUROC:\n", single.to_string(index=False))
    print("\nLeave-one-out:\n", loo.to_string(index=False))
    print(f"\nWrote: {de_path.relative_to(BASE)}")
    print(f"Wrote: {single_path.relative_to(BASE)}")
    print(f"Wrote: {loo_path.relative_to(BASE)}")


if __name__ == "__main__":
    main()
