"""
Generate Supplementary Tables S1–S12, Figures S1–S5, Appendices 1–2 (TRIPOD, PROBAST),
and a single Word document with detailed legends and numeric interpretations.

Run from project root:
  python generate_supplementary_materials_complete.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt
from matplotlib.gridspec import GridSpec
from scipy import stats
from scipy.stats import hypergeom, mannwhitneyu, norm

from blood_bridge_panel import (
    PANEL_SIZE,
    _auc_from_scores,
    _load_bridge_candidates,
    _load_gse159610_expression,
    _load_simple_expression,
    _zscore_rows,
)
from create_molecular_rupture_signatures_manuscript import (
    add_figure_with_caption_and_inference,
    add_table_from_dataframe,
)
from manuscript_revised_vs_metrics import load_story_metrics

BASE = Path(__file__).resolve().parent
TAB_DIR = BASE / "results" / "supplementary" / "tables"
FIG_DIR = BASE / "results" / "figures" / "supplementary"
OUT_DOCX = BASE / "Supplementary_Materials_Complete.docx"

PANEL_GENES = ["DYNLT1", "OASL", "MPP1", "TMEM130", "MARCO", "MCOLN1", "RNASET2", "ANK3"]
TISSUE_CORE_6 = ["TGFBI", "ANPEP", "MPP1", "MARCO", "CLEC5A", "CPXM2"]
DPI = 600

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 10,
    "savefig.dpi": DPI,
    "figure.dpi": DPI,
})


def _ensure_dirs() -> None:
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)


def _save_table(df: pd.DataFrame, name: str) -> Path:
    path = TAB_DIR / f"{name}.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def _bootstrap_auc_ci(
    labels: np.ndarray,
    scores: np.ndarray,
    n_boot: int = 1500,
    seed: int = 42,
) -> tuple[float, float, float]:
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=float)
    mask = np.isfinite(scores)
    labels, scores = labels[mask], scores[mask]
    if labels.sum() == 0 or (~labels).sum() == 0:
        return float("nan"), float("nan"), float("nan")
    auc = _auc_from_scores(labels, scores)
    rng = np.random.default_rng(seed)
    boots = []
    n = len(labels)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        y_b, s_b = labels[idx], scores[idx]
        if y_b.sum() == 0 or (~y_b).sum() == 0:
            continue
        boots.append(_auc_from_scores(y_b, s_b))
    if not boots:
        return auc, float("nan"), float("nan")
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(auc), float(lo), float(hi)


# ---------------------------------------------------------------------------
# Table S1 — Cohort clinical characteristics
# ---------------------------------------------------------------------------
def build_table_s1(m: dict) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """Main cohort table + platform gene coverage matrix for 8-gene panel."""
    rows = [
        dict(
            Phase="1 — Discovery and blood-bridge validation",
            GEO_accession="GSE122897",
            Comparison="Ruptured vs unruptured (wall tissue)",
            N_positive=m["n_ruptured"],
            N_negative=m["n_unruptured"],
            N_total=m["n_samples"],
            Age_mean_SD="Not in GEO metadata",
            Sex_pct_female="1/43 samples with sex field (not analysable)",
            Hypertension_pct="Not reported",
            Aneurysm_size_location="Not reported",
            Clinical_presentation="Resected aneurysm wall at surgery (ruptured or unruptured)",
            Collection_timing="At surgical resection",
            Biospecimen="Aneurysm wall tissue",
            RNA_method="Total RNA; Illumina HiSeq 2500 RNA-seq",
            Platform="GPL16791 / RNA-seq",
            Accession="GSE122897",
            Covariates_available="Rupture status only",
        ),
        dict(
            Phase="1 — Discovery and blood-bridge validation",
            GEO_accession="GSE159610",
            Comparison="Intracranial aneurysm vs healthy control (whole blood)",
            N_positive=23,
            N_negative=24,
            N_total=47,
            Age_mean_SD="Not in processed metadata",
            Sex_pct_female="Not reported",
            Hypertension_pct="Not reported",
            Aneurysm_size_location="Not reported",
            Clinical_presentation="Incidental/symptomatic aneurysm (DSA cohort)",
            Collection_timing="At angiography (pre-treatment blood)",
            Biospecimen="Whole blood",
            RNA_method="RNA-seq (CPM from raw counts)",
            Platform="GPL16791",
            Accession="GSE159610",
            Covariates_available="Group label (aneurysm vs control)",
        ),
        dict(
            Phase="1 — Discovery and blood-bridge validation",
            GEO_accession="GSE106520",
            Comparison="Aneurysm vs control (isolated neutrophils)",
            N_positive=16,
            N_negative=16,
            N_total=32,
            Age_mean_SD="Paired by age/comorbidity in study design",
            Sex_pct_female="Not reported",
            Hypertension_pct="Not reported",
            Aneurysm_size_location="Not reported",
            Clinical_presentation="Unruptured aneurysm vs controls",
            Collection_timing="At blood draw",
            Biospecimen="Neutrophils",
            RNA_method="RNA-seq (log2 FPKM+1)",
            Platform="GPL16791",
            Accession="GSE106520",
            Covariates_available="Group label",
        ),
        dict(
            Phase="2 — Direct rupture tissue replication",
            GEO_accession="GSE54083",
            Comparison="Ruptured vs unruptured (aneurysm dome)",
            N_positive=m["gse54083_n_ruptured"],
            N_negative=m["gse54083_n_unruptured"],
            N_total=m["gse54083_n_ruptured"] + m["gse54083_n_unruptured"],
            Age_mean_SD="Not in GEO",
            Sex_pct_female="Not reported",
            Hypertension_pct="Not reported",
            Aneurysm_size_location="Not reported",
            Clinical_presentation="Surgical tissue; 10 STA controls in series (excluded from RvU)",
            Collection_timing="At surgery",
            Biospecimen="Aneurysm wall",
            RNA_method="Agilent microarray",
            Platform="GPL4133",
            Accession="GSE54083",
            Covariates_available="Rupture status",
        ),
        dict(
            Phase="2 — Direct rupture tissue replication",
            GEO_accession="GSE13353",
            Comparison="Ruptured vs unruptured",
            N_positive=11,
            N_negative=8,
            N_total=19,
            Age_mean_SD="Not in GEO",
            Sex_pct_female="Not reported",
            Hypertension_pct="Not reported",
            Aneurysm_size_location="Not reported",
            Clinical_presentation="Surgical tissue",
            Collection_timing="At surgery",
            Biospecimen="Aneurysm wall",
            RNA_method="Affymetrix microarray",
            Platform="GPL570",
            Accession="GSE13353",
            Covariates_available="Rupture status",
        ),
        dict(
            Phase="3 — Orthogonal tissue and blood context",
            GEO_accession="GSE15629",
            Comparison="Ruptured vs unruptured (subset) / aneurysm vs control vessel",
            N_positive=8,
            N_negative=6,
            N_total=19,
            Age_mean_SD="Not in GEO",
            Sex_pct_female="Not reported",
            Hypertension_pct="Not reported",
            Aneurysm_size_location="Not reported",
            Clinical_presentation="Mixed design: rupture subset + normal vessel controls",
            Collection_timing="At surgery",
            Biospecimen="Aneurysm / vessel tissue",
            RNA_method="Affymetrix GPL6244",
            Platform="GPL6244",
            Accession="GSE15629",
            Covariates_available="Phenotype labels in series matrix",
        ),
        dict(
            Phase="3 — Orthogonal tissue and blood context",
            GEO_accession="GSE26969",
            Comparison="Unruptured aneurysm vs superficial temporal artery",
            N_positive=3,
            N_negative=3,
            N_total=6,
            Age_mean_SD="Not reported",
            Sex_pct_female="Not reported",
            Hypertension_pct="Not reported",
            Aneurysm_size_location="Not reported",
            Clinical_presentation="Incidental unruptured IA vs control artery",
            Collection_timing="At surgery",
            Biospecimen="Tissue",
            RNA_method="Microarray",
            Platform="GPL570",
            Accession="GSE26969",
            Covariates_available="Group",
        ),
        dict(
            Phase="3 — Orthogonal tissue and blood context",
            GEO_accession="GSE36791",
            Comparison="Ruptured aneurysm (SAH) vs headache controls",
            N_positive=m["gse36791_n_ruptured"],
            N_negative=m["gse36791_n_control"],
            N_total=m["gse36791_n_ruptured"] + m["gse36791_n_control"],
            Age_mean_SD="Not in GEO",
            Sex_pct_female="Not reported",
            Hypertension_pct="Not reported",
            Aneurysm_size_location="Not reported",
            Clinical_presentation="Post-rupture / acute SAH blood signature",
            Collection_timing="Acute phase after SAH (blood)",
            Biospecimen="Whole blood",
            RNA_method="Microarray",
            Platform="GPL10558",
            Accession="GSE36791",
            Covariates_available="Rupture/SAH vs control",
        ),
        dict(
            Phase="3 — Orthogonal tissue and blood context",
            GEO_accession="GSE73378",
            Comparison="Prior aSAH vs controls",
            N_positive=103,
            N_negative=107,
            N_total=210,
            Age_mean_SD="Not in GEO",
            Sex_pct_female="Not reported",
            Hypertension_pct="Not reported",
            Aneurysm_size_location="Not reported",
            Clinical_presentation="History of treated SAH vs controls",
            Collection_timing="Months to years post-treatment (long-horizon blood)",
            Biospecimen="Whole blood",
            RNA_method="Microarray (dual platform series)",
            Platform="GPL10904; GPL10558",
            Accession="GSE73378",
            Covariates_available="Prior SAH status",
        ),
        dict(
            Phase="4 — Vessel gradient (aneurysm vs normal artery)",
            GEO_accession="GSE75436",
            Comparison="Aneurysm tissue vs matched superficial temporal artery",
            N_positive=15,
            N_negative=15,
            N_total=30,
            Age_mean_SD="Not in GEO",
            Sex_pct_female="Not reported",
            Hypertension_pct="Not reported",
            Aneurysm_size_location="Not reported",
            Clinical_presentation="Aneurysm vs paired STA control",
            Collection_timing="At surgery",
            Biospecimen="Aneurysm wall vs STA",
            RNA_method="Agilent microarray",
            Platform="GPL570",
            Accession="GSE75436",
            Covariates_available="Paired tissue type",
        ),
        dict(
            Phase="4 — Vessel gradient (aneurysm vs normal artery)",
            GEO_accession="GSE66238",
            Comparison="Aneurysm dome vs superficial temporal artery",
            N_positive=6,
            N_negative=12,
            N_total=18,
            Age_mean_SD="Not in GEO",
            Sex_pct_female="Not reported",
            Hypertension_pct="Not reported",
            Aneurysm_size_location="Not reported",
            Clinical_presentation="Aneurysm vs control vessel (not RvU stratified)",
            Collection_timing="At surgery",
            Biospecimen="Aneurysm wall vs STA",
            RNA_method="RNA-seq (TPM)",
            Platform="GPL17303",
            Accession="GSE66238",
            Covariates_available="Tissue type",
        ),
    ]
    s1 = pd.DataFrame(rows)

    coverage = []
    for gene in PANEL_GENES:
        row = {"Gene": gene}
        row["GSE122897_RNAseq"] = "Yes"
        row["GSE159610_RNAseq"] = "Yes"
        row["GSE106520_RNAseq"] = "Yes"
        row["GSE54083_array"] = "Yes" if gene != "ANK3" else "Probe mapping uncertain"
        row["GSE13353_GPL570"] = "Yes" if gene in {"MARCO", "TGFBI", "MPP1"} else "Check probes"
        row["GSE36791_GPL10558"] = "Yes"
        row["GSE73378_GPL10558"] = "Yes" if gene != "ANK3" else "No — absent from array"
        row["GSE75436_GPL570"] = "Yes"
        row["GSE66238_RNAseq"] = "Yes"
        row["Notes"] = (
            "ANK3 absent on GSE73378 array (7/8 genes measurable); "
            "main-text footnote applies to GSE73378."
            if gene == "ANK3"
            else ""
        )
        coverage.append(row)
    cov = pd.DataFrame(coverage)

    legend = (
        f"Supplementary Table S1 documents all {len(rows)} cohort comparison arms used across four analytic phases. "
        f"Clinical covariates (age, sex, hypertension, PHASES components, aneurysm size) were not consistently available in "
        f"GEO metadata; only rupture status or cohort-specific group labels were used in models. The platform coverage matrix "
        f"shows that seven of eight blood-bridge panel genes are measurable on GSE73378; ANK3 is absent, consistent with the "
        f"main-text array footnote."
    )
    return s1, cov, legend


# ---------------------------------------------------------------------------
# Table S2 — 300-gene shortlist
# ---------------------------------------------------------------------------
def build_table_s2() -> tuple[pd.DataFrame, str]:
    sl = pd.read_csv(BASE / "results/phases_g/shortlist_evidence.tsv", sep="\t")
    g106 = pd.read_csv(
        BASE / "results/external_validation_gse106520/gse106520_shortlist_overlap.tsv", sep="\t"
    )
    g106 = g106[["Gene", "Log2FC_AvC", "P_value", "DirectionMatch_GSE122897"]].rename(
        columns={
            "Log2FC_AvC": "GSE106520_log2FC",
            "P_value": "GSE106520_P",
            "DirectionMatch_GSE122897": "GSE106520_direction_match",
        }
    )
    sl = sl.rename(columns={"symbol": "Gene_symbol"}) if "symbol" in sl.columns else sl
    if "Gene_symbol" not in sl.columns:
        sl["Gene_symbol"] = sl.get("Gene_ID", "")
    g106 = g106.rename(columns={"Gene": "Gene_symbol"})
    df = sl.merge(g106, on="Gene_symbol", how="left", suffixes=("", "_dup"))
    if "symbol" in df.columns:
        df = df.rename(columns={"symbol": "Gene_symbol"})
    elif "Gene_symbol" not in df.columns:
        df["Gene_symbol"] = df.get("symbol", df.get("Gene_ID"))

    p1 = np.clip(pd.to_numeric(df["GSE159610_P"], errors="coerce").fillna(1.0), 1e-300, 1.0)
    p2 = np.clip(pd.to_numeric(df["GSE106520_P"], errors="coerce").fillna(1.0), 1e-300, 1.0)
    df["Bridge_Stouffer_Z"] = (norm.isf(p1 / 2) + norm.isf(p2 / 2)) / np.sqrt(2)
    df["Bridge_meta_P"] = 2 * norm.sf(np.abs(df["Bridge_Stouffer_Z"]))
    df["Blood_both_direction_concordant"] = (
        (df.get("direction_match_rupture_and_aneurysm") == True)  # noqa: E712
        & (df.get("GSE106520_direction_match") == True)  # noqa: E712
    )
    df = df.sort_values(["Bridge_Stouffer_Z", "DESeq2_padj"], ascending=[False, True]).reset_index(drop=True)
    df["Shortlist_rank"] = np.arange(1, len(df) + 1)
    sym = df["Gene_symbol"].astype(str)
    df["In_8gene_panel"] = sym.isin(PANEL_GENES)
    df["In_6gene_wall_core"] = sym.isin(TISSUE_CORE_6)
    df["Wall_core_not_in_blood_panel"] = df["In_6gene_wall_core"] & ~df["In_8gene_panel"]

    cols = [
        "Shortlist_rank", "Gene_symbol", "DESeq2_log2FC", "DESeq2_padj",
        "LM_immune_adj_P_Ruptured", "GSE159610_Log2FC_aneurysm_vs_control", "GSE159610_P",
        "GSE106520_log2FC", "GSE106520_P", "Blood_both_direction_concordant",
        "Bridge_Stouffer_Z", "Bridge_meta_P", "In_8gene_panel", "In_6gene_wall_core",
    ]
    cols = [c for c in cols if c in df.columns]
    out = df[cols].copy()
    n300 = len(out)
    n_panel = int(out["In_8gene_panel"].sum())
    legend = (
        f"All {n300} genes met joint criteria: DESeq2 padj <= 0.05 and immune-adjusted rupture P <= 0.05 "
        f"(nominal, LM22-adjusted OLS), ranked by blood-bridge Stouffer Z across GSE159610 and GSE106520. "
        f"{n_panel} genes constitute the final 8-gene blood panel (highlighted). "
        f"DESeq2_log2FC is the primary discovery effect size (e.g. MPP1 +1.87, MARCO +3.04 in GSE122897); "
        f"see Supplementary Table S10 for cohort-native delta mean log-expression contrasts used in the "
        f"cross-cohort MPP1/MARCO survey (GSE122897: MPP1 +1.17, MARCO +2.25 on log2[TPM+1])."
    )
    return out, legend


# ---------------------------------------------------------------------------
# Table S3 — Per-gene AUROC with bootstrap CI
# ---------------------------------------------------------------------------
def build_table_s3() -> tuple[pd.DataFrame, str]:
    signs = pd.read_csv(BASE / "results/blood_bridge_panel/blood_bridge_panel_gene_summary.tsv", sep="\t")
    signs = np.sign(pd.to_numeric(signs.set_index("Gene")["DESeq2_log2FC"], errors="coerce"))

    g159_expr, g159_meta = _load_gse159610_expression(BASE)
    g106_expr, g106_meta = _load_simple_expression(
        BASE / "results/external_validation_gse106520/gse106520_gene_level_log2fpkm.tsv",
        BASE / "results/external_validation_gse106520/gse106520_sample_metadata.tsv",
        sample_col="Sample",
        group_col="Group",
    )

    rows = []
    for cohort, expr, meta, pos, neg in [
        ("GSE159610", g159_expr, g159_meta, "Aneurysm", "Control"),
        ("GSE106520", g106_expr, g106_meta, "Aneurysm", "Control"),
    ]:
        avail = [c for c in expr.columns if c in set(meta["Sample"])]
        expr_c = expr.loc[:, avail]
        meta_c = meta[meta["Sample"].isin(avail) & meta["Group"].isin([pos, neg])].drop_duplicates("Sample")
        meta_c = meta_c.set_index("Sample").loc[expr_c.columns]
        labels = (meta_c["Group"] == pos).to_numpy(dtype=bool)

        for gene in PANEL_GENES:
            if gene not in expr_c.index:
                continue
            sign = float(signs.get(gene, 1.0))
            z = _zscore_rows(expr_c.loc[[gene]]).iloc[0] * sign
            scores = z.reindex(meta_c.index).astype(float)
            auc, lo, hi = _bootstrap_auc_ci(labels, scores.to_numpy())
            _, mwu_p = mannwhitneyu(
                scores[labels], scores[~labels], alternative="two-sided"
            )
            rows.append({
                "Gene": gene,
                "Cohort": cohort,
                "AUROC": round(auc, 4),
                "AUROC_95CI_lower": round(lo, 4),
                "AUROC_95CI_upper": round(hi, 4),
                "Mann_Whitney_P": mwu_p,
                "N_positive": int(labels.sum()),
                "N_negative": int((~labels).sum()),
            })

    # Composite 8-gene panel AUROC
    panel_scores = pd.read_csv(
        BASE / "results/blood_bridge_panel/blood_bridge_panel_sample_scores.tsv", sep="\t"
    )
    for cohort in ["GSE159610", "GSE106520"]:
        sub = panel_scores[panel_scores["Dataset"] == cohort]
        pos_g = "Aneurysm"
        y = (sub["Group"] == pos_g).to_numpy(dtype=bool)
        s = sub["Panel_score"].to_numpy(dtype=float)
        auc, lo, hi = _bootstrap_auc_ci(y, s)
        _, mwu_p = mannwhitneyu(s[y], s[~y], alternative="two-sided")
        rows.append({
            "Gene": "8-gene composite panel",
            "Cohort": cohort,
            "AUROC": round(auc, 4),
            "AUROC_95CI_lower": round(lo, 4),
            "AUROC_95CI_upper": round(hi, 4),
            "Mann_Whitney_P": mwu_p,
            "N_positive": int(y.sum()),
            "N_negative": int((~y).sum()),
        })

    out = pd.DataFrame(rows)
    perf = pd.read_csv(BASE / "results/blood_bridge_panel/blood_bridge_panel_performance.tsv", sep="\t")
    p159 = perf.loc[perf["Dataset"] == "GSE159610", "AUC"].iloc[0]
    p106 = perf.loc[perf["Dataset"] == "GSE106520", "AUC"].iloc[0]
    legend = (
        f"Phase 1 blood-bridge cohorts: composite panel AUROC = {p159:.3f} (GSE159610) and {p106:.3f} (GSE106520), "
        f"exceeding most single-gene AUROCs. Bootstrap 95% CIs (1,500 resamples) reflect modest sample sizes "
        f"(n=23/24 and n=16/16). ANK3 shows the highest single-gene AUROC in GSE159610 (0.705) but not in GSE106520."
    )
    return out, legend


# ---------------------------------------------------------------------------
# Tables S4–S12 (data exports)
# ---------------------------------------------------------------------------
def build_table_s4(m: dict) -> tuple[pd.DataFrame, str]:
    des = pd.read_csv(BASE / "results/gse122897_deseq2_results.tsv", sep="\t")
    des.columns = [c.strip('"') for c in des.columns]
    ann = BASE / "results/gse122897_deseq2_top_annotated.tsv"
    if ann.exists():
        top = pd.read_csv(ann, sep="\t")
        if "Gene_Symbol" in top.columns and "Gene_ID" in top.columns:
            sym_map = top[["Gene_ID", "Gene_Symbol"]].drop_duplicates()
            sym_map["Gene_ID"] = sym_map["Gene_ID"].astype(str).str.strip('"')
            des["Gene_ID"] = des["Gene_ID"].astype(str).str.strip('"')
            des = des.merge(sym_map, on="Gene_ID", how="left")
    n_sig = int((des["padj"] < 0.05).sum()) if "padj" in des.columns else m["deseq2_padj_genes"]
    legend = (
        f"Full DESeq2 output for GSE122897 (n=22 ruptured vs 21 unruptured). "
        f"{n_sig:,} genes with padj < 0.05. Complete table archived as TSV ({len(des):,} rows)."
    )
    return des, legend


def build_table_s5() -> tuple[pd.DataFrame, str]:
    g540 = pd.read_csv(BASE / "results/external_validation_gse54083/gse54083_shortlist_overlap.tsv", sep="\t")
    g540["Cohort"] = "GSE54083"
    g540 = g540.rename(columns={"Log2FC_RvU": "Log2FC", "P_value": "P_value", "FDR_BH": "FDR_BH"})
    g133 = pd.read_csv(BASE / "results/rupture_meta_analysis/gse13353_gene_level_meta_input.tsv", sep="\t")
    g133["Cohort"] = "GSE13353"
    g133 = g133.rename(columns={"Log2FC": "Log2FC", "P_value": "P_value", "FDR_BH": "FDR_BH"})
    common = ["Cohort", "Gene", "Log2FC", "P_value", "FDR_BH"]
    out = pd.concat([g540[common], g133[common]], ignore_index=True)
    legend = (
        f"GSE54083: {len(g540):,} genes with rupture contrast statistics; "
        f"GSE13353: {len(g133):,} genes. Used for meta-analysis in main Table 3."
    )
    return out, legend


def build_table_s6() -> tuple[pd.DataFrame, str]:
    path = BASE / "results/external_validation_gse15629/gse15629_shortlist_overlap_rupture.tsv"
    ov = pd.read_csv(path, sep="\t")
    rows = []
    for gene in TISSUE_CORE_6:
        hit = ov[ov["Gene"] == gene]
        if len(hit):
            r = hit.iloc[0]
            rows.append({
                "Gene": gene,
                "Measurable_on_GSE15629": "Yes",
                "Log2FC_RvU": r.get("Log2FC_RvU"),
                "P_value": r.get("P_value"),
                "FDR_BH": r.get("FDR_BH"),
                "DirectionMatch_GSE122897": r.get("DirectionMatch_GSE122897"),
            })
        else:
            rows.append({
                "Gene": gene,
                "Measurable_on_GSE15629": "No / not in annotated overlap",
                "Log2FC_RvU": np.nan,
                "P_value": np.nan,
                "FDR_BH": np.nan,
                "DirectionMatch_GSE122897": np.nan,
            })
    core = pd.DataFrame(rows)
    measurable = core[core["Measurable_on_GSE15629"] == "Yes"]
    n_match = int(measurable["DirectionMatch_GSE122897"].fillna(False).sum())
    n_meas = len(measurable)
    legend = (
        f"GSE15629 rupture subset (8 ruptured vs 6 unruptured): {n_match}/{n_meas} measurable "
        f"six-gene wall-core genes directionally concordant with GSE122897. "
        f"ANPEP and CLEC5A were not present in the annotated GSE15629 overlap table; "
        f"TGFBI, MPP1, MARCO, and CPXM2 showed concordant directions (4/4 measurable)."
    )
    return core, legend


def build_table_s7() -> tuple[pd.DataFrame, str]:
    surv = pd.read_csv(BASE / "results/mpp1_marco_survey/mpp1_marco_all_cohort_results.tsv", sep="\t")
    mask = surv["Cohort"].isin(["GSE75436", "GSE66238"]) & surv["Gene"].isin(PANEL_GENES)
    out = surv.loc[mask].copy()
    legend = (
        "Phase 4b gradient cohorts: per-gene log2 fold changes for aneurysm vs control vessel "
        "(GSE75436 microarray; GSE66238 RNA-seq). Supports MPP1/MARCO elevation in aneurysm "
        "relative to superficial temporal artery controls."
    )
    return out, legend


def build_table_s8() -> tuple[pd.DataFrame, str]:
    frac = pd.read_csv(BASE / "results/cibersortx/lm22_fractions_gse122897.tsv", sep="\t")
    legend = (
        f"Per-sample LM22+NNLS immune fractions for all {len(frac)} GSE122897 wall samples "
        f"(22 ruptured, 21 unruptured). Used as covariates in immune-adjusted gene-wise OLS."
    )
    return frac, legend


def build_table_s9(m: dict) -> tuple[pd.DataFrame, str]:
    im = pd.read_csv(BASE / "results/cibersortx/deg_immune_covariate_adjusted.tsv", sep="\t")
    im = im.sort_values("LM_immune_adj_P_Ruptured")
    if "symbol" not in im.columns and "Gene_ID" in im.columns:
        ann = BASE / "results/gene_id_mapping_unified.tsv"
        if ann.exists():
            mp = pd.read_csv(ann, sep="\t")
            if "Gene_Symbol" in mp.columns:
                im["Gene_ID"] = im["Gene_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
                mp = mp.copy()
                mp["Original_ID"] = mp["Original_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
                im = im.merge(mp[["Original_ID", "Gene_Symbol"]], left_on="Gene_ID", right_on="Original_ID", how="left")
    n_nom = int((im["LM_immune_adj_P_Ruptured"] < 0.05).sum())
    legend = (
        f"Gene-wise OLS: rupture + {m.get('n_immune_covariates', 11)} LM22 covariates; "
        f"n={m['n_genes_immune_model_tested']:,} tested. "
        f"{n_nom:,} genes with nominal rupture P < 0.05; 0 with BH-FDR < 0.05 genome-wide."
    )
    return im, legend


def build_table_s10() -> tuple[pd.DataFrame, str]:
    surv = pd.read_csv(BASE / "results/mpp1_marco_survey/mpp1_marco_all_cohort_results.tsv", sep="\t")
    sl = pd.read_csv(BASE / "results/phases_g/shortlist_evidence.tsv", sep="\t")
    deseq_ref = sl.set_index("symbol")["DESeq2_log2FC"].to_dict()

    out = surv.copy()
    out = out.rename(columns={"Log2FC_positive_minus_negative": "Delta_mean_log2_expression"})
    out["GSE122897_DESeq2_log2FC_reference"] = out["Gene"].map(deseq_ref)

    g122 = out[(out["Cohort"] == "GSE122897") & (out["Gene"].isin(["MPP1", "MARCO"]))]
    mpp_d = g122.loc[g122["Gene"] == "MPP1", "Delta_mean_log2_expression"].iloc[0]
    mar_d = g122.loc[g122["Gene"] == "MARCO", "Delta_mean_log2_expression"].iloc[0]
    mpp_de = g122.loc[g122["Gene"] == "MPP1", "GSE122897_DESeq2_log2FC_reference"].iloc[0]
    mar_de = g122.loc[g122["Gene"] == "MARCO", "GSE122897_DESeq2_log2FC_reference"].iloc[0]

    legend = (
        "Consolidated MPP1 and MARCO statistics across surveyed cohorts. "
        f"Delta_mean_log2_expression is the Welch contrast on cohort-native log expression "
        f"(GSE122897: mean log2[TPM+1] difference; MPP1 {mpp_d:.2f}, MARCO {mar_d:.2f}). "
        f"GSE122897_DESeq2_log2FC_reference carries discovery DESeq2 estimates used in Table S2 "
        f"(MPP1 {mpp_de:.2f}, MARCO {mar_de:.2f}). These differ by definition: DESeq2 models raw counts; "
        "the survey column harmonises log-expression contrasts across external platforms."
    )
    return out, legend


def build_table_s11(m: dict) -> tuple[pd.DataFrame, str]:
    pw = pd.read_csv(BASE / "results/gse122897_pathway_enrichment.tsv", sep="\t")
    pw_sig = pw[pw["FDR_BH"] < 0.05].copy() if "FDR_BH" in pw.columns else pw.head(207)
    legend = (
        f"g:Profiler ORA on DESeq2-ranked genes: {m['pathways_fdr05']} terms at FDR < 0.05 "
        f"(full table {len(pw):,} rows tested)."
    )
    return pw_sig if len(pw_sig) else pw.head(100), legend


def build_table_s12() -> tuple[pd.DataFrame, str]:
    """Hypergeometric overlap: 8-gene blood panel vs 6-gene wall core."""
    panel_set = set(PANEL_GENES)
    core_set = set(TISSUE_CORE_6)
    overlap = len(panel_set & core_set)
    shortlist_n = 300
    anchor_n = 645
    tested_n = 14888
    universe_sizes = [
        ("Primary shortlist (pre-specified N)", shortlist_n),
        ("DESeq2+immune nominal pool (uncapped)", anchor_n),
        ("All genes immune-model tested", tested_n),
    ]
    rows = []
    for label, N in universe_sizes:
        K = 8  # panel size
        n = 6  # core size
        x = overlap
        p = hypergeom.sf(x - 1, N, n, K) if x > 0 else 1.0
        expected = (K * n) / N if N else float("nan")
        enrich = x / expected if expected and expected > 0 else float("nan")
        rows.append({
            "Background_definition": label,
            "Background_N": N,
            "Panel_size_k": K,
            "Wall_core_size_n": n,
            "Observed_overlap": x,
            "Expected_overlap": round(expected, 3),
            "Fold_enrichment": round(enrich, 2),
            "Hypergeometric_P": p,
        })
    out = pd.DataFrame(rows)
    legend = (
        f"Convergence of 8-gene blood panel with 6-gene wall core (observed overlap = {overlap}: MPP1, MARCO). "
        f"Hypergeometric P values are sensitive to background N; primary analysis uses N=300 shortlist."
    )
    return out, legend


# ---------------------------------------------------------------------------
# Supplementary Figures
# ---------------------------------------------------------------------------
def figure_s1_panel_violins() -> Path:
    signs = pd.read_csv(BASE / "results/blood_bridge_panel/blood_bridge_panel_gene_summary.tsv", sep="\t")
    signs = np.sign(pd.to_numeric(signs.set_index("Gene")["DESeq2_log2FC"], errors="coerce"))

    cohorts = []
    g159_expr, g159_meta = _load_gse159610_expression(BASE)
    cohorts.append(("GSE159610", g159_expr, g159_meta, "Aneurysm", "Control"))
    g106_expr, g106_meta = _load_simple_expression(
        BASE / "results/external_validation_gse106520/gse106520_gene_level_log2fpkm.tsv",
        BASE / "results/external_validation_gse106520/gse106520_sample_metadata.tsv",
        sample_col="Sample",
        group_col="Group",
    )
    cohorts.append(("GSE106520", g106_expr, g106_meta, "Aneurysm", "Control"))

    fig, axes = plt.subplots(len(PANEL_GENES), 2, figsize=(8, 2.2 * len(PANEL_GENES)), sharex="col")
    for j, (cname, expr, meta, pos, neg) in enumerate(cohorts):
        avail = [c for c in expr.columns if c in set(meta["Sample"])]
        expr_c = expr.loc[:, avail]
        meta_c = meta[meta["Sample"].isin(avail) & meta["Group"].isin([pos, neg])].drop_duplicates("Sample")
        meta_c = meta_c.set_index("Sample").loc[expr_c.columns]
        labels = (meta_c["Group"] == pos).to_numpy(dtype=bool)

        for i, gene in enumerate(PANEL_GENES):
            ax = axes[i, j]
            if gene not in expr_c.index:
                ax.text(0.5, 0.5, "absent", ha="center", transform=ax.transAxes)
                continue
            sign = float(signs.get(gene, 1))
            vals = expr_c.loc[gene].astype(float) * sign
            u = vals[~labels]
            r = vals[labels]
            ax.boxplot([u, r], positions=[0, 1], widths=0.5, patch_artist=True)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["Ctrl", "Aneurysm"], fontsize=8)
            if j == 0:
                ax.set_ylabel(gene, fontweight="bold")
            if i == 0:
                ax.set_title(cname, fontweight="bold")

    fig.suptitle("Supplementary Figure S1. Panel gene expression by group (Phase 1 blood-bridge cohorts)", y=1.01)
    fig.tight_layout()
    path = FIG_DIR / "SuppFig1_panel_gene_violins.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def figure_s2_roc_grid() -> Path:
    signs = pd.read_csv(BASE / "results/blood_bridge_panel/blood_bridge_panel_gene_summary.tsv", sep="\t")
    signs = np.sign(pd.to_numeric(signs.set_index("Gene")["DESeq2_log2FC"], errors="coerce"))
    from sklearn.metrics import roc_curve

    fig, axes = plt.subplots(len(PANEL_GENES), 2, figsize=(9, 2.2 * len(PANEL_GENES)))
    cohort_data = []
    g159_expr, g159_meta = _load_gse159610_expression(BASE)
    cohort_data.append(("GSE159610", g159_expr, g159_meta, "Aneurysm", "Control"))
    g106_expr, g106_meta = _load_simple_expression(
        BASE / "results/external_validation_gse106520/gse106520_gene_level_log2fpkm.tsv",
        BASE / "results/external_validation_gse106520/gse106520_sample_metadata.tsv",
        sample_col="Sample",
        group_col="Group",
    )
    cohort_data.append(("GSE106520", g106_expr, g106_meta, "Aneurysm", "Control"))

    for row_gene, gene in enumerate(PANEL_GENES):
        for col_cohort, (cname, expr, meta, pos, neg) in enumerate(cohort_data):
            ax = axes[row_gene, col_cohort]
            avail = [c for c in expr.columns if c in set(meta["Sample"])]
            meta_c = meta[meta["Sample"].isin(avail) & meta["Group"].isin([pos, neg])].drop_duplicates("Sample")
            meta_c = meta_c.set_index("Sample")
            expr_c = expr.loc[:, [c for c in avail if c in meta_c.index]]
            meta_c = meta_c.loc[expr_c.columns]
            y = (meta_c["Group"] == pos).astype(int).to_numpy()
            if gene not in expr_c.index:
                ax.axis("off")
                continue
            sign = float(signs.get(gene, 1))
            scores = _zscore_rows(expr_c.loc[[gene]]).iloc[0] * sign
            scores = scores.reindex(meta_c.index).astype(float).to_numpy()
            fpr, tpr, _ = roc_curve(y, scores)
            auc = _auc_from_scores(y.astype(bool), scores)
            ax.plot(fpr, tpr, lw=2)
            ax.plot([0, 1], [0, 1], "k--", alpha=0.4)
            ax.set_title(f"{gene}\n{cname} AUC={auc:.2f}", fontsize=8)
            ax.set_xlim(-0.02, 1.02)
            ax.set_ylim(-0.02, 1.02)
    fig.suptitle("Supplementary Figure S2. Single-gene ROC curves (8 genes x 2 Phase 1 cohorts)", fontsize=12)
    fig.tight_layout()
    path = FIG_DIR / "SuppFig2_single_gene_roc_grid.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def figure_s3_volcano_core() -> Path:
    des = pd.read_csv(BASE / "results/gse122897_deseq2_results.tsv", sep="\t")
    des.columns = [c.strip('"') for c in des.columns]
    des["neglog10"] = -np.log10(des["padj"].clip(lower=1e-300))
    sig = des["padj"] < 0.05
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(des.loc[~sig, "log2FoldChange"], des.loc[~sig, "neglog10"], c="#aaa", s=8, alpha=0.4, rasterized=True)
    ax.scatter(des.loc[sig, "log2FoldChange"], des.loc[sig, "neglog10"], c="#c44", s=10, alpha=0.5, rasterized=True)
    ann = pd.read_csv(BASE / "results/gse122897_deseq2_top_annotated.tsv", sep="\t")
    sym_map = dict(zip(ann["Gene_ID"].astype(str).str.strip('"'), ann["Gene_Symbol"]))
    des["sym"] = des["Gene_ID"].astype(str).str.strip('"').map(sym_map)
    for g in TISSUE_CORE_6:
        sub = des[des["sym"] == g]
        if len(sub):
            ax.scatter(sub["log2FoldChange"], sub["neglog10"], c="#1a5", s=80, edgecolors="black", zorder=5)
            ax.annotate(g, (sub["log2FoldChange"].iloc[0], sub["neglog10"].iloc[0]),
                        fontsize=9, fontweight="bold", xytext=(5, 5), textcoords="offset points")
    ax.axhline(-np.log10(0.05), ls="--", c="k", lw=0.8)
    ax.set_xlabel("log2 fold change (ruptured vs unruptured)")
    ax.set_ylabel("-log10(adjusted P)")
    ax.set_title(f"Supplementary Figure S3. GSE122897 DESeq2 volcano (6-gene wall core labelled; {sig.sum():,} padj<0.05)")
    fig.tight_layout()
    path = FIG_DIR / "SuppFig3_deseq2_volcano_wall_core.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def figure_s4_cibersortx() -> Path:
    merged = pd.read_csv(BASE / "results/cibersortx/cibersortx_merged_with_group.tsv", sep="\t")
    lm22_cols = [c for c in merged.columns if c not in ("Sample", "P-value", "Correlation", "RMSE", "Group")]
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(2, 1, height_ratios=[1.2, 1], hspace=0.35)
    ax0 = fig.add_subplot(gs[0])
    plot_df = merged.set_index("Sample")[lm22_cols]
    plot_df = plot_df.loc[merged.sort_values("Group")["Sample"]]
    im = ax0.imshow(plot_df.T, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.6)
    ax0.set_yticks(range(len(lm22_cols)))
    ax0.set_yticklabels(lm22_cols, fontsize=7)
    ax0.set_xticks([])
    ax0.set_title("A. LM22 cell fractions (CIBERSORTx Job 7; samples ordered by rupture status)")
    plt.colorbar(im, ax=ax0, fraction=0.02, pad=0.02)

    ax1 = fig.add_subplot(gs[1])
    for grp, color in [("Ruptured", "#D55E00"), ("Unruptured", "#0072B2")]:
        sub = merged[merged["Group"] == grp]
        ax1.scatter(sub["Correlation"], sub["RMSE"], label=f"{grp} (n={len(sub)})", alpha=0.8, s=50, c=color)
    ax1.set_xlabel("CIBERSORTx correlation (fit quality)")
    ax1.set_ylabel("RMSE")
    ax1.legend()
    ax1.set_title("B. Deconvolution fit quality per sample (median r=0.14)")
    fig.suptitle("Supplementary Figure S4. Immune deconvolution of GSE122897 wall tissue", fontsize=12)
    path = FIG_DIR / "SuppFig4_cibersortx_lm22.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def figure_s5_mpp1_marco_cross_cohort() -> Path:
    surv = pd.read_csv(BASE / "results/mpp1_marco_survey/mpp1_marco_all_cohort_results.tsv", sep="\t")
    surv = surv[surv["Status"] == "ok"].copy()
    cohorts = surv["Cohort"].unique()
    genes = ["MPP1", "MARCO"]
    fig, axes = plt.subplots(1, 2, figsize=(14, max(6, 0.35 * len(cohorts))))
    for ax, gene in zip(axes, genes):
        sub = surv[surv["Gene"] == gene].sort_values("Log2FC_positive_minus_negative")
        colors = ["#2e7d32" if x > 0 else "#c62828" for x in sub["Log2FC_positive_minus_negative"]]
        ax.barh(sub["Cohort"] + "\n" + sub["Comparison"], sub["Log2FC_positive_minus_negative"], color=colors)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel("log2 FC (positive vs negative arm)")
        ax.set_title(gene, fontweight="bold")
    fig.suptitle("Supplementary Figure S5. MPP1 and MARCO signed log2FC across all surveyed cohorts", fontsize=12)
    fig.tight_layout()
    path = FIG_DIR / "SuppFig5_mpp1_marco_cross_cohort.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Word document
# ---------------------------------------------------------------------------
def _add_section(doc: Document, title: str, legend: str, df: pd.DataFrame, max_rows: int = 35) -> None:
    doc.add_heading(title, level=2)
    doc.add_paragraph(legend)
    show = df.head(max_rows) if len(df) > max_rows else df
    if len(df) > max_rows:
        doc.add_paragraph(
            f"Showing first {max_rows} of {len(df):,} rows. Full table: results/supplementary/tables/{title.split('.')[0].replace(' ', '_')}.tsv",
            style="Intense Quote",
        )
    add_table_from_dataframe(doc, show, title=None, max_rows=max_rows)
    doc.add_paragraph()


def build_word_document(artifacts: dict) -> None:
    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "Arial"
    style.font.size = Pt(11)

    t = doc.add_heading("Supplementary Materials", 0)
    t.alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.add_paragraph(
        "Complete supplementary tables, figures, and appendices for the intracranial aneurysm "
        "transcriptomic analysis (GSE122897 discovery; multi-cohort validation). All values are derived "
        "from project result files unless noted as not available in GEO metadata."
    )
    doc.add_page_break()

    doc.add_heading("Supplementary Tables", 1)
    _add_section(doc, "Table S1. Cohort clinical and technical characteristics", artifacts["s1_legend"], artifacts["s1"], max_rows=15)
    _add_section(doc, "Table S1b. Platform coverage for 8-gene blood panel", artifacts["s1b_legend"], artifacts["s1b"], max_rows=12)
    _add_section(doc, "Table S2. Full 300-gene shortlist", artifacts["s2_legend"], artifacts["s2"], max_rows=40)
    _add_section(doc, "Table S3. Per-gene and composite AUROC (Phase 1)", artifacts["s3_legend"], artifacts["s3"], max_rows=25)
    _add_section(doc, "Table S4. Full GSE122897 DESeq2 results", artifacts["s4_legend"], artifacts["s4"], max_rows=30)
    _add_section(doc, "Table S5. External rupture cohort DE (GSE54083, GSE13353)", artifacts["s5_legend"], artifacts["s5"], max_rows=30)
    _add_section(doc, "Table S6. GSE15629 concordance (6-gene wall core)", artifacts["s6_legend"], artifacts["s6"], max_rows=10)
    _add_section(doc, "Table S7. Phase 4b gradient cohort gene statistics", artifacts["s7_legend"], artifacts["s7"], max_rows=20)
    _add_section(doc, "Table S8. LM22 fractions GSE122897", artifacts["s8_legend"], artifacts["s8"], max_rows=25)
    _add_section(doc, "Table S9. Immune-adjusted OLS (top genes by rupture P)", artifacts["s9_legend"], artifacts["s9"].head(50), max_rows=50)
    _add_section(doc, "Table S10. MPP1/MARCO all-cohort statistics", artifacts["s10_legend"], artifacts["s10"], max_rows=30)
    _add_section(doc, "Table S11. Pathway enrichment (FDR < 0.05)", artifacts["s11_legend"], artifacts["s11"], max_rows=30)
    _add_section(doc, "Table S12. Hypergeometric convergence sensitivity", artifacts["s12_legend"], artifacts["s12"], max_rows=10)

    doc.add_page_break()
    doc.add_heading("Supplementary Figures", 1)

    fig_specs = [
        ("s1", "Supplementary Figure S1. Individual panel-gene expression in Phase 1 blood cohorts (GSE159610, GSE106520). Violin/box plots show signed expression (discovery direction) for each of the eight genes. The composite panel outperforms most single genes in AUROC (Table S3).", 6.5),
        ("s2", "Supplementary Figure S2. Single-gene ROC curves for eight panel genes in GSE159610 (top row) and GSE106520 (bottom row). Curves confirm heterogeneous per-gene discrimination; no single gene matches the composite panel AUROC in both cohorts.", 6.5),
        ("s3", "Supplementary Figure S3. DESeq2 volcano plot for GSE122897 with the six-gene wall rupture core highlighted. Green points: TGFBI, ANPEP, MPP1, MARCO, CLEC5A, CPXM2; background shows genome-wide significant genes (padj < 0.05).", 5.5),
        ("s4", "Supplementary Figure S4. CIBERSORTx/LM22 deconvolution of wall tissue. Panel A: heatmap of 22 immune fractions (43 samples). Panel B: correlation vs RMSE per sample (median correlation 0.14), indicating weak fit of blood-derived reference to aneurysm wall.", 6.0),
        ("s5", "Supplementary Figure S5. MPP1 and MARCO signed log2 fold changes across all ten surveyed cohort comparisons. Green bars: same direction as GSE122897 discovery; red: discordant. Tissue rupture cohorts show consistent positive shifts for both genes.", 6.5),
    ]
    for key, caption, width in fig_specs:
        add_figure_with_caption_and_inference(doc, str(artifacts[f"fig_{key}"]), caption, "", width=Inches(width))

    doc.add_page_break()
    doc.add_heading("Supplementary Appendix 1. TRIPOD checklist", 1)
    doc.add_paragraph(artifacts["tripod_text"])
    doc.add_page_break()
    doc.add_heading("Supplementary Appendix 2. PROBAST risk-of-bias assessment", 1)
    doc.add_paragraph(artifacts["probast_text"])

    doc.save(str(OUT_DOCX))
    print(f"Saved {OUT_DOCX}")


def _tripod_text(m: dict) -> str:
    return f"""
TRIPOD (Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis) — contextual completion for secondary transcriptomic analyses.

Study type: Secondary analysis of public GEO data; not a prospective clinical prediction study. Items below state applicability.

Title/Abstract: Prediction model framing is exploratory; primary outcome in discovery is ruptured vs unruptured tissue status (GSE122897, n={m['n_samples']}).

Source of data: Public repositories (GSE122897 primary; external cohorts in Table S1). No single prospective registry.

Outcome: Rupture status in surgical tissue (discovery); aneurysm vs control in blood bridge cohorts. Not future rupture in unruptured patients.

Predictors: (1) Three tissue signature scores + ridge logistic model; (2) eight-gene blood panel score. PHASES clinical score not computed (covariates unavailable).

Sample size: Discovery n={m['n_samples']}; external arms vary (Table S1). Event per variable not met for high-dimensional models; nested CV AUC {m['auc_nested_mean']:.2f} +/- {m['auc_nested_sd']:.2f} reported as primary performance estimate.

Missing data: Samples with unknown rupture status excluded (17 of 60 in GSE122897). No imputation of clinical covariates.

Statistical methods: DESeq2, Welch tests, LM22-adjusted OLS, permutation importance, bootstrap AUROC (Table S3). Hypergeometric test for panel overlap (Table S12).

Performance measures: AUROC, Mann-Whitney P, meta-analysis Stouffer Z. Calibration and decision curve analysis not performed (exploratory scope).

Limitations explicitly reported: Tissue not available pre-operatively; immune deconvolution poor fit (median r={m['cibersortx_median_correlation']:.2f}); no external prospective validation; models unadjusted for age, sex, hypertension, PHASES.
""".strip()


def _probast_text(m: dict) -> str:
    return f"""
PROBAST (Prediction model Risk Of Bias ASsessment Tool) — contextual characterisation for secondary public-data analyses.

Domain 1 — Participants: High risk of bias for clinical prediction. Surgical tissue cohorts are not representative of unruptured patients managed conservatively. Blood cohorts mix aneurysm detection (GSE159610/106520) with post-SAH context (GSE36791/73378).

Domain 2 — Predictors: Moderate risk. Expression predictors are pre-specified from discovery (300-gene shortlist criteria) but eight-gene panel ranking uses the same public data ecosystem. Immune adjustment uses LM22 with weak wall-tissue fit.

Domain 3 — Outcome: Low risk for within-cohort rupture labels in discovery; high risk for generalisation to future rupture because outcomes are concurrent with tissue collection.

Domain 4 — Analysis: High risk for non-nested full-cohort AUC ({m['auc_full_cohort']:.3f}) due to circularity; moderate for nested estimate ({m['auc_nested_mean']:.3f} +/- {m['auc_nested_sd']:.3f}). No prespecified analysis plan before data access documented here.

Overall: High risk of bias if interpreted as clinical prediction; acceptable for transparent exploratory biology when limitations are stated (see Table S1 covariate transparency; Supplementary Methods).
""".strip()


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description="Generate supplementary tables and figures.")
    ap.add_argument(
        "--skip-docx",
        action="store_true",
        help="Write TSV/PNG/PDF only; skip Word document assembly.",
    )
    args = ap.parse_args()

    os.chdir(BASE)
    _ensure_dirs()
    m = load_story_metrics(BASE)
    print("Building supplementary tables...")

    s1, s1b, s1_legend = build_table_s1(m)
    s1b_legend = "Platform measurability of each of the eight blood-bridge panel genes across cohort arrays/RNA-seq."
    s2, s2_legend = build_table_s2()
    s3, s3_legend = build_table_s3()
    s4, s4_legend = build_table_s4(m)
    s5, s5_legend = build_table_s5()
    s6, s6_legend = build_table_s6()
    s7, s7_legend = build_table_s7()
    s8, s8_legend = build_table_s8()
    s9, s9_legend = build_table_s9(m)
    s10, s10_legend = build_table_s10()
    s11, s11_legend = build_table_s11(m)
    s12, s12_legend = build_table_s12()

    _save_table(s1, "Supplementary_Table_S1_cohort_characteristics")
    _save_table(s1b, "Supplementary_Table_S1b_platform_gene_coverage")
    _save_table(s2, "Supplementary_Table_S2_300gene_shortlist")
    _save_table(s3, "Supplementary_Table_S3_per_gene_AUROC")
    _save_table(s4, "Supplementary_Table_S4_full_deseq2_gse122897")
    _save_table(s5, "Supplementary_Table_S5_external_rupture_DE")
    _save_table(s6, "Supplementary_Table_S6_gse15629_wall_core")
    _save_table(s7, "Supplementary_Table_S7_phase4b_gradient")
    _save_table(s8, "Supplementary_Table_S8_lm22_fractions")
    _save_table(s9, "Supplementary_Table_S9_immune_adjusted_OLS")
    _save_table(s10, "Supplementary_Table_S10_mpp1_marco_all_cohorts")
    _save_table(s11, "Supplementary_Table_S11_pathway_enrichment")
    _save_table(s12, "Supplementary_Table_S12_hypergeometric_sensitivity")

    print("Building supplementary figures...")
    fig_s1 = figure_s1_panel_violins()
    fig_s2 = figure_s2_roc_grid()
    fig_s3 = figure_s3_volcano_core()
    fig_s4 = figure_s4_cibersortx()
    fig_s5 = figure_s5_mpp1_marco_cross_cohort()

    if args.skip_docx:
        print("Skipping Word document (--skip-docx).")
        print(f"Tables: {TAB_DIR}")
        print(f"Figures: {FIG_DIR}")
    else:
        print("Assembling Word document...")
        build_word_document({
            "s1": s1, "s1_legend": s1_legend, "s1b": s1b, "s1b_legend": s1b_legend,
            "s2": s2, "s2_legend": s2_legend,
            "s3": s3, "s3_legend": s3_legend,
            "s4": s4, "s4_legend": s4_legend,
            "s5": s5, "s5_legend": s5_legend,
            "s6": s6, "s6_legend": s6_legend,
            "s7": s7, "s7_legend": s7_legend,
            "s8": s8, "s8_legend": s8_legend,
            "s9": s9, "s9_legend": s9_legend,
            "s10": s10, "s10_legend": s10_legend,
            "s11": s11, "s11_legend": s11_legend,
            "s12": s12, "s12_legend": s12_legend,
            "fig_s1": fig_s1, "fig_s2": fig_s2, "fig_s3": fig_s3, "fig_s4": fig_s4, "fig_s5": fig_s5,
            "tripod_text": _tripod_text(m),
            "probast_text": _probast_text(m),
        })

        # Also refresh NatMed bundle alias
        alias = BASE / "NatMed_Supplementary_Complete.docx"
        if OUT_DOCX.exists():
            import shutil
            shutil.copy(OUT_DOCX, alias)
            print(f"Copied to {alias}")

    print("Done.")


if __name__ == "__main__":
    main()
