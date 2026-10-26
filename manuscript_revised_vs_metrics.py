"""
Load numeric facts for create_molecular_signatures_revised_vs.py narrative text.
Falls back to documented defaults if JSON/TSV missing.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def load_story_metrics(base: Path) -> dict:
    m: dict = {
        "n_samples": 42,
        "n_ruptured": 21,
        "n_unruptured": 21,
        "deseq2_padj_genes": 3094,
        "welch_fdr_genes": 4,
        "pathways_fdr05": 207,
        "auc_full_cohort": 0.827,
        "auc_cv_mean": 0.798,
        "auc_cv_sd": 0.156,
        "auc_nested_mean": 0.64,
        "auc_nested_sd": 0.17,
        "de_feature_importance_pct": 56.5,
        "ec_feature_importance_pct": 35.6,
        "ecm_feature_importance_pct": 7.8,
        "n_shortlist": 300,
        "n_gse159610_matched": 263,
        "n_blood_eqtl_catalogue": 0,
        "n_eqtlgen_cis": 0,
        "n_gwas_overlap_catalogue": 0,
        "n_gwas_overlap_eqtlgen": 0,
        "n_genes_immune_model_tested": 14888,
        "n_genes_immune_nominal_rupture": 1350,
        "n_genes_deseq2_and_nominal_rupture_immune": 645,
        "n_exprs_genes": 205,
        "n_prs_weights_snps": 204,
        "lm22_cell_types": 22,
        "lm22_fdr_sig_cell_types_005": 0,
        "n_genes_lm22_intersect": 516,
        "n_deseq2_fdr_among_immune_tested": 2187,
        "skip_eqtl": True,
        "n_cibersortx_job7_mixtures": 42,
        "n_cibersortx_gsm_parsed": 42,
        "cibersortx_median_correlation": 0.142,
        "cibersortx_median_rmse": 1.055,
        "shortlist_deseq2_padj_max": 0.05,
        "shortlist_lm_p_max": 0.05,
        "gse54083_n_ruptured": 8,
        "gse54083_n_unruptured": 5,
        "gse54083_n_shortlist_present": 0,
        "gse54083_n_shortlist_direction_match": 0,
        "gse54083_n_shortlist_p_lt_0_05": 0,
        "gse54083_n_shortlist_fdr_lt_0_05": 0,
        "gse13353_direction_matches": 6,
        "gse13353_genes_tested": 7,
        "gse13353_nominal_sig_genes": 3,
        "gse26969_n_aneurysm": 3,
        "gse26969_n_control": 3,
        "gse26969_n_shortlist_present": 0,
        "gse26969_n_shortlist_direction_match": 0,
        "gse26969_n_shortlist_p_lt_0_05": 0,
        "gse26969_n_shortlist_fdr_lt_0_05": 0,
        "gse36791_n_ruptured": 43,
        "gse36791_n_control": 18,
        "gse36791_n_shortlist_present": 0,
        "gse36791_n_shortlist_direction_match": 0,
        "gse36791_n_shortlist_p_lt_0_05": 0,
        "gse36791_n_shortlist_fdr_lt_0_05": 0,
    }

    job7 = base / "results/cibersortx/CIBERSORTx_Job7_Results.txt"
    if not job7.is_file():
        job7 = base / "CIBERSORTx_Job7_Results.txt"
    if job7.is_file():
        try:
            raw = job7.read_text(encoding="utf-8", errors="replace").splitlines()
            m["n_cibersortx_job7_mixtures"] = max(0, len(raw) - 1) if raw else 0
        except Exception:
            pass

    cx = base / "results/cibersortx/cibersortx_parse_summary.json"
    if cx.is_file():
        try:
            c = json.loads(cx.read_text(encoding="utf-8"))
            m["n_cibersortx_gsm_parsed"] = int(c.get("n_samples_parsed", m["n_cibersortx_gsm_parsed"]))
            qc = c.get("qc_summary") or {}
            cor = (qc.get("CIBERSORTx_Correlation") or {}).get("median")
            rmse = (qc.get("CIBERSORTx_RMSE") or {}).get("median")
            if cor is not None:
                m["cibersortx_median_correlation"] = float(cor)
            if rmse is not None:
                m["cibersortx_median_rmse"] = float(rmse)
        except Exception:
            pass

    ps = base / "results/phases_g/pipeline_summary.json"
    if ps.is_file():
        try:
            p = json.loads(ps.read_text(encoding="utf-8"))
            m["n_shortlist"] = int(p.get("n_shortlist_genes", m["n_shortlist"]))
            m["n_gse159610_matched"] = int(p.get("n_replication_gse159610_matched", m["n_gse159610_matched"]))
            m["n_blood_eqtl_catalogue"] = int(p.get("n_with_blood_eqtl_hit", 0))
            m["n_eqtlgen_cis"] = int(p.get("n_with_eqtlgen_cis_hit", 0))
            m["n_gwas_overlap_catalogue"] = int(p.get("n_eqtl_snps_in_IA_GWAS_file", 0))
            m["n_gwas_overlap_eqtlgen"] = int(p.get("n_eqtlgen_snps_in_IA_GWAS_file", 0))
            m["skip_eqtl"] = bool(p.get("skip_eqtl", True))
            filt = p.get("filters") or {}
            if "DESeq2_padj_max" in filt:
                m["shortlist_deseq2_padj_max"] = float(filt["DESeq2_padj_max"])
            if "LM_immune_adj_P_Ruptured_max" in filt:
                m["shortlist_lm_p_max"] = float(filt["LM_immune_adj_P_Ruptured_max"])
        except Exception:
            pass

    g54083 = base / "results/external_validation_gse54083/gse54083_validation_summary.json"
    if g54083.is_file():
        try:
            g = json.loads(g54083.read_text(encoding="utf-8"))
            m["gse54083_n_ruptured"] = int(g.get("n_ruptured", m["gse54083_n_ruptured"]))
            m["gse54083_n_unruptured"] = int(g.get("n_unruptured", m["gse54083_n_unruptured"]))
            m["gse54083_n_shortlist_present"] = int(g.get("n_shortlist_genes_present", 0))
            m["gse54083_n_shortlist_direction_match"] = int(g.get("n_shortlist_direction_match", 0))
            m["gse54083_n_shortlist_p_lt_0_05"] = int(g.get("n_shortlist_p_lt_0_05", 0))
            m["gse54083_n_shortlist_fdr_lt_0_05"] = int(g.get("n_shortlist_fdr_lt_0_05", 0))
        except Exception:
            pass

    g13353 = base / "results/external_validation_gse13353/gse13353_validation_summary.tsv"
    if g13353.is_file():
        try:
            s = pd.read_csv(g13353, sep="\t", header=None, names=["key", "value"])
            d = dict(zip(s["key"], s["value"]))
            if "direction_matches" in d:
                m["gse13353_direction_matches"] = int(float(d["direction_matches"]))
            if "genes_tested" in d:
                m["gse13353_genes_tested"] = int(float(d["genes_tested"]))
        except Exception:
            pass

    g13353_tbl = base / "results/external_validation_gse13353/gse13353_key_gene_validation.tsv"
    if g13353_tbl.is_file():
        try:
            t = pd.read_csv(g13353_tbl, sep="\t")
            if "P_value" in t.columns:
                m["gse13353_nominal_sig_genes"] = int((pd.to_numeric(t["P_value"], errors="coerce") < 0.05).sum())
        except Exception:
            pass

    g26969 = base / "results/external_validation_gse26969/gse26969_context_summary.json"
    if g26969.is_file():
        try:
            g = json.loads(g26969.read_text(encoding="utf-8"))
            m["gse26969_n_aneurysm"] = int(g.get("n_aneurysm", m["gse26969_n_aneurysm"]))
            m["gse26969_n_control"] = int(g.get("n_control", m["gse26969_n_control"]))
            m["gse26969_n_shortlist_present"] = int(g.get("n_shortlist_genes_present", 0))
            m["gse26969_n_shortlist_direction_match"] = int(g.get("n_shortlist_direction_match", 0))
            m["gse26969_n_shortlist_p_lt_0_05"] = int(g.get("n_shortlist_p_lt_0_05", 0))
            m["gse26969_n_shortlist_fdr_lt_0_05"] = int(g.get("n_shortlist_fdr_lt_0_05", 0))
        except Exception:
            pass

    g36791 = base / "results/external_validation_gse36791/gse36791_context_summary.json"
    if g36791.is_file():
        try:
            g = json.loads(g36791.read_text(encoding="utf-8"))
            m["gse36791_n_ruptured"] = int(g.get("n_ruptured", m["gse36791_n_ruptured"]))
            m["gse36791_n_control"] = int(g.get("n_control", m["gse36791_n_control"]))
            m["gse36791_n_shortlist_present"] = int(g.get("n_shortlist_genes_present", 0))
            m["gse36791_n_shortlist_direction_match"] = int(g.get("n_shortlist_direction_match", 0))
            m["gse36791_n_shortlist_p_lt_0_05"] = int(g.get("n_shortlist_p_lt_0_05", 0))
            m["gse36791_n_shortlist_fdr_lt_0_05"] = int(g.get("n_shortlist_fdr_lt_0_05", 0))
        except Exception:
            pass

    im = base / "results/cibersortx/deg_immune_adjustment_summary.json"
    if im.is_file():
        try:
            d = json.loads(im.read_text(encoding="utf-8"))
            m["n_genes_immune_model_tested"] = int(d.get("n_genes_tested", m["n_genes_immune_model_tested"]))
            m["n_genes_immune_nominal_rupture"] = int(
                d.get("n_rupture_coef_nominal_P_lt_0.05", m["n_genes_immune_nominal_rupture"])
            )
            m["n_genes_deseq2_and_nominal_rupture_immune"] = int(
                d.get(
                    "n_DESeq2_sig_and_rupture_coef_nominal_P_lt_0.05",
                    m["n_genes_deseq2_and_nominal_rupture_immune"],
                )
            )
            m["n_deseq2_fdr_among_immune_tested"] = int(
                d.get("n_DESeq2_FDR_lt_0.05_among_tested", m["n_deseq2_fdr_among_immune_tested"])
            )
            cov = d.get("immune_covariates_in_model") or []
            m["n_immune_covariates"] = len(cov)
        except Exception:
            m["n_immune_covariates"] = 11
    else:
        m["n_immune_covariates"] = 11

    lm22 = base / "results/cibersortx/lm22_deconvolution_summary.json"
    if lm22.is_file():
        try:
            L = json.loads(lm22.read_text(encoding="utf-8"))
            m["lm22_fdr_sig_cell_types_005"] = int(L.get("fdr_significant_cell_types_005", 0))
            m["n_genes_lm22_intersect"] = int(L.get("n_genes_intersect_lm22", m["n_genes_lm22_intersect"]))
            m["lm22_cell_types"] = int(L.get("cell_types", m["lm22_cell_types"]))
        except Exception:
            pass

    pa = base / "results/phases_g/PHASES_G_PRACTICAL_ANALYSES_SUMMARY.json"
    if pa.is_file():
        try:
            a = json.loads(pa.read_text(encoding="utf-8"))
            m["n_exprs_genes"] = int(a.get("n_genes_in_expr_weights", m["n_exprs_genes"]))
        except Exception:
            pass

    nested = base / "results/gse122897_nested_de_cv_summary.tsv"
    if nested.is_file():
        try:
            nd = pd.read_csv(nested, sep="\t")
            d = dict(zip(nd["Metric"], nd["Value"]))
            m["auc_nested_mean"] = float(d.get("Nested_CV_AUC_Mean", m["auc_nested_mean"]))
            m["auc_nested_sd"] = float(d.get("Nested_CV_AUC_Std", m["auc_nested_sd"]))
        except Exception:
            pass

    best = base / "results/phase5_best_model_summary.json"
    if best.is_file():
        try:
            b = json.loads(best.read_text(encoding="utf-8"))
            m["auc_cv_mean"] = float(b.get("cv_auc_mean", m["auc_cv_mean"]))
            m["auc_cv_sd"] = float(b.get("cv_auc_std", m["auc_cv_sd"]))
            m["auc_full_cohort"] = float(b.get("train_auc", m["auc_full_cohort"]))
        except Exception:
            pass

    prs_w = base / "results/phases_g/prs_rsid_weights_from_panel.tsv"
    if prs_w.is_file():
        try:
            lines = prs_w.read_text(encoding="utf-8").splitlines()
            m["n_prs_weights_snps"] = max(0, len(lines) - 1) if lines else 0
        except Exception:
            pass

    return m
