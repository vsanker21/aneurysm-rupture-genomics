"""
Empirical analyses for PHASES-G items 5–8 using available GSE122897 data.

Data reality (documented in outputs):
  - No per-patient germline genotypes in this workspace → cannot compute true SNP-based PRS.
  - GEO metadata for GSE122897 provides Group (ruptured/unruptured) and Sex only — not full PHASES
    (Population, Hypertension, Age, Size, Site, Earlier SAH).

What we compute instead (clearly labeled):
  5) Candidate SNP panel + eQTL weights from the unified shortlist, and an
     **expression-based eQTL-weighted score (EXPRS)** using log2(TPM+1) as a
     blood-tissue molecular proxy (NOT allele dosages).
  6) Logistic models: Sex vs Sex+EXPRS (+ optional Tier-1 science EXPRS); report AUC (CV),
     calibration-style split; define intended full PHASES+genotype architecture as missing data.
  7) Event counts, heuristic sample-size commentary for prospective validation.
  8) Consortium positioning block (structured JSON for grants/manuscript).

Inputs:
  - data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
  - results/gse122897_sample_metadata_supplementary.tsv
  - results/phases_g/shortlist_evidence_with_mouse_gse193533_stopgap.tsv (or --shortlist)

Outputs (results/phases_g/):
  - eqtl_snp_candidate_panel.tsv
  - exprs_gene_weights.tsv
  - exprs_sample_scores.tsv
  - phases_g_logistic_comparison.json
  - prospective_validation_scoping.json
  - consortium_positioning.json
  - PHASES_G_PRACTICAL_ANALYSES_SUMMARY.json
"""

from __future__ import annotations

import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

BASE = Path(r"e:\Aneurysms Genomics Project")
TPM_DEFAULT = BASE / "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
META_DEFAULT = BASE / "results/gse122897_sample_metadata_supplementary.tsv"
SHORTLIST_DEFAULT = BASE / "results/phases_g/shortlist_evidence_with_mouse_gse193533_stopgap.tsv"
TIER1_SCIENCE_DEFAULT = BASE / "results/phases_g/tier1_science_supplement.tsv"
OUT_DIR = BASE / "results/phases_g"


def _safe_float(x) -> float | None:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return None
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def load_tpm_samples(tpm_path: Path, samples: list[str]) -> pd.DataFrame:
    """Rows = Entrez GeneID, columns = GSM sample IDs (TPM)."""
    df = pd.read_csv(tpm_path, sep="\t", low_memory=False)
    if "GeneID" not in df.columns:
        raise ValueError("TPM matrix missing GeneID column")
    df = df.rename(columns={"GeneID": "Gene_ID"})
    df["Gene_ID"] = df["Gene_ID"].astype(str).str.replace(".0", "", regex=False)
    miss = [s for s in samples if s not in df.columns]
    if miss:
        raise ValueError(f"Samples missing from TPM: {miss[:5]}...")
    return df[["Gene_ID"] + samples].copy()


def build_expr_score(
    tpm_sub: pd.DataFrame,
    weights: pd.DataFrame,
    sample_cols: list[str],
    gene_col: str = "Gene_ID",
    wcol: str = "exprs_weight",
) -> pd.Series:
    """
    EXPRS_i = sum_g w_g * z_gi where z_gi is z-score of log2(TPM+1) across samples for gene g.
    """
    gids = weights[gene_col].astype(str).tolist()
    mat = tpm_sub[tpm_sub[gene_col].astype(str).isin(gids)].set_index(gene_col)
    mat = mat.reindex(gids)
    mat = mat[sample_cols].astype(float)
    X = np.log2(mat.to_numpy() + 1.0)
    mu = X.mean(axis=1, keepdims=True)
    sd = X.std(axis=1, keepdims=True)
    sd = np.where(sd < 1e-9, 1.0, sd)
    Z = (X - mu) / sd
    w = weights.set_index(gene_col)[wcol].astype(float).reindex(gids).fillna(0.0).to_numpy()
    scores = (w[:, None] * Z).sum(axis=0)
    return pd.Series(scores, index=sample_cols)


def cv_auc(X: np.ndarray, y: np.ndarray, n_splits: int = 5, seed: int = 42) -> dict:
    if len(np.unique(y)) < 2:
        return {"auc_mean": float("nan"), "auc_std": float("nan"), "folds": []}
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    aucs = []
    for tr, te in skf.split(X, y):
        pipe = Pipeline(
            [
                ("sc", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(
                        max_iter=2000,
                        class_weight="balanced",
                        solver="lbfgs",
                    ),
                ),
            ]
        )
        pipe.fit(X[tr], y[tr])
        pr = pipe.predict_proba(X[te])[:, 1]
        aucs.append(float(roc_auc_score(y[te], pr)))
    return {"auc_mean": float(np.mean(aucs)), "auc_std": float(np.std(aucs)), "folds": aucs}


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--tpm", type=Path, default=TPM_DEFAULT)
    ap.add_argument("--meta", type=Path, default=META_DEFAULT)
    ap.add_argument("--shortlist", type=Path, default=SHORTLIST_DEFAULT)
    ap.add_argument("--tier1-science", type=Path, default=TIER1_SCIENCE_DEFAULT)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    warnings.filterwarnings(
        "ignore",
        category=pd.errors.ParserWarning,
    )

    if not args.tpm.is_file():
        raise FileNotFoundError(f"TPM matrix not found: {args.tpm}")
    if not args.meta.is_file():
        raise FileNotFoundError(f"Metadata not found: {args.meta}")
    if not args.shortlist.is_file():
        raise FileNotFoundError(f"Shortlist not found: {args.shortlist}")

    meta = pd.read_csv(args.meta, sep="\t")
    meta["Group"] = meta["Group"].astype(str).str.strip()
    meta = meta[meta["Group"].isin({"Ruptured", "Unruptured"})].copy()
    samples = meta["Sample"].astype(str).tolist()
    meta["female"] = (meta["Sex"].astype(str).str.lower().eq("female")).astype(float)
    meta["unknown_sex"] = meta["Sex"].astype(str).str.lower().isin(["unknown", "nan", ""]).astype(float)
    y = meta["Group"].eq("Ruptured").astype(int).to_numpy()

    sl = pd.read_csv(args.shortlist, sep="\t")
    sl["Gene_ID"] = sl["Gene_ID"].astype(str).str.replace(".0", "", regex=False)

    # --- (5) SNP candidate panel (no genotypes) + gene weights for EXPRS ---
    cand = sl[
        sl["eQTLGen_top_SNP"].astype(str).str.startswith("rs", na=False)
        & sl["Gene_ID"].notna()
    ].copy()

    base_cols = [
        "Gene_ID",
        "symbol",
        "eQTLGen_top_SNP",
        "eQTLGen_top_Pvalue",
        "eQTLGen_top_Zscore",
        "DESeq2_log2FC",
        "LM_immune_adj_log2_coef_Ruptured",
    ]
    extra = [c for c in ("tier1_rank",) if c in cand.columns]
    snp_panel = cand[[c for c in base_cols + extra if c in cand.columns]].drop_duplicates(
        subset=["Gene_ID"]
    )

    rows_w = []
    for _, r in snp_panel.iterrows():
        gid = str(r["Gene_ID"])
        z = _safe_float(r.get("eQTLGen_top_Zscore"))
        if z is None:
            p = _safe_float(r.get("eQTLGen_top_Pvalue"))
            z = -math.log10(max(p, 1e-300)) if p is not None and p > 0 else 1.0
            z = min(10.0, z)
        eff = abs(float(z))
        eff = min(eff, 10.0)
        lm = _safe_float(r.get("LM_immune_adj_log2_coef_Ruptured"))
        d2 = _safe_float(r.get("DESeq2_log2FC"))
        sign = 0.0
        if lm is not None and lm != 0:
            sign = 1.0 if lm > 0 else -1.0
        elif d2 is not None and d2 != 0:
            sign = 1.0 if d2 > 0 else -1.0
        if sign == 0:
            continue
        w_expr = sign * eff
        rows_w.append({"Gene_ID": gid, "symbol": r.get("symbol"), "exprs_weight": w_expr, "eqtl_abs_z": eff, "direction_sign": sign})

    weights = pd.DataFrame(rows_w)
    snp_panel_path = args.out_dir / "eqtl_snp_candidate_panel.tsv"
    cand.to_csv(snp_panel_path, sep="\t", index=False)
    weights_path = args.out_dir / "exprs_gene_weights.tsv"
    weights.to_csv(weights_path, sep="\t", index=False)

    tpm_wide = load_tpm_samples(args.tpm, samples)
    tpm_weights = tpm_wide.merge(weights[["Gene_ID", "exprs_weight"]], on="Gene_ID", how="inner")
    if tpm_weights.empty:
        raise RuntimeError("No overlapping Gene_ID between TPM and weighted shortlist.")

    exprs_full = build_expr_score(
        tpm_weights,
        tpm_weights[["Gene_ID", "exprs_weight"]],
        samples,
    )

    tier1_sci = None
    if args.tier1_science.is_file():
        t1 = pd.read_csv(args.tier1_science, sep="\t")
        t1["Gene_ID"] = t1["Gene_ID"].astype(str).str.replace(".0", "", regex=False)
        tw = weights[weights["Gene_ID"].isin(set(t1["Gene_ID"].tolist()))]
        if len(tw) >= 3:
            tpm_t1 = tpm_wide.merge(tw[["Gene_ID", "exprs_weight"]], on="Gene_ID", how="inner")
            tier1_sci = build_expr_score(
                tpm_t1,
                tpm_t1[["Gene_ID", "exprs_weight"]],
                samples,
            )

    out_samples = meta[["Sample", "Group", "Sex"]].copy()
    out_samples["ruptured"] = y
    out_samples["sex_female"] = meta["female"].to_numpy()
    out_samples["sex_unknown"] = meta["unknown_sex"].to_numpy()
    ef_map = out_samples["Sample"].map(exprs_full)
    out_samples["EXPRS_eqtl_weighted_raw"] = ef_map
    efz = (ef_map - ef_map.mean()) / (ef_map.std() if ef_map.std() > 1e-9 else 1.0)
    out_samples["EXPRS_eqtl_weighted_z"] = efz
    if tier1_sci is not None:
        out_samples["EXPRS_tier1_science_raw"] = out_samples["Sample"].map(tier1_sci)
        tmap = out_samples["EXPRS_tier1_science_raw"]
        out_samples["EXPRS_tier1_science_z"] = (tmap - tmap.mean()) / (
            tmap.std() if tmap.std() > 1e-9 else 1.0
        )

    scores_path = args.out_dir / "exprs_sample_scores.tsv"
    out_samples.to_csv(scores_path, sep="\t", index=False)

    # --- (6) Logistic comparisons (minimal clinical covariates) ---
    X_sex = meta[["female"]].to_numpy(dtype=float)

    ef = out_samples["EXPRS_eqtl_weighted_z"].to_numpy(dtype=float)
    valid = np.isfinite(ef)
    X_combo = np.column_stack([X_sex[valid], ef[valid]])
    y_v = y[valid]

    res_models = {}

    def run(name: str, X: np.ndarray, y_: np.ndarray):
        cv = cv_auc(X.astype(float), y_.astype(int), n_splits=5)
        res_models[name] = {
            **cv,
            "n": int(len(y_)),
            "n_events_ruptured": int(y_.sum()),
        }

    run("Sex_only_female_indicator", X_sex, y)
    run("EXPRS_shortlist_eqtl_weighted_only", ef[valid].reshape(-1, 1), y_v)
    run("Sex_plus_EXPRS_shortlist", X_combo, y_v)

    if tier1_sci is not None:
        e1 = out_samples["EXPRS_tier1_science_z"].to_numpy(dtype=float)
        v1 = valid & np.isfinite(e1)
        X_t1 = np.column_stack([X_sex[v1], e1[v1]])
        run("Sex_plus_EXPRS_tier1_science_subset", X_t1, y[v1])

    comparison_path = args.out_dir / "phases_g_logistic_comparison.json"
    comparison_path.write_text(json.dumps(res_models, indent=2), encoding="utf-8")

    # --- (7) Prospective scoping (illustrative; not formal power calculation) ---
    n = int(len(y))
    n_r = int(y.sum())
    n_u = int(n - n_r)
    # Simple binomial SE for AUC ~ 0.5 baseline (very rough); focus is grant sizing not exact power
    validation = {
        "current_gse122897_discovery_subset": {"n_samples": n, "n_ruptured": n_r, "n_unruptured": n_u},
        "critical_limitation": (
            "True PHASES (+ radiographic/size/site/SAH hx) and germline SNP dosages absent; "
            "models here use Sex + transcriptome-derived EXPRS only."
        ),
        "illustrative_targets_for_grant_planning_only": {
            "independent_guias_ucas_style_uia_bank": (
                "Target ≥800–2500 untreated UIA patients with ≥24–60 month imaging follow-up; "
                "genome-wide genotyping or imputation for PRS; adjudicated rupture or treatment censoring."
            ),
            "rupture_or_high_risk_event_target_for_genetic_add_on": (
                "Order-of-magnitude ~200–500 pooled rupture events across sites for stable interaction "
                "with PHASES strata (rule-of-thumb; formal power via simulation with expected effect size)."
            ),
            "registry_nesting_examples_not_commitments": [
                "ISGC: harmonized genotype + IA phenotype working groups",
                "UCAS / prospective UIA imaging cohorts: phenotype density; requires DNA consent",
                "ISUIA legacy: outcome definitions; typically requires new biospecimen collection",
            ],
        },
    }
    val_path = args.out_dir / "prospective_validation_scoping.json"
    val_path.write_text(json.dumps(validation, indent=2), encoding="utf-8")

    # --- (8) Consortium positioning ---
    cons = {
        "positioning_statement": (
            "PHASES-G is a reproducible evidence-integration framework and candidate variable set; "
            "it is not a validated clinical calculator until multisite replication, calibration, and utility trials."
        ),
        "proposed_governance_sketch": [
            "Steering committee (clinical + genomics + biostatistics + ethics)",
            "Harmonization WG: PHASES phenotype dictionary, imaging QC, outcome ontology",
            "Genomics WG: SNP list versioning, QC, ancestry, federated options",
            "Clinical translation board: thresholds, decision-curve targets, regulatory path",
        ],
        "deliverables_for_community": [
            "Frozen SNP list + allele alignment spec (when genotypes available)",
            "Locked SAP for incremental validity vs PHASES alone",
            "Open analysis code and result templates (this repository direction)",
        ],
    }
    cons_path = args.out_dir / "consortium_positioning.json"
    cons_path.write_text(json.dumps(cons, indent=2), encoding="utf-8")

    master = {
        "inputs": {
            "tpm": str(args.tpm),
            "meta": str(args.meta),
            "shortlist": str(args.shortlist),
        },
        "severe_limitations_and_circular_logic_warning": (
            "SHORTLIST GENES ARE SELECTED USING RUPTURE ASSOCIATION IN GSE122897 (plus immune adjustment). "
            "Therefore EXPRS is built from outcomes-informed features — cross-validated AUC is expected to "
            "be vastly optimistic versus any independent patient set. Interpret as a PIPELINE DEMO / "
            "score construction illustration, NOT external validity. TRUE germline PRS requires SNP profiles "
            "(absent here). TRUE PHASES-G requires PHASES clinical fields + genotype (both largely absent)."
        ),
        "n_genes_in_expr_weights": int(len(weights)),
        "n_genes_overlapping_tpm": int(tpm_weights["Gene_ID"].nunique()),
        "outputs": {
            "eqtl_snp_candidate_panel": str(snp_panel_path),
            "exprs_gene_weights": str(weights_path),
            "exprs_sample_scores": str(scores_path),
            "logistic_comparison": str(comparison_path),
            "prospective_scoping": str(val_path),
            "consortium_positioning": str(cons_path),
        },
        "interpretation": (
            "EXPRS uses eQTL effect magnitudes as gene weights times z-scored log-expression — "
            "it is a molecular summary aligned to the shortlist, not a germline PRS."
        ),
    }
    (args.out_dir / "PHASES_G_PRACTICAL_ANALYSES_SUMMARY.json").write_text(
        json.dumps(master, indent=2), encoding="utf-8"
    )

    print(json.dumps(master, indent=2))


if __name__ == "__main__":
    main()
