"""
Define a science-filtered tier-1 gene set from results/phases_g/shortlist_evidence.tsv.

Goal: produce a small, publication-facing supplement (top 20–50) where rows are
explicitly supported by the key evidence columns the manuscript describes:
  - immune-adjusted rupture signal (LM coefficient + P)
  - blood cis-eQTL sentinel rsID + P (prefer eQTLGen when present; else Catalogue)
  - IA GWAS P at that sentinel SNP (hypothesis layer only; not colocalization; prefer GWAS-at-eQTLGen-SNP)
  - GSE159610 direction and P (aneurysm vs control; exploratory)
  - marker_celltype_hint (EC/SMC vs immune vs unclassified)

Outputs:
  results/phases_g/tier1_science_supplement.tsv
  results/phases_g/tier1_science_note.json
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
DEFAULT_SHORTLIST = BASE / "results/phases_g/shortlist_evidence.tsv"
DEFAULT_OUT_TSV = BASE / "results/phases_g/tier1_science_supplement.tsv"
DEFAULT_OUT_NOTE = BASE / "results/phases_g/tier1_science_note.json"


EXPORT_COLS = [
    "tier1_science_rank",
    "tier1_science_score",
    "tier1_science_rule",
    "symbol",
    "ensembl_gene",
    "Gene_ID",
    "LM_immune_adj_log2_coef_Ruptured",
    "LM_immune_adj_P_Ruptured",
    "DESeq2_log2FC",
    "DESeq2_padj",
    # Blood genetics (prefer eQTLGen columns when present)
    "eQTLGen_top_SNP",
    "eQTLGen_top_Pvalue",
    "blood_eqtl_top_rsid",
    "blood_eqtl_top_pvalue",
    "blood_eqtl_beta",
    "blood_eqtl_dataset_id",
    # GWAS at sentinel SNP (prefer eQTLGen sentinel when present)
    "IA_GWAS_P_at_eQTLGen_SNP",
    "IA_GWAS_BETA_at_eQTLGen_SNP",
    "IA_GWAS_P_at_eqtl_snp",
    "IA_GWAS_BETA_at_eqtl_snp",
    "GSE159610_Log2FC_aneurysm_vs_control",
    "GSE159610_P",
    "direction_match_rupture_and_aneurysm",
    "marker_celltype_hint",
    "is_known_lineage_marker",
]


def _as_num(df: pd.DataFrame, col: str) -> None:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")


def _marker_bucket(hint: str) -> str:
    h = (hint or "").strip().upper()
    if h in ("EC", "SMC"):
        return "mural_prior"
    if h == "IMMUNE":
        return "immune_prior"
    return "unclassified_prior"


def _score_rows(df: pd.DataFrame) -> pd.Series:
    # Evidence-strength components (bounded, transparent).
    # - rupture signal: stronger weight to immune-adjusted P + effect size
    p_lm = df["LM_immune_adj_P_Ruptured"].clip(lower=1e-300)
    s_lm = np.clip(-np.log10(p_lm) / 6.0, 0, 1)  # ~1 by P=1e-6
    b_lm = np.clip(np.abs(df["LM_immune_adj_log2_coef_Ruptured"]) / 1.25, 0, 1)

    # - blood eQTL: prefer eQTLGen when present; require rsID, then scale -log10(P)
    if "eQTLGen_top_Pvalue" in df.columns:
        p_eq = pd.to_numeric(df["eQTLGen_top_Pvalue"], errors="coerce").clip(lower=1e-300)
        has_rsid = df["eQTLGen_top_SNP"].astype(str).str.startswith("rs", na=False)
    else:
        p_eq = df["blood_eqtl_top_pvalue"].clip(lower=1e-300)
        has_rsid = df["blood_eqtl_top_rsid"].astype(str).str.startswith("rs", na=False)
    s_eq = np.clip(-np.log10(p_eq) / 12.0, 0, 1)
    s_eq = s_eq * (0.6 + 0.4 * has_rsid.astype(float))

    # - IA GWAS at sentinel SNP: hypothesis layer only (prefer GWAS-at-eQTLGen-SNP)
    if "IA_GWAS_P_at_eQTLGen_SNP" in df.columns:
        p_gwas = df["IA_GWAS_P_at_eQTLGen_SNP"]
    else:
        p_gwas = df["IA_GWAS_P_at_eqtl_snp"]
    has_gwas = p_gwas.notna() & np.isfinite(p_gwas) & (p_gwas > 0)
    s_gwas = np.where(
        has_gwas,
        np.clip(-np.log10(p_gwas.clip(lower=1e-300)) / 5.0, 0, 1),
        0.0,
    )
    s_gwas = np.where(has_gwas & (p_gwas < 0.05), np.maximum(s_gwas, 0.85), s_gwas)

    # - replication (exploratory): direction match + nominal P
    dm = df["direction_match_rupture_and_aneurysm"]
    s_rep = np.zeros(len(df))
    s_rep = np.where(dm == True, s_rep + 0.6, s_rep)  # noqa: E712
    s_rep = np.where(df["GSE159610_P"].notna() & (df["GSE159610_P"] < 0.05), s_rep + 0.4, s_rep)
    s_rep = np.clip(s_rep, 0, 1)

    # - marker prior bucket
    bucket = df["marker_celltype_hint"].astype(str).map(_marker_bucket)
    s_prior = bucket.map({"mural_prior": 1.0, "unclassified_prior": 0.7, "immune_prior": 0.35}).astype(float)

    score = (
        0.26 * pd.Series(s_lm, index=df.index)
        + 0.14 * pd.Series(b_lm, index=df.index)
        + 0.22 * pd.Series(s_eq, index=df.index)
        + 0.16 * pd.Series(s_gwas, index=df.index)
        + 0.12 * pd.Series(s_rep, index=df.index)
        + 0.10 * pd.Series(s_prior, index=df.index)
    )
    return score


def build_tier1_science(
    shortlist_tsv: Path | None = None,
    out_tsv: Path | None = None,
    out_note: Path | None = None,
    top_n: int = 50,
    require_blood_eqtl_rsid: bool = True,
    require_any_of: tuple[str, ...] = ("gwas_at_eqtl_snp", "replication_direction", "mural_marker_prior"),
) -> pd.DataFrame:
    path = Path(shortlist_tsv or DEFAULT_SHORTLIST)
    if not path.is_file():
        raise FileNotFoundError(f"Shortlist not found: {path}")

    df = pd.read_csv(path, sep="\t")
    for c in (
        "LM_immune_adj_log2_coef_Ruptured",
        "LM_immune_adj_P_Ruptured",
        "DESeq2_log2FC",
        "DESeq2_padj",
        "blood_eqtl_top_pvalue",
        "blood_eqtl_beta",
        "IA_GWAS_P_at_eqtl_snp",
        "IA_GWAS_BETA_at_eqtl_snp",
        "eQTLGen_top_Pvalue",
        "IA_GWAS_P_at_eQTLGen_SNP",
        "IA_GWAS_BETA_at_eQTLGen_SNP",
        "GSE159610_Log2FC_aneurysm_vs_control",
        "GSE159610_P",
    ):
        _as_num(df, c)

    if "eQTLGen_top_Pvalue" in df.columns:
        has_eqtl_rsid = df["eQTLGen_top_SNP"].astype(str).str.startswith("rs", na=False)
        has_eqtl_p = df["eQTLGen_top_Pvalue"].notna()
        gwas_p_col = "IA_GWAS_P_at_eQTLGen_SNP" if "IA_GWAS_P_at_eQTLGen_SNP" in df.columns else "IA_GWAS_P_at_eqtl_snp"
    else:
        has_eqtl_rsid = df["blood_eqtl_top_rsid"].astype(str).str.startswith("rs", na=False)
        has_eqtl_p = df["blood_eqtl_top_pvalue"].notna()
        gwas_p_col = "IA_GWAS_P_at_eqtl_snp"
    eqtl_ok = has_eqtl_p & has_eqtl_rsid if require_blood_eqtl_rsid else has_eqtl_p

    has_gwas = df[gwas_p_col].notna() & np.isfinite(df[gwas_p_col])
    has_rep_dir = df["direction_match_rupture_and_aneurysm"] == True  # noqa: E712
    has_mural_prior = df["marker_celltype_hint"].astype(str).str.upper().isin(["EC", "SMC"])

    anyof_flags = {
        "gwas_at_eqtl_snp": has_gwas,
        "replication_direction": has_rep_dir,
        "mural_marker_prior": has_mural_prior,
    }
    if not require_any_of:
        anyof_ok = pd.Series(True, index=df.index)
        anyof_rule = "eqtl_only"
    else:
        sel = [anyof_flags[k] for k in require_any_of if k in anyof_flags]
        anyof_ok = sel[0].copy()
        for s in sel[1:]:
            anyof_ok = anyof_ok | s
        anyof_rule = "any_of:" + ",".join(require_any_of)

    df["_tier1_science_rule"] = np.where(eqtl_ok & anyof_ok, anyof_rule, "excluded")
    keep = df[df["_tier1_science_rule"] != "excluded"].copy()
    if keep.empty:
        raise RuntimeError("Tier-1 science filter selected 0 genes; relax criteria.")

    keep["_tier1_science_score"] = _score_rows(keep)
    keep = keep.sort_values("_tier1_science_score", ascending=False).reset_index(drop=True)
    keep["tier1_science_rank"] = np.arange(1, len(keep) + 1)
    keep["tier1_science_score"] = keep["_tier1_science_score"].round(4)
    keep["tier1_science_rule"] = keep["_tier1_science_rule"]
    keep = keep.drop(columns=["_tier1_science_score", "_tier1_science_rule"])

    top = keep.head(int(top_n)).copy()

    out_path = Path(out_tsv or DEFAULT_OUT_TSV)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    top[EXPORT_COLS].to_csv(out_path, sep="\t", index=False)

    note = {
        "definition": {
            "require_blood_eqtl_rsid": bool(require_blood_eqtl_rsid),
            "require_any_of": list(require_any_of),
            "interpretation": "Tier-1 science filter = shortlist rows with blood cis-eQTL rsID (eQTLGen preferred when present; else Catalogue) + at least one additional support (GWAS-at-sentinel, replication direction, or mural marker prior).",
        },
        "weights": {
            "LM_immune_adj_neglog10_P_scaled": 0.26,
            "LM_immune_adj_abs_coef_scaled": 0.14,
            "blood_cis_eQTL_neglog10_P_scaled": 0.22,
            "IA_GWAS_at_eQTL_SNP_hypothesis_layer": 0.16,
            "GSE159610_replication_direction_and_P_exploratory": 0.12,
            "marker_celltype_prior_bucket": 0.10,
        },
        "caveats": [
            "This tier-1 science filter is for manuscript prioritization; not a validated clinical index.",
            "IA GWAS overlap remains rsID lookup at the eQTL sentinel SNP; formal genetics requires COLOC/SMR/HEIDI and allele harmonization.",
            "GSE159610 contrast is aneurysm vs control, not ruptured vs unruptured; direction match is exploratory.",
        ],
        "columns_in_tier1_science_supplement_tsv": EXPORT_COLS,
        "n_genes_in_shortlist": int(len(df)),
        "n_genes_passing_filter": int(len(keep)),
        "n_genes_exported": int(len(top)),
    }
    note_path = Path(out_note or DEFAULT_OUT_NOTE)
    note_path.write_text(json.dumps(note, indent=2), encoding="utf-8")

    return top


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description="Define tier-1 science-filtered supplement from shortlist.")
    ap.add_argument("--shortlist", type=Path, default=DEFAULT_SHORTLIST)
    ap.add_argument("--top-n", type=int, default=50)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT_TSV)
    ap.add_argument("--no-require-eqtl-rsid", action="store_true")
    ap.add_argument(
        "--require-any-of",
        type=str,
        default="gwas_at_eqtl_snp,replication_direction,mural_marker_prior",
        help="Comma list of flags; require at least one: gwas_at_eqtl_snp, replication_direction, mural_marker_prior. Empty string disables.",
    )
    args = ap.parse_args()

    req_any = tuple([s.strip() for s in args.require_any_of.split(",") if s.strip()])
    build_tier1_science(
        shortlist_tsv=args.shortlist,
        out_tsv=args.out,
        top_n=args.top_n,
        require_blood_eqtl_rsid=not args.no_require_eqtl_rsid,
        require_any_of=req_any,
    )
    print(f"Wrote {args.out}")
    print(f"Wrote {DEFAULT_OUT_NOTE}")


if __name__ == "__main__":
    main()

