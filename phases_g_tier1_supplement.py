"""
Build tier-1 prioritization from results/phases_g/shortlist_evidence.tsv.

Composite score (hypothesis prioritization, not a clinical index):
  - Tissue rupture signal: DESeq2 FDR and immune-fraction-adjusted rupture P (LM).
  - Blood cis-eQTL: Catalogue hit with rsID; stronger weight for lower eQTL P.
  - IA GWAS: rsID overlap with Bakker sumstats (hypothesis only; not colocalization).
  - Replication: GSE159610 aneurysm vs control direction vs rupture DE direction; nominal P.
  - scRNA-style prior: EC/SMC marker genes weighted higher than immune/unclassified.

Outputs:
  results/phases_g/tier1_supplement.tsv
  results/phases_g/tier1_scoring_note.json  (methods text for supplement / manuscript)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
DEFAULT_SHORTLIST = BASE / "results/phases_g/shortlist_evidence.tsv"
DEFAULT_OUT_TSV = BASE / "results/phases_g/tier1_supplement.tsv"
DEFAULT_OUT_TSV_WALL = BASE / "results/phases_g/tier1_supplement_vessel_wall.tsv"
DEFAULT_OUT_NOTE = BASE / "results/phases_g/tier1_scoring_note.json"
DEFAULT_MOUSE_EXPR = BASE / "results/phases_g/gse193533_mouse_celltype_gene_expression.tsv"

EXPORT_COLS = [
    "tier1_rank",
    "tier1_score",
    "symbol",
    "Gene_ID",
    "ensembl_gene",
    "DESeq2_log2FC",
    "DESeq2_padj",
    "LM_immune_adj_P_Ruptured",
    "blood_eqtl_top_rsid",
    "blood_eqtl_top_pvalue",
    "IA_GWAS_P_at_eqtl_snp",
    "eQTLGen_top_SNP",
    "eQTLGen_top_Pvalue",
    "IA_GWAS_P_at_eQTLGen_SNP",
    "GSE159610_Log2FC_aneurysm_vs_control",
    "GSE159610_P",
    "direction_match_rupture_and_aneurysm",
    "marker_celltype_hint",
    # Quantitative scRNA specificity (marker-table derived; human CVAE186 / CellMarker pipelines)
    "scrna_best_cell_type",
    "scrna_best_score",
    "scrna_wall_score",
    "scrna_immune_score",
    "scrna_wall_minus_immune",
    # Stopgap: mouse aneurysm-model scRNA (GSE193533) broad compartment evidence
    "mouse_gse193533_best_cell_type",
    "mouse_gse193533_best_mean_log1p",
    "mouse_gse193533_wall_minus_immune",
]


def _norm01(s: pd.Series) -> pd.Series:
    lo, hi = float(s.min()), float(s.max())
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return pd.Series(0.5, index=s.index, dtype=float)
    return ((s - lo) / (hi - lo)).clip(0, 1)


def _marker_score(hint: str) -> float:
    h = (hint or "").strip().upper()
    # Aggressive vessel-wall prioritization:
    #   - EC/SMC highest priority
    #   - pericyte + fibroblast/ECM are wall-adjacent/supporting
    #   - immune is penalized (reactive infiltration risk)
    #   - unclassified gets low prior (insufficient cell-type support)
    if h in ("EC", "SMC"):
        return 1.0
    if h in ("PERICYTE/PERIVASCULAR", "PERICYTE", "PERIVASCULAR"):
        return 0.85
    if h in ("FIBROBLAST/ECM", "FIBROBLAST", "ECM"):
        return 0.75
    if h.startswith("IMMUNE"):
        return 0.15
    if h == "UNCLASSIFIED":
        return 0.40
    return 0.50


def build_tier1_table(
    shortlist_tsv: Path | None = None,
    out_tsv: Path | None = None,
    out_tsv_wall: Path | None = None,
    out_note: Path | None = None,
    top_n: int = 40,
) -> pd.DataFrame:
    path = Path(shortlist_tsv or DEFAULT_SHORTLIST)
    if not path.is_file():
        raise FileNotFoundError(f"Shortlist not found: {path}")

    df = pd.read_csv(path, sep="\t")
    df["DESeq2_padj"] = pd.to_numeric(df["DESeq2_padj"], errors="coerce")
    df["LM_immune_adj_P_Ruptured"] = pd.to_numeric(
        df["LM_immune_adj_P_Ruptured"], errors="coerce"
    )
    df["blood_eqtl_top_pvalue"] = pd.to_numeric(
        df["blood_eqtl_top_pvalue"], errors="coerce"
    )
    if "eQTLGen_top_Pvalue" in df.columns:
        df["eQTLGen_top_Pvalue"] = pd.to_numeric(df["eQTLGen_top_Pvalue"], errors="coerce")
    df["IA_GWAS_P_at_eqtl_snp"] = pd.to_numeric(
        df["IA_GWAS_P_at_eqtl_snp"], errors="coerce"
    )
    if "IA_GWAS_P_at_eQTLGen_SNP" in df.columns:
        df["IA_GWAS_P_at_eQTLGen_SNP"] = pd.to_numeric(df["IA_GWAS_P_at_eQTLGen_SNP"], errors="coerce")
    df["GSE159610_P"] = pd.to_numeric(df["GSE159610_P"], errors="coerce")

    # Optional stopgap mouse scRNA layer (GSE193533): merge per-gene broad cell compartment evidence
    mouse_expr_path = DEFAULT_MOUSE_EXPR
    if mouse_expr_path.is_file():
        mex = pd.read_csv(mouse_expr_path, sep="\t")
        mex["mean_log1p"] = pd.to_numeric(mex["mean_log1p"], errors="coerce")

        # gene-level best compartment + best mean
        best = (
            mex.dropna(subset=["mean_log1p"])
            .sort_values(["gene_human", "mean_log1p"], ascending=[True, False])
            .drop_duplicates("gene_human")
            .rename(
                columns={
                    "gene_human": "symbol",
                    "cell_type": "mouse_gse193533_best_cell_type",
                    "mean_log1p": "mouse_gse193533_best_mean_log1p",
                }
            )[["symbol", "mouse_gse193533_best_cell_type", "mouse_gse193533_best_mean_log1p"]]
        )

        # wall minus immune (max wall - immune)
        wall = {"EC", "SMC", "FIBROBLAST"}
        rows = []
        for gene, sub in mex.groupby("gene_human"):
            sub = sub.dropna(subset=["mean_log1p"])
            if sub.empty:
                rows.append((gene, np.nan))
                continue
            immune = (
                float(sub.loc[sub["cell_type"] == "IMMUNE", "mean_log1p"].iloc[0])
                if (sub["cell_type"] == "IMMUNE").any()
                else np.nan
            )
            wall_max = (
                float(sub.loc[sub["cell_type"].isin(wall), "mean_log1p"].max())
                if (sub["cell_type"].isin(wall)).any()
                else np.nan
            )
            rows.append((gene, wall_max - immune))
        wmi = pd.DataFrame(rows, columns=["symbol", "mouse_gse193533_wall_minus_immune"])

        df = df.merge(best, on="symbol", how="left").merge(wmi, on="symbol", how="left")

    nld = -np.log10(df["DESeq2_padj"].clip(lower=1e-300))
    nll = -np.log10(df["LM_immune_adj_P_Ruptured"].clip(lower=1e-300))

    # Blood eQTL strength: prefer eQTLGen when present; otherwise fall back to Catalogue.
    if "eQTLGen_top_Pvalue" in df.columns:
        eqtl_p = df["eQTLGen_top_Pvalue"]
        has_rsid = df.get("eQTLGen_top_SNP", pd.Series(False, index=df.index)).astype(str).str.startswith(
            "rs", na=False
        )
        eqtl_strength_note = "eQTLGen_cis_FDR005_top_hit"
    else:
        eqtl_p = df["blood_eqtl_top_pvalue"]
        has_rsid = df["blood_eqtl_top_rsid"].notna() & (
            df["blood_eqtl_top_rsid"].astype(str).str.startswith("rs", na=False)
        )
        eqtl_strength_note = "eQTL_Catalogue_blood_top_hit"
    has_eqtl_p = pd.to_numeric(eqtl_p, errors="coerce").notna()
    nleq = np.where(has_eqtl_p, -np.log10(pd.to_numeric(eqtl_p, errors="coerce").clip(lower=1e-300)), 0.0)
    eqtl_strength = np.clip(nleq / 12.0, 0, 1) * (0.5 + 0.5 * has_rsid.astype(float))

    # GWAS overlap: prefer GWAS-at-eQTLGen-SNP when available, else fall back to Catalogue rsID overlap
    if "IA_GWAS_P_at_eQTLGen_SNP" in df.columns:
        gwas_p = df["IA_GWAS_P_at_eQTLGen_SNP"]
        gwas_note = "Bakker2020_at_eQTLGen_top_SNP"
    else:
        gwas_p = df["IA_GWAS_P_at_eqtl_snp"]
        gwas_note = "Bakker2020_at_Catalogue_top_rsid"
    has_gwas = gwas_p.notna() & np.isfinite(gwas_p) & (gwas_p > 0)
    gwas_score = np.zeros(len(df))
    gwas_score = np.where(has_gwas, np.clip(-np.log10(gwas_p.clip(lower=1e-300)) / 5.0, 0, 1), gwas_score)
    gwas_score = np.where(has_gwas & (gwas_p < 0.05), np.maximum(gwas_score, 0.85), gwas_score)

    dm = df["direction_match_rupture_and_aneurysm"]
    repl = np.zeros(len(df))
    repl = np.where(dm == True, repl + 0.55, repl)  # noqa: E712
    repl = np.where(df["GSE159610_P"].notna() & (df["GSE159610_P"] < 0.05), repl + 0.45, repl)
    repl = np.clip(repl, 0, 1)

    mscore = df["marker_celltype_hint"].astype(str).map(_marker_score)

    # Optional scRNA specificity layer (quantitative, marker-table derived).
    # If present, this can dominate the vessel-wall prioritization more defensibly than
    # a strict marker lookup. Range is approximately [-1, 1] given marker_strength in [0,1].
    if "scrna_wall_minus_immune" in df.columns:
        scrna_spec = pd.to_numeric(df["scrna_wall_minus_immune"], errors="coerce").fillna(0.0)
        # Map to [0,1] with a conservative clamp. Positive favors vessel-wall.
        scrna_spec01 = ((scrna_spec + 1.0) / 2.0).clip(0, 1)
    else:
        scrna_spec01 = pd.Series(0.5, index=df.index, dtype=float)

    df["_tier1_score"] = (
        # Weight vector tuned to favor EC/SMC-origin candidates more strongly.
        0.16 * _norm01(nld)
        + 0.16 * _norm01(nll)
        + 0.20 * pd.Series(eqtl_strength, index=df.index)
        + 0.12 * pd.Series(gwas_score, index=df.index)
        + 0.10 * pd.Series(repl, index=df.index)
        + 0.18 * mscore
        + 0.08 * scrna_spec01
    )

    df = df.sort_values("_tier1_score", ascending=False).reset_index(drop=True)
    df["tier1_rank"] = np.arange(1, len(df) + 1)
    df["tier1_score"] = df["_tier1_score"].round(4)
    df = df.drop(columns=["_tier1_score"])

    top = df.head(int(top_n)).copy()
    out_path = Path(out_tsv or DEFAULT_OUT_TSV)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # Backwards compatibility: only export columns that exist in this run's shortlist.
    export_cols = [c for c in EXPORT_COLS if c in top.columns]
    top[export_cols].to_csv(out_path, sep="\t", index=False)

    # Vessel-wall–prioritized view (more defensible “vulnerability” candidates):
    # keep EC/SMC/pericyte/fibroblast/ECM hints; exclude immune/unclassified by default.
    wall_hints = {
        "EC",
        "SMC",
        "PERICYTE/PERIVASCULAR",
        "PERICYTE",
        "PERIVASCULAR",
        "FIBROBLAST/ECM",
        "FIBROBLAST",
        "ECM",
    }
    wall = df[df["marker_celltype_hint"].astype(str).str.upper().isin(wall_hints)].copy()
    wall = wall.head(int(top_n))
    out_path_wall = Path(out_tsv_wall or DEFAULT_OUT_TSV_WALL)
    out_path_wall.parent.mkdir(parents=True, exist_ok=True)
    export_cols_wall = [c for c in EXPORT_COLS if c in wall.columns]
    wall[export_cols_wall].to_csv(out_path_wall, sep="\t", index=False)

    note = {
        "weights": {
            "DESeq2_neglog10_padj_normalized": 0.16,
            "LM_immune_adj_rupture_neglog10_P_normalized": 0.16,
            "blood_cis_eQTL_strength_with_rsid": 0.20,
            "IA_GWAS_at_eQTL_SNP_hypothesis_layer": 0.12,
            "GSE159610_replication_direction_and_P": 0.10,
            "marker_celltype_prior_EC_SMC_pericyte_fibroblast_immune_unclassified": 0.18,
            "scrna_wall_minus_immune_specificity_scaled": 0.08,
        },
        "signals_used": {
            "blood_cis_eQTL_strength_source": eqtl_strength_note,
            "IA_GWAS_overlap_source": gwas_note,
        },
        "caveats": [
            "Tier-1 score is for prioritization only; not validated for clinical use.",
            "IA GWAS overlap is rsID lookup, not colocalization or SMR.",
            "GSE159610 contrast is aneurysm vs control, not ruptured vs unruptured.",
            "Blood eQTL source depends on columns present in shortlist TSV (eQTLGen preferred when provided).",
            "Marker priors use a small curated panel; extend with IA scRNA gene sets when available.",
        ],
        "columns_in_tier1_supplement_tsv": EXPORT_COLS,
        "vessel_wall_view": {
            "output_tsv": str(out_path_wall),
            "included_marker_celltype_hint": sorted(list(wall_hints)),
            "n_exported": int(len(wall)),
            "note": "This view excludes immune/unclassified hints to emphasize vessel-wall origin candidates.",
        },
        "n_genes_in_shortlist": int(len(df)),
        "n_genes_exported": int(len(top)),
    }
    note_path = Path(out_note or DEFAULT_OUT_NOTE)
    note_path.write_text(json.dumps(note, indent=2), encoding="utf-8")

    return top


def main():
    import argparse

    ap = argparse.ArgumentParser(description="Build PHASES-G tier-1 supplement table.")
    ap.add_argument("--shortlist", type=Path, default=DEFAULT_SHORTLIST)
    ap.add_argument("--top-n", type=int, default=40)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT_TSV)
    args = ap.parse_args()
    build_tier1_table(
        shortlist_tsv=args.shortlist,
        out_tsv=args.out,
        top_n=args.top_n,
    )
    print(f"Wrote {args.out}")
    print(f"Wrote {DEFAULT_OUT_NOTE}")


if __name__ == "__main__":
    main()
