"""
Join an external ruptured-vs-unruptured (or other) tissue DE gene list to the
PHASES-G shortlist (GSE122897 + blood eQTL + Bakker lookup + GSE159610 + markers).

Outputs a bridge table for validation / manuscript: which narrative genes are in
the integrative shortlist, and whether optional tissue log2FC signs agree with
GSE122897 DESeq2 and immune-adjusted LM coefficients.

Usage:
  python link_tissue_de_to_phases_g.py --example
  python link_tissue_de_to_phases_g.py --de-tsv my_cohort.tsv --out results/phases_g/tissue_signature_phases_g_bridge.tsv

Also writes <stem>_summary.json next to the bridge TSV (counts, concordance, tier-1 hits).

Input TSV columns (header required):
  symbol  — HGNC
  signature_group — optional free text (e.g. upregulated_in_ruptured)
  tissue_log2FC, tissue_padj — optional; used only for concordance flags vs GSE122897
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
SHORTLIST = BASE / "results/phases_g/shortlist_evidence.tsv"
TIER1 = BASE / "results/phases_g/tier1_science_supplement.tsv"
EXAMPLE = BASE / "data/reference/tissue_rupture_signature_example.tsv"
DEFAULT_OUT = BASE / "results/phases_g/tissue_signature_phases_g_bridge.tsv"

SHORTLIST_MERGE_COLS = [
    "Gene_ID",
    "ensembl_gene",
    "DESeq2_log2FC",
    "DESeq2_padj",
    "LM_immune_adj_log2_coef_Ruptured",
    "LM_immune_adj_P_Ruptured",
    "blood_eqtl_top_rsid",
    "blood_eqtl_top_pvalue",
    "blood_eqtl_beta",
    "blood_eqtl_dataset_id",
    "IA_GWAS_P_at_eqtl_snp",
    "IA_GWAS_BETA_at_eqtl_snp",
    "GSE159610_Log2FC_aneurysm_vs_control",
    "GSE159610_P",
    "direction_match_rupture_and_aneurysm",
    "marker_celltype_hint",
    "is_known_lineage_marker",
]


def _sign_series(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    return np.sign(x)


def load_de(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#")
    if "symbol" not in df.columns:
        raise SystemExit("Input DE TSV must have a 'symbol' column.")
    df = df.copy()
    df["symbol"] = df["symbol"].astype(str).str.strip().str.upper()
    return df.dropna(subset=["symbol"])


def main() -> int:
    ap = argparse.ArgumentParser(description="Bridge tissue DE genes to PHASES-G shortlist")
    ap.add_argument("--de-tsv", type=Path, default=None, help="Tissue DE table with at least 'symbol'.")
    ap.add_argument("--example", action="store_true", help=f"Use {EXAMPLE} as input.")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    if args.example:
        de_path = EXAMPLE
    elif args.de_tsv and args.de_tsv.is_file():
        de_path = args.de_tsv
    else:
        print("Provide --de-tsv PATH or --example.", file=sys.stderr)
        return 1

    if not SHORTLIST.is_file():
        print(f"Missing shortlist: {SHORTLIST}", file=sys.stderr)
        return 1

    tissue = load_de(de_path)
    sl = pd.read_csv(SHORTLIST, sep="\t")
    sl["symbol"] = sl["symbol"].astype(str).str.strip().str.upper()
    sl_u = sl.drop_duplicates(subset=["symbol"], keep="first")

    merge_cols = ["symbol"] + [c for c in SHORTLIST_MERGE_COLS if c in sl_u.columns]
    out = tissue.merge(sl_u[merge_cols], on="symbol", how="left", suffixes=("", "_phases_g"))

    out["in_phases_g_shortlist"] = out["DESeq2_log2FC"].notna() | out["LM_immune_adj_log2_coef_Ruptured"].notna()

    if "tissue_log2FC" in out.columns:
        t = _sign_series(out["tissue_log2FC"])
        d = _sign_series(out["DESeq2_log2FC"])
        l = _sign_series(out["LM_immune_adj_log2_coef_Ruptured"])
        out["concordant_DESeq2_vs_tissue_fc"] = (t == d) & t.notna() & d.notna()
        out["concordant_LM_vs_tissue_fc"] = (t == l) & t.notna() & l.notna()
    else:
        out["concordant_DESeq2_vs_tissue_fc"] = pd.NA
        out["concordant_LM_vs_tissue_fc"] = pd.NA

    if TIER1.is_file() and "symbol" in pd.read_csv(TIER1, sep="\t", nrows=0).columns:
        t1 = pd.read_csv(TIER1, sep="\t")
        t1["symbol"] = t1["symbol"].astype(str).str.strip().str.upper()
        keep = ["symbol"] + [c for c in ("tier1_science_rank", "tier1_science_score", "tier1_science_rule") if c in t1.columns]
        out = out.merge(t1[keep], on="symbol", how="left")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    n_in = int(out["in_phases_g_shortlist"].sum())
    print(f"Wrote {args.out} ({len(out)} tissue rows; {n_in} matched PHASES-G shortlist)")

    summary_path = args.out.with_name(args.out.stem + "_summary.json")
    summ: dict = {
        "input_de_tsv": str(de_path.resolve()),
        "bridge_tsv": str(args.out.resolve()),
        "n_tissue_genes": int(len(out)),
        "n_in_phases_g_shortlist": n_in,
        "fraction_in_shortlist": round(n_in / max(len(out), 1), 6),
    }
    if "tier1_science_rank" in out.columns:
        tr = pd.to_numeric(out["tier1_science_rank"], errors="coerce")
        summ["n_tier1_science_hits"] = int(tr.notna().sum())
    if "concordant_DESeq2_vs_tissue_fc" in out.columns:
        s = out["concordant_DESeq2_vs_tissue_fc"]
        tested = s.notna()
        summ["n_tested_DESeq2_vs_tissue_direction"] = int(tested.sum())
        summ["n_concordant_DESeq2_vs_tissue"] = int((s.astype("boolean") == True).sum())
    if "concordant_LM_vs_tissue_fc" in out.columns:
        s = out["concordant_LM_vs_tissue_fc"]
        tested = s.notna()
        summ["n_tested_LM_vs_tissue_direction"] = int(tested.sum())
        summ["n_concordant_LM_vs_tissue"] = int((s.astype("boolean") == True).sum())
    summary_path.write_text(json.dumps(summ, indent=2), encoding="utf-8")
    print(f"Wrote {summary_path}")

    miss = out.loc[~out["in_phases_g_shortlist"], "symbol"].tolist()
    if miss:
        cap = 40
        tail = "" if len(miss) <= cap else f" … and {len(miss) - cap} more"
        print(
            "Not in shortlist (run pipeline / expand filters if expected):",
            ", ".join(miss[:cap]) + tail,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
