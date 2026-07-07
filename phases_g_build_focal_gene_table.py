"""
Build a manuscript-ready table linking focal genes to PHASES-G shortlist columns
and (if present) exploratory COLOC outputs.

Default: top N symbols from tier1_coloc_candidates.tsv (Bakker-present tier-1 science).

Outputs:
  results/phases_g/focal_genes_phases_g_table.tsv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent
SHORTLIST = BASE / "results/phases_g/shortlist_evidence.tsv"
COLCAND = BASE / "results/phases_g/tier1_coloc_candidates.tsv"
COLOC_SUM = BASE / "results/phases_g/stronger_genetics/coloc/coloc_summary.tsv"
OUT = BASE / "results/phases_g/focal_genes_phases_g_table.tsv"

SHORTLIST_COLS = [
    "symbol",
    "ensembl_gene",
    "Gene_ID",
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
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--symbols",
        type=str,
        default="",
        help="Comma-separated HGNC symbols (overrides --from-coloc-candidates).",
    )
    ap.add_argument(
        "--symbols-file",
        type=Path,
        default=None,
        help="Text file: one symbol per line.",
    )
    ap.add_argument(
        "--from-coloc-candidates",
        action="store_true",
        help="Take top N symbols from tier1_coloc_candidates.tsv in file order.",
    )
    ap.add_argument("--top-n", type=int, default=6)
    args = ap.parse_args()

    if not SHORTLIST.is_file():
        raise FileNotFoundError(SHORTLIST)

    if args.symbols.strip():
        syms = [s.strip().upper() for s in args.symbols.split(",") if s.strip()]
    elif args.symbols_file and args.symbols_file.is_file():
        syms = [
            ln.strip().upper()
            for ln in args.symbols_file.read_text(encoding="utf-8").splitlines()
            if ln.strip()
        ]
    elif args.from_coloc_candidates and COLCAND.is_file():
        cc = pd.read_csv(COLCAND, sep="\t")
        syms = [str(s).strip().upper() for s in cc["symbol"].head(int(args.top_n))]
    else:
        raise SystemExit(
            "Provide --symbols, --symbols-file, or --from-coloc-candidates (and tier1_coloc_candidates.tsv)."
        )

    sl = pd.read_csv(SHORTLIST, sep="\t")
    sl["symbol"] = sl["symbol"].astype(str).str.strip().str.upper()
    want = set(syms)
    sub = sl[sl["symbol"].isin(want)].drop_duplicates(subset=["symbol"], keep="first")
    missing = sorted(want - set(sub["symbol"].tolist()))
    cols = [c for c in SHORTLIST_COLS if c in sub.columns]
    out = sub[cols].copy()
    out["locus_id"] = out["symbol"] + "_" + out["blood_eqtl_top_rsid"].astype(str)

    if COLOC_SUM.is_file():
        co = pd.read_csv(COLOC_SUM, sep="\t")
        # tolerate extended columns (status, error)
        merge_cols = ["locus_id"] + [c for c in co.columns if c != "locus_id"]
        out = out.merge(co, on="locus_id", how="left", suffixes=("", "_coloc"))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    meta = {"symbols_requested": syms, "missing_from_shortlist": missing}
    (BASE / "results/phases_g/focal_genes_phases_g_meta.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8"
    )
    out.to_csv(OUT, sep="\t", index=False)
    print(f"Wrote {OUT} ({len(out)} rows)")
    if missing:
        print("Missing from shortlist:", ", ".join(missing))


if __name__ == "__main__":
    main()
