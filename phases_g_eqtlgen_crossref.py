"""
PHASES-G: eQTLGen blood cis-eQTL cross-reference for the anchor shortlist.

Uses eQTLGen cis-eQTL significant file (FDR<0.05) and extracts, for each gene
in the anchor shortlist, the top association by Pvalue (smallest).

Input:
  - results/phases_g/anchor_shortlist_ranked.tsv
  - data/external/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz

Outputs:
  - results/phases_g/eqtlgen_top_hit_per_gene.tsv
  - results/phases_g/anchor_shortlist_ranked_with_eqtlgen.tsv

Notes:
  - This uses the FDR<0.05 file (not the 4GB full file). For instrument discovery it is usually sufficient.
  - Variant overlap with IA GWAS is handled separately (next step) so we can use rsID joins consistently.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(r"e:\Aneurysms Genomics Project")
ANCHOR = BASE / r"results\phases_g\anchor_shortlist_ranked.tsv"
EQTLGEN = BASE / r"data\external\eqtlgen\2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
OUT_DIR = BASE / r"results\phases_g"


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--anchor", type=Path, default=ANCHOR)
    ap.add_argument("--eqtlgen", type=Path, default=EQTLGEN)
    ap.add_argument("--out-prefix", type=str, default="anchor_shortlist_ranked")
    args = ap.parse_args()

    if not args.anchor.is_file():
        raise FileNotFoundError(f"Missing anchor shortlist: {args.anchor}")
    if not args.eqtlgen.is_file():
        raise FileNotFoundError(f"Missing eQTLGen file: {args.eqtlgen}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    anchor = pd.read_csv(args.anchor, sep="\t", low_memory=False)
    # Prefer gene symbol if present; fall back to Ensembl if needed
    want_symbols = (
        anchor["symbol"].dropna().astype(str).str.strip().str.upper().unique().tolist()
        if "symbol" in anchor.columns
        else []
    )
    want_ensg = (
        anchor["ensembl_gene"].dropna().astype(str).str.strip().str.upper().unique().tolist()
        if "ensembl_gene" in anchor.columns
        else []
    )
    want_symbols_set = set(want_symbols)
    want_ensg_set = set(want_ensg)

    # Track best hit per gene symbol (and per ENSG as fallback)
    best_by_symbol: dict[str, dict] = {}
    best_by_ensg: dict[str, dict] = {}

    with gzip.open(args.eqtlgen, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        required = ["Pvalue", "SNP", "Zscore", "AssessedAllele", "OtherAllele", "Gene", "GeneSymbol", "FDR"]
        missing = [c for c in required if c not in col]
        if missing:
            raise ValueError(f"eQTLGen header missing columns: {missing}\nGot: {header[:30]}")

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            sym = parts[col["GeneSymbol"]].strip().upper()
            ensg = parts[col["Gene"]].strip().upper()

            keep = False
            if sym and sym in want_symbols_set:
                keep = True
            elif ensg and ensg in want_ensg_set:
                keep = True
            if not keep:
                continue

            try:
                p = float(parts[col["Pvalue"]])
            except Exception:
                continue

            row = {
                "GeneSymbol": parts[col["GeneSymbol"]].strip(),
                "Gene": parts[col["Gene"]].strip(),
                "eQTLGen_top_SNP": parts[col["SNP"]].strip(),
                "eQTLGen_top_Pvalue": p,
                "eQTLGen_top_Zscore": float(parts[col["Zscore"]]),
                "eQTLGen_top_AssessedAllele": parts[col["AssessedAllele"]].strip(),
                "eQTLGen_top_OtherAllele": parts[col["OtherAllele"]].strip(),
                "eQTLGen_top_FDR": float(parts[col["FDR"]]),
            }

            if sym and sym in want_symbols_set:
                prev = best_by_symbol.get(sym)
                if prev is None or p < float(prev["eQTLGen_top_Pvalue"]):
                    best_by_symbol[sym] = row
            if ensg and ensg in want_ensg_set:
                prev = best_by_ensg.get(ensg)
                if prev is None or p < float(prev["eQTLGen_top_Pvalue"]):
                    best_by_ensg[ensg] = row

    # Build a per-gene table
    hits = pd.DataFrame(list(best_by_symbol.values()))
    if len(hits) == 0 and len(best_by_ensg):
        hits = pd.DataFrame(list(best_by_ensg.values()))

    out_hits = OUT_DIR / "eqtlgen_top_hit_per_gene.tsv"
    hits.to_csv(out_hits, sep="\t", index=False)

    # Join back to anchor:
    out = anchor.copy()
    # Prefer symbol join; if empty symbols, attempt ENSG join
    if "symbol" in out.columns and len(best_by_symbol):
        out["_SYM_U"] = out["symbol"].astype(str).str.upper().str.strip()
        map_df = pd.DataFrame([{"_SYM_U": k, **v} for k, v in best_by_symbol.items()])
        out = out.merge(map_df, on="_SYM_U", how="left").drop(columns=["_SYM_U"])
    elif "ensembl_gene" in out.columns and len(best_by_ensg):
        out["_ENSG_U"] = out["ensembl_gene"].astype(str).str.upper().str.strip()
        map_df = pd.DataFrame([{"_ENSG_U": k, **v} for k, v in best_by_ensg.items()])
        out = out.merge(map_df, on="_ENSG_U", how="left").drop(columns=["_ENSG_U"])
    else:
        # still write anchor unchanged, but with empty columns
        for c in [
            "eQTLGen_top_SNP",
            "eQTLGen_top_Pvalue",
            "eQTLGen_top_Zscore",
            "eQTLGen_top_AssessedAllele",
            "eQTLGen_top_OtherAllele",
            "eQTLGen_top_FDR",
        ]:
            out[c] = np.nan

    out_anchor = OUT_DIR / f"{args.out_prefix}_with_eqtlgen.tsv"
    out.to_csv(out_anchor, sep="\t", index=False)

    n_hit = int(pd.to_numeric(out.get("eQTLGen_top_Pvalue"), errors="coerce").notna().sum())
    print(f"Wrote {out_hits}")
    print(f"Wrote {out_anchor}")
    print(f"Anchor genes with an eQTLGen cis hit: {n_hit}/{len(out)}")


if __name__ == "__main__":
    main()

