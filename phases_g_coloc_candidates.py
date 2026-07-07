"""
Build a COLOC-ready tier-1 subset by intersecting tier-1 science table with
stronger-genetics loci_manifest.json (requires sentinel rsID present in Bakker GWAS).

Outputs:
  results/phases_g/tier1_coloc_candidates.tsv
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent
TIER1 = BASE / "results/phases_g/tier1_science_supplement.tsv"
MANIFEST = BASE / "results/phases_g/stronger_genetics/loci_manifest.json"
OUT = BASE / "results/phases_g/tier1_coloc_candidates.tsv"


def main() -> None:
    if not TIER1.is_file():
        raise FileNotFoundError(TIER1)
    if not MANIFEST.is_file():
        raise FileNotFoundError(MANIFEST)

    t1 = pd.read_csv(TIER1, sep="\t")
    loci = json.loads(MANIFEST.read_text(encoding="utf-8"))
    man = pd.DataFrame(loci)
    man = man[["locus_id", "gwas_region_file"]].dropna()

    # stronger-genetics locus_id convention: SYMBOL_rsID
    t1["locus_id"] = t1["symbol"].astype(str).str.strip() + "_" + t1["blood_eqtl_top_rsid"].astype(str).str.strip()
    out = t1.merge(man, on="locus_id", how="inner")
    out = out.sort_values("tier1_science_rank", ascending=True)

    cols = [
        "tier1_science_rank",
        "tier1_science_score",
        "symbol",
        "ensembl_gene",
        "blood_eqtl_top_rsid",
        "blood_eqtl_top_pvalue",
        "blood_eqtl_dataset_id",
        "IA_GWAS_P_at_eqtl_snp",
        "LM_immune_adj_P_Ruptured",
        "GSE159610_P",
        "marker_celltype_hint",
        "locus_id",
        "gwas_region_file",
    ]
    cols = [c for c in cols if c in out.columns]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out[cols].to_csv(OUT, sep="\t", index=False)
    print(f"Wrote {OUT} ({len(out)} loci)")


if __name__ == "__main__":
    main()

