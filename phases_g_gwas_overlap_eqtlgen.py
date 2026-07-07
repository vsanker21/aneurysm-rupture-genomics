"""
PHASES-G: overlap eQTLGen sentinel cis-eQTL SNPs with Bakker 2020 IA GWAS.

Input:
  - results/phases_g/anchor_shortlist_ranked_with_eqtlgen.tsv
  - data/gwas/bakker_2020/downloaded/IA.GWAS.BakkerMK.2020.sumstats.Stage_2.txt.gz

Output:
  - results/phases_g/anchor_shortlist_ranked_with_eqtlgen_and_gwas.tsv

Method:
  - rsID intersection only (not colocalization).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(r"e:\Aneurysms Genomics Project")
IN_TSV = BASE / r"results\phases_g\anchor_shortlist_ranked_with_eqtlgen.tsv"
GWAS = BASE / r"data\gwas\bakker_2020\downloaded\IA.GWAS.BakkerMK.2020.sumstats.Stage_2.txt.gz"
OUT_TSV = BASE / r"results\phases_g\anchor_shortlist_ranked_with_eqtlgen_and_gwas.tsv"


def load_gwas_rows_for_rsids(path: Path, rsids: set[str]) -> pd.DataFrame:
    need = {r for r in rsids if r and str(r).startswith("rs")}
    hits = []
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", chunksize=400_000, low_memory=False):
        if "SNP" not in chunk.columns:
            raise ValueError(f"Unexpected GWAS columns: {chunk.columns.tolist()[:30]}")
        sub = chunk[chunk["SNP"].astype(str).isin(need)]
        if len(sub):
            hits.append(sub)
    if not hits:
        return pd.DataFrame()
    return pd.concat(hits, ignore_index=True).drop_duplicates(subset=["SNP"])


def main() -> None:
    if not IN_TSV.is_file():
        raise FileNotFoundError(f"Missing input: {IN_TSV}")
    if not GWAS.is_file():
        raise FileNotFoundError(f"Missing GWAS: {GWAS}")

    df = pd.read_csv(IN_TSV, sep="\t", low_memory=False)
    if "eQTLGen_top_SNP" not in df.columns:
        raise ValueError("Input missing column eQTLGen_top_SNP")

    rsids = set(df["eQTLGen_top_SNP"].dropna().astype(str).unique().tolist())
    rsids = {r for r in rsids if r.startswith("rs")}

    gwas_hits = load_gwas_rows_for_rsids(GWAS, rsids) if rsids else pd.DataFrame()
    gmap = {}
    if len(gwas_hits):
        for _, gr in gwas_hits.iterrows():
            gmap[str(gr["SNP"])] = {
                "IA_GWAS_P": gr.get("P", gr.get("PVAL", np.nan)),
                "IA_GWAS_BETA": gr.get("BETA", np.nan),
                "IA_GWAS_SE": gr.get("SE", np.nan),
                "IA_GWAS_OR": gr.get("OR", np.nan),
            }

    df["IA_GWAS_P_at_eQTLGen_SNP"] = df["eQTLGen_top_SNP"].map(
        lambda x: gmap.get(str(x), {}).get("IA_GWAS_P", np.nan) if pd.notna(x) else np.nan
    )
    df["IA_GWAS_BETA_at_eQTLGen_SNP"] = df["eQTLGen_top_SNP"].map(
        lambda x: gmap.get(str(x), {}).get("IA_GWAS_BETA", np.nan) if pd.notna(x) else np.nan
    )
    df["IA_GWAS_SE_at_eQTLGen_SNP"] = df["eQTLGen_top_SNP"].map(
        lambda x: gmap.get(str(x), {}).get("IA_GWAS_SE", np.nan) if pd.notna(x) else np.nan
    )
    df["IA_GWAS_OR_at_eQTLGen_SNP"] = df["eQTLGen_top_SNP"].map(
        lambda x: gmap.get(str(x), {}).get("IA_GWAS_OR", np.nan) if pd.notna(x) else np.nan
    )

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TSV, sep="\t", index=False)

    n_overlap = int(pd.to_numeric(df["IA_GWAS_P_at_eQTLGen_SNP"], errors="coerce").notna().sum())
    print(f"Wrote {OUT_TSV}")
    print(f"eQTLGen SNPs overlapped in IA GWAS: {n_overlap}/{len(df)}")


if __name__ == "__main__":
    main()

