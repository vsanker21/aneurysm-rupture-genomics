"""Weighted Stouffer meta on primary GSE122897 DESeq2 (21 vs 21 explicit labels).

GSE54083 and GSE13353 inputs are unchanged. Genome-wide BH-FDR is recomputed
on the overlapping gene set so FDR for the six-gene core is not taken from the
22 vs 21 run.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

from generate_supplementary_materials_complete import CORE_ENTREZ, TISSUE_CORE_6
from manuscript_results_detail import annotate_deseq2_with_symbols, load_gene_symbol_map

BASE = Path(r"E:\Aneurysms Genomics Project")
OUT = BASE / "results" / "deseq2_gse122897_21v21_explicit"
CORE_OUT = OUT / "six_gene_stouffer_meta_21v21.tsv"


def _effective_n(n_case: int, n_control: int) -> float:
    return float(np.sqrt((n_case * n_control) / max(n_case + n_control, 1)))


def _signed_z(lfc: pd.Series, p: pd.Series) -> pd.Series:
    eps = 1e-300
    return np.sign(lfc).replace(0, np.nan) * norm.isf(p.clip(lower=eps) / 2.0)


def main() -> None:
    de = pd.read_csv(BASE / "results/gse122897_deseq2_results.tsv", sep="\t")
    de["Gene_ID"] = de["Gene_ID"].astype(str).str.strip('"').str.replace(r"\.0$", "", regex=True)
    sym = load_gene_symbol_map(BASE)
    de = annotate_deseq2_with_symbols(de, sym)
    entrez_to_core = {v: k for k, v in CORE_ENTREZ.items()}
    de["display_symbol"] = de["Gene_ID"].map(entrez_to_core).fillna(de["display_symbol"])

    g122 = de[["display_symbol", "log2FoldChange", "pvalue", "padj"]].rename(
        columns={
            "display_symbol": "Gene",
            "log2FoldChange": "Log2FC_GSE122897",
            "pvalue": "P_GSE122897",
            "padj": "FDR_GSE122897",
        }
    ).drop_duplicates("Gene")

    g540 = pd.read_csv(
        BASE / "results/external_validation_gse54083/gse54083_gene_level_validation.tsv",
        sep="\t",
    )[["Gene", "Log2FC_RvU", "P_value", "FDR_BH"]].rename(
        columns={"Log2FC_RvU": "Log2FC_GSE54083", "P_value": "P_GSE54083", "FDR_BH": "FDR_GSE54083"}
    ).drop_duplicates("Gene")

    g133 = pd.read_csv(
        BASE / "results/rupture_meta_analysis/gse13353_gene_level_meta_input.tsv",
        sep="\t",
    )[["Gene", "Log2FC", "P_value", "FDR_BH"]].rename(
        columns={"Log2FC": "Log2FC_GSE13353", "P_value": "P_GSE13353", "FDR_BH": "FDR_GSE13353"}
    ).drop_duplicates("Gene")

    explicit = pd.read_csv(OUT / "sample_metadata_explicit_21v21.tsv", sep="\t")
    n122_r = int((explicit["Group"] == "Ruptured").sum())
    n122_u = int((explicit["Group"] == "Unruptured").sum())
    if n122_r != 21 or n122_u != 21:
        raise SystemExit(f"Expected 21 vs 21, got {n122_r} vs {n122_u}")

    meta540 = pd.read_csv(
        BASE / "results/external_validation_gse54083/gse54083_sample_metadata.tsv", sep="\t"
    )
    n540_r = int((meta540["Group"] == "Ruptured").sum())
    n540_u = int((meta540["Group"] == "Unruptured").sum())
    n133_r, n133_u = 11, 8

    merged = g122.merge(g540, on="Gene", how="inner").merge(g133, on="Gene", how="inner")
    merged["Z_GSE122897"] = _signed_z(merged["Log2FC_GSE122897"], merged["P_GSE122897"])
    merged["Z_GSE54083"] = _signed_z(merged["Log2FC_GSE54083"], merged["P_GSE54083"])
    merged["Z_GSE13353"] = _signed_z(merged["Log2FC_GSE13353"], merged["P_GSE13353"])
    merged["Weight_GSE122897"] = _effective_n(n122_r, n122_u)
    merged["Weight_GSE54083"] = _effective_n(n540_r, n540_u)
    merged["Weight_GSE13353"] = _effective_n(n133_r, n133_u)

    w = merged[["Weight_GSE122897", "Weight_GSE54083", "Weight_GSE13353"]].to_numpy(dtype=float)
    z = merged[["Z_GSE122897", "Z_GSE54083", "Z_GSE13353"]].to_numpy(dtype=float)
    merged["Meta_Z"] = np.nansum(w * z, axis=1) / np.sqrt(np.nansum(w**2, axis=1))
    merged["Meta_P"] = 2.0 * norm.sf(np.abs(merged["Meta_Z"]))
    merged["Meta_FDR_BH"] = multipletests(merged["Meta_P"].fillna(1.0).values, method="fdr_bh")[1]
    merged["All_three_directionally_concordant"] = (
        (np.sign(merged["Log2FC_GSE122897"]) == np.sign(merged["Log2FC_GSE54083"]))
        & (np.sign(merged["Log2FC_GSE122897"]) == np.sign(merged["Log2FC_GSE13353"]))
    )

    core = merged[merged["Gene"].isin(TISSUE_CORE_6)].copy()
    core["Gene"] = pd.Categorical(core["Gene"], TISSUE_CORE_6, ordered=True)
    core = core.sort_values("Gene")
    keep = [
        "Gene",
        "Log2FC_GSE122897",
        "P_GSE122897",
        "Z_GSE122897",
        "Weight_GSE122897",
        "Log2FC_GSE54083",
        "P_GSE54083",
        "Z_GSE54083",
        "Weight_GSE54083",
        "Log2FC_GSE13353",
        "P_GSE13353",
        "Z_GSE13353",
        "Weight_GSE13353",
        "Meta_Z",
        "Meta_P",
        "Meta_FDR_BH",
        "All_three_directionally_concordant",
    ]
    core[keep].to_csv(CORE_OUT, sep="\t", index=False)
    merged.to_csv(OUT / "rupture_consensus_meta_all_21v21.tsv", sep="\t", index=False)

    n_conc = int(
        (
            merged["All_three_directionally_concordant"]
            & (merged["Meta_FDR_BH"] < 0.05)
        ).sum()
    )
    summary = {
        "GSE122897_n": [n122_r, n122_u],
        "GSE122897_weight": _effective_n(n122_r, n122_u),
        "n_overlap": int(len(merged)),
        "n_concordant_fdr_lt_0.05": n_conc,
        "core": core[keep].to_dict(orient="records"),
    }
    (OUT / "six_gene_stouffer_meta_21v21.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )
    print(core[keep].to_string(index=False))
    print("n_overlap", len(merged), "concordant FDR<0.05", n_conc)
    print("weight 21v21", _effective_n(21, 21), "weight 22v21", _effective_n(22, 21))


if __name__ == "__main__":
    main()
