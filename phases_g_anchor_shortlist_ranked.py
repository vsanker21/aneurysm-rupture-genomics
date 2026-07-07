"""
PHASES-G: build an "anchor shortlist" with explicit rationale.

Goal:
  - Require persistence after immune-fraction adjustment in the rupture cohort (GSE122897)
  - Add orthogonal support from an external bulk cohort (currently GSE159610 aneurysm vs control)
  - Rank genes by a transparent evidence score

Outputs:
  - results/phases_g/anchor_shortlist_ranked.tsv
  - results/phases_g/anchor_shortlist_ranked_top50.tsv

Notes:
  - Replication cohort (GSE159610) is NOT rupture vs unruptured; treated as orthogonal support only.
  - This step does not call eQTL/GWAS. Those remain in run_phases_g_evidence_pipeline.py.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import requests


BASE = Path(r"e:\Aneurysms Genomics Project")
IMMUNE_ADJ = BASE / r"results\cibersortx\deg_immune_covariate_adjusted.tsv"
REPL = BASE / r"results\gse159610_full_de_results.tsv"
MARKERS = BASE / r"data\reference\scrna_aneurysm_celltype_markers.tsv"
OUT_DIR = BASE / r"results\phases_g"

MYGENE_BATCH = "https://mygene.info/v3/gene"


def _safe_neglog10(p: pd.Series) -> pd.Series:
    x = pd.to_numeric(p, errors="coerce").astype(float)
    x = np.clip(x, 1e-300, 1.0)
    return -np.log10(x)


def entrez_to_ensembl_symbol(entrez_ids: list[str], batch_size: int = 400) -> pd.DataFrame:
    rows = []
    for i in range(0, len(entrez_ids), batch_size):
        batch = entrez_ids[i : i + batch_size]
        r = requests.post(
            MYGENE_BATCH,
            data={
                "ids": ",".join(batch),
                "fields": "symbol,ensembl.gene,entrezgene",
                "species": "human",
            },
            timeout=120,
        )
        r.raise_for_status()
        data = r.json()
        if isinstance(data, dict):
            data = [data]
        for row in data:
            q = str(row.get("query", "")).strip()
            sym = row.get("symbol")
            eg = row.get("ensembl")
            ens = None
            if isinstance(eg, dict):
                ens = eg.get("gene")
            elif isinstance(eg, list) and eg:
                ens = eg[0].get("gene") if isinstance(eg[0], dict) else None
            if q and ens:
                rows.append({"Gene_ID": q, "symbol": sym, "ensembl_gene": ens})
    return pd.DataFrame(rows)


def marker_hint(symbol: str, markers: pd.DataFrame) -> str:
    if not isinstance(symbol, str) or not symbol:
        return "unclassified"
    u = symbol.upper()
    hit = markers[markers["gene_symbol"].astype(str).str.upper() == u]
    if len(hit):
        return str(hit.iloc[0]["cell_type_priority"])
    return "unclassified"


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--deseq2-padj-max", type=float, default=0.05)
    ap.add_argument("--lm-p-max", type=float, default=0.05)
    ap.add_argument("--max-genes", type=int, default=1000)
    ap.add_argument(
        "--require-direction-match",
        action="store_true",
        help="Require DESeq2 rupture direction to match GSE159610 aneurysm-vs-control Log2FC sign.",
    )
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    adj = pd.read_csv(IMMUNE_ADJ, sep="\t")
    adj["Gene_ID"] = adj["Gene_ID"].astype(str).str.replace(".0", "", regex=False)
    for c in [
        "DESeq2_log2FC",
        "DESeq2_padj",
        "LM_immune_adj_log2_coef_Ruptured",
        "LM_immune_adj_P_Ruptured",
        "LM_immune_adj_Padj_BH",
        "baseMean",
    ]:
        if c in adj.columns:
            adj[c] = pd.to_numeric(adj[c], errors="coerce")

    # Anchor filter: pre-DESeq2 significant AND persists after immune adjustment (nominal)
    short = adj[
        (adj["DESeq2_padj"] <= args.deseq2_padj_max) & (adj["LM_immune_adj_P_Ruptured"] <= args.lm_p_max)
    ].copy()
    short = short.sort_values("DESeq2_padj").head(args.max_genes)

    if len(short) == 0:
        raise SystemExit("No genes passed anchor filters. Relax thresholds.")

    # Map Entrez -> symbol + Ensembl for replication join
    idmap = entrez_to_ensembl_symbol(short["Gene_ID"].astype(str).tolist())
    m = short.merge(idmap, on="Gene_ID", how="left")

    # Replication cohort (orthogonal): join by Ensembl base ID
    repl = None
    if REPL.is_file():
        repl = pd.read_csv(REPL, sep="\t")
        repl["ensembl_base"] = repl["Gene_ID"].astype(str).str.split(".").str[0].str.strip()
        m["ensembl_base"] = m["ensembl_gene"].astype(str).str.split(".").str[0].str.strip()
        m = m.merge(
            repl[["ensembl_base", "Log2FC", "P_value", "P_adj_BH"]],
            on="ensembl_base",
            how="left",
            suffixes=("", "_GSE159610"),
        )
        m = m.rename(
            columns={
                "Log2FC": "GSE159610_Log2FC_aneurysm_vs_control",
                "P_value": "GSE159610_P",
                "P_adj_BH": "GSE159610_Padj_BH",
            }
        )
    else:
        m["GSE159610_Log2FC_aneurysm_vs_control"] = np.nan
        m["GSE159610_P"] = np.nan
        m["GSE159610_Padj_BH"] = np.nan

    # Direction match: rupture (DESeq2 log2FC) vs aneurysm/control (orthogonal)
    rdir = np.sign(pd.to_numeric(m["DESeq2_log2FC"], errors="coerce").astype(float))
    odir = np.sign(pd.to_numeric(m["GSE159610_Log2FC_aneurysm_vs_control"], errors="coerce").astype(float))
    m["direction_match_rupture_and_aneurysm"] = (rdir != 0) & (odir != 0) & (rdir == odir)

    if args.require_direction_match:
        m = m[m["direction_match_rupture_and_aneurysm"].fillna(False)].copy()

    # Marker hint
    markers = pd.read_csv(MARKERS, sep="\t")
    m["marker_celltype_hint"] = m["symbol"].map(lambda s: marker_hint(s, markers))

    # Transparent scoring
    #   - Core: immune-adjusted rupture p and DESeq2 padj
    #   - Support: orthogonal cohort p if direction matches
    m["score_core"] = _safe_neglog10(m["LM_immune_adj_P_Ruptured"]) + _safe_neglog10(m["DESeq2_padj"])
    repl_support = _safe_neglog10(m["GSE159610_P"]) * m["direction_match_rupture_and_aneurysm"].astype(int)
    m["score_orthogonal_support"] = repl_support.fillna(0.0)
    m["score_total"] = m["score_core"] + 0.5 * m["score_orthogonal_support"]

    # Rank (highest score_total first; tie-breaker by LM p then DESeq2 padj)
    m = m.sort_values(
        ["score_total", "LM_immune_adj_P_Ruptured", "DESeq2_padj"],
        ascending=[False, True, True],
    )
    m["rank"] = np.arange(1, len(m) + 1)

    cols = [
        "rank",
        "symbol",
        "Gene_ID",
        "ensembl_gene",
        "DESeq2_log2FC",
        "DESeq2_padj",
        "LM_immune_adj_log2_coef_Ruptured",
        "LM_immune_adj_P_Ruptured",
        "GSE159610_Log2FC_aneurysm_vs_control",
        "GSE159610_P",
        "direction_match_rupture_and_aneurysm",
        "marker_celltype_hint",
        "score_core",
        "score_orthogonal_support",
        "score_total",
        "baseMean",
    ]
    cols = [c for c in cols if c in m.columns]

    out_all = OUT_DIR / "anchor_shortlist_ranked.tsv"
    out_top = OUT_DIR / "anchor_shortlist_ranked_top50.tsv"
    m[cols].to_csv(out_all, sep="\t", index=False)
    m[cols].head(50).to_csv(out_top, sep="\t", index=False)

    print(f"Wrote {out_all}")
    print(f"Wrote {out_top}")
    print(
        f"Anchor shortlist size={len(m)} (filters: DESeq2_padj<={args.deseq2_padj_max}, "
        f"LM_p<={args.lm_p_max}, require_direction_match={bool(args.require_direction_match)})"
    )
    print("\nMarker hints (top):")
    print(m["marker_celltype_hint"].value_counts().head(10).to_string())


if __name__ == "__main__":
    main()

