"""
Compute scRNA-derived (marker-table based) cell-type specificity scores for genes.

This ingests a marker-strength table with columns:
  gene_symbol, cell_type, marker_strength

and produces per-gene:
  - scrna_best_cell_type
  - scrna_best_score
  - scrna_wall_score (max across EC/SMC/pericyte/fibroblast)
  - scrna_immune_score (max across immune)
  - scrna_wall_minus_immune (specificity proxy)

This is *not* a substitute for full scRNA differential expression; it's a quantitative,
reproducible improvement over strict marker lookup.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


WALL_TYPES = {"EC", "SMC", "PERICYTE", "FIBROBLAST"}


def load_marker_strengths(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#")
    df["gene_symbol"] = df["gene_symbol"].astype(str).str.strip()
    df["cell_type"] = df["cell_type"].astype(str).str.strip().str.upper()
    df["marker_strength"] = pd.to_numeric(df["marker_strength"], errors="coerce").fillna(0.0)

    # normalize cell type labels
    df["cell_type"] = df["cell_type"].replace(
        {
            "VSMC": "SMC",
            "SMOOTH_MUSCLE": "SMC",
            "ENDOTHELIAL": "EC",
            "MYELOID": "IMMUNE",
            "LYMPHOID": "IMMUNE",
        }
    )
    # unify immune
    df.loc[df["cell_type"].str.startswith("IMMUNE"), "cell_type"] = "IMMUNE"
    # unify pericyte/fibroblast
    df.loc[df["cell_type"].str.startswith("PERICYTE"), "cell_type"] = "PERICYTE"
    df.loc[df["cell_type"].str.startswith("FIBRO"), "cell_type"] = "FIBROBLAST"
    return df


def compute_specificity(marker_df: pd.DataFrame) -> pd.DataFrame:
    d = marker_df.copy()
    d["gene_symbol_u"] = d["gene_symbol"].astype(str).str.upper().str.strip()

    # best per gene/cell_type
    gct = (
        d.groupby(["gene_symbol_u", "cell_type"], as_index=False)["marker_strength"]
        .max()
        .rename(columns={"marker_strength": "score"})
    )

    # best cell type overall
    idx = gct.groupby("gene_symbol_u")["score"].idxmax()
    best = gct.loc[idx].rename(columns={"cell_type": "scrna_best_cell_type", "score": "scrna_best_score"})

    # wall score
    wall = gct[gct["cell_type"].isin(WALL_TYPES)].groupby("gene_symbol_u", as_index=False)["score"].max()
    wall = wall.rename(columns={"score": "scrna_wall_score"})

    # immune score
    immune = gct[gct["cell_type"] == "IMMUNE"].groupby("gene_symbol_u", as_index=False)["score"].max()
    immune = immune.rename(columns={"score": "scrna_immune_score"})

    out = best.merge(wall, on="gene_symbol_u", how="left").merge(immune, on="gene_symbol_u", how="left")
    out["scrna_wall_score"] = out["scrna_wall_score"].fillna(0.0)
    out["scrna_immune_score"] = out["scrna_immune_score"].fillna(0.0)
    out["scrna_wall_minus_immune"] = out["scrna_wall_score"] - out["scrna_immune_score"]
    return out


def attach_to_shortlist(shortlist: pd.DataFrame, spec: pd.DataFrame) -> pd.DataFrame:
    out = shortlist.copy()
    out["_SYM_U"] = out["symbol"].astype(str).str.upper().str.strip()
    spec2 = spec.rename(columns={"gene_symbol_u": "_SYM_U"})
    out = out.merge(spec2, on="_SYM_U", how="left").drop(columns=["_SYM_U"])
    return out


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--marker-strengths", type=Path, required=True)
    ap.add_argument("--shortlist-tsv", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, required=True)
    args = ap.parse_args()

    markers = load_marker_strengths(args.marker_strengths)
    spec = compute_specificity(markers)
    short = pd.read_csv(args.shortlist_tsv, sep="\t", low_memory=False)
    if "symbol" not in short.columns:
        raise ValueError("shortlist TSV must contain a 'symbol' column")
    out = attach_to_shortlist(short, spec)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"Wrote {args.out_tsv}")


if __name__ == "__main__":
    main()

