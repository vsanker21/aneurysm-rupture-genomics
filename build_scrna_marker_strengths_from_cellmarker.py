"""
Build a large marker-strength table from CellMarker downloads.

Inputs (downloaded earlier):
  - data/external/cellmarker/Cell_marker_Human.xlsx
  - data/external/cellmarker/Cell_marker_Seq.xlsx   (optional; scRNA-derived)

Outputs:
  - data/reference/scrna_celltype_marker_strengths_cellmarker.tsv

Approach:
  - CellMarker 2.0 exports one marker per row (Symbol column).
  - Map fine-grained cell names/types to broad buckets:
      EC / SMC / PERICYTE / FIBROBLAST / IMMUNE
  - For each (bucket, gene), compute how many independent rows support it across tissues/publications.
  - Convert support counts to a [0,1] marker_strength using a log-scaled normalization within bucket.

This is a pragmatic, publication-defensible reference-atlas approach when a single
aneurysm-wall scRNA marker supplement is not directly downloadable.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(r"e:\Aneurysms Genomics Project")
HUMAN_XLSX = BASE / r"data\external\cellmarker\Cell_marker_Human.xlsx"
SEQ_XLSX = BASE / r"data\external\cellmarker\Cell_marker_Seq.xlsx"
OUT = BASE / r"data\reference\scrna_celltype_marker_strengths_cellmarker.tsv"


CELLTYPE_MAP = {
    "endothelial cell": "EC",
    "endothelial cells": "EC",
    "vascular endothelial cell": "EC",
    "vascular endothelial cells": "EC",
    "smooth muscle cell": "SMC",
    "smooth muscle cells": "SMC",
    "vascular smooth muscle cell": "SMC",
    "vascular smooth muscle cells": "SMC",
    "pericyte": "PERICYTE",
    "pericytes": "PERICYTE",
    "fibroblast": "FIBROBLAST",
    "fibroblasts": "FIBROBLAST",
    # Immune buckets (many names will map here)
    "t cell": "IMMUNE",
    "t cells": "IMMUNE",
    "b cell": "IMMUNE",
    "b cells": "IMMUNE",
    "nk cell": "IMMUNE",
    "macrophage": "IMMUNE",
    "monocyte": "IMMUNE",
    "neutrophil": "IMMUNE",
    "dendritic cell": "IMMUNE",
    "mast cell": "IMMUNE",
    "immune cell": "IMMUNE",
    "immune cells": "IMMUNE",
    "leukocyte": "IMMUNE",
}


def _norm_celltype(x: str) -> str:
    s = str(x or "").strip().lower()
    if s in CELLTYPE_MAP:
        return CELLTYPE_MAP[s]
    # fuzzy fallbacks
    if "endothelial" in s:
        return "EC"
    if "smooth muscle" in s:
        return "SMC"
    if "pericyte" in s:
        return "PERICYTE"
    if "fibroblast" in s:
        return "FIBROBLAST"
    # immune: common tokens
    if any(t in s for t in ["t cell", "b cell", "nk", "macroph", "mono", "neutro", "mast", "dendritic", "immune", "leuk"]):
        return "IMMUNE"
    return ""


def _split_markers(val: str) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    s = str(val)
    # Common delimiters in CellMarker exports
    for ch in [";", ",", "|", "/", "\n", "\t"]:
        s = s.replace(ch, " ")
    toks = [t.strip().upper() for t in s.split() if t.strip()]
    # filter out obvious non-gene junk
    toks = [t for t in toks if t.isalnum() or ("-" in t)]
    return toks


def load_cellmarker(path: Path) -> pd.DataFrame:
    if not path.is_file():
        return pd.DataFrame()
    # Most CellMarker sheets load as first sheet; be permissive.
    df = pd.read_excel(path)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def main() -> None:
    if not HUMAN_XLSX.is_file():
        raise FileNotFoundError(f"Missing {HUMAN_XLSX}")

    df = load_cellmarker(HUMAN_XLSX)
    df_seq = load_cellmarker(SEQ_XLSX) if SEQ_XLSX.is_file() else pd.DataFrame()

    # Columns in CellMarker2.0: cell_name, cell_type, tissue_type, Symbol, PMID, etc.
    if "Symbol" not in df.columns:
        raise ValueError(f"Unexpected CellMarker columns (missing Symbol): {df.columns.tolist()[:40]}")
    tissue_col = "tissue_type" if "tissue_type" in df.columns else ""
    cell_name_col = "cell_name" if "cell_name" in df.columns else "cell_type"

    df["_bucket"] = df[cell_name_col].map(_norm_celltype)
    df["_gene"] = df["Symbol"].astype(str).str.upper().str.strip()
    df["_tissue"] = df[tissue_col].astype(str).str.strip() if tissue_col else ""

    df = df[df["_bucket"].astype(bool) & df["_gene"].astype(bool)].copy()

    # Support counts per (bucket, gene)
    support = (
        df.groupby(["_bucket", "_gene"], as_index=False)
        .size()
        .rename(columns={"size": "support_n"})
        .sort_values(["_bucket", "support_n"], ascending=[True, False])
    )

    # Convert to log-scaled [0,1] within bucket
    rows = []
    for bucket, sub in support.groupby("_bucket"):
        mx = int(sub["support_n"].max())
        denom = math.log1p(mx) if mx > 0 else 1.0
        for _, r in sub.iterrows():
            strength = float(math.log1p(int(r["support_n"])) / denom) if denom > 0 else 0.0
            rows.append(
                {
                    "gene_symbol": r["_gene"],
                    "cell_type": bucket,
                    "marker_strength": strength,
                    "source": "CellMarker_Human_aggregated_support",
                }
            )

    # Optionally add scRNA-derived markers, giving them a small boost.
    if len(df_seq):
        if "Symbol" in df_seq.columns:
            tissue_col2 = "tissue_type" if "tissue_type" in df_seq.columns else ""
            cell_name_col2 = "cell_name" if "cell_name" in df_seq.columns else "cell_type"
            df_seq["_bucket"] = df_seq[cell_name_col2].map(_norm_celltype)
            df_seq["_gene"] = df_seq["Symbol"].astype(str).str.upper().str.strip()
            df_seq = df_seq[df_seq["_bucket"].astype(bool) & df_seq["_gene"].astype(bool)].copy()
            support2 = (
                df_seq.groupby(["_bucket", "_gene"], as_index=False)
                .size()
                .rename(columns={"size": "support_n"})
            )
            for bucket, sub in support2.groupby("_bucket"):
                mx = int(sub["support_n"].max())
                denom = math.log1p(mx) if mx > 0 else 1.0
                for _, r in sub.iterrows():
                    strength = float(math.log1p(int(r["support_n"])) / denom) if denom > 0 else 0.0
                    rows.append(
                        {
                            "gene_symbol": r["_gene"],
                            "cell_type": bucket,
                            "marker_strength": min(1.0, 1.10 * strength),
                            "source": "CellMarker_Seq_aggregated_support",
                        }
                    )

    out = pd.DataFrame(rows)
    if len(out) == 0:
        raise RuntimeError("No markers extracted from CellMarker.")

    # Aggregate: keep max strength per gene/cell_type across tissues/sources
    agg = (
        out.groupby(["gene_symbol", "cell_type"], as_index=False)["marker_strength"]
        .max()
        .sort_values(["cell_type", "marker_strength"], ascending=[True, False])
    )
    agg["source"] = "CellMarker_aggregated"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    agg.to_csv(OUT, sep="\t", index=False)
    print(f"Wrote {OUT} rows={len(agg)} unique_gene_celltype={len(agg)}")
    print("Cell types:", agg["cell_type"].value_counts().to_string())


if __name__ == "__main__":
    main()

