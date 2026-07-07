"""
Build quantitative cell-type marker-strengths from CVAE186 (Cardiovascular Research 2024) supplemental tables.

Input:
  - results/cvae186/Supplemental_Tables.xlsx

Uses marker gene tables (as described in data/supplementary/cvae186/README.md):
  - Table S5: Marker genes for each identified cell type
  - Table S7: Marker genes for each endothelial cell subcluster
  - Table S8: Marker genes for each vascular smooth muscle cell subcluster
  - Table S9: Marker genes for each fibroblast subcluster

Output:
  - data/reference/scrna_celltype_marker_strengths_cvae186.tsv

The goal is not to perfectly reproduce the authors' statistics, but to extract
a large, quantitative marker resource that can be used to assign bulk shortlist
genes to likely cell origins (EC/SMC/fibroblast/pericyte/immune) without the
severe sparsity of canonical-marker lookup.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(r"e:\Aneurysms Genomics Project")
XLSX = BASE / r"results\cvae186\Supplemental_Tables.xlsx"
OUT = BASE / r"data\reference\scrna_celltype_marker_strengths_cvae186.tsv"


def _norm_celltype(name: str) -> str:
    s = str(name or "").strip().upper()
    # Normalize obvious families
    if "ENDOTHEL" in s or s.startswith("EC") or "EC_" in s:
        return "EC"
    if "SMC" in s or "VSMC" in s or "SMOOTH" in s:
        return "SMC"
    if "FIBRO" in s:
        return "FIBROBLAST"
    if "PERICY" in s:
        return "PERICYTE"
    # immune buckets
    if any(k in s for k in ["MAC", "MONO", "MYELOID", "NEUT", "DC", "DENDR", "T CELL", "B CELL", "NK", "MAST", "IMMUNE"]):
        return "IMMUNE"
    return "OTHER"


def _find_gene_col(df: pd.DataFrame) -> str:
    cols = [c for c in df.columns]
    # common possibilities
    for c in cols:
        cu = str(c).strip().upper()
        if cu in ("GENE", "GENE_SYMBOL", "GENESYMBOL", "SYMBOL"):
            return c
    # fallback: first column containing 'gene'
    for c in cols:
        if "gene" in str(c).lower():
            return c
    return cols[0]


def _find_cluster_col(df: pd.DataFrame) -> str | None:
    for c in df.columns:
        cu = str(c).strip().lower()
        if cu in ("cluster", "celltype", "cell_type", "cell type", "cell_type_name"):
            return c
        if "cluster" in cu or "cell" in cu and "type" in cu:
            return c
    return None


def _find_strength_cols(df: pd.DataFrame) -> tuple[str | None, str | None]:
    # log2FC column
    lfc = None
    padj = None
    for c in df.columns:
        cu = str(c).strip().lower()
        if lfc is None and ("log2fc" in cu or "avg_log2fc" in cu or "avg_logfc" in cu):
            lfc = c
        if padj is None and ("p_val_adj" in cu or "padj" in cu or "p_adj" in cu or "adj_p" in cu):
            padj = c
    return lfc, padj


def _strength_from_stats(log2fc: float | None, padj: float | None) -> float:
    # Conservative, bounded: combines effect + significance when available.
    s = 0.6
    if log2fc is not None and np.isfinite(log2fc):
        s = max(s, min(1.0, float(abs(log2fc)) / 2.0))  # abs(log2fc)=2 -> 1.0
    if padj is not None and np.isfinite(padj) and padj > 0:
        s = max(s, min(1.0, (-np.log10(padj)) / 20.0))  # padj 1e-20 -> 1.0
    return float(np.clip(s, 0.05, 1.0))


def read_sheet(xlsx: Path, sheet: str) -> pd.DataFrame:
    # CVAE186 sheets are formatted with a title row, then a header row, then data.
    # The first data-like row often repeats the header labels; header=2 yields columns:
    # gene, p_val, avg_logFC/avg_log2FC, pct.1, pct.2, p_val_adj, cluster
    return pd.read_excel(xlsx, sheet_name=sheet, header=2)


def extract_markers(sheet_name: str, df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    gene_col = _find_gene_col(df)
    cluster_col = _find_cluster_col(df)
    lfc_col, padj_col = _find_strength_cols(df)

    # Standardize gene symbols
    df["gene_symbol"] = df[gene_col].astype(str).str.upper().str.strip()
    df = df[df["gene_symbol"].str.fullmatch(r"[A-Z0-9][A-Z0-9\\-\\.]*", na=False)].copy()

    if cluster_col is None:
        df["cluster_name"] = ""
    else:
        df["cluster_name"] = df[cluster_col].astype(str).str.strip()

    if lfc_col:
        df["log2fc"] = pd.to_numeric(df[lfc_col], errors="coerce")
    else:
        df["log2fc"] = np.nan
    if padj_col:
        df["padj"] = pd.to_numeric(df[padj_col], errors="coerce")
    else:
        df["padj"] = np.nan

    # Use sheet context to assign broad lineage where possible
    su = sheet_name.strip().upper()
    if su.startswith("TABLE S7"):
        df["cell_type"] = "EC"
    elif su.startswith("TABLE S8"):
        df["cell_type"] = "SMC"
    elif su.startswith("TABLE S9"):
        df["cell_type"] = "FIBROBLAST"
    else:
        # Table S5: mixed cell types; infer from cluster label
        df["cell_type"] = df["cluster_name"].map(_norm_celltype)
    df["marker_strength"] = [
        _strength_from_stats(l2, pj) for l2, pj in zip(df["log2fc"].tolist(), df["padj"].tolist())
    ]
    df["source"] = f"CVAE186::{sheet_name}"

    return df[["gene_symbol", "cell_type", "marker_strength", "source", "cluster_name", "log2fc", "padj"]]


def main() -> None:
    if not XLSX.is_file():
        raise FileNotFoundError(f"Missing {XLSX}")

    xl = pd.ExcelFile(XLSX)
    sheets = xl.sheet_names

    # Select relevant sheets by name heuristics
    wanted = []
    for s in sheets:
        su = s.strip().upper()
        if su.startswith("TABLE S5") or su.startswith("TABLE S7") or su.startswith("TABLE S8") or su.startswith("TABLE S9"):
            wanted.append(s)
    if not wanted:
        raise RuntimeError(f"No expected sheets found. Sheets: {sheets}")

    rows = []
    for s in wanted:
        df = read_sheet(XLSX, s)
        rows.append(extract_markers(s, df))

    allm = pd.concat(rows, ignore_index=True)
    # Drop OTHER if present (should be minimal; indicates unmapped cluster labels in Table S5)
    allm = allm[allm["cell_type"] != "OTHER"].copy()

    # Aggregate to per-gene per broad cell_type, keep max strength
    agg = (
        allm.groupby(["gene_symbol", "cell_type"], as_index=False)["marker_strength"]
        .max()
        .sort_values(["cell_type", "marker_strength"], ascending=[True, False])
    )
    agg["source"] = "CVAE186_aggregated"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    agg.to_csv(OUT, sep="\t", index=False)
    print(f"Wrote {OUT} rows={len(agg)}")
    print("Cell types:", agg["cell_type"].value_counts().to_string())


if __name__ == "__main__":
    main()

