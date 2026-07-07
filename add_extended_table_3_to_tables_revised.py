"""
Add or replace Extended Data Table 3 (full 300-gene shortlist) in Tables_Revised.docx.

Embeds all 300 genes in the Word document; also writes the matching TSV to
Latest Submission Files for deposit alongside the manuscript.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from docx import Document
from docx.shared import Pt

BASE = Path(r"e:\Aneurysms Genomics Project")
SUB = BASE / "Latest Submission Files"
TABLES_REVISED = SUB / "Tables_Revised.docx"
SOURCE_TSV = BASE / "results/supplementary/tables/Supplementary_Table_S2_300gene_shortlist.tsv"
OUTPUT_TSV = SUB / "Extended_Data_Table_3_300gene_shortlist.tsv"

HEADING = "Extended Data Table 3 | Full 300-gene candidate shortlist"

LEGEND = (
    "Ranked shortlist of 300 genes from GSE122897 meeting joint discovery criteria: "
    "DESeq2 adjusted P ≤ 0.05 and nominal immune-adjusted rupture P ≤ 0.05 "
    "(LM22-adjusted OLS on log2[TPM+1]). Genes are ranked by blood-bridge Stouffer Z "
    "across GSE159610 and GSE106520. The eight-gene blood panel (In_8gene_panel) and "
    "six-gene wall rupture core (In_6gene_wall_core) are flagged. NR = gene not represented "
    "on the blood array/RNA-seq platform for that cohort. Primary discovery effect sizes "
    "are DESeq2_log2FC (e.g. MPP1 +1.87, MARCO +3.04). An identical tabular export is "
    f"provided as {OUTPUT_TSV.name}."
)

DISPLAY_COLUMNS = [
    "Shortlist_rank",
    "Gene_symbol",
    "DESeq2_log2FC",
    "DESeq2_padj",
    "LM_immune_adj_P_Ruptured",
    "GSE159610_Log2FC_aneurysm_vs_control",
    "GSE159610_P",
    "GSE106520_log2FC",
    "GSE106520_P",
    "Blood_both_direction_concordant",
    "Bridge_Stouffer_Z",
    "Bridge_meta_P",
    "In_8gene_panel",
    "In_6gene_wall_core",
]


def _ranked_evidence_lookup() -> pd.DataFrame:
    """Rebuild Shortlist_rank → symbol / Gene_ID mapping (same sort as Table S2)."""
    from scipy.stats import norm

    sl = pd.read_csv(BASE / "results/phases_g/shortlist_evidence.tsv", sep="\t")
    g106 = pd.read_csv(
        BASE / "results/external_validation_gse106520/gse106520_shortlist_overlap.tsv",
        sep="\t",
    )
    g106 = g106[["Gene", "P_value"]].rename(columns={"Gene": "Gene_symbol", "P_value": "GSE106520_P"})
    sl = sl.rename(columns={"symbol": "Gene_symbol"}) if "symbol" in sl.columns else sl
    symbol_col = "symbol" if "symbol" in sl.columns else "Gene_symbol"
    df = sl.merge(g106, on="Gene_symbol", how="left")
    p1 = pd.to_numeric(df.get("GSE159610_P"), errors="coerce").fillna(1.0).clip(1e-300, 1.0)
    p2 = pd.to_numeric(df["GSE106520_P"], errors="coerce").fillna(1.0).clip(1e-300, 1.0)
    df["Bridge_Stouffer_Z"] = (norm.isf(p1 / 2) + norm.isf(p2 / 2)) / (2**0.5)
    df = df.sort_values(["Bridge_Stouffer_Z", "DESeq2_padj"], ascending=[False, True]).reset_index(drop=True)
    df["Shortlist_rank"] = range(1, len(df) + 1)
    return df[["Shortlist_rank", "Gene_ID", symbol_col]].rename(columns={symbol_col: "symbol"})


def _prepare_display_df(df: pd.DataFrame) -> pd.DataFrame:
    """Backfill missing gene symbols using rank-matched evidence IDs."""
    out = df.copy()
    ref = _ranked_evidence_lookup()
    merged = out.merge(ref, on="Shortlist_rank", how="left", suffixes=("", "_ref"))
    out["Gene_symbol"] = merged["Gene_symbol"].fillna(merged["symbol"])
    entrez = "ENTREZ:" + merged["Gene_ID"].astype(str)
    out["Gene_symbol"] = out["Gene_symbol"].fillna(entrez)
    return out


def _ensure_source() -> pd.DataFrame:
    if not SOURCE_TSV.exists():
        raise FileNotFoundError(
            f"Missing shortlist TSV. Run: python generate_supplementary_materials_complete.py\n  {SOURCE_TSV}"
        )
    return _prepare_display_df(pd.read_csv(SOURCE_TSV, sep="\t"))


def _copy_supplementary_file(df: pd.DataFrame) -> None:
    SUB.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_TSV, sep="\t", index=False)
    print(f"Wrote {OUTPUT_TSV} ({len(df)} genes)")


def _delete_paragraph(paragraph) -> None:
    element = paragraph._element
    element.getparent().remove(element)


def _delete_table(table) -> None:
    element = table._element
    element.getparent().remove(element)


def _remove_existing_edt3(doc: Document) -> bool:
    """Remove prior Extended Data Table 3 section (heading, legend, manifest, preview)."""
    start_idx = None
    for i, para in enumerate(doc.paragraphs):
        if para.text.strip().startswith("Extended Data Table 3"):
            start_idx = i
            break
    if start_idx is None:
        return False

    # Remove paragraphs from EDT3 heading through end of document.
    for para in list(doc.paragraphs[start_idx:]):
        _delete_paragraph(para)

    # Remove manifest + preview tables (last two tables added for EDT3).
    while len(doc.tables) > 3:
        _delete_table(doc.tables[-1])

    return True


def _format_cell(val, col: str = "") -> str:
    bridge_cols = {
        "GSE159610_Log2FC_aneurysm_vs_control",
        "GSE159610_P",
        "GSE106520_log2FC",
        "GSE106520_P",
    }
    if col in bridge_cols and (pd.isna(val) or val == ""):
        return "NR"
    if pd.isna(val):
        return ""
    if isinstance(val, bool):
        return str(val)
    if isinstance(val, (int,)):
        return str(val)
    if isinstance(val, float):
        if abs(val) < 1e-3 or abs(val) >= 1000:
            return f"{val:.4g}"
        return f"{val:.4g}"
    return str(val)


def _add_full_gene_table(doc: Document, df: pd.DataFrame) -> None:
    cols = [c for c in DISPLAY_COLUMNS if c in df.columns]
    table = doc.add_table(rows=1, cols=len(cols))
    table.style = "Table"

    for j, col in enumerate(cols):
        cell = table.rows[0].cells[j]
        cell.text = col
        if cell.paragraphs[0].runs:
            cell.paragraphs[0].runs[0].bold = True

    for _, row in df[cols].iterrows():
        cells = table.add_row().cells
        for j, col in enumerate(cols):
            cells[j].text = _format_cell(row[col], col)


def patch_tables_revised(df: pd.DataFrame) -> None:
    if not TABLES_REVISED.exists():
        raise FileNotFoundError(TABLES_REVISED)

    doc = Document(str(TABLES_REVISED))
    replaced = _remove_existing_edt3(doc)
    if replaced:
        print(f"Removed prior Extended Data Table 3 section from {TABLES_REVISED.name}")

    doc.add_paragraph()
    heading = doc.add_paragraph(HEADING)
    if heading.runs:
        heading.runs[0].bold = True

    legend = doc.add_paragraph(LEGEND)
    for run in legend.runs:
        run.font.size = Pt(10)

    _add_full_gene_table(doc, df)

    doc.save(str(TABLES_REVISED))
    print(f"Updated {TABLES_REVISED} with full {len(df)}-gene table")


def main() -> None:
    df = _ensure_source()
    if len(df) != 300:
        print(f"Warning: expected 300 genes, found {len(df)}")
    _copy_supplementary_file(df)
    patch_tables_revised(df)
    print("Done.")


if __name__ == "__main__":
    main()
