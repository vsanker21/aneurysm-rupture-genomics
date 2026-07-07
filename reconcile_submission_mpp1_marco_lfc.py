"""Reconcile MPP1/MARCO log2FC definitions and update submission supplementary tables."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from docx import Document
from docx.shared import Pt

BASE = Path(r"e:\Aneurysms Genomics Project")
SUB = BASE / "Latest Submission Files"
TAB_DIR = BASE / "results" / "supplementary" / "tables"

FOOTNOTE_S2 = (
    "Footnote (effect-size definition). For GSE122897 discovery, DESeq2_log2FC is the primary rupture "
    "effect estimate (negative-binomial model on raw counts; n = 22 ruptured vs 21 unruptured). "
    "Examples: MPP1 +1.87, MARCO +3.04. These values are used for shortlisting, panel signing, and "
    "Table 3 meta-analysis inputs. They differ numerically from the cohort-native Welch contrasts in "
    "Supplementary Table S10 (see that table's footnote)."
)

FOOTNOTE_S10 = (
    "Footnote (effect-size definition). Log2FC_positive_minus_negative (renamed in the archived TSV as "
    "Delta_mean_log2_expression) is the Welch-group contrast on already log-transformed expression "
    "native to each cohort (GSE122897: difference of mean log2[TPM+1] between ruptured and unruptured; "
    "MPP1 +1.17, MARCO +2.25). This is not DESeq2 log2FoldChange and is not log2 of the ratio of "
    "arithmetic TPM means. The column GSE122897_DESeq2_log2FC_reference repeats the discovery DESeq2 "
    "estimates (MPP1 +1.87, MARCO +3.04) for cross-cohort direction matching. Use DESeq2 values when "
    "citing discovery magnitude; use Delta_mean_log2_expression when comparing effect sizes on a common "
    "log-expression scale across external cohorts."
)

RECONCILIATION = """MPP1 / MARCO GSE122897 fold-change reconciliation (confirmed)
==========================================================

Question
--------
Table 3 / Supplementary Table S2 report MPP1 log2FC +1.87 and MARCO +3.04, whereas
Supplementary Table S10 shows +1.17 and +2.25 for the GSE122897 row. Are these inconsistent?

Answer: No — different estimators, both correct.
-------------------------------------------------

| Gene  | DESeq2 log2FC (Tables S2, S3 meta, shortlist) | Delta mean log2 expression (Table S10 GSE122897 row) |
|-------|-----------------------------------------------|------------------------------------------------------|
| MPP1  | +1.873 (padj = 4.3e-05)                       | +1.174                                               |
| MARCO | +3.044 (padj = 8.9e-07)                       | +2.246                                               |

1. DESeq2 log2FoldChange (Tables S2, S3, shortlist_evidence.tsv)
   - Source: results/gse122897_deseq2_results.tsv
   - Model: DESeq2 negative binomial on raw counts
   - Use for: discovery ranking, 300-gene shortlist, 8-gene panel signing, three-cohort meta-analysis

2. Delta_mean_log2_expression (Supplementary Table S10, column formerly Log2FC_positive_minus_negative)
   - Source: mpp1_marco_multi_cohort_expression_survey.py via compute_gene_level_welch
   - Formula: mean(log2[TPM+1])_ruptured − mean(log2[TPM+1])_unruptured
   - Use for: cross-cohort survey with harmonised log-expression contrasts (arrays, blood, tissue)

Why they differ
---------------
DESeq2 estimates log2 fold change on counts with shrinkage and dispersion modelling.
The survey statistic is a simple difference of group means on the log2(TPM+1) scale.
For low-to-moderate expression genes these estimators need not coincide; both show
higher expression in ruptured tissue for MPP1 and MARCO.

Action taken
------------
- Supplementary Table S10 TSV regenerated with explicit column names and DESeq2 reference column.
- Footnotes added to Supplementary_Tables_Revised.docx (Tables S2 and S10).
- This reconciliation note saved alongside final submission files.
"""


def build_table_s10_reconciled() -> pd.DataFrame:
    raw = pd.read_csv(BASE / "results/mpp1_marco_survey/mpp1_marco_all_cohort_results.tsv", sep="\t")
    sl = pd.read_csv(BASE / "results/phases_g/shortlist_evidence.tsv", sep="\t")
    deseq_ref = sl.set_index("symbol")["DESeq2_log2FC"].to_dict()

    out = raw.copy()
    out = out.rename(columns={
        "Log2FC_positive_minus_negative": "Delta_mean_log2_expression",
    })
    out["GSE122897_DESeq2_log2FC_reference"] = out["Gene"].map(deseq_ref)
    front = [
        "Gene", "Cohort", "Comparison", "Biospecimen", "Matrix", "Evidence_tier", "Status",
        "N_positive", "N_negative", "Delta_mean_log2_expression",
        "GSE122897_DESeq2_log2FC_reference", "P_value", "T_statistic",
        "Direction_match_discovery", "Interpretation_note",
    ]
    front = [c for c in front if c in out.columns]
    rest = [c for c in out.columns if c not in front]
    out = out[front + rest]
    return out


def _insert_after_heading(doc: Document, heading_prefix: str, footnote: str) -> bool:
    for i, para in enumerate(doc.paragraphs):
        if para.text.strip().startswith(heading_prefix):
            # insert after next non-empty paragraph (legend line)
            j = i + 1
            while j < len(doc.paragraphs) and not doc.paragraphs[j].text.strip():
                j += 1
            if j < len(doc.paragraphs) and doc.paragraphs[j].text.strip().startswith("Footnote (effect-size"):
                doc.paragraphs[j].text = footnote
                return True
            new_p = doc.paragraphs[j]._element
            p = doc.paragraphs[i]._p
            from docx.oxml import OxmlElement
            from docx.text.paragraph import Paragraph

            new_el = OxmlElement("w:p")
            p.addnext(new_el)
            np = Paragraph(new_el, doc.paragraphs[i]._parent)
            np.text = footnote
            for run in np.runs:
                run.font.size = Pt(10)
                run.italic = True
            return True
    return False


def patch_supplementary_tables_docx() -> None:
    path = SUB / "Supplementary_Tables_Revised.docx"
    doc = Document(str(path))
    ok2 = _insert_after_heading(doc, "Supplementary Table S2", FOOTNOTE_S2)
    ok10 = _insert_after_heading(doc, "Supplementary Table S10", FOOTNOTE_S10)
    doc.save(str(path))
    print(f"Patched {path.name}: S2 footnote={ok2}, S10 footnote={ok10}")


FOOTNOTE_TABLE3 = (
    "Footnote (GSE122897 effect sizes in Table 3). Log2FC_GSE122897 columns for the discovery cohort "
    "are DESeq2 estimates from raw counts (e.g. MPP1 +1.87, MARCO +3.04). Supplementary Table S10 "
    "reports separate delta mean log2(TPM+1) contrasts (+1.17 and +2.25) for cross-cohort survey only."
)


def patch_tables_revised_docx() -> None:
    path = SUB / "Tables_Revised.docx"
    if not path.exists():
        return
    doc = Document(str(path))
    for para in doc.paragraphs:
        if para.text.strip().startswith("Table 3"):
            if not any("Footnote (GSE122897 effect sizes" in p.text for p in doc.paragraphs):
                doc.add_paragraph(FOOTNOTE_TABLE3)
            break
    doc.save(str(path))
    print(f"Patched {path.name} (Table 3 footnote if applicable)")


def patch_checklist() -> None:
    path = SUB / "Submission_Readiness_Checklist.docx"
    if not path.exists():
        return
    doc = Document(str(path))
    resolved = (
        "RESOLVED — MPP1/MARCO GSE122897 log2FC: Tables S2/S3 use DESeq2 LFC (MPP1 +1.87, MARCO +3.04); "
        "Table S10 uses delta mean log2(TPM+1) for cross-cohort survey (MPP1 +1.17, MARCO +2.25). "
        "Footnotes added; see MPP1_MARCO_LFC_Reconciliation.txt."
    )
    doc.add_paragraph(resolved)
    doc.save(str(path))
    print(f"Updated {path.name}")


def main() -> None:
    out = build_table_s10_reconciled()
    tsv = TAB_DIR / "Supplementary_Table_S10_mpp1_marco_all_cohorts.tsv"
    out.to_csv(tsv, sep="\t", index=False)
    print(f"Wrote reconciled {tsv}")

    (SUB / "MPP1_MARCO_LFC_Reconciliation.txt").write_text(RECONCILIATION, encoding="utf-8")
    print(f"Wrote {SUB / 'MPP1_MARCO_LFC_Reconciliation.txt'}")

    patch_supplementary_tables_docx()
    patch_tables_revised_docx()
    patch_checklist()

    # refresh generate_supplementary_materials_complete legend if present
    print("\nConfirmed values:")
    g = out[(out["Cohort"] == "GSE122897") & (out["Gene"].isin(["MPP1", "MARCO"]))]
    print(g[["Gene", "Delta_mean_log2_expression", "GSE122897_DESeq2_log2FC_reference"]].to_string(index=False))


if __name__ == "__main__":
    main()
