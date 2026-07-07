"""
Generate Nature Medicine correspondence Word document with embedded figures and tables.
Run: python create_nature_medicine_correspondence.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt, RGBColor
from docx.enum.table import WD_TABLE_ALIGNMENT

BASE = Path(__file__).resolve().parent
OUT_DOCX = BASE / "manuscripts" / "Nature_Medicine_Correspondence_Concordance_First_Panel_Screening.docx"
FIG_DIR = BASE / "results" / "figures" / "correspondence_nm"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_DOCX.parent.mkdir(parents=True, exist_ok=True)

TISSUE_CORE = ["TGFBI", "ANPEP", "MPP1", "MARCO", "CLEC5A", "CPXM2"]
BLOOD_BRIDGE = ["DYNLT1", "OASL", "MPP1", "TMEM130", "MARCO", "MCOLN1", "RNASET2", "ANK3"]
OVERLAP = sorted(set(TISSUE_CORE) & set(BLOOD_BRIDGE))


def _load_data() -> dict:
    transfer = pd.read_csv(
        BASE / "results/advanced_validation/external_signature_transfer_summary.tsv", sep="\t"
    )
    meta_conc = pd.read_csv(
        BASE / "results/rupture_meta_analysis/rupture_consensus_meta_concordant.tsv", sep="\t"
    )
    meta_summary = json.loads(
        (BASE / "results/rupture_meta_analysis/rupture_consensus_meta_summary.json").read_text(
            encoding="utf-8"
        )
    )
    blood_panel = pd.read_csv(
        BASE / "results/blood_bridge_panel/blood_bridge_panel_gene_summary.tsv", sep="\t"
    )
    phase7 = pd.read_csv(BASE / "results/phase7_optimized_model_metrics.tsv", sep="\t")
    nested = pd.read_csv(BASE / "results/gse122897_nested_de_cv_summary.tsv", sep="\t")
    mpp_marco = pd.read_csv(BASE / "results/mpp1_marco_survey/mpp1_marco_stouffer_meta.tsv", sep="\t")
    return {
        "transfer": transfer,
        "meta_conc": meta_conc,
        "meta_summary": meta_summary,
        "blood_panel": blood_panel,
        "phase7": phase7,
        "nested": nested,
        "mpp_marco": mpp_marco,
    }


def _fig1_workflow(out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    boxes = [
        (0.3, 4.2, 2.2, 1.1, "GSE122897\nDiscovery tissue\nn=43 (22R/21U)", "#e8f4fc"),
        (3.0, 4.2, 2.2, 1.1, "3-score model\n(DE + EC + ECM)", "#fff3cd"),
        (5.7, 4.2, 2.0, 1.1, "External transfer\nGSE13353 / GSE54083", "#fde8e8"),
        (8.0, 4.2, 1.7, 1.1, "Unstable\nAUROC", "#f5c6cb"),
        (0.3, 2.3, 2.8, 1.2, "Triple-cohort concordance\n16 genes overlap → 6 pass", "#d4edda"),
        (3.5, 2.3, 2.6, 1.2, "6-gene tissue core\nTGFBI ANPEP MPP1\nMARCO CLEC5A CPXM2", "#c3e6cb"),
        (6.5, 2.3, 2.8, 1.2, "8-gene blood bridge\n(aneurysm vs control)", "#cce5ff"),
        (3.5, 0.5, 3.2, 1.0, "No public blood\nRvU → CIAP / ROAR", "#e2e3e5"),
    ]
    for x, y, w, h, text, color in boxes:
        rect = FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
            linewidth=1.2, edgecolor="#333", facecolor=color, alpha=0.95,
        )
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9, wrap=True)

    arrows = [
        ((2.5, 4.75), (3.0, 4.75)),
        ((5.2, 4.75), (5.7, 4.75)),
        ((7.7, 4.75), (8.0, 4.75)),
        ((1.7, 4.2), (1.7, 3.5)),
        ((4.5, 4.2), (4.8, 3.5)),
        ((5.0, 2.3), (5.0, 1.5)),
        ((7.5, 2.9), (6.5, 2.9)),
    ]
    for (x0, y0), (x1, y1) in arrows:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", color="#444", lw=1.5))
    ax.text(5.0, 5.65, "Concordance-first aneurysm gene-panel screening framework",
            ha="center", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _fig2_transfer_auroc(transfer: pd.DataFrame, out_path: Path) -> None:
    direct = transfer[transfer["Role"].str.contains("direct", case=False, na=False)].copy()
    labels = direct["Dataset"].tolist()
    aucs = direct["AUC"].astype(float).tolist()
    colors = ["#2e7d32" if a >= 0.5 else "#c62828" for a in aucs]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axes[0]
    y = np.arange(len(labels))
    ax.barh(y, aucs, color=colors, edgecolor="#333", height=0.55)
    ax.axvline(0.5, color="gray", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("AUROC (fixed 3-score model, no retraining)")
    ax.set_xlim(0, 1)
    ax.set_title("A. Cohort-level classifier transfer")
    for i, v in enumerate(aucs):
        ax.text(v + 0.02, i, f"{v:.3f}", va="center", fontsize=9)

    ax2 = axes[1]
    genes = TISSUE_CORE
    meta = pd.read_csv(BASE / "results/rupture_meta_analysis/rupture_consensus_meta_concordant.tsv", sep="\t")
    lfc = [meta.loc[meta["Gene"] == g, "Log2FC_GSE122897"].iloc[0] for g in genes]
    dirs = ["+" if x > 0 else "−" for x in lfc]
    colors2 = ["#b54a4a" if d == "+" else "#4b7bb8" for d in dirs]
    ax2.barh(genes, [abs(x) for x in lfc], color=colors2, edgecolor="#333")
    ax2.set_xlabel("|log2FC| GSE122897 (direction concordant in all 3 cohorts)")
    ax2.set_title("B. Six-gene tissue core (stable direction)")
    fig.suptitle("Scores fail portability; sparse genes replicate", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _fig3_consensus_meta(meta_conc: pd.DataFrame, out_path: Path) -> None:
    df = meta_conc.copy()
    df = df.sort_values("Meta_Z", ascending=True)
    colors = np.where(df["Consensus_direction"] == "Up in ruptured", "#b54a4a", "#4b7bb8")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.barh(df["Gene"], df["Meta_Z"], color=colors, edgecolor="#333")
    ax.axvline(0, color="black", lw=1)
    ax.set_xlabel("Weighted Stouffer meta-Z (GSE122897 + GSE54083 + GSE13353)")
    ax.set_title("Six-gene rupture tissue core — complete triple-concordance set")
    for i, row in df.iterrows():
        ax.text(row["Meta_Z"], row["Gene"], f"  P={row['Meta_P']:.1e}", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _fig4_blood_bridge(out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    cohorts = ["GSE159610", "GSE106520"]
    aucs = [0.737, 0.852]
    axes[0].bar(cohorts, aucs, color=["#5b8db8", "#3d6e9e"], edgecolor="#333")
    axes[0].axhline(0.5, color="gray", ls="--")
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("AUROC")
    axes[0].set_title("A. 8-gene blood bridge panel\n(aneurysm vs control; discovery-signed)")
    for i, v in enumerate(aucs):
        axes[0].text(i, v + 0.03, f"{v:.3f}", ha="center", fontweight="bold")

    y = np.arange(len(BLOOD_BRIDGE))
    in_core = [1 if g in TISSUE_CORE else 0 for g in BLOOD_BRIDGE]
    colors = ["#2e7d32" if c else "#9e9e9e" for c in in_core]
    axes[1].barh(y, [1] * len(BLOOD_BRIDGE), color=colors, edgecolor="#333")
    axes[1].set_yticks(y)
    axes[1].set_yticklabels(BLOOD_BRIDGE)
    axes[1].set_xlim(0, 1.2)
    axes[1].set_xticks([])
    axes[1].set_title("B. Blood bridge genes (green = tissue core overlap)")
    axes[1].text(0.5, -0.8, "Overlap: MPP1, MARCO", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _fig5_panel_overlap(out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    c1 = plt.Circle((3.2, 3), 1.8, color="#c3e6cb", ec="#333", lw=1.5, alpha=0.85)
    c2 = plt.Circle((6.8, 3), 1.8, color="#cce5ff", ec="#333", lw=1.5, alpha=0.85)
    ax.add_patch(c1)
    ax.add_patch(c2)
    ax.text(2.0, 3, "6-gene\ntissue core", ha="center", va="center", fontsize=10)
    ax.text(8.0, 3, "8-gene\nblood bridge", ha="center", va="center", fontsize=10)
    ax.text(5.0, 3, "MPP1\nMARCO", ha="center", va="center", fontsize=11, fontweight="bold")
    ax.set_title("Dual-panel architecture for aneurysm screening", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def add_table(doc: Document, df: pd.DataFrame, title: str, note: str | None = None, max_rows: int | None = None) -> None:
    doc.add_heading(title, level=3)
    if note:
        p = doc.add_paragraph(note)
        p.runs[0].italic = True if p.runs else None
    show = df.head(max_rows) if max_rows else df
    table = doc.add_table(rows=1, cols=len(show.columns))
    table.style = "Light Grid Accent 1"
    hdr = table.rows[0].cells
    for i, col in enumerate(show.columns):
        hdr[i].text = str(col)
        if hdr[i].paragraphs[0].runs:
            hdr[i].paragraphs[0].runs[0].bold = True
    for _, row in show.iterrows():
        cells = table.add_row().cells
        for i, val in enumerate(row):
            if isinstance(val, float):
                if abs(val) < 0.001 or abs(val) > 1e4:
                    cells[i].text = f"{val:.2e}"
                else:
                    cells[i].text = f"{val:.4f}"
            else:
                cells[i].text = str(val) if pd.notna(val) else ""
    doc.add_paragraph()


def add_figure(doc: Document, path: Path, caption: str, width: float = 6.2) -> None:
    if path.is_file():
        doc.add_picture(str(path), width=Inches(width))
        cap = doc.add_paragraph(caption)
        cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
        if cap.runs:
            cap.runs[0].italic = True
            cap.runs[0].font.size = Pt(10)
    else:
        doc.add_paragraph(f"[Figure missing: {path}]")
    doc.add_paragraph()


def build_document(data: dict, figures: dict[str, Path]) -> Document:
    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "Calibri"
    style.font.size = Pt(11)

    # Title page
    t = doc.add_heading(
        "Concordance-first molecular screening for intracranial aneurysm rupture: "
        "when gene panels should replace portable risk scores",
        0,
    )
    t.alignment = WD_ALIGN_PARAGRAPH.CENTER
    sub = doc.add_paragraph("Correspondence — draft for Nature Medicine")
    sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
    sub.runs[0].italic = True
    sub.runs[0].font.size = Pt(12)
    doc.add_paragraph()
    for line in [
        "[Author list — to be completed]",
        "[Affiliations — to be completed]",
        "Corresponding author: [email]",
        "Word count (main text): ~950",
    ]:
        p = doc.add_paragraph(line)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.add_page_break()

    doc.add_heading("Main text", 1)
    ms = data["meta_summary"]
    n_conc = int(ms.get("n_directionally_concordant_all_three", 6))
    n_overlap = int(ms.get("n_overlap_all_three", 16))
    auc_disc = float(data["phase7"].loc[data["phase7"]["Metric"] == "ROC_AUC", "Optimized_Model"].iloc[0])
    auc_nested = float(data["nested"].loc[data["nested"]["Metric"] == "Nested_CV_AUC_Mean", "Value"].iloc[0])
    auc_nested_sd = float(data["nested"].loc[data["nested"]["Metric"] == "Nested_CV_AUC_Std", "Value"].iloc[0])

    paragraphs = [
        (
            "Intracranial aneurysm management still hinges on size, location, and morphology—proxies that "
            "explain only part of rupture risk. Molecular screening has been proposed for decades, yet clinical "
            "adoption remains sparse. A systematic reanalysis of public rupture transcriptomics suggests why: "
            "the field has optimized the wrong object. Portable multi-gene risk scores trained on small tissue "
            "cohorts do not travel reliably across studies, whereas sparse, directionally concordant gene panels—"
            "anchored in vessel-wall biology and selectively readable in blood—may be the unit of translation "
            "worth prospectively testing."
        ),
        (
            f"We integrated the three largest directly comparable human datasets of ruptured versus unruptured "
            f"aneurysm wall transcriptomes (GSE122897, n=43; GSE54083, n=13; GSE13353, n=19) with orthogonal "
            f"blood and genetics layers. Discovery in GSE122897 yielded a composite rupture model combining "
            f"differential-expression, endothelial-dysfunction, and extracellular-matrix program scores (same-cohort "
            f"AUROC {auc_disc:.2f}; nested five-fold cross-validation AUROC {auc_nested:.2f} ± {auc_nested_sd:.2f}). "
            f"When this fixed three-score classifier was projected externally without retraining, performance was "
            f"cohort-dependent: supportive in GSE13353 (AUROC 0.86) but inverted in GSE54083 (AUROC 0.18)—despite "
            f"overlapping biology (Figure 1; Table 2). This instability is not a technical nuisance; it defines what "
            f"public omics can and cannot claim about rupture prediction."
        ),
        (
            f"A gene-level analysis was more informative. Among only {n_overlap} genes measurable across all three "
            f"rupture tissue cohorts, {n_conc} showed identical directional effects: TGFBI, ANPEP, MPP1, MARCO, and "
            f"CLEC5A higher in ruptured wall, and CPXM2 lower (Stouffer meta P from ~10⁻⁸ to ~10⁻¹²; Figure 2; "
            f"Table 3). These six are not arbitrary top hits from a single differential-expression table; they are "
            f"the complete set of cross-platform rupture signals reproducible under the strictest public-data test "
            f"available today. Rupture biology in bulk tissue is low-dimensional and sparse, not a high-dimensional "
            f"classifier problem waiting for more parameters."
        ),
        (
            "Two genes—MPP1 and MARCO—emerge as bridge molecules. In tier-1 tissue meta-analysis across direct "
            "rupture cohorts, both remain directionally concordant (meta P ~10⁻⁶–10⁻⁷; Table 5). A separate eight-gene "
            "blood bridge panel (DYNLT1, OASL, MPP1, TMEM130, MARCO, MCOLN1, RNASET2, ANK3) was locked using "
            "discovery tissue sign and tested only in aneurysm-versus-control blood RNA-seq (GSE159610, GSE106520), "
            "where composite AUROCs reached ~0.74 and ~0.85 (Figure 3; Table 4). MPP1 and MARCO are the overlap "
            "between wall rupture core and blood accessibility. An exhaustive survey of public repositories identified "
            "no human blood mRNA cohort with ruptured versus unruptured intracranial aneurysm at baseline; blood omics "
            "therefore cannot validate prospective rupture risk today (Table 1)."
        ),
        (
            "Parallel integrative work retains 645 rupture-associated genes after adjustment for inferred immune "
            "infiltration in GSE122897, with blood cis-eQTL and intracranial aneurysm GWAS rsID overlap for genetic "
            "plausibility—composition-adjusted discovery, not yet transcriptome-wide FDR-robust and not formal "
            "colocalization. The convergent lesson is structural: composition-adjusted discovery → cross-cohort "
            "concordance filter → minimal blood bridge → genetics-informed prioritization → prospective outcome "
            "validation."
        ),
        (
            "We propose concordance-first panel screening for trials: (1) separate three intents—aneurysm detection, "
            "rupture biology in tissue, and prospective rupture risk; (2) pre-register minimal gene sets and evaluate "
            "directional replication before optimizing composite scores; (3) require incremental utility over PHASES on "
            "prospective pre-rupture blood (CIAP, ROAR-class cohorts; Table 7); (4) invest in human aneurysm single-cell "
            "and spatial profiling to assign genes to wall compartments and test causality."
        ),
        (
            "If validated prospectively, concordance-first screening would identify a biologically defined "
            "high-intrinsic-risk wall state partially legible in blood—an interpretable complement to diameter-based "
            "rules. The opportunity is not another machine-learning benchmark but a design principle: in rupture omics, "
            "replication of direction beats optimization of discrimination. Gene-panel screening should be built on "
            "what reproduces across tissues and platforms—and tested where public science ends, in blood drawn before "
            "the bleed."
        ),
    ]
    for text in paragraphs:
        doc.add_paragraph(text)

    doc.add_page_break()
    doc.add_heading("Figures", 1)
    add_figure(
        doc, figures["workflow"],
        "Figure 1 | Concordance-first framework for intracranial aneurysm gene-panel screening. "
        "Discovery tissue omics feed a three-score classifier with unstable external transfer; "
        "triple-cohort directional concordance yields a six-gene tissue core and eight-gene blood bridge; "
        "prospective rupture validation requires collaboration cohorts (CIAP, ROAR).",
        width=6.5,
    )
    add_figure(
        doc, figures["transfer"],
        "Figure 2 | Cohort-level classifier portability versus six-gene directional stability. "
        "(A) Fixed three-score model AUROC by direct rupture cohort. "
        "(B) Six-gene tissue core log2 fold-changes in GSE122897 (concordant across all three cohorts).",
        width=6.5,
    )
    add_figure(
        doc, figures["consensus"],
        "Figure 3 | Six-gene rupture tissue core — weighted Stouffer meta-Z across GSE122897, GSE54083, and GSE13353.",
        width=5.5,
    )
    add_figure(
        doc, figures["blood"],
        "Figure 4 | Blood bridge panel for aneurysm versus control context. "
        "(A) Composite panel AUROC in GSE159610 and GSE106520 (discovery-signed; not rupture-status validation). "
        "(B) Eight bridge genes; MPP1 and MARCO overlap the tissue core.",
        width=6.5,
    )
    add_figure(
        doc, figures["overlap"],
        "Figure 5 | Dual-panel architecture: six-gene tissue rupture core and eight-gene blood bridge with MPP1/MARCO overlap.",
        width=5.0,
    )

    doc.add_page_break()
    doc.add_heading("Tables", 1)

    tbl1 = pd.DataFrame([
        ["Aneurysm detection (blood vs control)", "Supported (bridge panel)", "GSE159610, GSE106520", "Not rupture-specific"],
        ["Rupture biology (tissue R vs U)", "Supported (six-gene core)", "GSE122897, GSE54083, GSE13353", "Retrospective labels"],
        ["Prospective rupture risk (pre-bleed blood)", "Not in public data", "None identified", "CIAP-1/2, ROAR required"],
    ], columns=["Screening intent", "Public-data status", "Example cohorts", "Interpretation"])
    add_table(doc, tbl1, "Table 1 | Three screening intents and public-data support")

    t2 = data["transfer"][
        ["Dataset", "Role", "Samples", "AUC", "AUC_bootstrap_ci_low", "AUC_bootstrap_ci_high", "AUC_permutation_pvalue"]
    ].copy()
    t2["AUC"] = t2["AUC"].map(lambda x: f"{float(x):.3f}")
    t2["AUC_bootstrap_ci_low"] = t2["AUC_bootstrap_ci_low"].map(lambda x: f"{float(x):.3f}")
    t2["AUC_bootstrap_ci_high"] = t2["AUC_bootstrap_ci_high"].map(lambda x: f"{float(x):.3f}")
    t2["AUC_permutation_pvalue"] = t2["AUC_permutation_pvalue"].map(lambda x: f"{float(x):.4f}")
    add_table(doc, t2, "Table 2 | External transfer of fixed three-score rupture model (no retraining)")

    t3 = data["meta_conc"][
        ["Gene", "Consensus_direction", "Meta_Z", "Meta_P", "Meta_FDR_BH",
         "Log2FC_GSE122897", "Log2FC_GSE54083", "Log2FC_GSE13353"]
    ].copy()
    t3["Meta_Z"] = t3["Meta_Z"].map(lambda x: f"{float(x):.3f}")
    t3["Meta_P"] = t3["Meta_P"].map(lambda x: f"{float(x):.2e}")
    t3["Meta_FDR_BH"] = t3["Meta_FDR_BH"].map(lambda x: f"{float(x):.2e}")
    for c in ["Log2FC_GSE122897", "Log2FC_GSE54083", "Log2FC_GSE13353"]:
        t3[c] = t3[c].map(lambda x: f"{float(x):.3f}")
    add_table(doc, t3, "Table 3 | Six-gene tissue rupture core (complete triple-concordance set)")

    t4 = data["blood_panel"][
        ["Panel_rank", "Gene", "DESeq2_log2FC", "DESeq2_padj",
         "GSE159610_Log2FC_aneurysm_vs_control", "GSE106520_Log2FC_aneurysm_vs_control", "Bridge_meta_P"]
    ].copy()
    t4["DESeq2_log2FC"] = t4["DESeq2_log2FC"].map(lambda x: f"{float(x):.3f}")
    t4["DESeq2_padj"] = t4["DESeq2_padj"].map(lambda x: f"{float(x):.2e}")
    t4["Bridge_meta_P"] = t4["Bridge_meta_P"].map(lambda x: f"{float(x):.3f}")
    add_table(doc, t4, "Table 4 | Eight-gene blood bridge panel (discovery-signed; aneurysm vs control)")

    t5 = pd.DataFrame([
        ["Discovery cohort (GSE122897) same-sample AUROC", f"{auc_disc:.3f}"],
        ["Nested CV AUROC (mean ± SD)", f"{auc_nested:.3f} ± {auc_nested_sd:.3f}"],
        ["Genes overlapping three tissue cohorts", str(n_overlap)],
        ["Genes directionally concordant (all three)", str(n_conc)],
        ["Immune-adjusted anchored genes (PHASES-G)", "645"],
        ["Public human blood mRNA R vs U cohorts", "0"],
    ], columns=["Metric", "Value"])
    add_table(doc, t5, "Table 5 | Discovery performance versus concordance filter summary")

    mm = data["mpp_marco"]
    t6 = mm[mm["Gene"].isin(["MPP1", "MARCO"])][
        ["Gene", "Meta_label", "Meta_Z", "Meta_P", "N_cohorts", "Fraction_concordant"]
    ].copy()
    t6["Meta_Z"] = t6["Meta_Z"].map(lambda x: f"{float(x):.3f}")
    t6["Meta_P"] = t6["Meta_P"].map(lambda x: f"{float(x):.2e}")
    add_table(doc, t6, "Table 6 | MPP1 and MARCO — tissue and blood meta-analysis layers")

    t7_cohorts = pd.DataFrame([
        ["CIAP-1", "very_high", "Prospective UIA + blood banking", "Contact investigators; pre-specify locked panel"],
        ["CIAP-2", "very_high", "Biomarker-framed rupture prediction", "Priority collaboration for validation"],
        ["ROAR", "high", "Rupture outcomes; biospecimen TBD", "Confirm PAXgene/RNA availability"],
        ["ISUIA", "moderate", "Natural history; no public biospecimen", "Investigator confirmation only"],
    ], columns=["Cohort", "Suitability", "Design match", "Recommended action"])
    add_table(doc, t7_cohorts, "Table 7 | Prospective validation targets (clinical cohort scout)")

    doc.add_page_break()
    doc.add_heading("Author contributions", 2)
    doc.add_paragraph("[To be completed]")
    doc.add_heading("Competing interests", 2)
    doc.add_paragraph("The authors declare no competing interests.")
    doc.add_heading("Data availability", 2)
    doc.add_paragraph(
        "All analyses use public GEO datasets (GSE122897, GSE54083, GSE13353, GSE159610, GSE106520, and others). "
        "Processed tables and figure source data are in the project repository under results/rupture_meta_analysis, "
        "results/blood_bridge_panel, and results/advanced_validation."
    )
    doc.add_heading("References (selected)", 2)
    refs = [
        "Bakker, M.K. et al. Genome-wide association study of intracranial aneurysms identifies 17 risk loci. Nat. Genet. (2020).",
        "Rupture consensus meta-analysis and blood bridge panel scripts — Aneurysms Genomics Project (2025–2026).",
        "Clinical cohort scout: CIAP-1/2, ROAR, ISUIA feasibility summary — results/clinical_cohort_scout/.",
    ]
    for i, r in enumerate(refs, 1):
        doc.add_paragraph(f"{i}. {r}")

    doc.add_heading("Figure and table legends (combined)", 2)
    doc.add_paragraph(
        "This section duplicates captions for journal submission checklists. All figures are embedded above "
        "at 300 dpi PNG generated programmatically from project TSV outputs."
    )
    return doc


def main() -> None:
    print("Loading data...")
    data = _load_data()
    print("Generating figures...")
    fig_paths = {
        "workflow": FIG_DIR / "fig1_concordance_workflow.png",
        "transfer": FIG_DIR / "fig2_transfer_vs_genes.png",
        "consensus": FIG_DIR / "fig3_six_gene_meta_z.png",
        "blood": FIG_DIR / "fig4_blood_bridge_panel.png",
        "overlap": FIG_DIR / "fig5_panel_overlap.png",
    }
    _fig1_workflow(fig_paths["workflow"])
    _fig2_transfer_auroc(data["transfer"], fig_paths["transfer"])
    _fig3_consensus_meta(data["meta_conc"], fig_paths["consensus"])
    _fig4_blood_bridge(fig_paths["blood"])
    _fig5_panel_overlap(fig_paths["overlap"])

    print("Building Word document...")
    doc = build_document(data, fig_paths)
    doc.save(OUT_DOCX)
    print(f"Saved: {OUT_DOCX}")
    print(f"Figures: {FIG_DIR}")


if __name__ == "__main__":
    main()
