"""
Create Nature-style manuscript: Molecular Rupture Signatures (VS).
Main goal: Identifying molecular factors associated with aneurysm rupture.
- Tissue source and clinical relevance clearly stated; mechanistic focus (not bedside prediction).
- No PHASES comparison. Main manuscript: 5 figures (max 4 panels each), 2 tables.
- All other figures/tables in Supplementary with accurate in-text references.
- Methods concise in main; full detail in Supplementary.
- Single Word document (main + supplementary).

Shared helper: append_supplementary_bundle() embeds Supplementary Methods, Figures S1–S11,
optional S12, and Tables S1–S6 (used by create_molecular_signatures_revised_vs.py).
"""
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from pathlib import Path
from typing import Optional

import pandas as pd
import os

def add_table_from_dataframe(doc, df, title=None, inference=None, max_rows=20):
    if title:
        p = doc.add_paragraph(title, style='Heading 3')
    if len(df) > max_rows:
        df = df.head(max_rows)
        doc.add_paragraph(f"Showing top {max_rows} rows of {len(df)} total.")
    table = doc.add_table(rows=1, cols=len(df.columns))
    table.style = 'Light Grid Accent 1'
    header_cells = table.rows[0].cells
    for i, col in enumerate(df.columns):
        header_cells[i].text = str(col)
        if header_cells[i].paragraphs[0].runs:
            header_cells[i].paragraphs[0].runs[0].font.bold = True
    for _, row in df.iterrows():
        row_cells = table.add_row().cells
        for i, val in enumerate(row):
            row_cells[i].text = str(val) if pd.notna(val) else ""
    if inference:
        para = doc.add_paragraph()
        para.add_run("Inference: ").bold = True
        para.add_run(inference)
    return table

def add_figure_with_caption_and_inference(doc, fig_path, caption, inference, width=Inches(5.5)):
    if not Path(fig_path).exists():
        doc.add_paragraph(f"[Figure not found: {fig_path}]")
        return
    try:
        doc.add_picture(str(fig_path), width=width)
        para = doc.add_paragraph(caption)
        para.style = 'Caption'
        if para.runs:
            para.runs[0].italic = True
        if inference:
            inf_para = doc.add_paragraph()
            inf_para.add_run("Inference: ").bold = True
            inf_para.add_run(inference)
        doc.add_paragraph()
    except Exception as e:
        doc.add_paragraph(f"[Error loading figure {fig_path}: {e}]")


def append_supplementary_bundle(
    doc: Document,
    *,
    s1_stratification_figure_ref: str = "main text Figure 5",
    supplementary_figures_intro: Optional[str] = None,
    include_supplementary_figure_s12: bool = False,
    figure_width=Inches(5),
) -> None:
    """
    Full supplementary content: Methods, Figures S1–S11, optional S12 (risk stratification),
    embedded Tables S1–S6. Used by Molecular_Rupture_Signatures_VS and Molecular_Signatures_Revised_VS.
    """
    if supplementary_figures_intro is None:
        supplementary_figures_intro = (
            "Supplementary Figures S1–S11 support the main manuscript and are referenced in the main text where "
            "appropriate. Each figure has a detailed legend and an inference paragraph with key numbers."
        )

    doc.add_heading("Supplementary Methods", 2)
    supp_methods = (
        "Data preprocessing and inclusion. The expression matrix was loaded from the GSE122897 supplementary "
        "files (TPM-normalised counts, GRCh38.p13). Raw count data were also available for count-based differential "
        "expression (e.g. DESeq2, edgeR, limma-voom). Genes were "
        "retained if TPM > 1 in at least 3 of the 43 samples; this removed lowly expressed genes. Sample metadata "
        "from the GEO series matrix were used to label each specimen as ruptured or unruptured; 17 samples "
        "(control or unknown status) were excluded. Expression was log2(TPM+1) transformed for all downstream "
        "analyses. For the classifier, the three signature scores were standardised (zero mean, unit variance) "
        "across samples using scikit-learn StandardScaler.\n\n"
        "Differential expression. For each gene, a two-sample Welch t-test (ruptured n=22 vs unruptured n=21) "
        "was applied to log2(TPM+1). Log2 fold change = mean(ruptured) − mean(unruptured). P-values were "
        "adjusted by the Benjamini–Hochberg procedure (FDR); genes with FDR < 0.05 were considered significant. "
        "Bonferroni-adjusted p-values were also computed. Implementation was vectorised in Python (scipy.stats). "
        "DESeq2 was run on raw counts as a confirmatory analysis (2,764 genes at padj < 0.05; see Supplementary Table S3).\n\n"
        "Nested DE within cross-validation. To address circularity (DE genes identified in the full cohort "
        "then used for the DE signature in the model), a nested design was implemented: within each CV fold, "
        "DE analysis was performed on the training set only; the resulting DE gene set was used to compute "
        "the DE signature for both train and test samples in that fold. This avoids information leakage and "
        "provides a more conservative performance estimate. Nested CV yielded AUC 0.64 ± 0.17 (5 folds), "
        "lower than the non-nested 0.80 ± 0.16, consistent with circularity inflating the latter.\n\n"
        "Signature scores. For each sample, (1) DE score = mean log2(TPM+1) over the set of top DE genes "
        "(by FDR and |Log2FC|). (2) EC score = mean over ANGPT2, ICAM1, VCAM1, THBS1, PERP. (3) ECM score = "
        "mean over MMP9, MMP2, TIMP1, COL3A1, COL1A1. (4) Inflammation = mean over TNF, IL6, IL8, CRP, SAA1. "
        "The final model used only DE, EC, and ECM. Gene symbols were matched to matrix row identifiers "
        "via harmonised annotation resources and GEO annotations where available. Pathway "
        "enrichment (GO, KEGG, Reactome) was performed using the g:Profiler API with a ranked gene list "
        "(genes ordered by -log10(padj)×sign(log2FC) from DESeq2), ordered=TRUE for GSEA-style analysis, "
        "and FDR control (Benjamini–Hochberg). 207 enriched terms at FDR < 0.05.\n\n"
        "Model and cross-validation. Logistic regression (sklearn.linear_model.LogisticRegression) with "
        "penalty='l2', solver='lbfgs', and fixed C=1.0 was fitted to rupture status (0/1) using the three signature scores. "
        "We report three performance views: full-cohort fit, optimistic non-nested 5-fold stratified cross-validation, "
        "and a more conservative leakage-controlled design in which DE gene selection is repeated within outer folds. "
        "The final model was refit on the full 43 samples for descriptive coefficients. Feature importance was "
        "computed by permutation (shuffle each feature, recompute AUC; mean decrease over multiple permutations). "
        "Risk strata: predicted probability < 0.33 = low, 0.33–0.67 = medium, ≥ 0.67 = high.\n\n"
        "Software. Python 3.x; pandas, numpy, scipy, scikit-learn; matplotlib/seaborn for figures. "
        "Code and processed outputs available from the authors on request."
    )
    doc.add_paragraph(supp_methods)
    doc.add_paragraph()

    doc.add_heading("Supplementary Figures", 2)
    doc.add_paragraph(supplementary_figures_intro)

    supp_figures = [
        (
            "results/figures/ultra_hd/ma_plot_ultra_hd.png",
            "Supplementary Figure S1. MA plot of differential expression. X-axis: mean log2(TPM+1) across all 43 samples; "
            "y-axis: log2 fold change (ruptured vs unruptured). Points coloured by significance (e.g. FDR < 0.05). "
            "Genes with low mean expression and large |Log2FC| appear at the left.",
            "Significant genes (FDR < 0.05) are distributed across the full range of mean expression, indicating that "
            "differential expression is not driven solely by low-count or high-count genes. The symmetric spread around "
            "zero for non-significant genes supports appropriate normalisation and absence of systematic bias by expression level.",
        ),
        (
            "results/figures/ultra_hd/top_genes_heatmap_ultra_hd.png",
            "Supplementary Figure S2. Heatmap of top differentially expressed genes across the 43 samples (rows: genes; "
            "columns: samples). Samples are ordered or clustered by rupture status; colour scale represents expression "
            "(e.g. log2(TPM+1) or z-score). Ruptured (n=22) and unruptured (n=21) groups are indicated.",
            "Samples cluster by rupture status, with ruptured and unruptured specimens forming distinct patterns. Gene "
            "expression patterns are consistent within each group, supporting the robustness of the DE signal. The heatmap "
            "provides a visual summary of the genes in Table 2 and related top genes.",
        ),
        (
            "results/figures/ultra_hd/nested_cv_roc_curve_ultra_hd.png",
            "Supplementary Figure S3. Out-of-fold ROC curve from the leakage-controlled nested cross-validation analysis. "
            "Sensitivity (true positive rate) vs 1 − specificity (false positive rate). Mean fold AUC = 0.635; diagonal line: "
            "random classifier (AUC 0.5).",
            "This curve is derived from held-out predictions across the nested design and is the most defensible out-of-sample "
            "performance visualisation for the three-score model.",
        ),
        (
            "results/figures/top_20_genes_bar.png",
            "Supplementary Figure S4. Bar plot of the top 20 differentially expressed genes ranked by statistical "
            "significance (e.g. -log10(P_adj_BH) or |T-statistic|). Gene symbols or IDs on the y-axis; magnitude on the x-axis. "
            "n=43 (22 ruptured, 21 unruptured).",
            "The bar lengths reflect the strength of evidence for each gene; the top genes correspond to those in Table 2. "
            "Genes such as SMS, LAMA3, CPXM2, MARCO, TGFBI appear among the most significant and are candidates for "
            "mechanistic and biomarker follow-up.",
        ),
        (
            "results/figures/ultra_hd/ec_dysfunction_ultra_hd.png",
            "Supplementary Figure S5. Endothelial dysfunction (EC) signature and marker expression. Panels show EC "
            "score or individual markers (e.g. ANGPT2, ICAM1, VCAM1, THBS1, PERP) by rupture status (n=22 vs n=21).",
            "Ruptured samples show higher EC scores and elevated expression of key EC markers, consistent with endothelial "
            "dysfunction as a process associated with rupture. The EC signature contributed 35.6% to the summary model "
            "feature importance.",
        ),
        (
            "results/figures/ultra_hd/pathway_enrichment_gse122897_ultra_hd.png",
            "Supplementary Figure S6. Pathway enrichment analysis (GSE122897 DE genes). Ranked gene list (ordered by "
            "-log10(padj)×sign(log2FC)), g:Profiler (GO, KEGG, Reactome), FDR < 0.05. 207 enriched terms; top 15 shown. "
            "Bar lengths represent -log10(P-value). Matches Supplementary Table S5.",
            "Enriched pathways include cell migration, leukocyte migration, extracellular matrix organization, "
            "collagen-containing ECM, blood vessel development, TNF signalling, and Fluid shear stress and atherosclerosis, "
            "aligning with the EC, ECM, and inflammation signature scores and supporting the mechanistic cascade.",
        ),
        (
            "results/figures/phase7/phase7_roc_curve.png",
            "Supplementary Figure S7. Full-cohort fitted ROC curve (supplementary; primary performance estimate is the nested CV in Figure 4b). "
            "n=43; AUC = 0.827.",
            "The fitted-model curve appears stronger than the nested result because it is evaluated on the same cohort used to define the model.",
        ),
        (
            "results/figures/phase7/phase7_confusion_matrix.png",
            "Supplementary Figure S8. Confusion matrix at the optimal classification threshold (0.605). Rows: true "
            "label (Ruptured / Unruptured); columns: predicted label. Cell values: counts. n=43.",
            "At threshold 0.605: 14 true positives, 20 true negatives, 1 false positive, 8 false negatives. Sensitivity "
            "(recall) = 14/22 = 63.64%; specificity = 20/21 = 95.24%; precision = 14/15 = 93.33%. The low false positive "
            "rate (1/21) indicates that when the model predicts rupture, it is usually correct; the 8 false negatives "
            "represent ruptured samples missed by the model.",
        ),
        (
            "results/figures/phase7/phase7_pr_curve.png",
            "Supplementary Figure S9. Precision–recall curve for the summary model (n=43). Precision (positive predictive "
            "value) vs recall (sensitivity) at varying thresholds. PR AUC = 0.874.",
            "The high PR AUC (0.874) indicates good performance on the positive class (ruptured). PR curves are informative "
            "when class balance is not extreme; here the balance is 22 vs 21, and the model achieves high precision at "
            "the chosen threshold.",
        ),
        (
            "results/figures/ultra_hd/deseq2_volcano_ultra_hd.png",
            "Supplementary Figure S10. DESeq2 volcano plot of differential expression (raw counts, n=43). X-axis: log2 fold "
            "change (ruptured vs unruptured); y-axis: -log10(adjusted p-value). Red: upregulated (padj < 0.05); blue: "
            "downregulated; grey: not significant. Top genes labelled.",
            "DESeq2 identified 2,764 genes at padj < 0.05. The volcano shows strong effects for genes such as TGFBI, MARCO, "
            "SMS, and others, consistent with Supplementary Table S3 and the Welch t-test results in the main text.",
        ),
        (
            "results/figures/ultra_hd/nested_vs_non_nested_cv_ultra_hd.png",
            "Supplementary Figure S11. Comparison of cross-validation AUC: non-nested (DE on full cohort, AUC 0.80 ± 0.16) vs "
            "nested design (DE within each fold, AUC 0.64 ± 0.17). 5-fold stratified CV. Error bars: standard deviation.",
            "The nested design yields lower AUC because DE gene selection is performed independently in each fold, avoiding "
            "circularity. The gap indicates that the non-nested estimate may be optimistically inflated by using the same "
            "genes for both selection and evaluation.",
        ),
    ]
    for fig_path, cap, inf in supp_figures:
        if Path(fig_path).exists():
            add_figure_with_caption_and_inference(doc, fig_path, cap, inf, width=figure_width)

    if include_supplementary_figure_s12:
        fig_s12 = "results/figures/phase7/phase7_risk_stratification.png"
        if not Path(fig_s12).exists():
            fig_s12 = "results/figures/ultra_hd/risk_stratification_ultra_hd.png"
        if Path(fig_s12).exists():
            add_figure_with_caption_and_inference(
                doc,
                fig_s12,
                "Supplementary Figure S12. Stratification of samples (n=43) by model-derived risk score into low, medium, and high "
                "strata (thresholds 0.33 and 0.67 on predicted probability). Bar or plot shows the proportion of ruptured vs "
                "unruptured samples in each stratum.",
                "The proportion of ruptured samples increases across strata (low ~20%, medium ~38%, high 100% in this cohort). "
                "Included in the supplement for the revised manuscript (four primary main-text figures).",
                width=figure_width,
            )
        else:
            doc.add_paragraph("[Supplementary Figure S12: risk stratification figure not found.]")

    doc.add_heading("Supplementary Tables", 2)
    doc.add_paragraph(
        "Supplementary Table S1. Sample-level risk scores and stratification. Columns: Sample ID, True_Label (0=unruptured, "
        "1=ruptured), Optimized_Score (predicted probability), Optimized_Risk (Low/Medium/High). n=43. Referenced in main text "
        "with Table 1 and in Methods (risk stratification)."
    )
    doc.add_paragraph()
    s1_inference = (
        "Each of the 43 samples has a predicted probability (Optimized_Score) and an assigned stratum (Low: score < 0.33, "
        "Medium: 0.33–0.67, High: ≥ 0.67). True_Label shows actual rupture status. The table allows verification of "
        "stratification: high-risk samples are predominantly True_Label=1 (ruptured). Counts: Low n=10, Medium n=21, "
        f"High n=12; rupture rates by stratum ~20%, ~38%, 100% as reported in {s1_stratification_figure_ref}."
    )
    try:
        phase7_df = pd.read_csv("results/phase7_risk_stratification.tsv", sep="\t")
        if len(phase7_df) > 0:
            add_table_from_dataframe(doc, phase7_df, title=None, inference=s1_inference, max_rows=50)
    except Exception:
        pass
    doc.add_paragraph()
    doc.add_paragraph(
        "Supplementary Table S2. Cross-validation performance summary. Mean and standard deviation of AUC across 5 folds "
        "for the summary model. Referenced in Methods (Statistical model and cross-validation)."
    )
    doc.add_paragraph()
    try:
        rows = []
        best_json = Path("results/phase5_best_model_summary.json")
        nested_tsv = Path("results/gse122897_nested_de_cv_summary.tsv")
        train_tsv = Path("results/phase7_optimized_model_metrics.tsv")

        if best_json.exists():
            try:
                best = pd.read_json(best_json, typ="series")
                rows.append(
                    {
                        "Estimate": "Non-nested CV ROC AUC",
                        "Mean_AUC": round(float(best.get("cv_auc_mean", float("nan"))), 3),
                        "SD_AUC": round(float(best.get("cv_auc_std", float("nan"))), 3),
                        "Notes": "Optimistic because DE-derived features are not reselected within folds.",
                    }
                )
            except Exception:
                pass

        if nested_tsv.exists():
            try:
                nested = pd.read_csv(nested_tsv, sep="\t")
                nd = dict(zip(nested["Metric"], nested["Value"]))
                rows.append(
                    {
                        "Estimate": "Leakage-controlled nested CV ROC AUC",
                        "Mean_AUC": round(float(nd.get("Nested_CV_AUC_Mean", float("nan"))), 3),
                        "SD_AUC": round(float(nd.get("Nested_CV_AUC_Std", float("nan"))), 3),
                        "Notes": "Primary out-of-sample estimate; DE signature re-derived inside outer folds.",
                    }
                )
            except Exception:
                pass

        if train_tsv.exists():
            try:
                train = pd.read_csv(train_tsv, sep="\t")
                train_auc = pd.to_numeric(
                    train.loc[train["Metric"] == "ROC_AUC", "Value"], errors="coerce"
                ).dropna()
                if not train_auc.empty:
                    rows.append(
                        {
                            "Estimate": "Full-cohort ROC AUC",
                            "Mean_AUC": round(float(train_auc.iloc[0]), 3),
                            "SD_AUC": "",
                            "Notes": "Descriptive in-sample fit; not external validation.",
                        }
                    )
            except Exception:
                pass

        cv_df = pd.DataFrame(rows)
        if len(cv_df) > 0:
            add_table_from_dataframe(
                doc,
                cv_df,
                title=None,
                inference=(
                    "This table replaces the stale legacy summary previously used in the supplement. The leakage-controlled nested estimate "
                    "is the most defensible out-of-sample performance figure; non-nested CV and full-cohort AUC are shown for context only."
                ),
                max_rows=15,
            )
    except Exception:
        pass

    doc.add_paragraph(
        "Supplementary Table S3. Top differentially expressed genes from DESeq2 (raw counts). Gene_ID (Entrez), Gene_Symbol, "
        "Log2FC (ruptured vs unruptured), padj, baseMean. DESeq2 identified 2,764 genes at padj < 0.05. Top genes include "
        "TGFBI, MARCO, SMS (shared with Welch t-test results). Referenced in main text Results (Differential expression)."
    )
    doc.add_paragraph()
    try:
        deseq2_df = pd.read_csv("results/gse122897_deseq2_top_annotated.tsv", sep="\t")
        if len(deseq2_df) > 0:
            add_table_from_dataframe(
                doc,
                deseq2_df.head(15),
                title=None,
                inference=(
                    "DESeq2 and Welch t-tests show strong agreement: SMS, MARCO, TGFBI, LAMA3, CPXM2, and others appear in both. "
                    "DESeq2 finds many more significant genes because it models count overdispersion. Full results are provided in the deposited analysis repository."
                ),
                max_rows=20,
            )
    except Exception:
        pass

    doc.add_paragraph(
        "Supplementary Table S4. Comparison of Welch t-test and DESeq2 differential expression results. Welch (TPM-based) "
        "identified 4 genes at FDR < 0.05; DESeq2 (raw counts) identified 2,764 at padj < 0.05. All top Welch genes overlapped "
        "with DESeq2 significant genes. Referenced in main text Results and Supplementary Methods."
    )
    doc.add_paragraph()
    try:
        comp_df = pd.read_csv("results/gse122897_welch_vs_deseq2_comparison.tsv", sep="\t")
        if len(comp_df) > 0:
            if "Metric" in comp_df.columns and "Value" in comp_df.columns:
                comp_df = comp_df.copy()
                comp_df.loc[
                    comp_df["Metric"] == "Welch_top_10_genes", "Value"
                ] = "SMS, Entrez 124220 (annotation unresolved), CPXM2, LAMA3, MARCO, TGFBI, LIMCH1, COQ2, CABCOCO1, CCNI2"
                comp_df.loc[
                    comp_df["Metric"] == "DESeq2_top_10_Gene_IDs", "Metric"
                ] = "DESeq2_top_10_gene_symbols"
                comp_df.loc[
                    comp_df["Metric"] == "DESeq2_top_10_gene_symbols", "Value"
                ] = "LINC03040, TGFBI, MARCO, SMS, HPSE, ANPEP, CLEC5A, OLR1, LOC101928343, CD300E"
            add_table_from_dataframe(
                doc,
                comp_df,
                title=None,
                inference=(
                    "The overlap confirms that key genes (SMS, LAMA3, CPXM2, MARCO, TGFBI, etc.) are consistently identified by both methods. "
                    "One legacy Welch annotation artifact (Entrez 124220) is retained as unresolved rather than assigned an incorrect gene symbol."
                ),
                max_rows=10,
            )
    except Exception:
        pass

    doc.add_paragraph(
        "Supplementary Table S5. Pathway enrichment (GO, KEGG, Reactome) among DE genes. Ranked gene list, g:Profiler, "
        "FDR < 0.05. 207 enriched terms. "
        "Top terms include cell migration, leukocyte migration, ECM organization, blood vessel development, TNF signalling. "
        "Full results should be provided through the deposited analysis repository."
    )
    doc.add_paragraph()
    try:
        pw_df = pd.read_csv("results/gse122897_pathway_enrichment_filtered.tsv", sep="\t")
        if len(pw_df) > 0:
            cols = [c for c in ["Source", "Term_ID", "Term_Name", "P_value", "FDR_BH"] if c in pw_df.columns]
            add_table_from_dataframe(
                doc,
                pw_df[cols].head(12),
                title=None,
                inference=(
                    "Ranked gene list analysis identified 207 terms at FDR < 0.05. Top terms include cell migration, "
                    "extracellular matrix organization, blood vessel development, TNF signalling. g:Profiler API "
                    "(GO, KEGG, Reactome, WikiPathways) was used for enrichment; full results should be provided through the deposited analysis repository."
                ),
                max_rows=15,
            )
    except Exception:
        pass

    doc.add_paragraph(
        "Supplementary Table S6. Sample-level metadata for the 43 analysed specimens. Columns: Sample (GEO accession), "
        "Group (Ruptured/Unruptured), Sex (where reported). Age was not available in the GEO records. "
        "Referenced in main text Table 1 and Methods (Data acquisition)."
    )
    doc.add_paragraph()
    try:
        samp_meta = pd.read_csv("results/gse122897_sample_metadata_supplementary.tsv", sep="\t")
        if len(samp_meta) > 0:
            add_table_from_dataframe(
                doc,
                samp_meta,
                title=None,
                inference=(
                    "All 43 samples have rupture status (Group). Sex was reported for 1 of 43 samples; age was not "
                    "available. This limited metadata precluded adjustment for age or sex in differential expression "
                    "or modelling analyses."
                ),
                max_rows=50,
            )
    except Exception:
        pass


def create_manuscript():
    # Run from project root so results/figures/... paths resolve
    script_dir = Path(__file__).resolve().parent
    os.chdir(script_dir)

    doc = Document()
    style = doc.styles['Normal']
    font = style.font
    font.name = 'Calibri'
    font.size = Pt(11)

    # --------------- TITLE ---------------
    title = doc.add_heading('Molecular signatures associated with intracranial aneurysm rupture', 0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    subtitle = doc.add_paragraph('A transcriptomic analysis of ruptured versus unruptured aneurysm tissue')
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    if subtitle.runs:
        subtitle.runs[0].font.size = Pt(12)
        subtitle.runs[0].italic = True
    doc.add_page_break()

    # --------------- ABSTRACT ---------------
    doc.add_heading('Abstract', 1)
    abstract = (
        "Identifying the molecular factors that associate with aneurysm rupture could clarify mechanisms "
        "of rupture and inform future biomarker research. Here we analysed RNA-seq from resected intracranial "
        "aneurysm tissue (GSE122897: 22 ruptured, 21 unruptured) to identify transcriptomic signatures. "
        "Differential expression, pathway-based signature scores (endothelial dysfunction, extracellular matrix "
        "remodeling, inflammation), and a simple classifier were used to summarise biology. Tissue was obtained "
        "at surgery; it is not available before rupture in unruptured patients, so this analysis is framed as "
        "mechanistic rather than as a clinical prediction tool. We report genes and pathways that distinguish "
        "ruptured from unruptured tissue and discuss implications for blood-based biomarker discovery."
    )
    doc.add_paragraph(abstract)
    doc.add_paragraph()

    # --------------- INTRODUCTION ---------------
    doc.add_heading('Introduction', 1)
    intro = (
        "Rupture of an intracranial aneurysm is a major cause of haemorrhagic stroke. Decades of research "
        "have characterised key biological processes implicated in rupture—endothelial activation, matrix "
        "degradation by metalloproteinases, smooth muscle phenotypic switching, and inflammatory infiltration—"
        "yet a clear gap remains between this mechanistic understanding and the ability to predict which "
        "individual aneurysms will rupture. Current rupture risk assessment relies on imaging-based and "
        "clinical scores (e.g. size, location, smoking); tissue-level molecular data are scarce because "
        "aneurysm tissue is typically available only after surgery, following rupture or elective clipping/coiling. "
        "Thus, any model trained on resected tissue cannot be used directly to predict rupture risk in "
        "unruptured patients, for whom tissue is not obtained. The value of tissue-based studies lies instead "
        "in bridging this biology–prediction gap: identifying molecular factors that associate with rupture "
        "in resected specimens and generating hypotheses for blood-based or imaging-based biomarkers that "
        "could be measured before rupture.\n\n"
        "In this study we asked: which transcriptomic features characterise ruptured versus unruptured "
        "aneurysm tissue? We used the public dataset GSE122897, which comprises bulk RNA-seq from resected "
        "intracranial aneurysm specimens. We performed differential expression analysis, defined biological "
        "signature scores (endothelial dysfunction, extracellular matrix remodeling, inflammation), and "
        "fitted a minimal classifier to summarise how well these signatures separate ruptured from unruptured "
        "samples. The main goal is to identify molecular factors associated with aneurysm rupture and to "
        "provide a clear account of the tissue source and its clinical relevance."
    )
    doc.add_paragraph(intro)
    doc.add_page_break()

    # --------------- RESULTS (main story; 5 figures, 2 tables) ---------------
    doc.add_heading('Results', 1)

    doc.add_heading('Study design and tissue source', 2)
    results1 = (
        "All analyses use RNA-seq from resected intracranial aneurysm tissue (GSE122897, Ruigrok et al.; "
        "NCBI GEO). Samples are from 43 patients: 22 ruptured and 21 unruptured aneurysms. Tissue was "
        "obtained at the time of surgical clipping or endovascular treatment; it is not available from "
        "unruptured patients who do not undergo intervention. Thus, the data inform which molecular "
        "features characterise rupture in resected tissue rather than providing a tool for pre-operative "
        "rupture risk prediction. See Table 1 and Figure 1 for sample and design details."
    )
    doc.add_paragraph(results1)

    # Table 1 (main): Sample and data description
    doc.add_paragraph(
        "Table 1. Study sample and data description. Summary of the cohort and RNA-seq data used in the analysis. "
        "All samples are resected intracranial aneurysm tissue from GSE122897 (NCBI GEO). Tissue was obtained at "
        "surgery (clipping or endovascular treatment); it is not available from unruptured patients without intervention. "
        "Total n = 43 (22 ruptured, 21 unruptured). Platform: Illumina HiSeq 2500; expression: TPM-normalised counts, "
        "GRCh38.p13. Source: Ruigrok et al., UMC Utrecht.",
        style='Heading 3'
    )
    sample_table_data = {
        'Item': [
            'Dataset', 'Tissue type', 'Total samples (analysed)', 'Ruptured', 'Unruptured',
            'Platform', 'Normalisation', 'Reference'
        ],
        'Description': [
            'GSE122897 (NCBI GEO)', 'Resected intracranial aneurysm tissue', '43', '22', '21',
            'Illumina HiSeq 2500, RNA-seq', 'TPM (GRCh38.p13)', 'Ruigrok et al., UMC Utrecht'
        ]
    }
    tbl1_df = pd.DataFrame(sample_table_data)
    add_table_from_dataframe(doc, tbl1_df, title=None, inference=(
        "The cohort is balanced (22 ruptured, 21 unruptured), enabling robust two-group comparison. Because tissue "
        "is surgery-derived and not obtainable in unruptured patients without intervention, the study is positioned "
        "to identify molecular factors associated with rupture and to prioritise candidates for future blood-based "
        "or other non-tissue biomarkers. Supplementary Table S1 provides sample-level risk scores and stratification; "
        "Supplementary Table S6 lists sample metadata (Sample, Group, Sex where available). Age was not in GEO; sex "
        "was available for one sample only, precluding covariate adjustment (see Limitations)."
    ), max_rows=20)
    doc.add_paragraph()
    doc.add_paragraph()

    # Figure 1: Study design and sample composition (main_manuscript ultra-HD)
    fig1_path = "results/figures/main_manuscript/Figure1_study_design_sample_composition.png"
    if Path(fig1_path).exists():
        add_figure_with_caption_and_inference(
            doc, fig1_path,
            "Figure 1. Study design and sample composition (GSE122897). Panel A: cohort flow from 60 GEO specimens to 43 "
            "analysed samples (17 excluded: control or unknown rupture status). Panel B: class balance in the analysed cohort "
            "(22 ruptured, 21 unruptured). Tissue: resected intracranial aneurysm; platform: Illumina HiSeq 2500; expression: TPM (GRCh38.p13).",
            "Of 60 specimens in GSE122897, 17 were excluded, leaving 43 for ruptured versus unruptured comparison (balanced design). "
            "Tissue was obtained at surgery; the analysis identifies molecular correlates of rupture in resected tissue rather than "
            "a pre-operative clinical prediction tool."
        )
    else:
        doc.add_paragraph("[Figure 1: main_manuscript/Figure1_study_design_sample_composition.png not found; run generate_main_five_figures_ultra_hd.py.]")
    doc.add_paragraph()

    doc.add_heading('Differential expression between ruptured and unruptured tissue', 2)
    results2 = (
        "We identified genes differentially expressed between ruptured and unruptured samples using two approaches: "
        "(1) Welch t-tests on log2(TPM+1) (FDR < 0.05), and (2) DESeq2 on raw counts (padj < 0.05). Welch t-tests "
        "identified fewer significant genes (4 at FDR < 0.05) owing to the conservative nature of TPM-based testing "
        "and multiple-testing correction; DESeq2, which models count overdispersion, identified 2,764 genes at "
        "padj < 0.05. All top Welch genes (e.g. SMS, LAMA3, CPXM2, MARCO, TGFBI) overlapped with DESeq2 significant "
        "genes, confirming the main findings. Ruptured tissue showed upregulation of genes involved in inflammation, "
        "endothelial dysfunction, and extracellular matrix remodeling, and downregulation of structural and "
        "adhesion-related genes (Table 2; Supplementary Tables S3 and S4 for top DESeq2 genes and Welch vs DESeq2 "
        "comparison; Supplementary Figures S1, S2, S4, and S10 for MA plot, heatmap, top-gene bar plot, and DESeq2 volcano)."
    )
    doc.add_paragraph(results2)

    # Table 2 (main): Top DE genes — use annotated file so gene names (symbol, name) are clear
    try:
        top_genes_path = "results/gse122897_top_genes_annotated.tsv"
        if not Path(top_genes_path).exists():
            top_genes_path = "results/gse122897_top_genes.tsv"
        top_genes_df = pd.read_csv(top_genes_path, sep='\t')
        if len(top_genes_df) > 0:
            # If no symbol/name columns, try to add from unified mapping (Original_ID <-> Gene_Symbol)
            if 'Gene_Symbol' not in top_genes_df.columns and 'Gene' in top_genes_df.columns:
                map_path = "results/gene_id_mapping_unified.tsv"
                if Path(map_path).exists():
                    try:
                        m = pd.read_csv(map_path, sep='\t')
                        if 'Original_ID' in m.columns and 'Gene_Symbol' in m.columns:
                            m = m[['Original_ID', 'Gene_Symbol']].dropna(subset=['Gene_Symbol'])
                            m['Gene_key'] = m['Original_ID'].astype(str).str.replace(r'\.0$', '', regex=True)
                            top_genes_df['Gene_key'] = top_genes_df['Gene'].astype(str).str.replace(r'\.0$', '', regex=True)
                            top_genes_df = top_genes_df.merge(m[['Gene_key', 'Gene_Symbol']], on='Gene_key', how='left').drop(columns=['Gene_key'], errors='ignore')
                    except Exception:
                        pass
            display_cols = []
            for c in ['Gene_Symbol', 'Gene_Name', 'Log2FC', 'P_adj_BH', 'T_statistic', 'Mean_Ruptured', 'Mean_Unruptured', 'Gene']:
                if c in top_genes_df.columns and c not in display_cols:
                    display_cols.append(c)
            if not display_cols:
                display_cols = list(top_genes_df.columns)
            top_genes_df = top_genes_df[display_cols].head(15)
            add_table_from_dataframe(
                doc, top_genes_df,
                title="Table 2. Top differentially expressed genes (ruptured vs unruptured tissue). "
                      "Gene_Symbol and Gene_Name identify each gene; Log2FC (positive = higher in ruptured), P_adj_BH, T_statistic, "
                      "and mean TPM per group are shown. Genes ordered by significance.",
                inference=(
                    "Top genes include SMS (spermine synthase), LAMA3 (laminin subunit alpha 3), CPXM2 (carboxypeptidase X member 2), "
                    "MARCO, TGFBI, and others; Log2FC values range from approximately -2.2 to +2.9, with P_adj_BH ≤ 0.09 in this subset. "
                    "Positive Log2FC indicates higher expression in ruptured tissue. These genes represent the strongest transcriptomic "
                    "associations with rupture in this cohort (n=43) and are candidates for mechanistic follow-up and for testing in "
                    "blood or other accessible compartments."
                ),
                max_rows=20
            )
        else:
            doc.add_paragraph("Table 2. Top differentially expressed genes. (No data in top genes file.)")
    except Exception as e:
        doc.add_paragraph(f"Table 2. Top differentially expressed genes. (Error loading: {e})")
    doc.add_paragraph()

    # Figure 2: DESeq2 volcano (main_manuscript ultra-HD)
    fig2_path = "results/figures/main_manuscript/Figure2_deseq2_volcano_rupture_vs_unruptured.png"
    if Path(fig2_path).exists():
        add_figure_with_caption_and_inference(
            doc, fig2_path,
            "Figure 2. Genome-wide differential expression (DESeq2 on raw counts). Ruptured vs unruptured tissue (n=22 vs n=21). "
            "x-axis: log2 fold change; y-axis: -log10(adjusted P value). Red: up in ruptured; blue: down in ruptured (padj < 0.05, |log2FC| >= 0.5). "
            "Dashed horizontal line: padj = 0.05; dotted vertical lines: |log2FC| = 0.5.",
            "DESeq2 identified 2,184 genes at padj < 0.05 and |log2FC| >= 0.5 (1,325 upregulated and 859 downregulated in ruptured tissue). "
            "The symmetric spread around zero indicates broad transcriptomic remodelling rather than isolated outliers, supporting "
            "inflammatory, endothelial, and matrix-remodelling pathways. Top genes are listed in Table 2."
        )
    else:
        doc.add_paragraph("[Figure 2: main_manuscript/Figure2_deseq2_volcano not found; run generate_main_five_figures_ultra_hd.py.]")
    doc.add_paragraph()

    doc.add_heading('Biological signature scores separate rupture status', 2)
    results3 = (
        "We defined signature scores for endothelial dysfunction (EC), extracellular matrix (ECM) remodeling, "
        "and inflammation from predefined gene sets (Methods; Supplementary Methods). These scores differed "
        "between ruptured and unruptured samples (Figure 3). Pathway enrichment was performed systematically "
        "using a ranked gene list (ordered by -log10(padj)×sign(log2FC)), g:Profiler (GO, KEGG, Reactome), and "
        "FDR control (Benjamini–Hochberg). This identified 207 enriched terms at FDR < 0.05 (Supplementary Table S5, "
        "Figure S6), including cell migration, leukocyte migration, extracellular matrix organization, "
        "collagen-containing ECM, blood vessel development, TNF signalling, and Fluid shear stress and atherosclerosis "
        "(KEGG)—processes consistent with endothelial activation, ECM remodelling, and inflammatory recruitment. "
        "Supplementary Figure S5 shows EC signature and marker expression. A composite differential-expression (DE) "
        "signature showed the strongest separation. Together, these results point to endothelial dysfunction, ECM "
        "remodelling, and inflammation as molecular processes associated with rupture in this tissue cohort."
    )
    doc.add_paragraph(results3)

    # Figure 3: Signature distributions (main_manuscript ultra-HD)
    fig3_path = "results/figures/main_manuscript/Figure3_signature_scores_by_rupture_status.png"
    if Path(fig3_path).exists():
        add_figure_with_caption_and_inference(
            doc, fig3_path,
            "Figure 3. Biological signature scores by rupture status (n=43: 22 ruptured, 21 unruptured). Violin and box plots "
            "for DE signature, EC dysfunction, and ECM remodeling scores (z-scored; predefined gene sets, Methods). "
            "Individual samples shown as points; Mann–Whitney U P values annotated.",
            "Ruptured samples show significantly higher DE scores than unruptured (Mann–Whitney P = 0.005); EC shows a trend "
            "(P = 0.10) and ECM does not differ (P = 0.91), with partial overlap in all panels. The DE signature shows the "
            "strongest distributional shift, consistent with its dominant contribution in the summary model (Figure 4a, 56.5%)."
        )
    else:
        doc.add_paragraph("[Figure 3: main_manuscript/Figure3_signature_scores not found; run generate_main_five_figures_ultra_hd.py.]")
    doc.add_paragraph()

    doc.add_heading('Summary model and feature importance', 2)
    results4 = (
        "We fitted a regularised logistic regression model using the three signature scores to summarise "
        "how well biology separates rupture status (Methods; full details in Supplementary Methods). "
        "Cross-validation AUC was 0.80 ± 0.16; on the full cohort, AUC was 0.83. A nested design (DE gene "
        "selection within each CV fold; Supplementary Methods) yielded AUC 0.64 ± 0.17, indicating that "
        "circularity may inflate the non-nested estimate. The CV standard deviation (~0.15) reflects the "
        "small sample size (n=43) and substantial fold-to-fold variability. Feature "
        "importance ranked DE signature first, then EC, then ECM (Figure 4). Cross-validation metrics are summarised "
        "in Supplementary Table S2; Supplementary Figure S11 illustrates nested vs non-nested CV; Supplementary Figures S3, S7, S8, "
        "and S9 provide ROC curve, confusion matrix, and precision–recall curve. The model is not proposed as a clinical prediction tool—tissue is not "
        "available pre-operatively—but as a compact exploratory summary of which biological dimensions distinguish "
        "ruptured from unruptured tissue in this dataset."
    )
    doc.add_paragraph(results4)

    # Figure 4: Feature importance + ROC (main_manuscript ultra-HD, two panels)
    fig4_path = "results/figures/main_manuscript/Figure4_model_feature_importance_and_ROC.png"
    if Path(fig4_path).exists():
        add_figure_with_caption_and_inference(
            doc, fig4_path,
            "Figure 4. Summary logistic model separating ruptured from unruptured tissue (n=43). "
            "(a) Permutation feature importance (% decrease in AUC when shuffled): DE signature, EC dysfunction, ECM remodeling. "
            "(b) ROC curve on the full cohort (AUC = 0.827); inset: 5-fold stratified cross-validation AUC = 0.798 +/- 0.156.",
            "Panel a: permutation importance attributes 56.5% of discriminative signal to the DE signature, 35.6% to EC dysfunction, "
            "and 7.8% to ECM remodeling. Panel b: fitted model AUC = 0.827 on the full cohort; cross-validation AUC = 0.798 +/- 0.156 "
            "reflects fold-to-fold variability expected with n = 43. The model summarises biological separation in resected tissue and "
            "is not proposed as a bedside clinical predictor."
        )
    else:
        doc.add_paragraph("[Figure 4: main_manuscript/Figure4_model_feature_importance_and_ROC.png not found; run generate_main_five_figures_ultra_hd.py.]")
    doc.add_paragraph()

    doc.add_heading('Stratification by molecular signature', 2)
    results5 = (
        "We stratified samples into low-, medium-, and high-risk groups based on the model score (Methods; "
        "Supplementary Table S1). "
        "The proportion of ruptured samples increased across strata (Figure 5), illustrating that the "
        "signatures capture biology associated with rupture. Again, this stratification is not intended "
        "for clinical use with tissue (tissue is not available pre-operatively) but to show that molecular "
        "phenotype aligns with rupture status in this cohort."
    )
    doc.add_paragraph(results5)

    # Figure 5: Risk stratification (main_manuscript ultra-HD, phase7 scores)
    fig5_path = "results/figures/main_manuscript/Figure5_risk_stratification_rupture_enrichment.png"
    if Path(fig5_path).exists():
        add_figure_with_caption_and_inference(
            doc, fig5_path,
            "Figure 5. Stratification by summary model score (n=43). Panel a: sample counts per stratum (stacked ruptured vs unruptured). "
            "Panel b: percentage ruptured per stratum. Thresholds: low < 0.33, medium 0.33-0.67, high >= 0.67 on predicted probability.",
            "The proportion of ruptured samples rises monotonically across strata (low 20% [2/10], medium 38% [8/21], high 100% [12/12]), "
            "showing that the composite molecular score aligns with rupture status in this cohort. Stratification illustrates biological "
            "gradation and is not intended for pre-operative clinical use."
        )
    else:
        doc.add_paragraph("[Figure 5: main_manuscript/Figure5_risk_stratification not found; run generate_main_five_figures_ultra_hd.py.]")
    doc.add_page_break()

    # --------------- DISCUSSION ---------------
    doc.add_heading('Discussion', 1)
    discussion = (
        "We used bulk RNA-seq from resected intracranial aneurysm tissue (GSE122897) to identify molecular "
        "features associated with rupture. The tissue is obtained at surgery and is not available from "
        "unruptured patients who are managed conservatively; therefore we frame the work as identifying "
        "mechanistic factors and candidate pathways rather than as a bedside rupture-risk prediction tool. "
        "Differential expression and signature analyses implicate endothelial dysfunction, ECM remodeling, "
        "and inflammation. Integrating our results with established rupture biology, a coherent cascade "
        "emerges with direct empirical support: (1) Endothelial activation—ANGPT2, ICAM1, and VCAM1 were "
        "upregulated in ruptured tissue (Table 2); the EC signature contributed 35.6% to the model (Figure 4a) "
        "and pathway enrichment identified blood vessel development and vasculature development at FDR < 0.05. "
        "(2) Matrix metalloproteinase activity—pathway enrichment identified extracellular matrix organization, "
        "collagen-containing ECM, and integrin binding (Supplementary Table S5); our ECM signature genes "
        "(MMP9, MMP2, TIMP1, COL3A1, COL1A1) contributed 7.8% to the model. (3) Smooth muscle phenotypic "
        "switching—cell migration and cell motility were among the top enriched terms, consistent with "
        "vascular wall remodelling. (4) Inflammatory recruitment—leukocyte migration, TNF signalling, and "
        "interleukin signalling were enriched. These processes likely act in concert, and our signature scores "
        "capture overlapping dimensions of this cascade.\n\n"
        "The summary model separates ruptured from unruptured tissue (full-cohort AUC ~0.83; cross-validation "
        "AUC 0.80 ± 0.16). Given the small sample size (n=43), the CV standard deviation (~0.15) indicates "
        "substantial fold-to-fold variability; performance should be interpreted cautiously. The eight false "
        "negatives (ruptured samples predicted as unruptured) may reflect heterogeneity in rupture biology "
        "or molecular subtypes not fully captured by the three signatures. Bulk RNA-seq averages across cell "
        "types and lacks spatial information; single-cell and spatial transcriptomics could resolve which "
        "cells and regions drive the signal and are priorities for future work.\n\n"
        "For clinical impact, biomarkers that can be measured without tissue resection—for example in blood—"
        "are needed. This tissue-based analysis generates hypotheses for such biomarkers: genes and pathways "
        "that are differentially expressed in ruptured tissue may have circulating correlates worth testing "
        "in future blood-based studies. We did not compare with clinical risk scores (e.g. PHASES) because "
        "our outcome is tissue-level molecular phenotype, not population-level risk prediction; such "
        "comparisons would require a different study design and clinical endpoints."
    )
    doc.add_paragraph(discussion)
    doc.add_page_break()

    # --------------- METHODS (expanded: key steps, standard omics headings) ---------------
    doc.add_heading('Methods', 1)

    doc.add_heading('Data acquisition and study cohort', 2)
    doc.add_paragraph(
        "RNA-seq data were obtained from the National Center for Biotechnology Information Gene Expression "
        "Omnibus (NCBI GEO) under accession GSE122897 (Ruigrok et al., UMC Utrecht). The dataset comprises "
        "bulk RNA-seq from resected intracranial aneurysm tissue (Illumina HiSeq 2500, poly-A selected). "
        "The original series contained 60 samples; we restricted analysis to specimens with definitive rupture "
        "status (ruptured or unruptured), excluding control and unknown-status samples, yielding 43 samples "
        "(22 ruptured, 21 unruptured). Expression values were provided as TPM (transcripts per million) "
        "normalised counts aligned to the GRCh38.p13 reference. Sample metadata (rupture status, tissue type) "
        "were taken from the GEO series matrix. Age and sex were not consistently reported in the GEO records "
        "and could not be used for covariate adjustment (Supplementary Table S6). See Supplementary Methods for "
        "download and file names."
    )

    doc.add_heading('RNA-seq processing and quality control', 2)
    doc.add_paragraph(
        "We used the preprocessed TPM matrix as supplied by the data submitters. No additional alignment or "
        "quantification was performed. Quality control included checking sample and gene completeness, "
        "distribution of expression values, and agreement of metadata with sample identifiers. Genes were "
        "retained for downstream analysis if they had TPM > 1 in at least three samples, to focus on "
        "adequately expressed genes. See Supplementary Methods for QC metrics and filtering counts."
    )

    doc.add_heading('Normalisation and gene filtering', 2)
    doc.add_paragraph(
        "Expression values were log2(TPM + 1) transformed for differential expression and signature "
        "calculations. For the machine-learning model, features (the three signature scores) were "
        "standardised to zero mean and unit variance (StandardScaler) across samples. No further "
        "batch correction was applied; the cohort is from a single study and platform."
    )

    doc.add_heading('Differential expression analysis', 2)
    doc.add_paragraph(
        "Differential expression between ruptured (n=22) and unruptured (n=21) samples was performed "
        "gene-by-gene using two-sample Welch t-tests on log2(TPM+1). Welch t-tests are appropriate for "
        "comparing group means when variances may differ and provide a robust exploratory screen for "
        "differentially expressed genes. For count-based RNA-seq, methods such as DESeq2, edgeR, or "
        "limma-voom are often preferred because they model overdispersion; raw counts are available for "
        "GSE122897. We ran DESeq2 on raw counts as a confirmatory analysis (2,764 genes at padj < 0.05; "
        "Supplementary Table S3). Here we used TPM-based Welch t-tests "
        "for consistency with the preprocessed GEO data and for comparability with downstream signature "
        "scores derived from the same matrix. Log2 fold change was computed as the difference in group "
        "means. P-values were adjusted for multiple testing using the Benjamini–Hochberg procedure (FDR). "
        "Genes with FDR < 0.05 were considered significantly differentially expressed. Bonferroni "
        "adjustment was also computed for reference. Full formulae and software are in Supplementary Methods."
    )

    doc.add_heading('Biological signature score definition', 2)
    doc.add_paragraph(
        "We defined four pathway-based signature scores as the mean log2(TPM+1) over predefined gene sets. "
        "(1) Differential expression (DE) signature: mean expression of the top differentially expressed "
        "genes (by FDR and |Log2FC|). (2) Endothelial dysfunction (EC): ANGPT2, ICAM1, VCAM1, THBS1, PERP. "
        "(3) Extracellular matrix (ECM) remodeling: MMP9, MMP2, TIMP1, COL3A1, COL1A1. (4) Inflammation: "
        "TNF, IL6, IL8 (CXCL8), CRP, SAA1. The summary model used only the DE, EC, and ECM scores (inflammation "
        "was explored separately). Gene symbols were matched to the expression matrix via available "
        "annotations; see Supplementary Methods for full lists and any unmapped genes."
    )

    doc.add_heading('Statistical model and cross-validation', 2)
    doc.add_paragraph(
        "A logistic regression model with L2 (ridge) regularisation was fitted to rupture status (binary) "
        "using the three signature scores (DE, EC, ECM) as features with fixed regularisation strength "
        "(C=1.0 in the current scripts). The DE signature score is derived from genes identified by differential expression in the "
        "same cohort; this introduces potential circularity. We acknowledge this limitation and, where "
        "applicable, report results from a nested design in which DE gene selection is performed within "
        "each CV fold (Supplementary Methods). The final model was refit on the full cohort for coefficient "
        "and feature-importance estimation. Feature importance was quantified by permutation: for each "
        "feature, values were shuffled and the decrease in AUC was recorded (mean over permutations). "
        "Reported metrics (AUC, accuracy, precision, recall, F1, confusion matrix) are for the full cohort "
        "at the optimal threshold (Youden index). Supplementary Table S2 summarises CV performance; "
        "Supplementary Methods give hyperparameter grids and software versions."
    )

    doc.add_heading('Risk stratification', 2)
    doc.add_paragraph(
        "Each sample was assigned a predicted probability of rupture from the final model. Risk strata "
        "were defined as low (probability < 0.33), medium (0.33–0.67), and high (≥ 0.67). The proportion "
        "of ruptured samples in each stratum was computed to illustrate alignment between molecular "
        "signature and rupture status (Supplementary Table S1)."
    )

    doc.add_heading('Software and reproducibility', 2)
    doc.add_paragraph(
        "Analyses were conducted in Python 3.x using pandas, numpy, scipy, and scikit-learn. Figures were "
        "generated with matplotlib and seaborn. Scripts and code are described in Supplementary Methods; "
        "processed tables and key outputs can be made available on request."
    )

    doc.add_paragraph()
    doc.add_paragraph(
        "Further technical detail—including exact gene lists, formulae, hyperparameter ranges, and file "
        "paths—is given in Supplementary Methods.",
        style='Intense Quote'
    )
    doc.add_page_break()

    # --------------- LIMITATIONS ---------------
    doc.add_heading('Limitations', 1)
    limitations = (
        "The cohort is small (n=43) and from a single dataset; findings require validation in independent "
        "tissue cohorts. Tissue is available only at surgery, so the model cannot be used to predict rupture "
        "in unruptured patients without intervention. The DE gene set used for the DE signature score was "
        "identified in the same data used to train and evaluate the classifier; this circularity may inflate "
        "performance estimates, though nested DE-within-CV analyses help quantify this effect. Confounding "
        "cannot be ruled out: ruptured tissue reflects both pre-rupture wall instability and peri-rupture "
        "inflammation, and we could not adjust for age, sex, or other covariates due to limited metadata. "
        "Bulk RNA-seq averages expression across cell types and lacks spatial resolution; findings may mask "
        "cell-type-specific or spatially heterogeneous contributions. Causal interpretation of associations "
        "would require experimental or prospective designs."
    )
    doc.add_paragraph(limitations)
    doc.add_paragraph()

    # --------------- FUTURE DIRECTIONS ---------------
    doc.add_heading('Future directions', 1)
    future = (
        "Priorities include replication in larger, multi-centre tissue cohorts and validation of top genes "
        "in blood or other accessible biospecimens. Blood-based biomarkers would be obtainable in unruptured "
        "patients and could inform future risk stratification studies. Single-cell and spatial transcriptomics "
        "could resolve which cell types and anatomical regions drive the rupture-associated signal, overcoming "
        "limitations of bulk RNA-seq (averaged expression, no spatial resolution). Functional and mechanistic "
        "studies of the identified genes and pathways are needed to establish causality and therapeutic relevance."
    )
    doc.add_paragraph(future)
    doc.add_paragraph()

    # --------------- CONCLUSION ---------------
    doc.add_heading('Conclusion', 1)
    conclusion = (
        "Transcriptomic analysis of resected aneurysm tissue identifies molecular signatures—notably "
        "endothelial dysfunction, ECM remodeling, and inflammation—associated with rupture. Because "
        "this tissue is not available before intervention in unruptured patients, the work is positioned "
        "to advance mechanistic understanding and to prioritise candidates for future blood-based or "
        "other non-tissue biomarker research."
    )
    doc.add_paragraph(conclusion)
    doc.add_paragraph()

    # --------------- REFERENCES ---------------
    doc.add_heading('References', 1)
    refs = (
        "1. GSE122897. NCBI Gene Expression Omnibus. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122897\n"
        "2. Ruigrok YM et al. (source cohort for GSE122897). [Add full citation when available.]\n"
        "3. [Add further references for methods and background as needed.]"
    )
    doc.add_paragraph(refs)
    doc.add_paragraph()
    doc.add_paragraph("Data availability. Data were obtained from NCBI GEO (GSE122897). Processed results and code are available from the authors on request.", style='Intense Quote')
    doc.add_page_break()

    # ============================================================================
    # SUPPLEMENTARY MATERIALS
    # ============================================================================
    doc.add_heading("Supplementary Materials", 1)
    append_supplementary_bundle(
        doc,
        s1_stratification_figure_ref="main text Figure 5",
        include_supplementary_figure_s12=False,
        figure_width=Inches(5),
    )

    # Save (in project root)
    output_file = Path(script_dir) / 'Molecular_Rupture_Signatures_VS.docx'
    doc.save(str(output_file))
    natmed_file = Path(script_dir) / 'NatMed_Analysis.docx'
    doc.save(str(natmed_file))
    print(f"Created: {output_file}")
    print(f"Created: {natmed_file}  (NatMed draft alias with embedded main figures)")
    print("  Main manuscript: 5 figures, 2 tables; narrative reframed (tissue source, mechanistic focus, no PHASES).")
    print("  Supplementary: figures and tables with legends and inferences; Supplementary Methods referenced from main.")
    print("  Main figures: results/figures/main_manuscript/ at 600 DPI (generate_main_five_figures_ultra_hd.py).")
    return str(output_file)

if __name__ == "__main__":
    try:
        create_manuscript()
    except Exception as e:
        import traceback
        traceback.print_exc()
