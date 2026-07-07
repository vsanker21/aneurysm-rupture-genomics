"""
Molecular_Signatures_Revised_VS — Nature-style manuscript (revised story arc).

Main text: **7 figures** and **5 tables**. All other panels (including risk stratification)
live in **Supplementary Materials** with explicit in-text pointers.

Narrative (single arc): **initial bulk DE** → **CIBERSORTx Job 7 + LM22/NNLS +
immune-adjusted re-testing** → **pathway / signature scores + ridge tissue model** →
**combined external validation hierarchy** → **external sample-level signature transfer +
direct-rupture meta-analysis** → **PHASES–G** (eQTL, GWAS rsID overlap, EXPRS placeholder,
prespecified germline PRS when genotypes exist). Key numbers load from
`manuscript_revised_vs_metrics.py` (JSON under `results/` with safe fallbacks).

Output: Molecular_Signatures_Revised_VS.docx (project root).

Re-run after refreshing figures. Requires: python-docx, pandas.

  python create_molecular_signatures_revised_vs.py

Supplementary Methods, Supplementary Figures S1–S12, and embedded Supplementary
Tables S1–S6 are produced in one run via append_supplementary_bundle() from
create_molecular_rupture_signatures_manuscript.py.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

from combined_validation_summary import build_combined_validation_assets
from create_molecular_rupture_signatures_manuscript import (
    add_figure_with_caption_and_inference,
    add_table_from_dataframe,
    append_supplementary_bundle,
)
from manuscript_revised_vs_metrics import load_story_metrics
from manuscript_results_detail import (
    deseq2_top_prose_sentence,
    gse159610_shortlist_stats,
    load_catalogue_eqtl_examples,
    load_cibersortx_merged_exemplars,
    load_deseq2_top_table,
    load_immune_adjusted_exemplars,
    load_lm22_contrast_table,
    load_phases_g_genetics_table,
)


def create_revised_manuscript() -> str:
    script_dir = Path(__file__).resolve().parent
    os.chdir(script_dir)
    _base = script_dir
    _m = load_story_metrics(_base)
    _gse_stats = gse159610_shortlist_stats(_base)
    _combined_validation_df = build_combined_validation_assets(_base)
    _transfer_summary_path = _base / "results" / "advanced_validation" / "external_signature_transfer_summary.tsv"
    _transfer_summary_df = pd.read_csv(_transfer_summary_path, sep="\t") if _transfer_summary_path.is_file() else pd.DataFrame()
    _transfer_by_dataset = (
        {str(r["Dataset"]): r for r in _transfer_summary_df.to_dict(orient="records")}
        if not _transfer_summary_df.empty and "Dataset" in _transfer_summary_df.columns
        else {}
    )
    _meta_top_path = _base / "results" / "rupture_meta_analysis" / "rupture_consensus_meta_top25.tsv"
    _meta_top_df = pd.read_csv(_meta_top_path, sep="\t") if _meta_top_path.is_file() else pd.DataFrame()
    _meta_summary_path = _base / "results" / "rupture_meta_analysis" / "rupture_consensus_meta_summary.json"
    _meta_summary = json.loads(_meta_summary_path.read_text(encoding="utf-8")) if _meta_summary_path.is_file() else {}
    _t133 = _transfer_by_dataset.get("GSE13353", {})
    _t540 = _transfer_by_dataset.get("GSE54083", {})
    _t367 = _transfer_by_dataset.get("GSE36791", {})
    _t159 = _transfer_by_dataset.get("GSE159610", {})
    _t269 = _transfer_by_dataset.get("GSE26969", {})

    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "Calibri"
    style.font.size = Pt(11)

    # --- Title ---
    t = doc.add_heading(
        "Integrated tissue signatures and a genetics-informed framework for aneurysm rupture biology",
        0,
    )
    t.alignment = WD_ALIGN_PARAGRAPH.CENTER
    sub = doc.add_paragraph(
        "Transcriptomic analysis of resected aneurysm tissue linked to a prespecified PHASES–genetics integration path (PHASES-G)"
    )
    sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
    if sub.runs:
        sub.runs[0].font.size = Pt(12)
        sub.runs[0].italic = True
    doc.add_page_break()

    # --- Abstract ---
    doc.add_heading("Abstract", 1)
    doc.add_paragraph(
        f"Intracranial aneurysm rupture reflects arterial wall biology that bulk RNA-seq can quantify when tissue "
        f"is obtained at surgery. In GSE122897 (n={_m['n_samples']}; {_m['n_ruptured']} ruptured, {_m['n_unruptured']} unruptured) "
        f"we first mapped rupture-associated expression: DESeq2 identified {_m['deseq2_padj_genes']:,} genes with false-discovery rate "
        f"(FDR) <0.05, while Welch tests on log₂(TPM+1) yielded {_m['welch_fdr_genes']} genes at transcriptome-wide FDR in this cohort "
        f"(Fig. 2; Table 2; Supplementary Tables S3–S4). We then used LM22-based immune deconvolution: **CIBERSORTx** portal output "
        f"(**Job 7**, archived as `CIBERSORTx_Job7_Results.txt`, {_m['n_cibersortx_job7_mixtures']} mixtures) and, for the immune-fraction "
        f"covariates in the adjusted model, reproducible **LM22+NNLS** proportions on GSE122897 TPM covering all {_m['n_samples']} specimens "
        f"(`lm22_fractions_gse122897.tsv`; {_m['n_genes_lm22_intersect']} LM22 signature genes detected). We re-modelled expression with "
        f"{_m['n_immune_covariates']} selected fraction covariates; among {_m['n_genes_immune_model_tested']:,} tested genes "
        f"no rupture coefficient survived Benjamini–Hochberg FDR<0.05 across that tested set, yet {_m['n_deseq2_fdr_among_immune_tested']:,} "
        f"genes retained DESeq2 FDR<0.05 within it and {_m['n_genes_deseq2_and_nominal_rupture_immune']:,} showed joint DESeq2 significance "
        f"and nominal rupture association post-adjustment—filters that feed a capped mechanistic shortlist of {_m['n_shortlist']} genes. "
        f"Pathway enrichment (g:Profiler) yielded {_m['pathways_fdr05']} terms at FDR<0.05. A ridge logistic classifier on three "
        f"predefined signature scores achieved AUROC {_m['auc_full_cohort']:.3f} on the full cohort, {_m['auc_cv_mean']:.3f}±{_m['auc_cv_sd']:.3f} "
        f"under optimistic non-nested cross-validation, and a more conservative leakage-controlled nested estimate of "
        f"{_m['auc_nested_mean']:.3f}±{_m['auc_nested_sd']:.3f} (Supplementary Fig. S11). "
        f"In independent rupture microarray GSE54083 ({_m['gse54083_n_ruptured']} ruptured, {_m['gse54083_n_unruptured']} unruptured), "
        f"{_m['gse54083_n_shortlist_present']} shortlisted genes were measurable, {_m['gse54083_n_shortlist_direction_match']} matched the "
        f"GSE122897 rupture direction, and {_m['gse54083_n_shortlist_p_lt_0_05']} were nominally significant (none at FDR<0.05); "
        f"a targeted seven-gene check in GSE13353 showed {_m['gse13353_direction_matches']}/{_m['gse13353_genes_tested']} directional matches "
        f"({_m['gse13353_nominal_sig_genes']} nominally significant). GSE26969 and GSE36791 were added as orthogonal tissue and blood context, "
        f"and external sample-level transfer of the fixed three-score model was supportive in GSE13353 (AUROC {float(_t133.get('AUC', float('nan'))):.3f}) "
        f"but failed in GSE54083 (AUROC {float(_t540.get('AUC', float('nan'))):.3f}). A weighted direct-rupture meta-analysis identified "
        f"{int(_meta_summary.get('n_directionally_concordant_all_three', 0))} directionally concordant genes at meta-analysis FDR<0.05. GSE159610 is retained as orthogonal aneurysm-versus-control context, "
        f"not rupture validation. Finally, **PHASES–G** connects the shortlist to blood cis-eQTL evidence (eQTL Catalogue: {_m['n_blood_eqtl_catalogue']} genes; "
        f"eQTLGen cis: {_m['n_eqtlgen_cis']}) and rsID-level overlap with intracranial aneurysm GWAS ({_m['n_gwas_overlap_catalogue']} / "
        f"{_m['n_gwas_overlap_eqtlgen']} sentinel SNPs for Catalogue / eQTLGen layers), plus a prespecified germline polygenic-score path "
        f"({_m['n_prs_weights_snps']} SNP weights on the eQTL panel when dosages exist) and an expression-based eQTL-weighted score "
        f"(EXPRS; {_m['n_exprs_genes']} genes) where genotypes are absent. Tissue classifiers are not pre-operative bedside tools; they "
        f"prioritise mechanisms, immune-aware gene lists, orthogonal bulk checks, and a reproducible bridge to genetics and prospective "
        f"PHASES co-measurement (Supplementary Materials)."
    )
    doc.add_page_break()

    # --- Introduction ---
    doc.add_heading("Introduction", 1)
    doc.add_paragraph(
        "Ruptured intracranial aneurysms cause substantial morbidity and mortality. Clinical scores such as PHASES "
        "summarise anatomy and history but do not, by themselves, integrate germline or circulating molecular evidence. "
        "Resected aneurysm wall RNA offers a direct readout of biology at intervention, yet tissue is unavailable for "
        "the large population with unruptured lesions managed conservatively. The scientific payoff of bulk surgical "
        "transcriptomes is therefore **mechanistic prioritisation**: quantifying which genes and programmes track "
        "rupture in sampled wall, stress-testing those signals against immune admixture, asking whether they echo in "
        "orthogonal bulk cohorts, and only then layering genetics-informed bridges that can, in principle, travel "
        "with DNA and clinical covariates on the same individuals.\n\n"
        "Here we re-analyse GSE122897 as a single connected narrative. We begin with conventional differential "
        "expression to establish the magnitude and genome-wide structure of rupture-associated shifts (Table 2; Fig. 2). "
        "We next integrate **CIBERSORTx Job 7** mixture tables (archived under `results/cibersortx/`) with reproducible "
        "LM22+NNLS fractions for covariate adjustment, then apply immune-fraction-adjusted re-testing to separate wall-resident "
        "programmes from signals plausibly driven by leukocyte content, and we propagate the resulting immune-aware "
        f"gene prioritisation into a capped shortlist ({_m['n_shortlist']} genes) used throughout downstream integration. "
        "We then collapse expression into interpretable pathway signatures and a regularised tissue classifier with "
        "explicit cross-validation reporting (Figs. 3–4; Supplementary Fig. S11). Finally, we connect that shortlist "
        "to independent rupture replication in GSE54083 and GSE13353, add GSE26969 and GSE36791 as orthogonal context, "
        "test external score transfer and direct-rupture gene-level consensus, use GSE159610 only as "
        "aneurysm-versus-control orthogonal context, and then extend into **PHASES–G**, a prespecified "
        "stack of blood cis-eQTL queries, IA GWAS rsID overlap, optional germline PRS when genotypes exist, and EXPRS "
        "as an RNA-only placeholder—fully documented in Supplementary Note 1 so that genetics steps are not conflated "
        "with tissue labels."
    )
    doc.add_page_break()

    # --- Results ---
    doc.add_heading("Results", 1)
    doc.add_paragraph(
        "We organised results to follow the analysis pipeline: bulk differential expression in GSE122897; "
        "CIBERSORTx-based immune deconvolution and immune-fraction-adjusted re-testing; pathway- and signature-level "
        "summaries with multivariate tissue discrimination; combined external validation and orthogonal context; "
        "external sample-level transfer plus direct-rupture consensus analysis; and the PHASES–genetics "
        "(PHASES–G) evidence stack linking the immune-aware shortlist to blood eQTLs, GWAS, and optional germline scoring."
    )

    doc.add_heading("Cohort, data provenance, and analysis overview", 2)
    doc.add_paragraph(
        "We used 43 specimens (22 ruptured, 21 unruptured) from GSE122897 after excluding non-informative statuses "
        "(Methods). TPM-normalised expression (GRCh38.p13) supports differential expression, predefined signature "
        "scores, and a regularised logistic classifier. **Table 1** summarises the cohort; **Figure 1** illustrates "
        "design and sample composition. Extended sample metadata appear in **Supplementary Table S6**."
    )

    # Table 1 (same structure as VS manuscript)
    doc.add_paragraph(
        "Table 1. Study sample and RNA-seq data (GSE122897). Resected intracranial aneurysm tissue; surgery-derived; "
        "not available from unruptured patients without intervention.",
        style="Heading 3",
    )
    sample_table_data = {
        "Item": [
            "Dataset",
            "Tissue type",
            "Total samples (analysed)",
            "Ruptured",
            "Unruptured",
            "Platform",
            "Normalisation",
            "Reference",
        ],
        "Description": [
            "GSE122897 (NCBI GEO)",
            "Resected intracranial aneurysm tissue",
            "43",
            "22",
            "21",
            "Illumina HiSeq 2500, RNA-seq",
            "TPM (GRCh38.p13)",
            "Ruigrok et al., UMC Utrecht",
        ],
    }
    add_table_from_dataframe(
        doc,
        pd.DataFrame(sample_table_data),
        title=None,
        inference=(
            "Balanced rupture/unruptured groups support two-class comparisons. The surgical origin of tissue defines "
            "the scope of inference (wall biology at resection, not pre-operative risk in the general UIA population)."
        ),
        max_rows=20,
    )
    doc.add_paragraph()

    fig1 = "results/figures/dataset_coverage.png"
    if Path(fig1).exists():
        add_figure_with_caption_and_inference(
            doc,
            fig1,
            "Figure 1 | Study design and sample composition. Data source, tissue type, and class sizes (n = 43).",
            "Design figure orients readers to the surgical sampling frame and the distinction between tissue-level "
            "discovery and clinical prediction in unruptured patients.",
        )
    else:
        doc.add_paragraph("[Figure 1: results/figures/dataset_coverage.png not found.]")
    doc.add_paragraph()

    doc.add_heading("Initial bulk analysis: thousands of rupture-associated genes", 2)
    doc.add_paragraph(
        f"We first performed standard two-group comparisons on all specimens passing QC. DESeq2 (raw counts; "
        f"ruptured versus unruptured) identified {_m['deseq2_padj_genes']:,} genes with FDR-adjusted P (padj) <0.05—i.e., "
        f"a broad transcriptional remodelling rather than a handful of isolated hits. In parallel, Welch t-tests on "
        f"log₂(TPM+1) highlighted the strongest univariate shifts; {_m['welch_fdr_genes']} genes met a transcriptome-wide "
        f"Benjamini–Hochberg FDR threshold in this cohort, with extensive overlap of top-ranked genes between methods "
        f"(**Table 2**; **Supplementary Tables S3–S4**; **Supplementary Figures S1, S2, S4, S10** for MA plots, heatmaps, "
        f"method comparisons, and an alternate DESeq2 volcano). **Figure 2** displays the genome-wide volcano view used "
        f"in the main text. This initial pass establishes that rupture status is associated with a strong, coordinated "
        f"expression programme in sampled wall—setting up the question of how much is immune-infiltration-driven versus "
        f"intrinsic vascular programmes."
    )
    _dseq_prose = deseq2_top_prose_sentence(_base, n=8)
    if _dseq_prose:
        doc.add_paragraph(_dseq_prose)

    try:
        top_genes_path = "results/gse122897_deseq2_top_annotated.tsv"
        if not Path(top_genes_path).exists():
            top_genes_path = "results/gse122897_top_genes_annotated.tsv"
        if not Path(top_genes_path).exists():
            top_genes_path = "results/gse122897_top_genes.tsv"
        top_genes_df = pd.read_csv(top_genes_path, sep="\t")
        if len(top_genes_df) > 0:
            rename_map = {
                "Gene_ID": "Gene",
                "P_adj": "P_adj_BH",
                "log2FoldChange": "Log2FC",
                "padj": "P_adj_BH",
                "baseMean": "baseMean",
            }
            top_genes_df = top_genes_df.rename(columns=rename_map)
            if "Gene_Symbol" not in top_genes_df.columns and "Gene" in top_genes_df.columns:
                map_path = "results/gene_id_mapping_unified.tsv"
                if Path(map_path).exists():
                    try:
                        m = pd.read_csv(map_path, sep="\t")
                        if "Original_ID" in m.columns and "Gene_Symbol" in m.columns:
                            m = m[["Original_ID", "Gene_Symbol"]].dropna(subset=["Gene_Symbol"])
                            m["Gene_key"] = m["Original_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
                            top_genes_df["Gene_key"] = top_genes_df["Gene"].astype(str).str.replace(
                                r"\.0$", "", regex=True
                            )
                            top_genes_df = top_genes_df.merge(
                                m[["Gene_key", "Gene_Symbol"]], on="Gene_key", how="left"
                            ).drop(columns=["Gene_key"], errors="ignore")
                    except Exception:
                        pass
            display_cols = []
            for c in [
                "Gene_Symbol",
                "Gene_Name",
                "Log2FC",
                "P_adj_BH",
                "T_statistic",
                "Mean_Ruptured",
                "Mean_Unruptured",
                "Gene",
            ]:
                if c in top_genes_df.columns and c not in display_cols:
                    display_cols.append(c)
            if not display_cols:
                display_cols = list(top_genes_df.columns)
            top_genes_df = top_genes_df[display_cols].head(15)
            add_table_from_dataframe(
                doc,
                top_genes_df,
                title="Table 2 | Top differentially expressed genes (DESeq2; ruptured vs unruptured).",
                inference=(
                    "Representative strongest DESeq2 associations with resolved gene symbols where available; full DESeq2 lists in Supplementary Table S3."
                ),
                max_rows=20,
            )
    except Exception as e:
        doc.add_paragraph(f"Table 2. (Error loading top genes: {e})")
    doc.add_paragraph()

    _dseq_tbl = load_deseq2_top_table(_base, n=12)
    if not _dseq_tbl.empty:
        doc.add_paragraph(
            "**Numeric detail | DESeq2 genes with lowest padj (GSE122897, ruptured vs unruptured).** "
            "Welch-ranked genes in Table 2 emphasise univariate TPM contrasts; this panel lists the top DESeq2 hits by "
            "adjusted P-value for direct comparison with immune-adjusted models below.",
            style="Heading 3",
        )
        add_table_from_dataframe(
            doc,
            _dseq_tbl,
            title=None,
            inference="Full DESeq2 output: **Supplementary Table S3**. Gene symbols from `gene_id_mapping_unified.tsv` where available.",
            max_rows=20,
        )
        doc.add_paragraph()

    fig2 = "results/figures/ultra_hd/volcano_plot_ultra_hd.png"
    if Path(fig2).exists():
        add_figure_with_caption_and_inference(
            doc,
            fig2,
            "Figure 2 | Volcano plot of differential expression (n = 22 vs 21). Effect size vs significance (FDR).",
            "Highlights genome-wide separation; numerical detail in Table 2 and Supplementary Tables.",
        )
    else:
        doc.add_paragraph("[Figure 2: volcano_plot_ultra_hd.png not found.]")
    doc.add_paragraph()

    doc.add_heading("CIBERSORTx Job 7 mixtures, LM22+NNLS fractions, and immune-fraction-adjusted re-testing", 2)
    doc.add_paragraph(
        f"Because bulk aneurysm wall RNA integrates vascular cells and infiltrating leukocytes, we next asked how inferred "
        f"immune composition relates to rupture and whether gene–rupture associations persist after linear adjustment for "
        f"those fractions. **CIBERSORTx** (LM22 signature basis) was executed via the CIBERSORTx portal; the returned mixture matrix "
        f"for **Job 7** is archived as `results/cibersortx/CIBERSORTx_Job7_Results.txt` (a duplicate may also sit at the "
        f"repository root). That export contains {_m['n_cibersortx_job7_mixtures']} samples with per-sample CIBERSORTx P-value, "
        f"Pearson correlation, and RMSE. Pat-style mixture labels were rewritten to GEO GSM identifiers using "
        f"`rewrite_cibersortx_mixture_to_gsm.py`, yielding `CIBERSORTx_Job7_Results_GSM_mixture.txt`; parsed QC medians in this "
        f"workspace are correlation {_m['cibersortx_median_correlation']:.3f} and RMSE {_m['cibersortx_median_rmse']:.3f} "
        f"(`results/cibersortx/cibersortx_parse_summary.json`; n={_m['n_cibersortx_gsm_parsed']} GSM rows). "
        f"Separately, to feed the genome-wide OLS adjustment with proportions aligned **one-to-one** to all {_m['n_samples']} "
        f"RNA profiles, we used `lm22_fractions_gse122897.tsv` from `immune_deconvolution_lm22_gse122897.py`—the same {_m['lm22_cell_types']}-column "
        f"LM22 basis deconvolved by constrained non-negative least squares on log-transformed TPM ({_m['n_genes_lm22_intersect']} "
        f"LM22 genes detected in the mixture). Group contrasts on those NNLS fractions (Welch tests with FDR) are summarised in "
        f"`lm22_deconvolution_summary.json`: {_m['lm22_fdr_sig_cell_types_005']} cell types reached FDR<0.05, consistent with "
        f"subtle composition shifts rather than one dominant infiltrate class. We then regressed log₂(TPM+1) on rupture status "
        f"while including {_m['n_immune_covariates']} selected LM22 fraction covariates (`deg_immune_adjusted_gse122897.py`; "
        f"**Supplementary Methods**). Among {_m['n_genes_immune_model_tested']:,} adequately expressed genes tested in this "
        f"framework, {_m['n_deseq2_fdr_among_immune_tested']:,} retained DESeq2 FDR<0.05 from the primary analysis, "
        f"{_m['n_genes_immune_nominal_rupture']:,} showed nominal (P<0.05) rupture coefficients after immune adjustment, and "
        f"critically **zero** genes passed Benjamini–Hochberg FDR<0.05 for the rupture coefficient when FDR was computed across "
        f"this tested gene set. For mechanistic prioritisation we nevertheless retained genes that are both DESeq2-significant "
        f"in tissue **and** nominally rupture-associated after immune adjustment (n={_m['n_genes_deseq2_and_nominal_rupture_immune']}), "
        f"then applied directional replication filters and a hard cap to arrive at the {_m['n_shortlist']}-gene shortlist used "
        f"in all subsequent **PHASES–G** panels (**Supplementary Note 1**; `results/phases_g/shortlist_evidence*.tsv`). "
        f"Deconvolution diagnostics include `results/cibersortx/figures_ultra_hd/lm22_fractions_by_group_ultra_hd.png`."
    )
    doc.add_paragraph()

    _lm22_tbl = load_lm22_contrast_table(_base)
    if not _lm22_tbl.empty:
        doc.add_paragraph(
            "**Numeric detail | LM22 cell-type fractions (NNLS on TPM): ruptured vs unruptured means and contrasts.** "
            "Proportions sum to one per sample; P and FDR from Welch tests across n=22 vs 21.",
            style="Heading 3",
        )
        add_table_from_dataframe(
            doc,
            _lm22_tbl,
            title=None,
            inference="No compartment reached FDR<0.05; largest nominal shifts involve macrophage and T-cell subsets.",
            max_rows=25,
        )
        doc.add_paragraph()

    _cx_merged = load_cibersortx_merged_exemplars(_base, n=8)
    if not _cx_merged.empty:
        doc.add_paragraph(
            "**Numeric detail | CIBERSORTx Job 7 mixtures (GSM-merged) — representative samples spanning the correlation range.** "
            "Values are portal fractions for the listed LM22 columns; P-value, Correlation, and RMSE are CIBERSORTx fit diagnostics.",
            style="Heading 3",
        )
        add_table_from_dataframe(
            doc,
            _cx_merged,
            title=None,
            inference="Full Job 7 matrix: `results/cibersortx/CIBERSORTx_Job7_Results.txt`. Merged table: `cibersortx_merged_with_group.tsv`.",
            max_rows=12,
        )
        doc.add_paragraph()

    _imm_ex = load_immune_adjusted_exemplars(_base, n=10)
    if not _imm_ex.empty:
        doc.add_paragraph(
            "**Numeric detail | Immune-fraction–adjusted model: genes with strongest nominal rupture coefficients among DESeq2 hits.** "
            "OLS on log₂(TPM+1) with selected LM22 covariates (`deg_immune_covariate_adjusted.tsv`).",
            style="Heading 3",
        )
        add_table_from_dataframe(
            doc,
            _imm_ex,
            title=None,
            inference="BH-FDR for the rupture coefficient across the tested gene set was zero at 0.05; LM P values are nominal.",
            max_rows=15,
        )
        doc.add_paragraph()

    doc.add_heading("Pathway signatures and multivariate tissue discrimination", 2)
    doc.add_paragraph(
        f"To translate single-gene lists into interpretable biology, we scored each specimen with predefined differential-expression, "
        f"endothelial, and extracellular-matrix signatures (Methods). The three composite scores separate ruptured from unruptured "
        f"tissue (**Figure 3**). Systematic enrichment on the ranked DESeq2 gene list (g:Profiler) returned {_m['pathways_fdr05']} "
        f"terms at FDR<0.05 (**Supplementary Table S5**; **Supplementary Figure S6**), and endothelial marker panels supporting "
        f"interpretation appear in **Supplementary Figure S5**. We then trained a ridge logistic regression on the three signature "
        f"scores alone. **Figure 4a** shows permutation-based feature importance: the differential-expression score contributes "
        f"the largest share (~{_m['de_feature_importance_pct']:.1f}%), with endothelial (~{_m['ec_feature_importance_pct']:.1f}%) "
        f"and matrix (~{_m['ecm_feature_importance_pct']:.1f}%) scores adding complementary signal. **Figure 4b** reports the "
        f"receiver-operating-characteristic curve for the fitted model on the full cohort (AUROC {_m['auc_full_cohort']:.3f}; "
        f"n={_m['n_samples']}). Non-nested stratified cross-validation yielded mean AUROC {_m['auc_cv_mean']:.3f}±{_m['auc_cv_sd']:.3f}, "
        f"whereas the leakage-controlled nested estimate was {_m['auc_nested_mean']:.3f}±{_m['auc_nested_sd']:.3f} "
        f"(**Supplementary Table S2**; **Supplementary Figure S11**). We treat the nested estimate as the more defensible "
        f"out-of-sample summary because DE-derived features are re-selected within outer folds. Confusion matrices, precision–recall curves, and alternate ROC "
        f"layouts are in **Supplementary Figures S3, S7–S9**; risk stratification by predicted probability is reserved for "
        f"**Supplementary Figure S12** so the main text stays within four composite figures."
    )
    fig3 = "results/figures/ultra_hd/feature_distributions_ultra_hd.png"
    if Path(fig3).exists():
        add_figure_with_caption_and_inference(
            doc,
            fig3,
            "Figure 3 | Signature score distributions by rupture status (DE, EC, ECM).",
            "Connects single-gene signals to interpretable programs; DE contributes most in the multivariate model below.",
        )
    else:
        doc.add_paragraph("[Figure 3: feature_distributions_ultra_hd.png not found.]")
    doc.add_paragraph()

    fig4a = "results/figures/ultra_hd/feature_importance_ultra_hd.png"
    fig4b = "results/figures/ultra_hd/roc_curve_ultra_hd.png"
    if Path(fig4a).exists():
        add_figure_with_caption_and_inference(
            doc,
            fig4a,
            "Figure 4a | Relative contribution of signature scores (permutation importance; Methods).",
            "DE signature dominates; EC and ECM add complementary information.",
        )
    if Path(fig4b).exists():
        add_figure_with_caption_and_inference(
            doc,
            fig4b,
            "Figure 4b | ROC curve for the summary model (full cohort, n = 43).",
            "Illustrates discrimination in this surgical dataset; nested CV in Supplementary Figure S11 is more conservative.",
        )
    if not Path(fig4a).exists():
        doc.add_paragraph("[Figure 4: feature_importance / ROC assets not found.]")
    doc.add_paragraph()

    doc.add_heading("External context: direct rupture support and orthogonal tissue-blood validation hierarchy", 2)
    _n_sl, _n_gse, _n_dir = _gse_stats
    doc.add_paragraph(
        f"Independent of the tissue classifier, we sought support in external datasets with a hierarchy of relevance. "
        f"First, **GSE54083** provides an independent rupture-versus-unrupture aneurysm microarray dataset "
        f"({_m['gse54083_n_ruptured']} ruptured, {_m['gse54083_n_unruptured']} unruptured; total "
        f"n={_m['gse54083_n_ruptured'] + _m['gse54083_n_unruptured']}). Within that small cohort, "
        f"{_m['gse54083_n_shortlist_present']} immune-aware shortlist genes were measurable, "
        f"{_m['gse54083_n_shortlist_direction_match']} matched the GSE122897 rupture direction, and "
        f"{_m['gse54083_n_shortlist_p_lt_0_05']} were nominally significant, although none survived FDR<0.05. "
        f"Second, a targeted seven-gene check in **GSE13353** showed {_m['gse13353_direction_matches']}/"
        f"{_m['gse13353_genes_tested']} direction matches with {_m['gse13353_nominal_sig_genes']} nominally significant genes. "
        f"Third, **GSE26969** benchmarks unruptured aneurysm wall against non-aneurysmal superficial temporal artery "
        f"({_m['gse26969_n_aneurysm']} aneurysm, {_m['gse26969_n_control']} control), where "
        f"{_m['gse26969_n_shortlist_direction_match']}/{_m['gse26969_n_shortlist_present']} measurable shortlist genes matched the "
        f"GSE122897 direction and {_m['gse26969_n_shortlist_p_lt_0_05']} were nominally significant, but none survived FDR. "
        f"Fourth, **GSE36791** provides blood rupture/SAH context "
        f"({_m['gse36791_n_ruptured']} ruptured, {_m['gse36791_n_control']} control), where "
        f"{_m['gse36791_n_shortlist_direction_match']}/{_m['gse36791_n_shortlist_present']} measurable shortlist genes matched direction, "
        f"{_m['gse36791_n_shortlist_p_lt_0_05']} were nominally significant, and {_m['gse36791_n_shortlist_fdr_lt_0_05']} survived FDR, "
        f"although this blood signal is biologically confounded by acute post-rupture systemic response. Finally, in **GSE159610** "
        f"(aneurysm vs control blood RNA-seq), {_n_gse} shortlist genes had a non-missing contrast and {_n_dir} shared the same sign as "
        f"GSE122897 rupture effects. We therefore use GSE26969, GSE36791, and GSE159610 as orthogonal context rather than direct "
        f"rupture replication. **Figure 5** and **Table 3** summarise this hierarchy."
    )
    _combined_validation_fig = _base / "results" / "figures" / "ultra_hd" / "combined_validation_summary_ultra_hd.png"
    if _combined_validation_fig.exists():
        add_figure_with_caption_and_inference(
            doc,
            _combined_validation_fig,
            "Figure 5 | Combined hierarchy of external transcriptomic support across public datasets.",
            "Panel A summarises direction concordance of measurable or targeted genes with the GSE122897 shortlist; panel B summarises nominal and FDR-significant genes within each cohort-specific comparison. Direct rupture datasets remain the most relevant for replication, whereas GSE26969, GSE36791, and GSE159610 provide orthogonal tissue or blood context only.",
            width=Inches(6.5),
        )
    _combined_validation_tbl = _combined_validation_df[
        [
            "Dataset",
            "Role",
            "Samples",
            "Measured_or_tested_genes",
            "Direction_matched_genes",
            "Direction_match_pct",
            "Nominal_P_lt_0_05",
            "FDR_lt_0_05",
            "Interpretation",
        ]
    ].rename(
        columns={
            "Measured_or_tested_genes": "Genes measured/tested",
            "Direction_matched_genes": "Dir. matched",
            "Direction_match_pct": "Dir. match (%)",
            "Nominal_P_lt_0_05": "Nominal P<0.05",
            "FDR_lt_0_05": "FDR<0.05",
        }
    )
    add_table_from_dataframe(
        doc,
        _combined_validation_tbl,
        title="Table 3 | Combined external-validation hierarchy across public datasets.",
        inference="Direct rupture replication remains limited by small public cohorts, but the broader context layers help distinguish true rupture support from aneurysm-versus-control tissue shifts and systemic blood response after rupture.",
        max_rows=10,
    )
    doc.add_paragraph()

    doc.add_heading("External sample-level transfer and direct-rupture consensus", 2)
    doc.add_paragraph(
        f"To move beyond gene-level overlap, we projected the fixed GSE122897 three-score model into each external cohort without "
        f"retraining. Transfer was heterogeneous. In the smaller direct rupture cohort **GSE13353**, discrimination was supportive "
        f"(AUROC {float(_t133.get('AUC', float('nan'))):.3f}, bootstrap 95% CI {float(_t133.get('AUC_bootstrap_ci_low', float('nan'))):.3f}–"
        f"{float(_t133.get('AUC_bootstrap_ci_high', float('nan'))):.3f}, permutation P={float(_t133.get('AUC_permutation_pvalue', float('nan'))):.4f}), "
        f"whereas **GSE54083** showed inverse separation (AUROC {float(_t540.get('AUC', float('nan'))):.3f}; mean probability "
        f"{float(_t540.get('Mean_probability_positive', float('nan'))):.3f} in ruptured vs "
        f"{float(_t540.get('Mean_probability_negative', float('nan'))):.3f} in unruptured). This discordance argues against stable portable "
        f"rupture discrimination from the current tissue-derived score set. Orthogonal blood rupture/SAH context in **GSE36791** showed only "
        f"modest transfer (AUROC {float(_t367.get('AUC', float('nan'))):.3f}), while aneurysm-versus-control context in **GSE159610** and "
        f"**GSE26969** remained near-null or inverse (AUROC {float(_t159.get('AUC', float('nan'))):.3f} and {float(_t269.get('AUC', float('nan'))):.3f}, respectively).\n\n"
        f"We therefore asked whether gene-level consensus was more reproducible than cohort-level prediction. Among "
        f"{int(_meta_summary.get('n_overlap_all_three', 0))} genes measurable across **GSE122897**, **GSE54083**, and **GSE13353**, "
        f"{int(_meta_summary.get('n_directionally_concordant_all_three', 0))} were directionally concordant across all three cohorts and "
        f"remained significant at meta-analysis FDR<0.05. The strongest concordant genes were TGFBI, ANPEP, MPP1, MARCO, and CLEC5A "
        f"(higher in ruptured), and CPXM2 (lower in ruptured), suggesting that a smaller rupture-associated molecular core is more "
        f"reproducible than the full cohort-level classifier. **Figure 6**, **Figure 7**, **Table 4**, and **Table 5** summarise these analyses."
    )
    _transfer_fig = _base / "results" / "figures" / "ultra_hd" / "external_signature_transfer_ultra_hd.png"
    if _transfer_fig.exists():
        add_figure_with_caption_and_inference(
            doc,
            _transfer_fig,
            "Figure 6 | External sample-level transfer of the fixed three-score rupture model across public cohorts.",
            "Supportive discrimination in GSE13353 did not generalise to GSE54083, indicating that score portability is unstable across the currently available direct rupture cohorts; orthogonal blood and tissue contexts provide biology, not clinical validation.",
            width=Inches(6.5),
        )
    if not _transfer_summary_df.empty:
        _transfer_tbl = _transfer_summary_df[
            [
                "Dataset",
                "Role",
                "Positive_group",
                "Negative_group",
                "Samples",
                "AUC",
                "AUC_bootstrap_ci_low",
                "AUC_bootstrap_ci_high",
                "AUC_permutation_pvalue",
                "Mean_probability_positive",
                "Mean_probability_negative",
            ]
        ].copy()
        _transfer_tbl["Comparison"] = _transfer_tbl["Positive_group"].astype(str) + " vs " + _transfer_tbl["Negative_group"].astype(str)
        _transfer_tbl["AUROC (95% CI)"] = _transfer_tbl.apply(
            lambda r: f"{float(r['AUC']):.3f} ({float(r['AUC_bootstrap_ci_low']):.3f}-{float(r['AUC_bootstrap_ci_high']):.3f})",
            axis=1,
        )
        _transfer_tbl["Permutation P"] = _transfer_tbl["AUC_permutation_pvalue"].map(lambda x: f"{float(x):.4f}")
        _transfer_tbl["Mean prob. positive"] = _transfer_tbl["Mean_probability_positive"].map(lambda x: f"{float(x):.3f}")
        _transfer_tbl["Mean prob. comparator"] = _transfer_tbl["Mean_probability_negative"].map(lambda x: f"{float(x):.3f}")
        _transfer_tbl["Interpretation"] = _transfer_tbl["Dataset"].map(
            {
                "GSE13353": "Supportive direct transfer in a small rupture cohort.",
                "GSE54083": "Inverse direct transfer despite gene-level overlap; portability is unstable.",
                "GSE36791": "Modest blood rupture/SAH signal; likely acute systemic context.",
                "GSE159610": "Near-null blood aneurysm-control context only.",
                "GSE26969": "Inverse tissue aneurysm-control context; not a rupture cohort.",
            }
        )
        _transfer_tbl = _transfer_tbl[
            [
                "Dataset",
                "Role",
                "Comparison",
                "Samples",
                "AUROC (95% CI)",
                "Permutation P",
                "Mean prob. positive",
                "Mean prob. comparator",
                "Interpretation",
            ]
        ]
        add_table_from_dataframe(
            doc,
            _transfer_tbl,
            title="Table 4 | External sample-level transfer of the fixed three-score rupture model.",
            inference="A single supportive external transfer result is outweighed by failure in the other direct rupture cohort, so these projections strengthen biological interpretation but not clinical validation.",
            max_rows=10,
        )
        doc.add_paragraph()
    _meta_fig = _base / "results" / "figures" / "ultra_hd" / "rupture_consensus_meta_ultra_hd.png"
    if _meta_fig.exists():
        add_figure_with_caption_and_inference(
            doc,
            _meta_fig,
            "Figure 7 | Weighted gene-level rupture consensus across the three direct public rupture cohorts.",
            "Although cohort-level score transfer was unstable, a smaller gene-level core remained directionally concordant across GSE122897, GSE54083, and GSE13353, led by TGFBI, ANPEP, MPP1, MARCO, CLEC5A, and CPXM2.",
            width=Inches(6.2),
        )
    if not _meta_top_df.empty:
        _meta_tbl = _meta_top_df[
            [
                "Gene",
                "Consensus_direction",
                "Meta_Z",
                "Meta_P",
                "Meta_FDR_BH",
                "Log2FC_GSE122897",
                "Log2FC_GSE54083",
                "Log2FC_GSE13353",
            ]
        ].copy()
        _meta_tbl["Meta_Z"] = _meta_tbl["Meta_Z"].map(lambda x: f"{float(x):.3f}")
        _meta_tbl["Meta_P"] = _meta_tbl["Meta_P"].map(lambda x: f"{float(x):.2e}")
        _meta_tbl["Meta_FDR_BH"] = _meta_tbl["Meta_FDR_BH"].map(lambda x: f"{float(x):.2e}")
        _meta_tbl["GSE122897 log2FC"] = _meta_tbl["Log2FC_GSE122897"].map(lambda x: f"{float(x):.3f}")
        _meta_tbl["GSE54083 log2FC"] = _meta_tbl["Log2FC_GSE54083"].map(lambda x: f"{float(x):.3f}")
        _meta_tbl["GSE13353 log2FC"] = _meta_tbl["Log2FC_GSE13353"].map(lambda x: f"{float(x):.3f}")
        _meta_tbl = _meta_tbl.rename(columns={"Consensus_direction": "Consensus direction"})[
            [
                "Gene",
                "Consensus direction",
                "Meta_Z",
                "Meta_P",
                "Meta_FDR_BH",
                "GSE122897 log2FC",
                "GSE54083 log2FC",
                "GSE13353 log2FC",
            ]
        ]
        add_table_from_dataframe(
            doc,
            _meta_tbl,
            title="Table 5 | Top directionally concordant genes in the three-cohort rupture meta-analysis.",
            inference="These genes are the most reproducible transcript-level rupture signals across the public direct rupture cohorts and are more stable than the transferable cohort-level classifier.",
            max_rows=10,
        )
        doc.add_paragraph()

    doc.add_heading("PHASES–G: blood eQTLs, IA GWAS overlap, EXPRS, and prespecified germline PRS", 2)
    doc.add_paragraph(
        f"We finally layer **PHASES–G**, a reproducible genetics-facing stack built on the {_m['n_shortlist']}-gene immune-aware "
        f"shortlist (not on reusing tissue labels as genotypes). Querying blood cis-eQTL resources, {_m['n_blood_eqtl_catalogue']} "
        f"short-listed genes had a sentinel cis hit in the eQTL Catalogue whole-blood panel (per-gene top hit under the pipeline’s "
        f"pagination cap; **Supplementary Note 1**), and {_m['n_eqtlgen_cis']} genes had cis evidence in the eQTLGen layer we "
        f"attached. For IA GWAS overlap we intersected those sentinel SNPs (rsID-level) with the Bakker et al. summary-statistic "
        f"file: {_m['n_gwas_overlap_catalogue']} Catalogue-derived SNPs and {_m['n_gwas_overlap_eqtlgen']} eQTLGen-derived SNPs "
        f"appear in the GWAS rsID list—useful for hypothesis generation but **not** formal colocalisation (no COLOC/SMR in this "
        f"public RNA-only setting). Where individual-level genotypes and PHASES covariates exist on harmonised sample IDs, the "
        f"repository implements LD-aware clumping and a germline polygenic score using {_m['n_prs_weights_snps']} SNP weights "
        f"derived from the eQTL panel (`prs_rsid_weights_from_panel.tsv`). Public GEO submissions lack dosages; as an RNA-only "
        f"placeholder we also compute EXPRS—eQTL effect sizes used as weights on z-scored log-expression—for {_m['n_exprs_genes']} "
        f"short-listed genes with TPM support. EXPRS and nested-cohort AUROC demos are explicitly labelled as **pipeline "
        f"illustrations** in `PHASES_G_PRACTICAL_ANALYSES_SUMMARY.json` because the shortlist is outcome-informed in GSE122897; "
        f"they must not be read as independent external validity. Prospective validation requires DNA banking, PHASES fields, "
        f"and pre-registered analysis (`phases_g_cohort_validation.py`; **Supplementary Note 1**)."
    )
    _cat_tbl = load_catalogue_eqtl_examples(_base, n=8)
    if not _cat_tbl.empty:
        doc.add_paragraph(
            "**Numeric detail | eQTL Catalogue (whole-blood) sentinel hits on ranked shortlist genes.**",
            style="Heading 3",
        )
        add_table_from_dataframe(
            doc,
            _cat_tbl,
            title=None,
            inference="Pagination-limited scan per gene/dataset; see `pipeline_summary.json` limitations.",
            max_rows=12,
        )
        doc.add_paragraph()

    _ph_tbl = load_phases_g_genetics_table(_base, n=14)
    if not _ph_tbl.empty:
        doc.add_paragraph(
            "**Numeric detail | eQTLGen cis sentinel SNPs on shortlist genes with IA GWAS summary-statistics lookup at the same rsID.** "
            "P values are reported exactly as stored (very small eQTLGen P values reflect discrete cis-eQTL catalogue rounding).",
            style="Heading 3",
        )
        add_table_from_dataframe(
            doc,
            _ph_tbl,
            title=None,
            inference="rsID intersection only—not COLOC/SMR. Full panel: `eqtl_snp_candidate_panel.tsv` / `shortlist_evidence.tsv`.",
            max_rows=18,
        )
        doc.add_paragraph()

    doc.add_page_break()

    # --- Discussion ---
    doc.add_heading("Discussion", 1)
    doc.add_paragraph(
        f"We presented GSE122897 as one continuous story rather than disconnected modules. The initial bulk analysis "
        f"established a large rupture-associated programme ({_m['deseq2_padj_genes']:,} DESeq2 FDR hits), motivating an "
        f"immune-aware second pass: archived **CIBERSORTx Job 7** mixtures plus LM22+NNLS fractions fed a linear model in which strict FDR on rupture coefficients "
        f"vanished, yet thousands of genes retained partial evidence—culminating in a capped {_m['n_shortlist']}-gene list for "
        f"genetics. Pathway enrichment and a three-signature ridge model compress that biology into interpretable scores with "
        f"transparent performance reporting; the key out-of-sample estimate is the leakage-controlled nested AUROC "
        f"({_m['auc_nested_mean']:.3f}±{_m['auc_nested_sd']:.3f}), whereas the non-nested CV "
        f"({_m['auc_cv_mean']:.3f}±{_m['auc_cv_sd']:.3f}) and full-cohort fit ({_m['auc_full_cohort']:.3f}) are optimistic by comparison. "
        f"None of these numbers implies bedside risk prediction for unruptured patients.\n\n"
        f"External support is strongest at the level of a smaller reproducible gene set rather than a portable classifier. "
        f"GSE54083 and GSE13353 still provide limited direct rupture support at the gene level, while GSE26969, GSE36791, and "
        f"GSE159610 remain orthogonal context rather than rupture validation. When we projected the fixed three-score model across "
        f"cohorts, performance was markedly asymmetric: GSE13353 showed supportive direct transfer (AUROC {float(_t133.get('AUC', float('nan'))):.3f}), "
        f"whereas GSE54083 showed inverse transfer (AUROC {float(_t540.get('AUC', float('nan'))):.3f}), despite gene-level directional overlap. "
        f"This instability argues against interpreting the current tissue-derived classifier as externally validated. By contrast, the stricter "
        f"three-cohort gene-level meta-analysis yielded {int(_meta_summary.get('n_directionally_concordant_all_three', 0))} directionally concordant "
        f"genes at meta-analysis FDR<0.05, led by TGFBI, ANPEP, MPP1, MARCO, CLEC5A, and CPXM2, suggesting that a small rupture-associated "
        f"molecular core is more reproducible than cohort-level risk discrimination. **PHASES–G** then maps the "
        f"shortlist onto blood eQTLs ({_m['n_blood_eqtl_catalogue']} / {_m['n_eqtlgen_cis']} genes in the Catalogue / eQTLGen "
        f"layers we queried) and IA GWAS rsID overlap ({_m['n_gwas_overlap_catalogue']} / {_m['n_gwas_overlap_eqtlgen']} SNPs), "
        f"preserving a clean separation between tissue discovery and germline inference. EXPRS and related cohort demos should "
        f"be viewed as scaffolding until genotype–PHASES datasets exist.\n\n"
        f"Limitations include small n, single-centre bulk RNA, incomplete clinical metadata, confounding by surgical timing, "
        f"and the optimism inherent in outcome-informed feature selection for EXPRS; nested CV (**Supplementary Figure S11**) "
        f"and explicit caveats in `pipeline_summary.json` partially temper over-interpretation. Future work should pair "
        f"multi-centre tissue and single-cell spatial maps with pre-registered blood biomarker studies and PHASES+PRS validation "
        f"cohorts with DNA banking."
    )

    doc.add_heading("Methods", 1)
    doc.add_paragraph(
        "Below we summarise all analysis steps in the order used in the Results. Command-line invocations, exact file "
        "paths, software versions, and extended parameter lists are documented in **Supplementary Methods** (embedded "
        "in this document after the main text) and in **Supplementary Note 1** for the PHASES–genetics (PHASES–G) layer."
    )

    doc.add_heading("Study design, data acquisition, and cohort definition", 2)
    doc.add_paragraph(
        "We performed a secondary analysis of publicly deposited bulk RNA-seq from resected intracranial aneurysm wall "
        f"(GSE122897; NCBI GEO). The submitters provided TPM-normalised expression (GRCh38.p13) and raw counts for the same "
        f"samples. We restricted to specimens with definitive **ruptured** or **unruptured** wall status, excluding "
        f"non-informative labels in the series metadata, yielding n={_m['n_samples']} samples ({_m['n_ruptured']} ruptured, "
        f"{_m['n_unruptured']} unruptured). Age, sex, and detailed procedural covariates were not available consistently in "
        f"GEO-derived metadata and were not included as regression covariates (**Supplementary Table S6**). The study is "
        f"descriptive of molecular differences at resection and does not constitute a prospective clinical prediction study."
    )

    doc.add_heading("Expression preprocessing, filtering, and transformations", 2)
    doc.add_paragraph(
        "We used the submitter-provided TPM matrix without re-alignment or re-quantification. For univariate screens and "
        "signature construction, values were transformed as log₂(TPM+1). Unless stated otherwise, genes were retained when "
        "TPM > 1 in at least three samples (same threshold as the primary project QC). DESeq2 used gene-level raw counts "
        "with default filtering via DESeq2’s independent filtering where applicable. For the ridge logistic classifier, "
        "the three signature scores (below) were standardised to zero mean and unit variance across the n=43 samples "
        "(scikit-learn `StandardScaler`). No cross-study batch correction was applied within GSE122897."
    )

    doc.add_heading("Differential expression: Welch t-tests and DESeq2", 2)
    doc.add_paragraph(
        f"We computed two parallel differential-expression views. **(1) Welch t-tests** on log₂(TPM+1), comparing ruptured "
        f"vs unruptured for each gene; log₂ fold change was the difference in group means. Two-sided P-values were adjusted "
        f"with Benjamini–Hochberg FDR; genes with FDR < 0.05 were counted as significant in this TPM-based screen (four genes "
        f"at transcriptome-wide FDR in this cohort; full distribution in **Figure 2** and **Supplementary Figures**). "
        f"**(2) DESeq2** was run on raw counts as a count-based confirmatory analysis (`run_deseq2_gse122897.py`); we report "
        f"{_m['deseq2_padj_genes']:,} genes with adjusted P (padj) < 0.05 (**Supplementary Table S3**; **Supplementary Figure S10**). "
        f"Ranked lists for pathway analysis used DESeq2 statistics as described below. HGNC-style symbols for reporting were "
        f"resolved with `results/gene_id_mapping_unified.tsv` where necessary."
    )

    doc.add_heading("Immune deconvolution: CIBERSORTx (Job 7) and local LM22 with non-negative least squares", 2)
    doc.add_paragraph(
        f"To estimate leukocyte admixture in bulk wall RNA we used the LM22 immune reference framework. **CIBERSORTx** "
        f"analyses were executed through the CIBERSORTx portal; the Job 7 mixture output is archived as "
        f"`results/cibersortx/CIBERSORTx_Job7_Results.txt` (n={_m['n_cibersortx_job7_mixtures']} mixture rows in the archived "
        f"export), with per-sample **P-value**, Pearson **correlation**, and **RMSE** quality fields as returned by the portal. "
        f"Pat-style mixture identifiers were rewritten to GEO sample accessions with `rewrite_cibersortx_mixture_to_gsm.py` "
        f"using `cibersortx_pat_to_gsm_map.tsv`, producing `CIBERSORTx_Job7_Results_GSM_mixture.txt` and merged tables such as "
        f"`cibersortx_merged_with_group.tsv` for inspection alongside phenotype labels.\n\n"
        f"For analyses that require proportions aligned **one-to-one** to every RNA profile in GSE122897, we computed "
        f"`results/cibersortx/lm22_fractions_gse122897.tsv` with `immune_deconvolution_lm22_gse122897.py`: LM22 signature genes "
        f"were intersected with the TPM matrix, both mixture and signatures were transformed with log₂(x+1), and per-sample "
        f"proportions were obtained with **non-negative least squares** (scipy.optimize.nnls) and renormalised to sum to one. "
        f"This implements the same biological prior (LM22) in a fully local, reproducible form and is the fraction table "
        f"referenced by the immune-adjusted regression (`deg_immune_adjustment_summary.json` reports {_m['n_genes_lm22_intersect']} "
        f"intersecting LM22 genes). Group contrasts on NNLS proportions used Welch tests with Benjamini–Hochberg FDR across "
        f"the {_m['lm22_cell_types']} compartments (`lm22_ruptured_vs_unruptured.tsv`; summary JSON: "
        f"`lm22_deconvolution_summary.json`)."
    )

    doc.add_heading("Immune-fraction–adjusted, gene-level linear models", 2)
    doc.add_paragraph(
        f"We fitted ordinary least squares (OLS) models gene-wise on log₂(TPM+1) using statsmodels. The design included a "
        f"binary **ruptured** indicator and {_m['n_immune_covariates']} LM22 proportion covariates selected for the adjustment "
        f"set (exact column names listed in `deg_immune_adjustment_summary.json`). Models were fit for genes passing expression "
        f"filters (`n_genes_tested` = {_m['n_genes_immune_model_tested']:,} in the archived summary). For each gene we recorded "
        f"the DESeq2 padj from the primary analysis, the immune-adjusted coefficient for rupture, and two-sided P-values. "
        f"Benjamini–Hochberg FDR for the rupture coefficient was computed **across the tested gene set** in that run "
        f"(not across the full transcriptome). Outputs were written to `deg_immune_covariate_adjusted.tsv`. Mechanistic "
        f"shortlisting used joint criteria (DESeq2 padj ≤ {_m['shortlist_deseq2_padj_max']} and immune-adjusted rupture "
        f"P ≤ {_m['shortlist_lm_p_max']}; capped at {_m['n_shortlist']} genes unless otherwise noted in the build script) as "
        f"recorded in `pipeline_summary.json`."
    )

    doc.add_heading("Pathway and functional enrichment", 2)
    doc.add_paragraph(
        f"We queried **g:Profiler** with a DESeq2-derived ranked gene list (ordered statistic; ordered_query=TRUE) for "
        f"Gene Ontology, KEGG, and Reactome sources with g:SCS multiple-testing control. We report {_m['pathways_fdr05']} terms "
        f"at g:Profiler FDR < 0.05 (**Supplementary Table S5**; **Supplementary Figure S6**). This step is independent of the "
        f"three-score classifier and summarises genome-wide programme-level themes."
    )

    doc.add_heading("Predefined tissue signature scores (DE, EC, ECM)", 2)
    doc.add_paragraph(
        "Composite scores were defined as the mean of log₂(TPM+1) across predefined gene sets after symbol matching to the "
        "expression matrix. **DE score:** mean over top differential genes ranked by DESeq2 evidence (see primary pipeline "
        "and **Supplementary Methods** for the exact selection rule). **EC (endothelial dysfunction) score:** mean over "
        "ANGPT2, ICAM1, VCAM1, THBS1, PERP. **ECM (extracellular matrix remodelling) score:** mean over MMP9, MMP2, TIMP1, "
        "COL3A1, COL1A1. An **inflammation** score (TNF, IL6, IL8/CXCL8, CRP, SAA1) was computed for exploration but excluded "
        "from the primary three-feature ridge model. Unmapped symbols were handled as documented in the supplementary gene "
        "lists."
    )

    doc.add_heading("Supervised classification, cross-validation, and variable importance", 2)
    doc.add_paragraph(
        f"We trained a **binomial logistic regression with L2 (ridge) penalty** on rupture status using only the three "
        f"standardised scores (DE, EC, ECM) with fixed regularisation strength (**C=1.0** in the current scripts). "
        f"We report full-cohort fit, optimistic non-nested five-fold stratified cross-validation, and a leakage-controlled "
        f"nested design in which DE genes are re-selected within each outer fold. The final coefficients were refit on all {_m['n_samples']} "
        f"samples for reporting. **Permutation importance** "
        f"was estimated by shuffling each score, re-evaluating AUC, and averaging the decrease over repeated permutations "
        f"(values underlying **Figure 4a**). The primary ROC curve (**Figure 4b**) reflects full-cohort predictions at the "
        f"optimised decision rule (Youden threshold; see **Supplementary Methods**). Repeated CV summaries are in "
        f"**Supplementary Table S2**; **Supplementary Figure S11** contrasts non-nested CV with a **nested** design in which "
        f"DE genes for the DE score are re-selected within each outer fold to reduce selection bias. Risk stratification by "
        f"predicted probability (low / medium / high) is described in **Supplementary Methods** and **Supplementary Figure S12**."
    )

    doc.add_heading("External transcriptomic context, transfer testing, direct-rupture consensus, and PHASES–G", 2)
    doc.add_paragraph(
        f"The PHASES–G layer operates on the immune-aware capped shortlist (n={_m['n_shortlist']} in the current build). "
        f"**External transcriptomic hierarchy:** we first prioritised direct rupture support in GSE54083 and GSE13353, then "
        f"interpreted GSE26969, GSE36791, and GSE159610 as orthogonal tissue or blood context only. **Sample-level transfer:** "
        f"the fixed three-score classifier was projected unchanged into each public cohort and quantified with AUROC, bootstrap "
        f"confidence intervals, and label-permutation P values; portability was supportive in GSE13353 (AUROC {float(_t133.get('AUC', float('nan'))):.3f}) "
        f"but unstable in GSE54083 (AUROC {float(_t540.get('AUC', float('nan'))):.3f}). **Direct-rupture consensus:** genes measurable "
        f"across GSE122897, GSE54083, and GSE13353 were combined with a weighted Stouffer meta-analysis, yielding "
        f"{int(_meta_summary.get('n_directionally_concordant_all_three', 0))} directionally concordant genes at meta-analysis FDR<0.05. "
        f"**Orthogonal bulk check:** for each shortlist gene with a resolvable symbol, we compared the sign of the GSE122897 "
        f"rupture DESeq2 log₂FC to the sign of a published **GSE159610** aneurysm-vs-control contrast when available "
        f"({_m['n_gse159610_matched']} genes with non-missing GSE159610 estimates in `pipeline_summary.json`). "
        f"**eQTL Catalogue:** programmatic queries (HTTP API) against whole-blood expression datasets were used to attach a "
        f"sentinel cis-association per gene under project pagination limits (`n_with_blood_eqtl_hit` = {_m['n_blood_eqtl_catalogue']}). "
        f"**eQTLGen:** cis-eQTLs at FDR < 0.05 from the bundled probe-level file were cross-referenced per gene "
        f"(`n_with_eqtlgen_cis_hit` = {_m['n_eqtlgen_cis']}). **IA GWAS overlap:** sentinel rsIDs from each layer were intersected "
        f"with rsIDs present in Bakker et al. stage-2 intracranial aneurysm summary statistics "
        f"(Catalogue layer: {_m['n_gwas_overlap_catalogue']} SNPs; eQTLGen layer: {_m['n_gwas_overlap_eqtlgen']} SNPs). "
        f"This is **rsID intersection only**—not formal colocalisation; selected loci for COLOC/SMR-style follow-up are under "
        f"`results/phases_g/stronger_genetics/` where computed.\n\n"
        f"**EXPRS (expression-based score):** when individual genotypes are unavailable (typical for GEO RNA deposits), we "
        f"computed a weighted z-score of log-expression using eQTL effect sizes as weights for {_m['n_exprs_genes']} genes with "
        f"overlapping TPM support (`exprs_gene_weights.tsv`, `exprs_sample_scores.tsv`). **Germline PRS path:** when cohort "
        f"dosages and PHASES covariates exist on harmonised `sample_id`s, the repository provides LD clumping (PLINK 2), "
        f"construction of SNP weights (`prs_rsid_weights_from_panel.tsv`; {_m['n_prs_weights_snps']} SNPs in the current panel), "
        f"and validation scaffolds (`phases_g_prs_from_genotypes.py`, `validate_prs_genotype_inputs.py`, "
        f"`phases_g_cohort_validation.py`; see `data/genotypes/README_GENOTYPE_INPUT.txt`)."
    )

    doc.add_heading("Software, reporting, and reproducibility", 2)
    doc.add_paragraph(
        "Analyses were implemented in **Python 3** using pandas, numpy, scipy, statsmodels, and scikit-learn; figures used "
        "matplotlib and seaborn. DESeq2 was run via R/Bioconductor as referenced in the supplementary scripts. CIBERSORTx "
        "results were obtained from the Stanford/NIAID CIBERSORTx portal (Job 7 export archived locally). g:Profiler was accessed "
        "through its public API. Analysis-ready tables and JSON summaries are written under `results/` as cited above; "
        "`verify_phases_g_objectives.py` and `results/phases_g/PHASES_G_REQUIREMENTS_STATUS.md` track pipeline completeness."
    )

    doc.add_heading("Data availability", 1)
    doc.add_paragraph(
        "GSE122897 (NCBI GEO). Processed tables and analysis scripts: project repository / authors upon request."
    )

    doc.add_heading("References", 1)
    doc.add_paragraph(
        "[Numbered references as in prior manuscript; add Ruigrok/GSE122897, Bakker IA GWAS, eQTLGen/eQTL Catalogue as needed.]"
    )
    doc.add_page_break()

    # --- Supplementary (full Methods + figs S1–S12 + embedded tables S1–S6) ---
    doc.add_heading("Supplementary Materials", 1)

    doc.add_heading("Supplementary Note 1. Repository outputs, deconvolution caveats, and genetics prioritisation layer", 2)
    doc.add_paragraph(
        "Key output files including shortlist gene tables, eQTL overlap panels, and PRS scaffolds are deposited in the public code repository "
        "(DOI: [to be assigned at submission]). The same repository should also host the deconvolution summaries, immune-adjusted "
        "regression tables, external transfer outputs, direct-rupture meta-analysis tables, and figure assets used in this manuscript.\n\n"
        "Because the LM22 reference is leukocyte-derived and the observed CIBERSORTx fit was weak, the deconvolution-derived "
        "fractions and downstream immune-adjusted models should be interpreted as sensitivity analyses rather than validated "
        "cell-composition estimates. The downstream genetics-informed prioritisation layer links the tissue shortlist to blood cis-eQTL "
        "resources, rsID-level IA GWAS overlap, expression-weighted exploratory scores, and optional PRS scaffolds for future cohorts "
        "with same-subject genotype and clinical data. In the repository this layer is historically labelled PHASES-G, but the PHASES "
        "score itself was not computed in the present manuscript."
    )
    doc.add_paragraph()

    supp_fig_intro = (
        "Supplementary Figures S1–S11 support the main-text Results and Methods; Supplementary Figure S3 duplicates the ROC "
        "layout for the supplement (primary ROC: main-text Figure 4b). Supplementary Figure S12 shows risk stratification "
        "by model score (paired with Supplementary Table S1) and replaces the fifth main-text figure used in earlier drafts. "
        "Each figure includes a legend and inference paragraph."
    )
    append_supplementary_bundle(
        doc,
        s1_stratification_figure_ref="Supplementary Figure S12",
        supplementary_figures_intro=supp_fig_intro,
        include_supplementary_figure_s12=True,
        figure_width=Inches(5),
    )

    out = script_dir / "Molecular_Signatures_Revised_VS.docx"
    doc.save(str(out))
    print(f"Created: {out}")
    print(
        "Main: Figures 1-7 and Tables 1-5 plus embedded Results numeric panels; "
        "stratification -> Supp Fig S12; PHASES-G -> Supp Note 1. "
        "Full Supplementary Methods + S1-S12 + embedded Tables S1-S6 in same file."
    )
    return str(out)


if __name__ == "__main__":
    create_revised_manuscript()
