"""
The Lancet Digital Health — style draft manuscript (single Word file).

Produces one .docx with: journal-style title page, **Summary** (semistructured:
Background, Methods, Findings, Interpretation, Funding per journal convention),
**Research in context** panel, IMRAD order (Introduction, Methods, Results,
Discussion), declarations, Vancouver-style references, embedded main figures
and tables, and a full **Appendix** (Supplementary Methods, Figures S1–S12,
Tables S1–S6).

This is a **draft for author preparation**, not a substitute for the journal’s
live “instructions for authors” (word limits, structured forms, and file uploads
change over time).

Output: Manuscript_Lancet_Digital_Health.docx (project root).

  python create_manuscript_lancet_digital_health.py

Analytic content matches `create_molecular_signatures_revised_vs.py`; layout and
front matter follow Lancet Digital Health–style reporting.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

from blood_bridge_panel import build_blood_bridge_panel_assets
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


def create_lancet_digital_health_manuscript() -> str:
    script_dir = Path(__file__).resolve().parent
    os.chdir(script_dir)
    _base = script_dir
    _m = load_story_metrics(_base)
    _gse_stats = gse159610_shortlist_stats(_base)
    _combined_validation_df = build_combined_validation_assets(_base)
    _combined_validation_by_dataset = (
        {str(r["Dataset"]): r for r in _combined_validation_df.to_dict(orient="records")}
        if not _combined_validation_df.empty and "Dataset" in _combined_validation_df.columns
        else {}
    )
    _panel_summary = build_blood_bridge_panel_assets(_base)
    _panel_perf_path = _base / "results" / "blood_bridge_panel" / "blood_bridge_panel_performance.tsv"
    _panel_gene_path = _base / "results" / "blood_bridge_panel" / "blood_bridge_panel_gene_summary.tsv"
    _panel_perf_df = pd.read_csv(_panel_perf_path, sep="\t") if _panel_perf_path.is_file() else pd.DataFrame()
    _panel_gene_df = pd.read_csv(_panel_gene_path, sep="\t") if _panel_gene_path.is_file() else pd.DataFrame()
    _panel_perf_by_dataset = (
        {str(r["Dataset"]): r for r in _panel_perf_df.to_dict(orient="records")}
        if not _panel_perf_df.empty and "Dataset" in _panel_perf_df.columns
        else {}
    )
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
    _p106 = _panel_perf_by_dataset.get("GSE106520", {})
    _p159 = _panel_perf_by_dataset.get("GSE159610", {})
    _p367 = _panel_perf_by_dataset.get("GSE36791", {})
    _p733 = _panel_perf_by_dataset.get("GSE73378", {})
    _cv156r = _combined_validation_by_dataset.get("GSE15629 (rupture subset)", {})
    _cv156c = _combined_validation_by_dataset.get("GSE15629 (aneurysm vs control)", {})
    _cv106 = _combined_validation_by_dataset.get("GSE106520", {})
    _cv733 = _combined_validation_by_dataset.get("GSE73378", {})
    _panel_gene_list = _panel_summary.get("panel_genes", []) if isinstance(_panel_summary, dict) else []
    _panel_gene_text = ", ".join(str(x) for x in _panel_gene_list)

    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "Arial"
    style.font.size = Pt(11)

    # --- Journal-style cover (draft) ---
    cov = doc.add_paragraph()
    cov.add_run("Draft prepared in The Lancet Digital Health article style").italic = True
    cov.alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.add_paragraph("Article type: Original Research (secondary analysis)", style="Intense Quote")
    doc.add_paragraph()
    t = doc.add_heading(
        "Secondary analysis of aneurysm wall transcriptomes identifies exploratory rupture-associated programmes and a genetics-informed prioritisation framework",
        0,
    )
    t.alignment = WD_ALIGN_PARAGRAPH.CENTER
    auth = doc.add_paragraph("[Author names and affiliations to be inserted]")
    auth.alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.add_paragraph()
    doc.add_paragraph(
        "Correspondence: [institutional address and email to be inserted at submission]"
    )
    doc.add_page_break()

    # --- Summary (semistructured abstract) ---
    doc.add_heading("Summary", 1)
    doc.add_heading("Background", 2)
    doc.add_paragraph(
        "Public aneurysm-wall transcriptomic datasets offer a way to study rupture-associated biology, but they are small, "
        "surgically ascertained, and unsuitable for direct bedside prediction. We therefore aimed to perform a conservative "
        "secondary analysis focused on mechanistic signal detection, sensitivity to immune admixture, limited external "
        "gene-level support, and genetics-informed prioritisation rather than clinical validation."
    )
    doc.add_heading("Methods", 2)
    doc.add_paragraph(
        f"We analysed GSE122897 bulk RNA-seq from resected aneurysm wall (n={_m['n_samples']}; {_m['n_ruptured']} ruptured, "
        f"{_m['n_unruptured']} unruptured). We used DESeq2 on raw counts and Welch tests on log₂(TPM+1), estimated LM22 "
        f"immune fractions using archived CIBERSORTx output and local non-negative least squares, and fitted immune-adjusted "
        f"ordinary least squares models gene-wise. We summarised expression into three predefined tissue scores and assessed a "
        f"ridge logistic model with full-cohort fit, optimistic non-nested cross-validation, and a leakage-controlled nested "
        f"design. External transcriptomic support was sought in rupture microarrays GSE54083 and GSE13353, the mixed-design tissue study "
        f"GSE15629, aneurysm-versus-control tissue context in GSE26969, aneurysm-versus-control blood bridge cohorts GSE159610 and "
        f"GSE106520, and post-rupture / subarachnoid-haemorrhage blood context in GSE36791 and GSE73378. We additionally projected the "
        f"three-score model sample-wise into each public cohort, derived a reduced blood-testable panel by ranking genes that remained "
        f"directionally concordant in both GSE159610 and GSE106520, and meta-analysed gene-level rupture signals across the three direct "
        f"rupture datasets. We then attached blood cis-eQTL and intracranial aneurysm GWAS rsID support to shortlisted genes as a "
        f"genetics-informed prioritisation framework."
    )
    doc.add_heading("Findings", 2)
    doc.add_paragraph(
        f"DESeq2 identified {_m['deseq2_padj_genes']:,} genes with false-discovery rate (FDR) <0·05. After linear adjustment for "
        f"{_m['n_immune_covariates']} LM22-derived covariates across {_m['n_genes_immune_model_tested']:,} tested genes, no rupture "
        f"coefficient reached Benjamini–Hochberg FDR<0·05 within that tested set; {_m['n_genes_deseq2_and_nominal_rupture_immune']:,} "
        f"genes met joint DESeq2 significance and nominal rupture association and fed a capped shortlist of {_m['n_shortlist']} genes. "
        f"g:Profiler returned {_m['pathways_fdr05']} terms at FDR<0·05. The three-score ridge model achieved AUROC "
        f"{_m['auc_nested_mean']:.3f} (SD {_m['auc_nested_sd']:.3f}) in the leakage-controlled nested analysis, compared with "
        f"{_m['auc_cv_mean']:.3f} (SD {_m['auc_cv_sd']:.3f}) in non-nested cross-validation and {_m['auc_full_cohort']:.3f} on the "
        f"full cohort. In independent rupture microarray GSE54083 ({_m['gse54083_n_ruptured']} ruptured, {_m['gse54083_n_unruptured']} "
        f"unruptured), {_m['gse54083_n_shortlist_present']} shortlist genes were measurable, {_m['gse54083_n_shortlist_direction_match']} "
        f"matched GSE122897 direction, and {_m['gse54083_n_shortlist_p_lt_0_05']} were nominally significant; none survived FDR. "
        f"A targeted seven-gene check in GSE13353 showed {_m['gse13353_direction_matches']}/{_m['gse13353_genes_tested']} directional matches "
        f"({_m['gse13353_nominal_sig_genes']} nominally significant). In GSE15629, a secondary rupture subset yielded 138/251 directionally "
        f"matched shortlist genes with 10 nominal signals, whereas the aneurysm-versus-control contrast yielded 167/251 directionally matched "
        f"genes with 54 nominal signals and therefore serves as tissue context rather than rupture replication. The two main blood bridge cohorts "
        f"were GSE159610 (143/263 directionally matched shortlist genes) and GSE106520 (123/262). A reduced eight-gene blood-testable panel "
        f"({ _panel_gene_text }) achieved AUROC {float(_p159.get('AUC', float('nan'))):.3f} in GSE159610 and {float(_p106.get('AUC', float('nan'))):.3f} in GSE106520. "
        f"In blood microarray GSE36791 ({_m['gse36791_n_ruptured']} ruptured, {_m['gse36791_n_control']} control), {_m['gse36791_n_shortlist_direction_match']}/"
        f"{_m['gse36791_n_shortlist_present']} measurable shortlist genes matched direction; {_m['gse36791_n_shortlist_p_lt_0_05']} were nominally "
        f"significant and {_m['gse36791_n_shortlist_fdr_lt_0_05']} survived FDR, although this dataset mainly reflects post-rupture blood response; "
        f"the reduced panel scored AUROC {float(_p367.get('AUC', float('nan'))):.3f} in GSE36791 but only {float(_p733.get('AUC', float('nan'))):.3f} in GSE73378, "
        f"supporting its interpretation as bridge biology rather than long-horizon post-SAH signal. In sample-level transfer of the original three-score "
        f"tissue model, AUROC was {float(_t133.get('AUC', float('nan'))):.3f} in GSE13353 but only {float(_t540.get('AUC', float('nan'))):.3f} in GSE54083, "
        f"indicating unstable portability across direct rupture cohorts. "
        f"Weighted meta-analysis across GSE122897, GSE54083, and GSE13353 identified "
        f"{int(_meta_summary.get('n_directionally_concordant_all_three', 0))} directionally concordant genes with meta-analysis FDR<0·05. "
        f"The genetics-informed prioritisation layer added blood eQTL and IA GWAS rsID overlap as metadata, not as external validation."
    )
    doc.add_heading("Interpretation", 2)
    doc.add_paragraph(
        "This study supports a broad exploratory rupture-associated wall transcriptomic programme and identifies candidates that "
        "remain of interest after immune-aware filtering and limited external checking. The work does not validate a clinical "
        "prediction model, does not establish cell-intrinsic biology, and does not demonstrate added value over PHASES. "
        "Cross-cohort sample-level transfer remained unstable for the original tissue-derived classifier, but a smaller blood-accessible panel transported more consistently across the two bridge cohorts. "
        "That reduced panel remains exploratory because it was selected from outcome-informed public data and was only benchmarked in aneurysm-presence and post-rupture blood settings. "
        "The genetics-informed prioritisation layer should be interpreted as a transparent exploratory framework whose clinical assessment requires independent "
        "cohorts with accessible biomarkers, harmonised clinical covariates, and same-subject genotypes."
    )
    doc.add_heading("Funding", 2)
    doc.add_paragraph(
        "This secondary analysis used publicly deposited data. Funding, sponsor roles, and grant numbers are to be declared at "
        "submission per journal policy."
    )
    doc.add_page_break()

    # --- Research in context ---
    doc.add_heading("Research in context", 1)
    doc.add_heading("Evidence before this study", 2)
    doc.add_paragraph(
        "Small public surgical cohorts have suggested that ruptured and unruptured aneurysm walls differ substantially in bulk "
        "expression, but the same datasets also carry major limitations: small sample size, incomplete clinical metadata, immune "
        "admixture, and limited replication. Public resources additionally provide rupture microarrays, aneurysm-versus-control tissue "
        "and blood expression datasets, blood rupture/SAH transcriptomes, blood cis-eQTL catalogues, and intracranial aneurysm GWAS "
        "summary statistics, but these layers are often combined in ways that blur discovery, replication, and genetics-based hypothesis generation."
    )
    doc.add_heading("Added value of this study", 2)
    doc.add_paragraph(
        "We rewrote the analysis around a stricter evidentiary hierarchy. The study now foregrounds DESeq2-based rupture biology, "
        "shows that immune-adjusted rupture coefficients do not survive FDR, makes the leakage-controlled nested performance estimate "
        "the primary modelling result, uses GSE54083 as the main independent rupture replication dataset, keeps GSE13353 as targeted "
        "support, adds GSE15629 as a mixed-design tissue cohort, separates GSE159610 and GSE106520 as blood bridge cohorts from GSE36791 "
        "and GSE73378 as post-rupture / SAH context only, and adds both sample-level external signature transfer and a three-cohort "
        "rupture meta-analysis. Those added analyses show mixed external portability of the original tissue classifier but motivate a reduced "
        "eight-gene blood-testable panel that transports better across the two bridge cohorts. Genetics-linked outputs are retained only as candidate-gene prioritisation layers."
    )
    doc.add_heading("Implications of all the available evidence", 2)
    doc.add_paragraph(
        "Bulk surgical transcriptomics can still be useful when framed as exploratory biology rather than digital-health validation. "
        "Analyses in this area should report leakage-controlled performance, distinguish rupture replication from aneurysm-control context, "
        "separate blood bridge cohorts from post-rupture systemic context, test whether score transfer is stable across cohorts, and avoid presenting rsID overlap or expression-derived scores as clinical validation. "
        "Progress toward clinically relevant risk modelling will require accessible sample types, full PHASES variables, germline data, and independent validation cohorts."
    )
    doc.add_page_break()

    # --- Introduction ---
    doc.add_heading("Introduction", 1)
    doc.add_paragraph(
        "Ruptured intracranial aneurysms remain clinically devastating, but the molecular evidence most readily available in "
        "public repositories comes from resected wall tissue obtained at surgery. Such tissue provides a direct view of biology "
        "at intervention, yet it is not available in the larger population of patients with unruptured lesions and therefore "
        "cannot by itself support bedside risk stratification. The most defensible use of these datasets is mechanistic: to "
        "identify rupture-associated programmes, test how strongly those signals depend on inferred immune composition, and ask "
        "whether prioritised genes show any support in independent expression resources.\n\n"
        "Here we reframe GSE122897 accordingly. We treat it as a secondary exploratory tissue study, not as a clinical prediction "
        "cohort. We begin with differential expression, then apply immune-aware retesting, collapse the signal into interpretable "
        "pathway and tissue scores, report conservative as well as optimistic model-performance estimates, seek limited external "
        "rupture support in GSE54083, GSE13353, and a secondary GSE15629 rupture subset, add GSE26969 and the aneurysm-versus-control "
        "arm of GSE15629 as orthogonal tissue context, use GSE159610 and GSE106520 as blood bridge cohorts, keep GSE36791 and GSE73378 "
        "as post-rupture / SAH context only, test cross-cohort signature transfer and reduced blood-panel transportability, and finally "
        "attach blood eQTL and GWAS metadata through a genetics-informed prioritisation layer rather than a validated risk framework."
    )
    doc.add_page_break()

    # --- Methods (IMRAD: before Results) ---
    doc.add_heading("Methods", 1)
    doc.add_paragraph(
        "We report analyses in the order used for the Results. Additional implementation details, software versions, "
        "and extended parameters are in the **Appendix (Supplementary Methods)**. Repository and reproducibility notes are in "
        "**Supplementary Note 1** within the Appendix."
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
        f"**(2) DESeq2** was run on raw counts as a count-based confirmatory analysis; we report "
        f"{_m['deseq2_padj_genes']:,} genes with adjusted P (padj) < 0.05 (**Supplementary Table S3**; **Supplementary Figure S10**). "
        f"Ranked lists for pathway analysis used DESeq2 statistics as described below. HGNC-style symbols for reporting were "
        f"resolved through harmonised annotation tables, with unresolved legacy annotations left as IDs rather than assigned incorrect symbols."
    )

    doc.add_heading("Immune deconvolution: CIBERSORTx (Job 7) and local LM22 with non-negative least squares", 2)
    doc.add_paragraph(
        f"To estimate leukocyte admixture in bulk wall RNA we used the LM22 immune reference framework. **CIBERSORTx** "
        f"analyses were executed through the CIBERSORTx portal; the archived export contained "
        f"{_m['n_cibersortx_job7_mixtures']} mixture rows with per-sample **P-value**, Pearson **correlation**, and **RMSE** "
        f"quality fields. The median correlation in this workspace was only {_m['cibersortx_median_correlation']:.3f}, with "
        f"median RMSE {_m['cibersortx_median_rmse']:.3f}, indicating weak fit of the LM22 blood-derived reference to aneurysm-wall "
        f"tissue. We therefore interpret this layer as a coarse immune-sensitivity analysis rather than validated cell-type "
        f"quantification.\n\n"
        f"For analyses that require proportions aligned **one-to-one** to every RNA profile in GSE122897, we computed "
        f"local LM22+NNLS fractions from the TPM matrix: LM22 signature genes "
        f"were intersected with the TPM matrix, both mixture and signatures were transformed with log₂(x+1), and per-sample "
        f"proportions were obtained with **non-negative least squares** (scipy.optimize.nnls) and renormalised to sum to one. "
        f"This local implementation uses the same LM22 biological prior in reproducible form. Because LM22 is a leukocyte "
        f"reference rather than a vessel-wall reference and only {_m['n_genes_lm22_intersect']} signature genes intersected the "
        f"tissue matrix, the resulting fractions should likewise be treated cautiously. Group contrasts on NNLS proportions used "
        f"Welch tests with Benjamini–Hochberg FDR across the {_m['lm22_cell_types']} compartments."
    )

    doc.add_heading("Immune-fraction–adjusted, gene-level linear models", 2)
    doc.add_paragraph(
        f"We fitted ordinary least squares (OLS) models gene-wise on log₂(TPM+1) using statsmodels. The design included a "
        f"binary **ruptured** indicator and {_m['n_immune_covariates']} LM22 proportion covariates selected for the adjustment "
        f"set. Models were fit for genes passing expression filters (n tested = {_m['n_genes_immune_model_tested']:,}). For each gene we recorded "
        f"the DESeq2 padj from the primary analysis, the immune-adjusted coefficient for rupture, and two-sided P-values. "
        f"Benjamini–Hochberg FDR for the rupture coefficient was computed **across the tested gene set** in that run "
        f"(not across the full transcriptome). Mechanistic "
        f"shortlisting used joint criteria (DESeq2 padj ≤ {_m['shortlist_deseq2_padj_max']} and immune-adjusted rupture "
        f"P ≤ {_m['shortlist_lm_p_max']}; capped at {_m['n_shortlist']} genes unless otherwise noted in the repository build summary). "
        f"Given the weak LM22/CIBERSORTx fit, this adjustment should be interpreted as a robustness filter against coarse inferred "
        f"immune content, not as proof of cell-intrinsic biology."
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
        f"standardised scores (DE, EC, ECM) with fixed regularisation strength (**C=1.0** in the current scripts). We report "
        f"full-cohort fit, optimistic non-nested five-fold stratified cross-validation, and a leakage-controlled nested design "
        f"in which DE genes for the DE score are re-selected within each outer fold. The final coefficients were refit on all "
        f"{_m['n_samples']} samples for reporting. **Permutation importance** "
        f"was estimated by shuffling each score, re-evaluating AUC, and averaging the decrease over repeated permutations "
        f"(values underlying **Figure 4a**). The primary ROC curve (**Figure 4b**) is drawn from held-out predictions in the "
        f"leakage-controlled nested design. The full-cohort fitted ROC is retained only in the supplement. **Supplementary Table S2** and "
        f"**Supplementary Figure S11** show why the nested estimate is the more defensible out-of-sample summary. Risk stratification by "
        f"predicted probability (low / medium / high) is described in **Supplementary Methods** and **Supplementary Figure S12**."
    )

    doc.add_heading("External transcriptomic context and genetics-informed prioritisation", 2)
    doc.add_paragraph(
        f"The genetics-informed prioritisation layer operates on the immune-aware capped shortlist (n={_m['n_shortlist']} in the current build) and is "
        f"used here for prioritisation only. The repository historically labels this layer **PHASES-G**, but the clinical PHASES score "
        f"itself was **not computed** in this manuscript and was not used as a comparator. External transcriptomic support was sought first in rupture-specific datasets: "
        f"GSE54083 was processed gene-wise from raw Agilent files, GSE13353 contributes a targeted seven-gene rupture check, and "
        f"GSE15629 contributes both a secondary rupture subset and an aneurysm-versus-control tissue contrast. We then added GSE26969 as "
        f"aneurysm-wall versus superficial temporal artery context, designated GSE159610 and GSE106520 as the main blood bridge cohorts "
        f"for reduced-panel derivation, and kept GSE36791 and GSE73378 as post-rupture / SAH blood context only because they do not reflect "
        f"prospective aneurysm-wall rupture biology.\n\n"
        f"Blood cis-eQTL metadata were attached through the eQTL Catalogue ({_m['n_blood_eqtl_catalogue']} genes with a returned blood cis-eQTL hit) "
        f"and eQTLGen ({_m['n_eqtlgen_cis']} genes with cis evidence). Sentinel rsIDs from these layers were intersected with "
        f"Bakker et al. intracranial aneurysm summary statistics (Catalogue layer: {_m['n_gwas_overlap_catalogue']} SNPs; eQTLGen "
        f"layer: {_m['n_gwas_overlap_eqtlgen']} SNPs). This is **rsID intersection only**, not formal colocalisation. When "
        f"individual-level genotypes are unavailable, the repository can compute EXPRS for {_m['n_exprs_genes']} genes with "
        f"overlapping TPM support; when genotype and PHASES covariates exist on the same subjects, the repository also contains "
        f"optional PRS scaffolds using {_m['n_prs_weights_snps']} SNP weights. None of those downstream layers constitutes clinical validation "
        f"in the present public RNA-only setting."
    )

    doc.add_heading("External sample-level signature transfer, blood-bridge panel derivation, and direct-rupture meta-analysis", 2)
    doc.add_paragraph(
        "To test whether the three-score tissue classifier was portable beyond the discovery cohort, we recalculated the same "
        "DE, EC, and ECM scores in each external expression dataset using cohort-specific gene availability, applied the saved "
        "GSE122897 ridge-model scaler and coefficients without refitting, and quantified discrimination with AUROC, area under "
        "the precision-recall curve, Mann-Whitney tests, non-parametric bootstrap confidence intervals, and label-permutation "
        "P values. These projections were interpreted according to cohort role: direct rupture transfer in GSE54083 and GSE13353, "
        "mixed-design secondary rupture support in GSE15629, orthogonal tissue context in GSE26969 and the aneurysm-versus-control arm of GSE15629, "
        "and blood context in GSE159610, GSE106520, GSE36791, and GSE73378.\n\n"
        "Because those projections were unstable across the direct rupture cohorts, we derived a reduced blood-testable panel from the "
        "immune-aware shortlist by retaining only genes that remained directionally concordant in both GSE159610 and GSE106520, ranking them "
        "by combined unweighted Stouffer Z across the two bridge-cohort P values, and carrying forward the top eight genes. For each cohort "
        "we then z-scored each panel gene within cohort, multiplied by the discovery-cohort log2 fold-change sign, averaged across available "
        "genes, and evaluated the resulting panel score with AUROC and Mann-Whitney testing. GSE36791 and GSE73378 were deliberately treated "
        "as context-only readouts rather than as discovery or optimisation cohorts.\n\n"
        "Separately, we constructed a stricter cross-cohort gene-level synthesis restricted to the three direct rupture datasets "
        "(GSE122897, GSE54083, GSE13353). For genes measurable in all three cohorts, we converted two-sided P values to signed Z "
        "statistics using the direction of log2 fold change and combined them with a weighted Stouffer method using effective "
        "sample-size weights sqrt(nR*nU/(nR+nU)). Benjamini-Hochberg FDR was applied to the meta-analysis P values, and genes "
        "with concordant direction across all three cohorts were highlighted as the most reproducible rupture-associated core."
    )

    doc.add_heading("Software, reporting, and reproducibility", 2)
    doc.add_paragraph(
        "Analyses were implemented in **Python 3** using pandas, numpy, scipy, statsmodels, and scikit-learn; figures used "
        "matplotlib and seaborn. DESeq2 was run via R/Bioconductor as referenced in the supplementary scripts. CIBERSORTx "
        "results were obtained from the Stanford/NIAID CIBERSORTx portal. g:Profiler was accessed through its public API. "
        "Analysis-ready tables, figure assets, and code used for this draft should be deposited in a public repository and cited "
        "by DOI or permanent URL at submission."
    )
    doc.add_page_break()

    # --- Results ---
    doc.add_heading("Results", 1)
    doc.add_paragraph(
        "We organised results to follow the analysis pipeline: bulk differential expression in GSE122897; "
        "CIBERSORTx-based immune deconvolution and immune-fraction-adjusted re-testing; pathway- and signature-level "
        "summaries with multivariate tissue discrimination; independent rupture microarray support plus orthogonal aneurysm-control "
        "context; bridge-versus-context external validation with a reduced blood-panel transportability analysis; direct-rupture consensus analysis; and a genetics-informed prioritisation layer "
        "linking the immune-aware shortlist to blood eQTLs and GWAS metadata."
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
            inference="Full DESeq2 output appears in **Supplementary Table S3**; unresolved legacy annotations were left as IDs rather than assigned speculative symbols.",
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
        f"for **Job 7** contained {_m['n_cibersortx_job7_mixtures']} samples with per-sample CIBERSORTx P-value, "
        f"Pearson correlation, and RMSE. Parsed QC medians in this workspace were correlation {_m['cibersortx_median_correlation']:.3f} "
        f"and RMSE {_m['cibersortx_median_rmse']:.3f} across {_m['n_cibersortx_gsm_parsed']} mapped samples. Correlations this low "
        f"imply that LM22 is poorly matched to aneurysm-wall tissue, which includes vascular wall cell types not represented in a "
        f"blood-derived leukocyte reference. We therefore treat the immune-deconvolution layer as a coarse sensitivity analysis "
        f"rather than as high-confidence cellular quantification. "
        f"Separately, to feed the genome-wide OLS adjustment with proportions aligned **one-to-one** to all {_m['n_samples']} "
        f"RNA profiles, we used the same {_m['lm22_cell_types']}-column LM22 basis deconvolved by constrained non-negative least "
        f"squares on log-transformed TPM ({_m['n_genes_lm22_intersect']} LM22 genes detected in the mixture). Group contrasts on "
        f"those NNLS fractions showed {_m['lm22_fdr_sig_cell_types_005']} cell types at FDR<0.05, consistent with subtle or unstable "
        f"composition shifts rather than a robust deconvolution signal. We then regressed log₂(TPM+1) on rupture status "
        f"while including {_m['n_immune_covariates']} selected LM22 fraction covariates (**Supplementary Methods**). Among {_m['n_genes_immune_model_tested']:,} adequately expressed genes tested in this "
        f"framework, {_m['n_deseq2_fdr_among_immune_tested']:,} retained DESeq2 FDR<0.05 from the primary analysis, "
        f"{_m['n_genes_immune_nominal_rupture']:,} showed nominal (P<0.05) rupture coefficients after immune adjustment, and "
        f"critically **zero** genes passed Benjamini–Hochberg FDR<0.05 for the rupture coefficient when FDR was computed across "
        f"this tested gene set. For mechanistic prioritisation we nevertheless retained genes that are both DESeq2-significant "
        f"in tissue **and** nominally rupture-associated after immune adjustment (n={_m['n_genes_deseq2_and_nominal_rupture_immune']}), "
        f"then applied directional replication filters and a hard cap to arrive at the {_m['n_shortlist']}-gene shortlist used "
        f"in the downstream genetics-informed prioritisation panels. The low CIBERSORTx correlations should temper any strong "
        f"claim that this shortlist is fully immune-independent."
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
            inference="These examples illustrate the wide range of fit quality across samples; the generally low correlations support treating LM22-based deconvolution as a coarse sensitivity analysis only.",
            max_rows=12,
        )
        doc.add_paragraph()

    _imm_ex = load_immune_adjusted_exemplars(_base, n=10)
    if not _imm_ex.empty:
        doc.add_paragraph(
            "**Numeric detail | Immune-fraction–adjusted model: genes with strongest nominal rupture coefficients among DESeq2 hits.** "
            "OLS on log₂(TPM+1) with selected LM22 covariates.",
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
        f"receiver-operating-characteristic curve from the leakage-controlled nested analysis (mean fold AUROC {_m['auc_nested_mean']:.3f}; "
        f"n={_m['n_samples']}). Non-nested stratified cross-validation yielded mean AUROC {_m['auc_cv_mean']:.3f}±{_m['auc_cv_sd']:.3f}, "
        f"whereas the leakage-controlled nested estimate was {_m['auc_nested_mean']:.3f}±{_m['auc_nested_sd']:.3f} "
        f"(**Supplementary Table S2**; **Supplementary Figure S11**). We interpret the nested estimate as the most defensible "
        f"out-of-sample summary because DE-derived features are re-selected within outer folds. The stronger full-cohort ROC "
        f"({_m['auc_full_cohort']:.3f}) is retained only as a supplementary descriptive fit. Confusion matrices, precision–recall curves, and alternate ROC "
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
    fig4b = "results/figures/ultra_hd/nested_cv_roc_curve_ultra_hd.png"
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
            "Figure 4b | Out-of-fold ROC curve from the leakage-controlled nested cross-validation analysis (n = 43).",
            "This is the primary performance visualisation because it reflects held-out predictions rather than the fitted full-cohort model.",
        )
    if not Path(fig4a).exists():
        doc.add_paragraph("[Figure 4: feature_importance / ROC assets not found.]")
    doc.add_paragraph()

    doc.add_heading("Independent rupture support and orthogonal tissue-blood context", 2)
    _n_sl, _n_gse, _n_dir = _gse_stats
    doc.add_paragraph(
        f"External support was sought in a hierarchy matched to the biological question. First, **GSE54083** provides an "
        f"independent rupture-versus-unrupture aneurysm microarray dataset ({_m['gse54083_n_ruptured']} ruptured, "
        f"{_m['gse54083_n_unruptured']} unruptured; total n={_m['gse54083_n_ruptured'] + _m['gse54083_n_unruptured']}). Within that small dataset, {_m['gse54083_n_shortlist_present']} immune-aware "
        f"shortlist genes were measurable, {_m['gse54083_n_shortlist_direction_match']} matched the GSE122897 rupture direction, "
        f"and {_m['gse54083_n_shortlist_p_lt_0_05']} were nominally significant, although none survived FDR<0.05. This is modest "
        f"and should be interpreted cautiously given the external sample size. Second, a targeted seven-gene check in **GSE13353** showed "
        f"{_m['gse13353_direction_matches']}/{_m['gse13353_genes_tested']} direction matches with {_m['gse13353_nominal_sig_genes']} "
        f"nominally significant genes. Third, **GSE26969** benchmarks unruptured aneurysm wall against non-aneurysmal superficial temporal "
        f"artery ({_m['gse26969_n_aneurysm']} aneurysm, {_m['gse26969_n_control']} control): {_m['gse26969_n_shortlist_direction_match']}/"
        f"{_m['gse26969_n_shortlist_present']} measurable shortlist genes matched the GSE122897 rupture sign and "
        f"{_m['gse26969_n_shortlist_p_lt_0_05']} were nominally significant, but none survived FDR. Fourth, **GSE15629** was treated in "
        f"two ways because its parent study mixes ruptured, unruptured, and control tissue. The rupture subset "
        f"({_cv156r.get('Samples', '8 R / 6 U')}) yielded {int(_cv156r.get('Direction_matched_genes', 0))}/{int(_cv156r.get('Measured_or_tested_genes', 1))} directionally matched "
        f"shortlist genes with {int(_cv156r.get('Nominal_P_lt_0_05', 0))} nominally significant signals but no FDR survivors, whereas the aneurysm-versus-control contrast "
        f"({_cv156c.get('Samples', '14 A / 5 C')}) showed {int(_cv156c.get('Direction_matched_genes', 0))}/{int(_cv156c.get('Measured_or_tested_genes', 1))} directionally matched genes and "
        f"{int(_cv156c.get('Nominal_P_lt_0_05', 0))} nominal signals, again without being interpreted as rupture replication.\n\n"
        f"We then separated blood cohorts into bridge versus context-only roles. In **GSE159610** (aneurysm vs control blood RNA-seq), "
        f"{_n_gse} shortlist genes had a non-missing contrast and {_n_dir} shared the same sign as GSE122897 rupture effects. In **GSE106520** "
        f"(unruptured aneurysm vs control neutrophil RNA-seq), {int(_cv106.get('Direction_matched_genes', 0))}/{int(_cv106.get('Measured_or_tested_genes', 1))} measurable shortlist genes matched direction "
        f"with {int(_cv106.get('Nominal_P_lt_0_05', 0))} nominal signals. We therefore treated GSE159610 and GSE106520 as the main blood bridge cohorts for reduced-panel derivation. "
        f"By contrast, **GSE36791** provides acute blood rupture/SAH context ({_m['gse36791_n_ruptured']} ruptured, {_m['gse36791_n_control']} control): "
        f"{_m['gse36791_n_shortlist_direction_match']}/{_m['gse36791_n_shortlist_present']} measurable shortlist genes matched direction, "
        f"{_m['gse36791_n_shortlist_p_lt_0_05']} were nominally significant, and {_m['gse36791_n_shortlist_fdr_lt_0_05']} survived FDR, "
        f"but this signal is biologically confounded by acute post-rupture systemic response. Because GSE36791 compares blood RNA after "
        f"rupture and subarachnoid haemorrhage against controls, it more likely reflects coagulation, inflammatory, and stress responses "
        f"to established bleeding than aneurysm-wall biology predictive of future rupture. **GSE73378** adds long-horizon whole-blood "
        f"context after past aneurysmal SAH ({_cv733.get('Samples', '103 past SAH / 107 C')}), with {int(_cv733.get('Direction_matched_genes', 0))}/"
        f"{int(_cv733.get('Measured_or_tested_genes', 1))} measurable shortlist genes matching direction but no FDR survivors. We therefore "
        f"use GSE26969, both GSE15629 context layers, GSE36791, and GSE73378 as tissue or blood context rather than direct rupture replication. "
        f"**Figure 5** and **Table 3** summarise this expanded hierarchy in a single view."
    )
    _combined_validation_fig = _base / "results" / "figures" / "ultra_hd" / "combined_validation_summary_ultra_hd.png"
    if _combined_validation_fig.exists():
        add_figure_with_caption_and_inference(
            doc,
            _combined_validation_fig,
            "Figure 5 | Combined hierarchy of external transcriptomic support across public datasets.",
            "Panel A summarises direction concordance of measurable or targeted genes with the GSE122897 shortlist; panel B summarises nominal and FDR-significant genes within each cohort-specific comparison. Direct rupture datasets remain the most relevant for replication, GSE15629 and GSE26969 provide tissue context, GSE159610 and GSE106520 act as blood bridges, and GSE36791 plus GSE73378 are retained as post-rupture / SAH blood context only.",
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
            inference="Direct rupture replication remains limited by small public cohorts, but the expanded hierarchy now separates secondary tissue context from the two blood bridge cohorts used for reduced-panel derivation and from post-rupture / SAH blood context.",
        max_rows=10,
    )
    doc.add_paragraph()

    doc.add_heading("Transportability, reduced blood panel, and direct-rupture consensus", 2)
    doc.add_paragraph(
        f"To move beyond gene-level overlap, we projected the fixed GSE122897 three-score model into each external cohort without "
        f"retraining. Transfer was heterogeneous. In the smaller direct rupture cohort **GSE13353**, discrimination was supportive "
        f"(AUROC {float(_t133.get('AUC', float('nan'))):.3f}, bootstrap 95% CI {float(_t133.get('AUC_bootstrap_ci_low', float('nan'))):.3f}–"
        f"{float(_t133.get('AUC_bootstrap_ci_high', float('nan'))):.3f}, permutation P={float(_t133.get('AUC_permutation_pvalue', float('nan'))):.4f}), "
        f"with higher mean transferred rupture probability in ruptured than unruptured samples. By contrast, **GSE54083** showed inverse "
        f"separation (AUROC {float(_t540.get('AUC', float('nan'))):.3f}; mean probability {float(_t540.get('Mean_probability_positive', float('nan'))):.3f} in ruptured "
        f"vs {float(_t540.get('Mean_probability_negative', float('nan'))):.3f} in unruptured), despite its gene-level directional overlap. "
        f"This discordance argues against stable portable rupture discrimination from the present tissue-derived score set. Orthogonal blood "
        f"rupture/SAH context in **GSE36791** showed only modest transfer (AUROC {float(_t367.get('AUC', float('nan'))):.3f}), while aneurysm-versus-control "
        f"context in **GSE159610** and **GSE26969** remained near-null or inverse (AUROC {float(_t159.get('AUC', float('nan'))):.3f} and "
        f"{float(_t269.get('AUC', float('nan'))):.3f}, respectively).\n\n"
        f"We therefore shifted the transportability target from a tissue-only classifier to a smaller blood-accessible panel. Restricting the shortlist to genes "
        f"that remained directionally concordant in both bridge cohorts (**GSE159610** and **GSE106520**) produced a reduced eight-gene panel "
        f"({_panel_gene_text}). This panel separated aneurysm from control blood in both bridge cohorts, with AUROC "
        f"{float(_p159.get('AUC', float('nan'))):.3f} in GSE159610 and {float(_p106.get('AUC', float('nan'))):.3f} in GSE106520. When carried into "
        f"the context-only cohorts, the same panel retained moderate signal in acute post-rupture **GSE36791** (AUROC {float(_p367.get('AUC', float('nan'))):.3f}) "
        f"but was essentially null in long-horizon post-SAH **GSE73378** (AUROC {float(_p733.get('AUC', float('nan'))):.3f}), reinforcing the distinction between bridge biology "
        f"and long-term post-bleed state. **Figure 6** and **Table 4** summarise this reduced blood-panel transportability.\n\n"
        f"We then asked whether gene-level consensus was more reproducible than cohort-level prediction. Among "
        f"{int(_meta_summary.get('n_overlap_all_three', 0))} genes measurable across **GSE122897**, **GSE54083**, and **GSE13353**, "
        f"{int(_meta_summary.get('n_directionally_concordant_all_three', 0))} were directionally concordant across all three cohorts and "
        f"remained significant at meta-analysis FDR<0.05. The strongest concordant genes were "
        f"TGFBI, ANPEP, MPP1, MARCO, and CLEC5A (higher in ruptured), and CPXM2 (lower in ruptured), suggesting that a smaller "
        f"rupture-associated molecular core is more reproducible than the full cohort-level classifier. **Figure 7** and **Table 5** summarise the direct-rupture "
        f"meta-analysis."
    )
    _blood_panel_fig = _base / "results" / "figures" / "ultra_hd" / "blood_bridge_panel_transport_ultra_hd.png"
    if _blood_panel_fig.exists():
        add_figure_with_caption_and_inference(
            doc,
            _blood_panel_fig,
            "Figure 6 | Reduced blood-testable panel anchored on the two blood bridge cohorts.",
            "The eight-gene panel was selected only from genes that remained directionally concordant in GSE159610 and GSE106520. It transported across the bridge cohorts and acute post-rupture blood context, but not across long-horizon post-SAH whole blood, so it remains a biologically motivated bridge signature rather than a validated clinical test.",
            width=Inches(6.5),
        )
    if not _panel_perf_df.empty:
        _panel_perf_tbl = _panel_perf_df[
            [
                "Dataset",
                "Role",
                "Samples",
                "Panel_genes_measured",
                "AUC",
                "Mann_Whitney_P",
                "Mean_score_positive",
                "Mean_score_negative",
                "Score_difference",
            ]
        ].copy()
        _panel_perf_tbl["Genes measured"] = _panel_perf_tbl["Panel_genes_measured"].map(lambda x: f"{int(x)}")
        _panel_perf_tbl["AUROC"] = _panel_perf_tbl["AUC"].map(lambda x: f"{float(x):.3f}")
        _panel_perf_tbl["Mann-Whitney P"] = _panel_perf_tbl["Mann_Whitney_P"].map(lambda x: f"{float(x):.4f}")
        _panel_perf_tbl["Mean score positive"] = _panel_perf_tbl["Mean_score_positive"].map(lambda x: f"{float(x):.3f}")
        _panel_perf_tbl["Mean score comparator"] = _panel_perf_tbl["Mean_score_negative"].map(lambda x: f"{float(x):.3f}")
        _panel_perf_tbl["Delta"] = _panel_perf_tbl["Score_difference"].map(lambda x: f"{float(x):.3f}")
        _panel_perf_tbl["Interpretation"] = _panel_perf_tbl["Dataset"].map(
            {
                "GSE159610": "Primary bridge cohort; supportive aneurysm-versus-control separation.",
                "GSE106520": "Independent neutrophil bridge cohort; strongest bridge separation.",
                "GSE36791": "Moderate signal in acute post-rupture / SAH blood context only.",
                "GSE73378": "Near-null in long-horizon post-SAH context; not a bridge cohort.",
            }
        )
        _panel_perf_tbl = _panel_perf_tbl[
            [
                "Dataset",
                "Role",
                "Samples",
                "Genes measured",
                "AUROC",
                "Mann-Whitney P",
                "Mean score positive",
                "Mean score comparator",
                "Delta",
                "Interpretation",
            ]
        ]
        add_table_from_dataframe(
            doc,
            _panel_perf_tbl,
            title="Table 4 | Reduced blood-testable panel transport across bridge and context blood cohorts.",
            inference="The reduced panel transports across the two bridge cohorts used for derivation, shows only moderate signal in acute post-rupture blood, and is near-null in long-horizon post-SAH whole blood, supporting a transportability story rather than clinical validation.",
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
        _meta_core_n = int(_meta_summary.get("n_directionally_concordant_all_three", 0)) or 6
        _meta_tbl = _meta_top_df.head(_meta_core_n)[
            [
                "Gene",
                "Meta_Z",
                "Meta_P",
                "Meta_FDR_BH",
                "Log2FC_GSE122897",
                "Log2FC_GSE54083",
                "Log2FC_GSE13353",
            ]
        ].copy()
        _meta_tbl["Stouffer Z"] = _meta_tbl["Meta_Z"].map(lambda x: f"{float(x):.3f}")
        _meta_tbl["meta-analysis P"] = _meta_tbl["Meta_P"].map(lambda x: f"{float(x):.2e}")
        _meta_tbl["meta-analysis FDR"] = _meta_tbl["Meta_FDR_BH"].map(lambda x: f"{float(x):.2e}")
        _meta_tbl["GSE122897 log2FC"] = _meta_tbl["Log2FC_GSE122897"].map(lambda x: f"{float(x):.3f}")
        _meta_tbl["GSE54083 log2FC"] = _meta_tbl["Log2FC_GSE54083"].map(lambda x: f"{float(x):.3f}")
        _meta_tbl["GSE13353 log2FC"] = _meta_tbl["Log2FC_GSE13353"].map(lambda x: f"{float(x):.3f}")
        _meta_tbl = _meta_tbl[
            [
                "Gene",
                "GSE122897 log2FC",
                "GSE54083 log2FC",
                "GSE13353 log2FC",
                "Stouffer Z",
                "meta-analysis P",
                "meta-analysis FDR",
            ]
        ]
        add_table_from_dataframe(
            doc,
            _meta_tbl,
            title="Table 5 | Top directionally concordant genes in the three-cohort rupture meta-analysis.",
            inference="These genes are the most reproducible transcript-level rupture signals across the public direct rupture cohorts and are more stable than the transferable cohort-level classifier.",
            max_rows=max(6, _meta_core_n),
        )
        doc.add_paragraph()

    doc.add_heading("Genetics-informed prioritisation: blood eQTLs, IA GWAS overlap, EXPRS, and optional PRS scaffolds", 2)
    doc.add_paragraph(
        f"We finally layer a reproducible genetics-facing prioritisation stack onto the {_m['n_shortlist']}-gene immune-aware "
        f"shortlist (not on reusing tissue labels as genotypes). In the codebase this layer is historically labelled **PHASES-G**, "
        f"but no PHASES score was calculated here. Querying blood cis-eQTL resources, {_m['n_blood_eqtl_catalogue']} "
        f"short-listed genes had a sentinel cis hit in the eQTL Catalogue whole-blood panel (per-gene top hit under the pipeline’s "
        f"pagination cap; **Supplementary Note 1**), and {_m['n_eqtlgen_cis']} genes had cis evidence in the eQTLGen layer we "
        f"attached. For IA GWAS overlap we intersected those sentinel SNPs (rsID-level) with the Bakker et al. summary-statistic "
        f"file: {_m['n_gwas_overlap_catalogue']} Catalogue-derived SNPs and {_m['n_gwas_overlap_eqtlgen']} eQTLGen-derived SNPs "
        f"appear in the GWAS rsID list—useful for hypothesis generation but **not** formal colocalisation (no COLOC/SMR in this "
        f"public RNA-only setting). Where individual-level genotypes and PHASES covariates exist on harmonised sample IDs, the "
        f"repository implements LD-aware clumping and a germline polygenic score using {_m['n_prs_weights_snps']} SNP weights "
        f"derived from the eQTL panel. Public GEO submissions lack dosages; as an RNA-only "
        f"placeholder we also compute EXPRS—eQTL effect sizes used as weights on z-scored log-expression—for {_m['n_exprs_genes']} "
        f"short-listed genes with TPM support. EXPRS and related AUROC demonstrations are explicitly treated here as **pipeline "
        f"illustrations** because the shortlist is outcome-informed in GSE122897; "
        f"they must not be read as independent external validity. Prospective validation requires DNA banking, PHASES fields, "
        f"and pre-registered analysis in an independent cohort (**Supplementary Note 1**)."
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
            inference="Returned blood eQTL hits are exploratory prioritisation metadata and may be incomplete because the scan was programmatic rather than exhaustive.",
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
            inference="rsID intersection only, not formal COLOC or SMR evidence.",
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
        f"genetics. Pathway enrichment and a three-signature ridge model compress that biology into interpretable scores, but the "
        f"key out-of-sample estimate is the leakage-controlled nested AUROC ({_m['auc_nested_mean']:.3f}±{_m['auc_nested_sd']:.3f}), "
        f"with the non-nested CV ({_m['auc_cv_mean']:.3f}±{_m['auc_cv_sd']:.3f}) and full-cohort fit ({_m['auc_full_cohort']:.3f}) "
        f"best viewed as optimistic comparators rather than primary validation.\n\n"
        f"A central caution is that the immune-adjustment layer rests on a weakly fitting leukocyte reference: median CIBERSORTx "
        f"correlation was only {_m['cibersortx_median_correlation']:.3f}, consistent with mismatch between LM22 and aneurysm-wall tissue. "
        f"Accordingly, we interpret immune-adjusted persistence as a robustness screen against coarse inferred immune content, not as "
        f"reliable cell-resolved deconvolution.\n\n"
        f"External support is strongest in the small rupture microarray GSE54083, where "
        f"{_m['gse54083_n_shortlist_direction_match']}/{_m['gse54083_n_shortlist_present']} shortlist genes matched direction and "
        f"{_m['gse54083_n_shortlist_p_lt_0_05']} were nominally significant, but none survived FDR. GSE13353 adds a small targeted "
        f"seven-gene check. GSE15629 adds both a secondary rupture subset and a stronger aneurysm-versus-control tissue contrast, "
        f"underscoring how mixed-design public cohorts can still be informative when their roles are separated. GSE26969 shows that "
        f"{_m['gse26969_n_shortlist_p_lt_0_05']} shortlist genes are nominally dysregulated between aneurysm wall and non-aneurysmal artery, "
        f"while GSE36791 shows a much stronger blood-side rupture/SAH signature ({_m['gse36791_n_shortlist_fdr_lt_0_05']} shortlist genes at FDR<0.05) "
        f"that is biologically informative but not equivalent to direct rupture validation. The bridge cohorts GSE159610 and GSE106520 are more useful "
        f"for asking whether a reduced blood-accessible signal exists at all, whereas GSE36791 and GSE73378 are better interpreted as acute or long-horizon "
        f"post-bleed context.\n\n"
        f"When we projected the fixed three-score model across cohorts, performance was markedly asymmetric: GSE13353 showed supportive "
        f"direct transfer (AUROC {float(_t133.get('AUC', float('nan'))):.3f}), whereas GSE54083 showed inverse transfer "
        f"(AUROC {float(_t540.get('AUC', float('nan'))):.3f}), despite gene-level directional overlap. This instability argues that the "
        f"current tissue-derived classifier should not be interpreted as externally validated. By contrast, the reduced eight-gene blood panel "
        f"({ _panel_gene_text }) transported in both bridge cohorts (AUROC {float(_p159.get('AUC', float('nan'))):.3f} in GSE159610 and "
        f"{float(_p106.get('AUC', float('nan'))):.3f} in GSE106520), showed only moderate acute post-rupture signal in GSE36791 (AUROC "
        f"{float(_p367.get('AUC', float('nan'))):.3f}), and was null in long-horizon GSE73378 (AUROC {float(_p733.get('AUC', float('nan'))):.3f}). "
        f"That pattern is more consistent with a transportable aneurysm-biology bridge than with a fully generalisable rupture-risk test. The stricter three-cohort "
        f"gene-level meta-analysis yielded {int(_meta_summary.get('n_directionally_concordant_all_three', 0))} directionally concordant "
        f"genes at meta-analysis FDR<0.05, led by TGFBI, ANPEP, MPP1, MARCO, CLEC5A, and CPXM2, suggesting that a small "
        f"rupture-associated consensus core is more reproducible than cohort-level risk discrimination. "
        f"The genetics-informed prioritisation layer then maps the "
        f"shortlist onto blood eQTLs ({_m['n_blood_eqtl_catalogue']} / {_m['n_eqtlgen_cis']} genes in the Catalogue / eQTLGen "
        f"layers we queried) and IA GWAS rsID overlap ({_m['n_gwas_overlap_catalogue']} / {_m['n_gwas_overlap_eqtlgen']} SNPs), "
        f"preserving a clean separation between tissue discovery and germline inference. EXPRS and related cohort demos should "
        f"be viewed as scaffolding until genotype–PHASES datasets exist.\n\n"
        f"Limitations remain substantial: small n, single public surgical cohort, incomplete clinical metadata, likely confounding "
        f"by surgical timing and cellular admixture, no same-subject genotype data, and no true external validation of a rupture-risk "
        f"model in a clinically deployable sample type. Future work should pair accessible biomarker sampling with harmonised PHASES "
        f"covariates, germline genotypes, and pre-registered external validation rather than extending outcome-informed tissue models."
    )

    doc.add_heading("Contributors", 1)
    doc.add_paragraph(
        "[To be completed at submission: initials and full names with CRediT roles—Conceptualisation, Data curation, "
        "Formal analysis, Funding acquisition, Investigation, Methodology, Project administration, Resources, Software, "
        "Supervision, Validation, Visualisation, Writing—original draft, Writing—review and editing.]"
    )

    doc.add_heading("Declaration of interests", 1)
    doc.add_paragraph(
        "The authors declare no competing interests. [Update if any financial or non-financial relationships exist.]"
    )

    doc.add_heading("Role of the funding source", 1)
    doc.add_paragraph(
        "The funders of the study had no role in study design, data collection, data analysis, data interpretation, or writing "
        "of the report. [Replace with sponsor-specific text if applicable.]"
    )

    doc.add_heading("Acknowledgments", 1)
    doc.add_paragraph(
        "We thank the investigators and participants contributing to the public GSE122897, GSE54083, GSE13353, GSE15629, GSE26969, GSE106520, GSE36791, GSE73378, and GSE159610 "
        "datasets, the GEO curators, and the developers of open-source analysis tools cited in the Appendix."
    )

    doc.add_heading("Data sharing", 1)
    doc.add_paragraph(
        "Processed tables, figure assets, and analysis scripts reside in the project repository supporting this draft. "
        "This includes the external validation and context outputs generated from GSE54083, GSE13353, GSE26969, GSE36791, and GSE159610 alongside the primary GSE122897 "
        "results. Individual-level genotypes were not available for GSE122897; therefore germline PRS and full genetics-linked "
        "validation demonstrations are provided as reproducible scaffolds only. Any future sharing of additional data will comply "
        "with institutional and participant consent requirements."
    )

    doc.add_heading("Data availability", 1)
    doc.add_paragraph(
        "Primary dataset: GSE122897 (NCBI GEO). External expression datasets used for contextual support: GSE54083, GSE13353, GSE26969, GSE36791, and "
        "GSE159610 (NCBI GEO). Processed tables and analysis scripts: project repository / authors upon request."
    )

    doc.add_heading("References", 1)
    doc.add_paragraph(
        "References are numbered in order of appearance (Vancouver style). Replace bracketed placeholders with final "
        "citations at submission; journal limits may apply (e.g. maximum reference count—verify current instructions)."
    )
    refs = [
        "Greving JP, Wermer MJH, Brown RD Jr, et al. Development of the PHASES score for prediction of risk of rupture of intracranial aneurysms: a pooled analysis of six prospective cohort studies. Lancet Neurol 2014; 13: 59-66.",
        "Kleinloog R, Verweij BH, van der Vlies P, et al. RNA sequencing analysis of intracranial aneurysm walls reveals involvement of lysosomes and immunoglobulins in rupture. Stroke 2016; 47: 1286-93.",
        "Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. "
        "Genome Biol 2014; 15: 550.",
        "Newman AM, Steen CB, Liu CL, et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol 2019; 37: 773-82.",
        "Steen CB, Liu CL, Alizadeh AA, Newman AM. Profiling cell type abundance and expression in bulk tissues with CIBERSORTx. Nat Biotechnol 2019; 37: 1475-80.",
        "Bakker MK, et al. Genome-wide association study of intracranial aneurysms in European cohorts. Nat Genet 2020; "
        "52: 1303–13. (summary statistics as used in pipeline)",
        "The GTEx Consortium. The Genotype-Tissue Expression (GTEx) pilot analysis. Nature 2013; 496: 214–20. "
        "(context for eQTL interpretation)",
        "Kerimov N, Hayhurst JD, Peikova K, et al. A compendium of uniformly processed human gene expression and splicing quantitative trait loci. Nat Genet 2021; 53: 1290-99.",
        "Võsa U, et al. Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores "
        "for blood gene expression. Nat Commun 2021; 12: 2434. (eQTLGen resource)",
        "Raudvere U, et al. g:Profiler: a web server for functional enrichment analysis and conversions of gene lists "
        "(2019 update). Nucleic Acids Res 2019; 47: W191–98.",
        "Pedregosa F, et al. Scikit-learn: machine learning in Python. J Mach Learn Res 2011; 12: 2825–30.",
        "Kurki MI, et al. Gene expression profiling of intracranial aneurysms (GSE54083). NCBI GEO. "
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54083",
        "Pera J, et al. Comparison of gene expression between ruptured and unruptured human intracranial aneurysms (GSE13353). "
        "NCBI GEO. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13353",
        "Li L, Yang X, Jiang F, et al. Transcriptome-wide characterization of gene expression associated with unruptured intracranial aneurysms (GSE26969). NCBI GEO. "
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26969",
        "Pera J, Korostynski M, Piechota M, et al. Gene expression profiling of blood in ruptured intracranial aneurysms (GSE36791). NCBI GEO. "
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36791",
        "Hasan D, et al. RNA expression in human cerebral aneurysms and control arteries (GSE159610). NCBI GEO. "
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159610",
    ]
    for i, r in enumerate(refs, start=1):
        doc.add_paragraph(f"{i}. {r}")
    doc.add_page_break()

    # --- Appendix: supplementary materials embedded in the same file ---
    doc.add_heading("Appendix", 1)
    doc.add_paragraph(
        "The Appendix contains Supplementary Methods, Supplementary Figures S1–S12, embedded Supplementary Tables S1–S6, "
        "and Supplementary Note 1 (repository outputs and reproducibility notes). This satisfies a single-file submission draft; "
        "the journal may require separate upload of high-resolution figures at acceptance."
    )
    doc.add_paragraph()

    # --- Supplementary (full Methods + figs S1–S12 + embedded tables S1–S6) ---
    doc.add_heading("Supplementary Materials", 2)

    doc.add_heading("Supplementary Note 1. Repository outputs, deconvolution caveats, and genetics prioritisation layer", 2)
    doc.add_paragraph(
        "Key output files including shortlist gene tables, eQTL overlap panels, and PRS scaffolds are deposited in the public code repository "
        "(DOI: [to be assigned at submission]). The same repository should also host the external transfer summaries, "
        "direct-rupture meta-analysis tables, and figure assets used in this appendix. The deconvolution outputs comprise archived CIBERSORTx portal "
        "results, locally reproduced LM22+NNLS fractions, and immune-adjusted gene-level regression tables. Because the LM22 reference "
        "is leukocyte-derived and the observed CIBERSORTx fit was weak, these outputs should be interpreted as sensitivity analyses "
        "rather than validated cell-composition estimates.\n\n"
        "The downstream genetics-informed prioritisation layer links the tissue shortlist to blood cis-eQTL resources, rsID-level IA "
        "GWAS overlap, expression-weighted exploratory scores, and optional PRS scaffolds for future cohorts with same-subject genotype "
        "and clinical data. In the repository this layer is historically labelled PHASES-G, but the PHASES score itself was not computed "
        "in the present manuscript."
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

    out = script_dir / "Manuscript_Lancet_Digital_Health.docx"
    doc.save(str(out))
    print(f"Created: {out}")
    print(
        "Lancet Digital Health-style draft: Summary (structured), Research in context, IMRAD, declarations, "
        "Vancouver references, Appendix with full supplementary bundle. Main Figures 1-7 and Tables 1-5 embedded."
    )
    return str(out)


if __name__ == "__main__":
    create_lancet_digital_health_manuscript()
