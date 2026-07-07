"""
Create comprehensive Nature Journal Review Document
"""
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from pathlib import Path

def create_review_document():
    """Create comprehensive Nature journal review document"""
    
    doc = Document()
    
    # Set default font
    style = doc.styles['Normal']
    font = style.font
    font.name = 'Calibri'
    font.size = Pt(11)
    
    # Title
    title = doc.add_heading('Nature Journal Review: Critical Assessment', 0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    
    subtitle = doc.add_paragraph('Manuscript: Multi-Omics Integrative Analysis for Cerebral Aneurysm Rupture Risk Prediction')
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    subtitle.runs[0].italic = True
    
    doc.add_paragraph()
    doc.add_paragraph('Reviewer: Nature Journal Standards')
    doc.add_paragraph('Date: November 2025')
    doc.add_paragraph('Decision: MAJOR REVISION REQUIRED')
    
    doc.add_page_break()
    
    # Executive Summary
    doc.add_heading('Executive Summary', 1)
    
    summary = """
    This manuscript presents a computational study developing a multi-omics model for predicting 
    cerebral aneurysm rupture risk. While the work shows promise and addresses an important clinical 
    question, several critical issues prevent it from meeting Nature publication standards in its 
    current form.
    
    The study uses only 43 samples for model development and validation, lacks independent external 
    validation, has incomplete methodology descriptions, and makes premature claims about clinical 
    translation. The statistical rigor needs improvement, and experimental validation is absent.
    
    RECOMMENDATION: MAJOR REVISION REQUIRED
    
    The manuscript requires substantial additional work before it can be considered for Nature. 
    Alternatively, this work may be better suited for a specialized journal (e.g., Stroke, Journal 
    of Neurointerventional Surgery) as a computational discovery study, with clear acknowledgment 
    of limitations and future validation needs.
    """
    
    doc.add_paragraph(summary.strip())
    
    # Critical Issues
    doc.add_heading('CRITICAL ISSUES (Must Address for Publication)', 1)
    
    # Issue 1
    doc.add_heading('1. Sample Size Limitations', 2)
    doc.add_paragraph('Severity: MAJOR', style='Intense Quote')
    
    doc.add_paragraph("""
    The study uses only 43 samples (22 ruptured, 21 unruptured) for model development and validation. 
    This is critically insufficient for a Nature-level publication.
    """)
    
    doc.add_heading('Issues:', 3)
    issues1 = [
        "Power analysis not provided - cannot determine if sample size is adequate",
        "Cross-validation with n=43 and 5 folds means ~8-9 samples per fold - highly unstable",
        "CV AUC SD of ±0.156 indicates extreme variability across folds",
        "No independent test set - all data used for both training and validation",
        "Risk of overfitting despite regularization (model may still be overfit)"
    ]
    for issue in issues1:
        doc.add_paragraph(issue, style='List Bullet')
    
    doc.add_heading('Required Changes:', 3)
    changes1 = [
        "Provide power analysis showing required sample size for adequate statistical power",
        "Acknowledge sample size as a major limitation in abstract and discussion",
        "Use more conservative validation (leave-one-out CV or nested CV)",
        "Report confidence intervals for all performance metrics",
        "Reframe as a pilot/exploratory study, not definitive validation",
        "Consider this a proof-of-concept requiring larger validation"
    ]
    for change in changes1:
        doc.add_paragraph(change, style='List Bullet')
    
    # Issue 2
    doc.add_heading('2. No Independent External Validation', 2)
    doc.add_paragraph('Severity: MAJOR', style='Intense Quote')
    
    doc.add_paragraph("""
    The manuscript lacks independent external validation, which is essential for Nature. All validation 
    is internal (cross-validation on the same dataset).
    """)
    
    doc.add_heading('Issues:', 3)
    issues2 = [
        "All validation is internal (cross-validation on same dataset)",
        "No external dataset validation (different institution, platform, or time period)",
        "No temporal validation (data from different time periods)",
        "No prospective validation cohort",
        "Single dataset (GSE122897) used for everything - dataset-specific findings may not generalize"
    ]
    for issue in issues2:
        doc.add_paragraph(issue, style='List Bullet')
    
    doc.add_heading('Required Changes:', 3)
    changes2 = [
        "Validate on at least one independent external dataset (e.g., GSE159610, GSE26969)",
        "If external data unavailable, clearly state this as a critical limitation",
        "Report performance on external data separately from internal CV",
        "Acknowledge that internal CV alone is insufficient for Nature",
        "Consider using other GEO datasets for validation",
        "Clearly separate discovery cohort from validation cohort"
    ]
    for change in changes2:
        doc.add_paragraph(change, style='List Bullet')
    
    # Issue 3
    doc.add_heading('3. Statistical Methods and Reporting', 2)
    doc.add_paragraph('Severity: MAJOR', style='Intense Quote')
    
    doc.add_paragraph("""
    Several statistical issues need addressing to meet Nature standards.
    """)
    
    doc.add_heading('Issues:', 3)
    issues3 = [
        "Multiple testing correction: Multiple thresholds tested (0.01, 0.05, 0.1, 0.2) without justification",
        "Effect sizes: Log2FC reported but no confidence intervals",
        "Model performance: Only point estimates reported, no confidence intervals",
        "Cross-validation: 5-fold CV with n=43 is unstable - should use nested CV or LOOCV",
        "Overfitting assessment: Gap of 0.029 may still indicate overfitting given small n",
        "No permutation testing to assess if model performance is significantly above chance"
    ]
    for issue in issues3:
        doc.add_paragraph(issue, style='List Bullet')
    
    doc.add_heading('Required Changes:', 3)
    changes3 = [
        "Report 95% confidence intervals for all metrics (AUC, accuracy, precision, recall)",
        "Use nested cross-validation for hyperparameter tuning to avoid overfitting",
        "Perform permutation testing (shuffle labels 1000+ times) to show model performance is above chance",
        "Report effect sizes (Log2FC) with confidence intervals",
        "Justify statistical thresholds used (why FDR 0.05 vs 0.01?)",
        "Report bootstrap confidence intervals for model performance",
        "Use DeLong test or similar to compare ROC curves statistically"
    ]
    for change in changes3:
        doc.add_paragraph(change, style='List Bullet')
    
    # Issue 4
    doc.add_heading('4. Methodology Details Insufficient', 2)
    doc.add_paragraph('Severity: MAJOR', style='Intense Quote')
    
    doc.add_paragraph("""
    Several methodological details are missing or unclear, preventing reproducibility and full understanding.
    """)
    
    doc.add_heading('Issues:', 3)
    issues4 = [
        "Feature engineering: Exact formulas for signature scores not provided",
        "Gene selection criteria for signature scores unclear",
        "Weighting scheme for signature scores not described",
        "Model selection: Why Logistic Regression chosen? Need formal comparison with statistical tests",
        "Hyperparameter tuning process not detailed (search space, selection criteria)",
        "Risk stratification thresholds: How were 0.33 and 0.67 chosen? Need clinical or statistical justification",
        "Missing data handling: Not mentioned",
        "Batch effects: Not addressed (if samples from different batches)",
        "Quality control: Criteria for sample exclusion not detailed"
    ]
    for issue in issues4:
        doc.add_paragraph(issue, style='List Bullet')
    
    doc.add_heading('Required Changes:', 3)
    changes4 = [
        "Provide exact mathematical formulas for all signature scores in Methods section",
        "Detail gene selection criteria (how many genes, which genes, why?)",
        "Describe weighting scheme (equal weights? significance-weighted? fold-change weighted?)",
        "Provide formal model comparison with statistical tests (not just performance metrics)",
        "Detail hyperparameter search space and selection criteria",
        "Justify all thresholds (0.33, 0.67) with clinical rationale or data-driven optimization",
        "Describe missing data handling explicitly",
        "Address potential batch effects (PCA, batch correction if needed)",
        "Provide detailed QC criteria and sample exclusion rationale"
    ]
    for change in changes4:
        doc.add_paragraph(change, style='List Bullet')
    
    # Issue 5
    doc.add_heading('5. Insufficient Comparison with Existing Methods', 2)
    doc.add_paragraph('Severity: MAJOR', style='Intense Quote')
    
    doc.add_paragraph("""
    The manuscript claims superiority over PHASES score but comparison is inadequate.
    """)
    
    doc.add_heading('Issues:', 3)
    issues5 = [
        "PHASES score not actually calculated on this dataset",
        "Comparison is to literature-reported PHASES performance (AUC 0.65-0.72), not direct comparison",
        "No direct head-to-head comparison on same patients",
        "Other existing models not compared (e.g., other ML models from literature)",
        "No statistical test for superiority (e.g., DeLong test for ROC comparison)"
    ]
    for issue in issues5:
        doc.add_paragraph(issue, style='List Bullet')
    
    doc.add_heading('Required Changes:', 3)
    changes5 = [
        "Calculate PHASES score for all patients in dataset (requires clinical data: age, hypertension, size, etc.)",
        "Perform direct comparison on same patients (head-to-head)",
        "Use statistical test (DeLong test) to compare ROC curves and show if improvement is significant",
        "Compare with other published ML models for aneurysm rupture if applicable",
        "Report p-value for superiority claim",
        "If clinical data unavailable for PHASES, acknowledge this limitation and tone down claims"
    ]
    for change in changes5:
        doc.add_paragraph(change, style='List Bullet')
    
    # Issue 6
    doc.add_heading('6. Lack of Experimental Validation', 2)
    doc.add_paragraph('Severity: MAJOR', style='Intense Quote')
    
    doc.add_paragraph("""
    No experimental validation of findings, which is typically expected for Nature.
    """)
    
    doc.add_heading('Issues:', 3)
    issues6 = [
        "All analyses are computational - no wet lab validation",
        "Biomarker candidates not validated experimentally (qPCR, targeted sequencing)",
        "No functional studies to confirm mechanisms (knockdown, overexpression)",
        "No animal model validation",
        "No independent cohort validation with experimental methods"
    ]
    for issue in issues6:
        doc.add_paragraph(issue, style='List Bullet')
    
    doc.add_heading('Required Changes:', 3)
    changes6 = [
        "Validate top 5-10 biomarkers in independent samples using qPCR or targeted sequencing",
        "Perform functional studies on key genes (e.g., CRISPR knockdown, overexpression in cell culture)",
        "If experimental validation not possible, clearly state this as a major limitation",
        "Reframe as computational discovery study requiring experimental validation",
        "Consider animal model studies for top candidates",
        "Acknowledge that findings are computational predictions requiring validation"
    ]
    for change in changes6:
        doc.add_paragraph(change, style='List Bullet')
    
    # Issue 7
    doc.add_heading('7. Reproducibility Concerns', 2)
    doc.add_paragraph('Severity: MAJOR', style='Intense Quote')
    
    doc.add_paragraph("""
    Reproducibility information is insufficient for Nature standards.
    """)
    
    doc.add_heading('Issues:', 3)
    issues7 = [
        "Code availability: Not mentioned in manuscript",
        "Data availability: GSE122897 is public, but processed data (signature scores, predictions) not provided",
        "Random seeds: Not all mentioned (some use random_state=42, others unclear)",
        "Software versions: Not specified",
        "Exact parameters: Some model parameters not fully specified",
        "Environment: No requirements.txt or conda environment file"
    ]
    for issue in issues7:
        doc.add_paragraph(issue, style='List Bullet')
    
    doc.add_heading('Required Changes:', 3)
    changes7 = [
        "Provide code repository (GitHub) with all analysis scripts - make publicly available",
        "Provide processed data files (signature scores, model predictions) in repository",
        "Specify all random seeds used in all analyses",
        "Provide software versions (requirements.txt, conda environment.yml)",
        "Include detailed README with step-by-step reproduction instructions",
        "Add Data Availability and Code Availability statements to manuscript",
        "Ensure all analyses can be reproduced from raw data"
    ]
    for change in changes7:
        doc.add_paragraph(change, style='List Bullet')
    
    # Issue 8
    doc.add_heading('8. Clinical Translation Claims Unsupported', 2)
    doc.add_paragraph('Severity: MAJOR', style='Intense Quote')
    
    doc.add_paragraph("""
    Claims about clinical utility are premature given the study design.
    """)
    
    doc.add_heading('Issues:', 3)
    issues8 = [
        "No prospective validation",
        "No cost-effectiveness analysis",
        "No comparison with current clinical workflow",
        "No assessment of clinical decision impact",
        "Risk thresholds (0.33, 0.67) not clinically validated",
        "No discussion of implementation challenges"
    ]
    for issue in issues8:
        doc.add_paragraph(issue, style='List Bullet')
    
    doc.add_heading('Required Changes:', 3)
    changes8 = [
        "Tone down clinical translation claims throughout manuscript",
        "Clearly state this is a discovery/development study, not a validated clinical tool",
        "Discuss what would be needed for clinical implementation (prospective study, cost analysis, etc.)",
        "Acknowledge that clinical validation is future work",
        "Provide realistic timeline for translation (years, not months)",
        "Separate 'clinical potential' from 'clinical validation'",
        "Add section on 'Path to Clinical Translation' with specific requirements"
    ]
    for change in changes8:
        doc.add_paragraph(change, style='List Bullet')
    
    doc.add_page_break()
    
    # Moderate Issues
    doc.add_heading('MODERATE ISSUES (Should Address)', 1)
    
    doc.add_heading('1. Figure Quality and Presentation', 2)
    doc.add_paragraph('Severity: MODERATE', style='Intense Quote')
    
    doc.add_paragraph("""
    Some figures could be improved for publication quality.
    """)
    
    doc.add_heading('Issues:', 3)
    fig_issues = [
        "Some figures may be too complex (information overload)",
        "Color schemes may not be colorblind-friendly",
        "Figure legends could be more detailed",
        "Some figures lack statistical annotations (p-values, effect sizes, confidence intervals)"
    ]
    for issue in fig_issues:
        doc.add_paragraph(issue, style='List Bullet')
    
    doc.add_heading('Required Changes:', 3)
    fig_changes = [
        "Simplify complex figures or split into multiple panels",
        "Use colorblind-friendly palettes (viridis, cividis, or test with colorblind simulators)",
        "Add detailed figure legends explaining all elements",
        "Add statistical annotations to all relevant figures (p-values, effect sizes, CIs)",
        "Ensure all figures are publication-quality (300+ DPI, vector format preferred)",
        "Add sample sizes to all figures"
    ]
    for change in fig_changes:
        doc.add_paragraph(change, style='List Bullet')
    
    doc.add_heading('2. Literature Context', 2)
    doc.add_paragraph('Severity: MODERATE', style='Intense Quote')
    
    doc.add_paragraph("""
    Literature review and context could be strengthened.
    """)
    
    doc.add_heading('Issues:', 3)
    lit_issues = [
        "Limited discussion of previous ML models for aneurysm rupture",
        "Incomplete comparison with other omics studies",
        "Missing discussion of why this approach is novel",
        "Limited mechanistic discussion linking findings to biology"
    ]
    for issue in lit_issues:
        doc.add_paragraph(issue, style='List Bullet')
    
    doc.add_heading('Required Changes:', 3)
    lit_changes = [
        "Comprehensive literature review of existing prediction models (ML and clinical)",
        "Clear articulation of novelty - what makes this different?",
        "Better mechanistic discussion linking findings to known biology",
        "Discussion of how this fits into existing knowledge",
        "Compare with similar multi-omics studies in other diseases",
        "Discuss limitations of previous approaches that this addresses"
    ]
    for change in lit_changes:
        doc.add_paragraph(change, style='List Bullet')
    
    doc.add_page_break()
    
    # Minor Issues
    doc.add_heading('MINOR ISSUES (Nice to Address)', 1)
    
    doc.add_heading('1. Writing and Presentation', 2)
    doc.add_paragraph('Severity: MINOR', style='Intense Quote')
    
    doc.add_heading('Issues:', 3)
    writing_issues = [
        "Some technical terms could be better explained",
        "Abbreviations not always defined on first use",
        "Some sentences could be more concise",
        "Consistency in terminology throughout"
    ]
    for issue in writing_issues:
        doc.add_paragraph(issue, style='List Bullet')
    
    doc.add_heading('Required Changes:', 3)
    writing_changes = [
        "Define all abbreviations on first use",
        "Simplify technical language where possible (without losing accuracy)",
        "Ensure consistent terminology throughout (e.g., 'ruptured' vs 'rupture' consistently)",
        "Professional copy-editing recommended",
        "Ensure Nature style guide compliance"
    ]
    for change in writing_changes:
        doc.add_paragraph(change, style='List Bullet')
    
    doc.add_page_break()
    
    # Further Advancements
    doc.add_heading('FURTHER ADVANCEMENTS NEEDED FOR NATURE PUBLICATION', 1)
    
    advancements = [
        {
            'title': '1. Data Integration',
            'items': [
                "Integrate additional omics layers (genomics/GWAS, proteomics, metabolomics)",
                "Include imaging data (aneurysm morphology, hemodynamics, size, location)",
                "Integrate clinical variables more systematically (age, hypertension, smoking, etc.)",
                "Multi-center data integration to increase sample size and generalizability"
            ]
        },
        {
            'title': '2. Model Development',
            'items': [
                "Develop ensemble methods combining multiple algorithms",
                "Use deep learning for non-linear pattern discovery (with proper validation)",
                "Implement transfer learning from related diseases (vascular, cardiovascular)",
                "Develop patient-specific models (personalized medicine approach)"
            ]
        },
        {
            'title': '3. Validation',
            'items': [
                "Multi-center prospective validation study (n > 200 minimum)",
                "Longitudinal follow-up validation (time-to-event analysis)",
                "Cross-platform validation (different sequencing technologies)",
                "Validation in diverse populations (different ethnicities, geographies)"
            ]
        },
        {
            'title': '4. Mechanistic Understanding',
            'items': [
                "Functional genomics studies (CRISPR-Cas9, RNAi screens)",
                "Animal model validation (mouse/rat aneurysm models)",
                "Single-cell resolution studies to identify cell-type-specific changes",
                "Pathway perturbation experiments to validate mechanisms"
            ]
        },
        {
            'title': '5. Clinical Translation',
            'items': [
                "Develop clinical decision support tool (CDSS)",
                "Cost-effectiveness analysis comparing to current standard of care",
                "Implementation study in clinical setting",
                "Clinical trial design for intervention studies based on risk stratification"
            ]
        }
    ]
    
    for adv in advancements:
        doc.add_heading(adv['title'], 2)
        for item in adv['items']:
            doc.add_paragraph(item, style='List Bullet')
    
    doc.add_page_break()
    
    # Overall Assessment
    doc.add_heading('OVERALL ASSESSMENT AND RECOMMENDATION', 1)
    
    assessment = """
    CURRENT STATUS: MAJOR REVISION REQUIRED
    
    The manuscript presents interesting findings and addresses an important clinical question. 
    However, several critical issues must be addressed before it can be considered for Nature 
    publication.
    
    CRITICAL BARRIERS TO PUBLICATION:
    
    1. Sample Size: n=43 is insufficient for Nature. Minimum n > 200 recommended for validation.
    2. External Validation: No independent validation - essential for Nature.
    3. Statistical Rigor: Missing CIs, permutation testing, nested CV.
    4. Methodology: Incomplete descriptions prevent reproducibility.
    5. Experimental Validation: No wet lab validation of findings.
    6. Reproducibility: Code and data not provided.
    7. Clinical Claims: Premature claims about clinical utility.
    
    RECOMMENDATION:
    
    Option 1 (Preferred for Nature): 
    - Expand sample size to n > 200 through multi-center collaboration
    - Perform independent external validation on separate dataset
    - Validate top biomarkers experimentally (qPCR, functional studies)
    - Provide complete reproducibility package (code, data, environment)
    - Tone down clinical claims, reframe as discovery study
    - Add all required statistical analyses (CIs, permutation testing, nested CV)
    - Complete methodology descriptions
    
    Option 2 (Alternative):
    - Submit to specialized journal (Stroke, Journal of Neurointerventional Surgery, 
      Scientific Reports) as computational discovery study
    - Clearly acknowledge all limitations
    - Position as proof-of-concept requiring future validation
    - Focus on biological insights rather than clinical translation
    
    TIMELINE FOR NATURE PUBLICATION:
    If pursuing Nature, expect 12-24 months of additional work:
    - 3-6 months: Expand dataset, external validation
    - 6-12 months: Experimental validation of biomarkers
    - 3-6 months: Complete statistical analyses, reproducibility package
    - 3-6 months: Manuscript revision, resubmission
    
    The work shows promise but needs substantial additional validation and experimental 
    work before reaching Nature standards. The computational approach is sound, but 
    the validation framework needs significant strengthening.
    """
    
    doc.add_paragraph(assessment.strip())
    
    # Save
    output_file = 'Nature_Journal_Review_Critical_Assessment.docx'
    doc.save(output_file)
    print(f"\n✓ Review document created: {output_file}")
    
    return output_file

if __name__ == "__main__":
    create_review_document()

