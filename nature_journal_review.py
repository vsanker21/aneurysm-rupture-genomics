"""
Nature Journal Review: Critical Assessment of the Manuscript
Reviewing as a strict Nature journal reviewer
"""
from docx import Document
from pathlib import Path
import pandas as pd

def conduct_nature_review():
    """Conduct a comprehensive Nature journal-style review"""
    
    print("=" * 80)
    print("NATURE JOURNAL REVIEW: CRITICAL ASSESSMENT")
    print("=" * 80)
    print("\nReviewing: Complete_Analysis_Report_With_Inferences.docx")
    print("Reviewer Perspective: Strict Nature Journal Standards")
    print("=" * 80)
    
    # Load document
    doc = Document('Complete_Analysis_Report_With_Inferences.docx')
    
    # Load key data files for verification
    try:
        features_df = pd.read_csv('results/integrated_features.tsv', sep='\t')
        risk_df = pd.read_csv('results/phase7_risk_stratification.tsv', sep='\t')
        de_results = pd.read_csv('results/gse122897_differential_expression.tsv', sep='\t')
    except:
        features_df = None
        risk_df = None
        de_results = None
    
    review = []
    
    # ============================================================================
    # CRITICAL ISSUE 1: SAMPLE SIZE
    # ============================================================================
    review.append({
        'category': 'CRITICAL',
        'issue': 'Sample Size Limitations',
        'severity': 'Major',
        'description': """
        The study uses only 43 samples (22 ruptured, 21 unruptured) for model development and validation.
        This is critically insufficient for a Nature-level publication.
        
        Issues:
        1. Power analysis not provided - cannot determine if sample size is adequate
        2. Cross-validation with n=43 and 5 folds means ~8-9 samples per fold - highly unstable
        3. CV AUC SD of ±0.156 indicates extreme variability
        4. No independent test set - all data used for both training and validation
        5. Risk of overfitting despite regularization (model may still be overfit)
        
        Required Changes:
        - Provide power analysis showing required sample size
        - Acknowledge sample size as major limitation
        - Use more conservative validation (leave-one-out or nested CV)
        - Report confidence intervals for all metrics
        - Consider this a pilot/exploratory study, not definitive validation
        """,
        'recommendation': 'MAJOR REVISION REQUIRED - Sample size insufficient for Nature publication'
    })
    
    # ============================================================================
    # CRITICAL ISSUE 2: LACK OF INDEPENDENT VALIDATION
    # ============================================================================
    review.append({
        'category': 'CRITICAL',
        'issue': 'No Independent External Validation',
        'severity': 'Major',
        'description': """
        The manuscript lacks independent external validation, which is essential for Nature.
        
        Issues:
        1. All validation is internal (cross-validation on same dataset)
        2. No external dataset validation (e.g., different institution, different platform)
        3. No temporal validation (data from different time periods)
        4. No prospective validation cohort
        5. Single dataset (GSE122897) used for everything
        
        Required Changes:
        - Validate on at least one independent external dataset
        - If external data unavailable, clearly state this as a limitation
        - Consider using other GEO datasets (GSE159610, GSE26969) for validation
        - Report performance on external data separately
        - Acknowledge that internal CV alone is insufficient
        """,
        'recommendation': 'MAJOR REVISION REQUIRED - External validation essential'
    })
    
    # ============================================================================
    # CRITICAL ISSUE 3: STATISTICAL RIGOR
    # ============================================================================
    review.append({
        'category': 'CRITICAL',
        'issue': 'Statistical Methods and Reporting',
        'severity': 'Major',
        'description': """
        Several statistical issues need addressing:
        
        Issues:
        1. Multiple testing correction: Only FDR mentioned, but multiple thresholds tested (0.01, 0.05, 0.1, 0.2)
           - Need to justify threshold selection
           - Should report results at most conservative threshold
        2. Effect sizes: Log2FC reported but no confidence intervals
        3. Model performance: Only point estimates reported, no confidence intervals
        4. Cross-validation: 5-fold CV with n=43 is unstable - should use nested CV or LOOCV
        5. Overfitting assessment: Gap of 0.029 may still indicate overfitting given small n
        6. No permutation testing to assess significance of model performance
        
        Required Changes:
        - Report 95% confidence intervals for all metrics
        - Use nested cross-validation for hyperparameter tuning
        - Perform permutation testing (shuffle labels) to show model performance is above chance
        - Report effect sizes with confidence intervals
        - Justify statistical thresholds used
        """,
        'recommendation': 'MAJOR REVISION REQUIRED - Statistical rigor needs improvement'
    })
    
    # ============================================================================
    # CRITICAL ISSUE 4: METHODOLOGY CLARITY
    # ============================================================================
    review.append({
        'category': 'CRITICAL',
        'issue': 'Methodology Details Insufficient',
        'severity': 'Major',
        'description': """
        Several methodological details are missing or unclear:
        
        Issues:
        1. Feature engineering: How exactly are signature scores calculated?
           - Formula not provided
           - Gene selection criteria unclear
           - Weighting scheme not described
        2. Model selection: Why Logistic Regression chosen over others?
           - Need formal model comparison with statistical tests
           - Hyperparameter tuning process not detailed
        3. Risk stratification thresholds: How were 0.33 and 0.67 chosen?
           - Need clinical justification or data-driven optimization
        4. Missing data handling: Not mentioned
        5. Batch effects: Not addressed (if samples from different batches)
        6. Quality control: Criteria for sample exclusion not detailed
        
        Required Changes:
        - Provide exact formulas for all signature scores
        - Detail hyperparameter search space and selection criteria
        - Justify all thresholds with clinical or statistical rationale
        - Describe missing data handling
        - Address potential batch effects
        - Provide detailed QC criteria
        """,
        'recommendation': 'MAJOR REVISION REQUIRED - Methodology needs complete description'
    })
    
    # ============================================================================
    # CRITICAL ISSUE 5: COMPARISON WITH EXISTING METHODS
    # ============================================================================
    review.append({
        'category': 'CRITICAL',
        'issue': 'Insufficient Comparison with Existing Methods',
        'severity': 'Major',
        'description': """
        The manuscript claims superiority over PHASES score but comparison is inadequate.
        
        Issues:
        1. PHASES score not actually calculated on this dataset
        2. Comparison is to literature-reported PHASES performance (AUC 0.65-0.72)
        3. No direct head-to-head comparison on same patients
        4. Other existing models not compared (e.g., other ML models from literature)
        5. No statistical test for superiority
        
        Required Changes:
        - Calculate PHASES score for all patients in dataset
        - Perform direct comparison on same patients
        - Use statistical test (e.g., DeLong test) to compare ROC curves
        - Compare with other published models if applicable
        - Report if improvement is statistically significant
        """,
        'recommendation': 'MAJOR REVISION REQUIRED - Direct comparison needed'
    })
    
    # ============================================================================
    # CRITICAL ISSUE 6: BIOLOGICAL VALIDATION
    # ============================================================================
    review.append({
        'category': 'CRITICAL',
        'issue': 'Lack of Experimental Validation',
        'severity': 'Major',
        'description': """
        No experimental validation of findings, which is expected for Nature.
        
        Issues:
        1. All analyses are computational - no wet lab validation
        2. Biomarker candidates not validated experimentally
        3. No functional studies to confirm mechanisms
        4. No animal model validation
        5. No independent cohort validation with qPCR/sequencing
        
        Required Changes:
        - Validate top biomarkers in independent samples using qPCR or targeted sequencing
        - Perform functional studies on key genes (e.g., knockdown/overexpression)
        - If experimental validation not possible, clearly state this as limitation
        - Consider this a computational discovery study, not a validated biomarker study
        """,
        'recommendation': 'MAJOR REVISION REQUIRED - Experimental validation needed or clearly state limitation'
    })
    
    # ============================================================================
    # CRITICAL ISSUE 7: REPRODUCIBILITY
    # ============================================================================
    review.append({
        'category': 'CRITICAL',
        'issue': 'Reproducibility Concerns',
        'severity': 'Major',
        'description': """
        Reproducibility information is insufficient for Nature standards.
        
        Issues:
        1. Code availability: Not mentioned
        2. Data availability: GSE122897 is public, but processed data not provided
        3. Random seeds: Not all mentioned (some analyses use random_state=42, others unclear)
        4. Software versions: Not specified
        5. Exact parameters: Some model parameters not fully specified
        
        Required Changes:
        - Provide code repository (GitHub) with all analysis scripts
        - Provide processed data files (signature scores, model predictions)
        - Specify all random seeds used
        - Provide software versions and environment files (requirements.txt, conda env)
        - Include detailed README for reproducibility
        """,
        'recommendation': 'MAJOR REVISION REQUIRED - Reproducibility must be ensured'
    })
    
    # ============================================================================
    # CRITICAL ISSUE 8: CLINICAL TRANSLATION
    # ============================================================================
    review.append({
        'category': 'CRITICAL',
        'issue': 'Clinical Translation Claims Unsupported',
        'severity': 'Major',
        'description': """
        Claims about clinical utility are premature given the study design.
        
        Issues:
        1. No prospective validation
        2. No cost-effectiveness analysis
        3. No comparison with current clinical workflow
        4. No assessment of clinical decision impact
        5. Risk thresholds not clinically validated
        6. No discussion of implementation challenges
        
        Required Changes:
        - Tone down clinical translation claims
        - Clearly state this is a discovery/development study
        - Discuss what would be needed for clinical implementation
        - Acknowledge that clinical validation is future work
        - Provide realistic timeline for translation
        """,
        'recommendation': 'MAJOR REVISION REQUIRED - Clinical claims need moderation'
    })
    
    # ============================================================================
    # MODERATE ISSUES
    # ============================================================================
    
    review.append({
        'category': 'MODERATE',
        'issue': 'Figure Quality and Presentation',
        'severity': 'Moderate',
        'description': """
        Some figures could be improved:
        
        Issues:
        1. Some figures may be too complex (information overload)
        2. Color schemes may not be colorblind-friendly
        3. Figure legends could be more detailed
        4. Some figures lack statistical annotations (e.g., p-values, effect sizes)
        
        Required Changes:
        - Simplify complex figures or split into panels
        - Use colorblind-friendly palettes
        - Add statistical annotations to all relevant figures
        - Ensure all figures are publication-quality (300+ DPI)
        """,
        'recommendation': 'MODERATE REVISION - Improve figure presentation'
    })
    
    review.append({
        'category': 'MODERATE',
        'issue': 'Literature Context',
        'severity': 'Moderate',
        'description': """
        Literature review and context could be strengthened:
        
        Issues:
        1. Limited discussion of previous ML models for aneurysm rupture
        2. Incomplete comparison with other omics studies
        3. Missing discussion of why this approach is novel
        4. Limited mechanistic discussion
        
        Required Changes:
        - Comprehensive literature review of existing prediction models
        - Clear articulation of novelty
        - Better mechanistic discussion linking findings to biology
        - Discussion of how this fits into existing knowledge
        """,
        'recommendation': 'MODERATE REVISION - Strengthen literature context'
    })
    
    # ============================================================================
    # MINOR ISSUES
    # ============================================================================
    
    review.append({
        'category': 'MINOR',
        'issue': 'Writing and Presentation',
        'severity': 'Minor',
        'description': """
        Minor writing improvements needed:
        
        Issues:
        1. Some technical terms could be better explained
        2. Abbreviations not always defined on first use
        3. Some sentences could be more concise
        4. Consistency in terminology
        
        Required Changes:
        - Define all abbreviations on first use
        - Simplify technical language where possible
        - Ensure consistent terminology throughout
        - Professional copy-editing recommended
        """,
        'recommendation': 'MINOR REVISION - Improve writing clarity'
    })
    
    # ============================================================================
    # FURTHER ADVANCEMENTS NEEDED
    # ============================================================================
    
    advancements = [
        {
            'area': 'Data Integration',
            'suggestion': """
            - Integrate additional omics layers (genomics, proteomics, metabolomics)
            - Include imaging data (aneurysm morphology, hemodynamics)
            - Integrate clinical variables more systematically
            - Multi-center data integration
            """
        },
        {
            'area': 'Model Development',
            'suggestion': """
            - Develop ensemble methods combining multiple algorithms
            - Use deep learning for non-linear pattern discovery
            - Implement transfer learning from related diseases
            - Develop patient-specific models (personalized medicine)
            """
        },
        {
            'area': 'Validation',
            'suggestion': """
            - Multi-center prospective validation study
            - Longitudinal follow-up validation
            - Cross-platform validation (different sequencing technologies)
            - Validation in diverse populations (ethnicity, geography)
            """
        },
        {
            'area': 'Mechanistic Understanding',
            'suggestion': """
            - Functional genomics studies (CRISPR, RNAi)
            - Animal model validation
            - Single-cell resolution studies
            - Pathway perturbation experiments
            """
        },
        {
            'area': 'Clinical Translation',
            'suggestion': """
            - Develop clinical decision support tool
            - Cost-effectiveness analysis
            - Implementation study
            - Clinical trial design
            """
        }
    ]
    
    # Print review
    print("\n" + "=" * 80)
    print("CRITICAL ISSUES (Must Address for Publication)")
    print("=" * 80)
    
    for i, item in enumerate([r for r in review if r['category'] == 'CRITICAL'], 1):
        print(f"\n{i}. {item['issue']} ({item['severity']})")
        print("-" * 80)
        print(item['description'])
        print(f"\nRecommendation: {item['recommendation']}")
    
    print("\n" + "=" * 80)
    print("MODERATE ISSUES (Should Address)")
    print("=" * 80)
    
    for i, item in enumerate([r for r in review if r['category'] == 'MODERATE'], 1):
        print(f"\n{i}. {item['issue']} ({item['severity']})")
        print("-" * 80)
        print(item['description'])
        print(f"\nRecommendation: {item['recommendation']}")
    
    print("\n" + "=" * 80)
    print("MINOR ISSUES (Nice to Address)")
    print("=" * 80)
    
    for i, item in enumerate([r for r in review if r['category'] == 'MINOR'], 1):
        print(f"\n{i}. {item['issue']} ({item['severity']})")
        print("-" * 80)
        print(item['description'])
        print(f"\nRecommendation: {item['recommendation']}")
    
    print("\n" + "=" * 80)
    print("FURTHER ADVANCEMENTS NEEDED FOR NATURE PUBLICATION")
    print("=" * 80)
    
    for i, adv in enumerate(advancements, 1):
        print(f"\n{i}. {adv['area']}")
        print("-" * 80)
        print(adv['suggestion'])
    
    print("\n" + "=" * 80)
    print("OVERALL ASSESSMENT")
    print("=" * 80)
    print("""
    Current Status: MAJOR REVISION REQUIRED
    
    The manuscript presents interesting findings but has several critical issues that 
    must be addressed before it can be considered for Nature publication:
    
    1. Sample size is insufficient (n=43) - this is a major limitation
    2. No independent external validation - essential for Nature
    3. Statistical rigor needs improvement (CIs, permutation testing, nested CV)
    4. Methodology details incomplete (formulas, thresholds, QC)
    5. No experimental validation of findings
    6. Reproducibility concerns (code, data, parameters)
    7. Clinical translation claims premature
    
    Recommendation:
    - This work is better suited for a specialized journal (e.g., Stroke, 
      Journal of Neurointerventional Surgery) as a computational discovery study
    - For Nature consideration, would need:
      * Larger sample size (n > 200 minimum)
      * Independent external validation
      * Experimental validation of key findings
      * Prospective clinical validation
      * Complete reproducibility package
    
    The work shows promise but needs substantial additional validation and 
    experimental work before reaching Nature standards.
    """)
    
    print("\n" + "=" * 80)
    print("REVIEW COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    conduct_nature_review()

