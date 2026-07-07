"""
Comprehensive Analysis Pipeline
Runs all advanced analyses, ML models, multi-omics integration, and generates reports
"""
import os
import sys
import subprocess
import pandas as pd
import numpy as np
import time
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

def run_script(script_name, description, timeout=1800):
    """Run a Python script and report results"""
    print("\n" + "=" * 80)
    print(f"RUNNING: {description}")
    print("=" * 80)
    
    if not os.path.exists(script_name):
        print(f"[WARNING] Script not found: {script_name}")
        return False
    
    try:
        result = subprocess.run(
            [sys.executable, script_name],
            capture_output=True,
            text=True,
            timeout=timeout
        )
        
        if result.returncode == 0:
            print(f"[OK] {description} completed successfully")
            if result.stdout:
                # Print last 500 chars of output
                output_lines = result.stdout.split('\n')
                if len(output_lines) > 20:
                    print("\n... (output truncated) ...")
                    print("\n".join(output_lines[-20:]))
                else:
                    print(result.stdout)
            return True
        else:
            print(f"[ERROR] {description} failed with return code {result.returncode}")
            if result.stderr:
                print("Error output (last 500 chars):")
                print(result.stderr[-500:])
            return False
            
    except subprocess.TimeoutExpired:
        print(f"[ERROR] {description} timed out after {timeout} seconds")
        return False
    except Exception as e:
        print(f"[ERROR] {description} failed: {e}")
        return False

def check_prerequisites():
    """Check if required data files exist"""
    print("\n" + "=" * 80)
    print("CHECKING PREREQUISITES")
    print("=" * 80)
    
    required_files = [
        "results/gse122897_differential_expression.tsv",
        "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    ]
    
    missing = []
    for file in required_files:
        if os.path.exists(file):
            print(f"  [OK] {file}")
        else:
            print(f"  [MISSING] {file}")
            missing.append(file)
    
    if missing:
        print(f"\n[WARNING] {len(missing)} required file(s) missing")
        print("Some analyses may fail or use fallback methods")
    
    return len(missing) == 0

def run_advanced_analyses():
    """Run advanced analyses: gene ID conversion, pathway enrichment, literature annotation"""
    print("\n" + "=" * 80)
    print("PHASE 1: ADVANCED ANALYSES")
    print("=" * 80)
    
    results = {}
    
    # 1. Gene ID Conversion
    if os.path.exists("convert_and_integrate_gene_ids.py"):
        results['Gene_ID_Conversion'] = run_script(
            "convert_and_integrate_gene_ids.py",
            "Gene ID Conversion and Integration",
            timeout=1200
        )
    else:
        print("[WARNING] convert_and_integrate_gene_ids.py not found")
        results['Gene_ID_Conversion'] = False
    
    time.sleep(2)
    
    # 2. Pathway Enrichment
    if os.path.exists("pathway_enrichment_analysis.py"):
        results['Pathway_Enrichment'] = run_script(
            "pathway_enrichment_analysis.py",
            "Pathway Enrichment Analysis",
            timeout=1800
        )
    else:
        print("[WARNING] pathway_enrichment_analysis.py not found")
        results['Pathway_Enrichment'] = False
    
    time.sleep(2)
    
    # 3. Literature Annotation
    if os.path.exists("annotate_top_genes_literature.py"):
        results['Literature_Annotation'] = run_script(
            "annotate_top_genes_literature.py",
            "Literature Annotation",
            timeout=1800
        )
    else:
        print("[WARNING] annotate_top_genes_literature.py not found")
        results['Literature_Annotation'] = False
    
    return results

def run_ml_classification():
    """Run machine learning classification models"""
    print("\n" + "=" * 80)
    print("PHASE 2: MACHINE LEARNING CLASSIFICATION")
    print("=" * 80)
    
    results = {}
    
    if os.path.exists("ml_classification_blood_genes.py"):
        results['ML_Classification'] = run_script(
            "ml_classification_blood_genes.py",
            "Machine Learning Classification Models",
            timeout=1200
        )
    else:
        print("[WARNING] ml_classification_blood_genes.py not found")
        results['ML_Classification'] = False
    
    return results

def run_endothelial_dysfunction():
    """Run endothelial dysfunction analysis"""
    print("\n" + "=" * 80)
    print("PHASE 3: ENDOTHELIAL DYSFUNCTION ANALYSIS")
    print("=" * 80)
    
    results = {}
    
    if os.path.exists("endothelial_dysfunction_analysis.py"):
        results['EC_Dysfunction'] = run_script(
            "endothelial_dysfunction_analysis.py",
            "Endothelial Dysfunction Analysis",
            timeout=1200
        )
    else:
        print("[WARNING] endothelial_dysfunction_analysis.py not found")
        results['EC_Dysfunction'] = False
    
    return results

def run_multi_omics_integration():
    """Run multi-omics integration with ML models"""
    print("\n" + "=" * 80)
    print("PHASE 4: MULTI-OMICS INTEGRATION")
    print("=" * 80)
    
    results = {}
    
    if os.path.exists("multi_omics_integration.py"):
        results['Multi_Omics_Integration'] = run_script(
            "multi_omics_integration.py",
            "Multi-Omics Integration with ML Models",
            timeout=1800
        )
    else:
        print("[WARNING] multi_omics_integration.py not found")
        results['Multi_Omics_Integration'] = False
    
    return results

def run_validation_analyses():
    """Run validation and comparison analyses"""
    print("\n" + "=" * 80)
    print("PHASE 5: VALIDATION AND COMPARISON ANALYSES")
    print("=" * 80)
    
    results = {}
    
    # Run validation pipeline if available
    if os.path.exists("run_validation_pipeline.py"):
        results['Validation_Pipeline'] = run_script(
            "run_validation_pipeline.py",
            "Validation Pipeline",
            timeout=1800
        )
    else:
        print("[WARNING] run_validation_pipeline.py not found")
        results['Validation_Pipeline'] = False
    
    # Additional validation scripts
    validation_scripts = [
        ("validate_de_genes.py", "Gene Validation"),
        ("review_biomarker_candidates.py", "Biomarker Candidate Review")
    ]
    
    for script, description in validation_scripts:
        if os.path.exists(script):
            results[description] = run_script(script, description, timeout=600)
            time.sleep(2)
        else:
            print(f"[WARNING] {script} not found")
            results[description] = False
    
    return results

def generate_comprehensive_report(all_results):
    """Generate comprehensive report summarizing all analyses"""
    print("\n" + "=" * 80)
    print("GENERATING COMPREHENSIVE REPORT")
    print("=" * 80)
    
    os.makedirs("results", exist_ok=True)
    
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("COMPREHENSIVE ANALYSIS REPORT")
    report_lines.append("=" * 80)
    report_lines.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append("\n" + "=" * 80)
    
    # Summary of all analyses
    report_lines.append("\n## ANALYSIS SUMMARY")
    report_lines.append("=" * 80)
    
    total_analyses = len(all_results)
    successful = sum(1 for r in all_results.values() if r)
    failed = total_analyses - successful
    
    report_lines.append(f"\nTotal Analyses: {total_analyses}")
    report_lines.append(f"Successful: {successful}")
    report_lines.append(f"Failed: {failed}")
    report_lines.append(f"Success Rate: {successful/total_analyses*100:.1f}%")
    
    # Detailed results by phase
    report_lines.append("\n\n## DETAILED RESULTS BY PHASE")
    report_lines.append("=" * 80)
    
    phases = {
        "Phase 1: Advanced Analyses": ["Gene_ID_Conversion", "Pathway_Enrichment", "Literature_Annotation"],
        "Phase 2: ML Classification": ["ML_Classification"],
        "Phase 3: Endothelial Dysfunction": ["EC_Dysfunction"],
        "Phase 4: Multi-Omics Integration": ["Multi_Omics_Integration"],
        "Phase 5: Validation": ["Validation_Pipeline", "Gene Validation", "Biomarker Candidate Review"]
    }
    
    for phase_name, analysis_keys in phases.items():
        report_lines.append(f"\n### {phase_name}")
        report_lines.append("-" * 80)
        for key in analysis_keys:
            if key in all_results:
                status = "✓ COMPLETE" if all_results[key] else "✗ FAILED"
                report_lines.append(f"  {status}: {key}")
    
    # Key output files
    report_lines.append("\n\n## KEY OUTPUT FILES")
    report_lines.append("=" * 80)
    
    output_files = {
        "Gene ID Conversion": [
            "results/gene_id_mapping_unified.tsv",
            "results/integrated_cross_dataset_genes.tsv"
        ],
        "Pathway Enrichment": [
            "results/pathway_enrichment_fdr_0.05.tsv",
            "results/pathway_enrichment_summary.tsv"
        ],
        "Literature Annotation": [
            "results/gse159610_top_genes_literature_annotated.tsv",
            "results/gse159610_genes_with_pubmed_results.tsv"
        ],
        "ML Classification": [
            "results/blood_ml_feature_importance.tsv",
            "results/blood_ml_cv_results.tsv",
            "results/blood_ml_summary.tsv"
        ],
        "EC Dysfunction": [
            "results/ec_dysfunction_scores.tsv",
            "results/ec_dysfunction_markers.tsv"
        ],
        "Multi-Omics Integration": [
            "results/integrated_features.tsv",
            "results/risk_stratification.tsv",
            "results/feature_importance.tsv"
        ]
    }
    
    for category, files in output_files.items():
        report_lines.append(f"\n### {category}")
        for file in files:
            if os.path.exists(file):
                size = os.path.getsize(file) / 1024  # KB
                report_lines.append(f"  [OK] {file} ({size:.1f} KB)")
            else:
                report_lines.append(f"  [MISSING] {file}")
    
    # Statistics summary
    report_lines.append("\n\n## STATISTICS SUMMARY")
    report_lines.append("=" * 80)
    
    # Try to load and summarize key results
    try:
        # DE results
        de_file = "results/gse122897_differential_expression.tsv"
        if os.path.exists(de_file):
            de_df = pd.read_csv(de_file, sep='\t')
            sig_genes = (de_df['P_adj_BH'] < 0.05).sum() if 'P_adj_BH' in de_df.columns else 0
            report_lines.append(f"\nDifferential Expression:")
            report_lines.append(f"  Total genes tested: {len(de_df):,}")
            report_lines.append(f"  Significant genes (P_adj < 0.05): {sig_genes:,}")
        
        # Pathway enrichment
        pathway_file = "results/pathway_enrichment_fdr_0.05.tsv"
        if os.path.exists(pathway_file):
            pathway_df = pd.read_csv(pathway_file, sep='\t')
            report_lines.append(f"\nPathway Enrichment:")
            report_lines.append(f"  Enriched pathways (FDR < 0.05): {len(pathway_df):,}")
            if len(pathway_df) > 0:
                top_pathway = pathway_df.iloc[0]['Term'] if 'Term' in pathway_df.columns else "N/A"
                report_lines.append(f"  Top pathway: {top_pathway}")
        
        # ML results
        ml_file = "results/blood_ml_cv_results.tsv"
        if os.path.exists(ml_file):
            ml_df = pd.read_csv(ml_file, sep='\t')
            if 'Mean_AUC' in ml_df.columns:
                best_auc = ml_df['Mean_AUC'].max()
                best_model = ml_df.loc[ml_df['Mean_AUC'].idxmax(), 'Model'] if 'Model' in ml_df.columns else "N/A"
                report_lines.append(f"\nMachine Learning:")
                report_lines.append(f"  Best model: {best_model}")
                report_lines.append(f"  Best CV AUC: {best_auc:.3f}")
        
        # EC dysfunction
        ec_file = "results/ec_dysfunction_scores.tsv"
        if os.path.exists(ec_file):
            ec_df = pd.read_csv(ec_file, sep='\t')
            report_lines.append(f"\nEndothelial Dysfunction:")
            report_lines.append(f"  Samples analyzed: {len(ec_df):,}")
            if 'EC_Dysfunction_Score' in ec_df.columns:
                mean_score = ec_df['EC_Dysfunction_Score'].mean()
                report_lines.append(f"  Mean dysfunction score: {mean_score:.3f}")
        
        # Multi-omics integration
        risk_file = "results/risk_stratification.tsv"
        if os.path.exists(risk_file):
            risk_df = pd.read_csv(risk_file, sep='\t')
            report_lines.append(f"\nRisk Stratification:")
            report_lines.append(f"  Samples stratified: {len(risk_df):,}")
            if 'Risk_Category' in risk_df.columns:
                risk_counts = risk_df['Risk_Category'].value_counts()
                for category, count in risk_counts.items():
                    report_lines.append(f"    {category}: {count}")
    
    except Exception as e:
        report_lines.append(f"\n[WARNING] Could not generate statistics: {e}")
    
    # Next steps
    report_lines.append("\n\n## NEXT STEPS")
    report_lines.append("=" * 80)
    report_lines.append("\n1. Review pathway enrichment results for biological insights")
    report_lines.append("2. Examine literature annotations for top candidate genes")
    report_lines.append("3. Validate ML model performance on independent datasets")
    report_lines.append("4. Compare multi-omics model with clinical scores (PHASES)")
    report_lines.append("5. Prepare publication-quality figures and tables")
    
    # Write report
    report_text = "\n".join(report_lines)
    report_file = "results/COMPREHENSIVE_ANALYSIS_REPORT.txt"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    print(f"\n[OK] Comprehensive report saved to {report_file}")
    print("\n" + report_text)
    
    return report_file

def generate_summary_table(all_results):
    """Generate summary table of all results"""
    print("\n" + "=" * 80)
    print("GENERATING SUMMARY TABLE")
    print("=" * 80)
    
    summary_data = []
    
    for analysis_name, success in all_results.items():
        summary_data.append({
            'Analysis': analysis_name,
            'Status': 'Complete' if success else 'Failed',
            'Timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = "results/analysis_summary_table.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    
    print(f"\n[OK] Summary table saved to {summary_file}")
    print("\n" + summary_df.to_string(index=False))
    
    return summary_file

def main():
    """Main comprehensive analysis pipeline"""
    
    print("=" * 80)
    print("COMPREHENSIVE ANALYSIS PIPELINE")
    print("=" * 80)
    print("\nThis pipeline will run:")
    print("  1. Advanced analyses (gene ID conversion, pathway enrichment, literature annotation)")
    print("  2. Machine learning classification models")
    print("  3. Endothelial dysfunction analysis")
    print("  4. Multi-omics integration with ML models")
    print("  5. Validation and comparison analyses")
    print("  6. Comprehensive reports and summaries")
    print("\n[NOTE] This may take 1-3 hours depending on API rates and data size")
    
    start_time = time.time()
    
    # Check prerequisites
    check_prerequisites()
    
    # Ensure results directory exists
    os.makedirs("results", exist_ok=True)
    
    # Collect all results
    all_results = {}
    
    # Phase 1: Advanced Analyses
    phase1_results = run_advanced_analyses()
    all_results.update(phase1_results)
    
    time.sleep(5)
    
    # Phase 2: ML Classification
    phase2_results = run_ml_classification()
    all_results.update(phase2_results)
    
    time.sleep(5)
    
    # Phase 3: Endothelial Dysfunction
    phase3_results = run_endothelial_dysfunction()
    all_results.update(phase3_results)
    
    time.sleep(5)
    
    # Phase 4: Multi-Omics Integration
    phase4_results = run_multi_omics_integration()
    all_results.update(phase4_results)
    
    time.sleep(5)
    
    # Phase 5: Validation
    phase5_results = run_validation_analyses()
    all_results.update(phase5_results)
    
    # Generate reports
    report_file = generate_comprehensive_report(all_results)
    summary_file = generate_summary_table(all_results)
    
    # Final summary
    elapsed_time = time.time() - start_time
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    
    print("\n" + "=" * 80)
    print("PIPELINE COMPLETE")
    print("=" * 80)
    print(f"\nTotal time: {hours}h {minutes}m")
    print(f"\nSuccessful analyses: {sum(1 for r in all_results.values() if r)}/{len(all_results)}")
    print(f"\nReports generated:")
    print(f"  - {report_file}")
    print(f"  - {summary_file}")
    print("\n" + "=" * 80)

if __name__ == "__main__":
    main()






