"""
Generate Final Comprehensive Report
Summarizes all completed analyses and generates publication-ready summaries
"""
import pandas as pd
import numpy as np
import os
from datetime import datetime
import glob

def check_file_exists(filepath):
    """Check if file exists and return size"""
    if os.path.exists(filepath):
        size_kb = os.path.getsize(filepath) / 1024
        return True, size_kb
    return False, 0

def load_summary_statistics():
    """Load and summarize key statistics from results"""
    stats = {}
    
    # DE Results
    de_file = "results/gse122897_differential_expression.tsv"
    if os.path.exists(de_file):
        de_df = pd.read_csv(de_file, sep='\t')
        stats['DE'] = {
            'total_genes': len(de_df),
            'sig_genes_005': (de_df['P_adj_BH'] < 0.05).sum() if 'P_adj_BH' in de_df.columns else 0,
            'sig_genes_001': (de_df['P_adj_BH'] < 0.01).sum() if 'P_adj_BH' in de_df.columns else 0,
            'upregulated': (de_df['Log2FC'] > 0).sum() if 'Log2FC' in de_df.columns else 0,
            'downregulated': (de_df['Log2FC'] < 0).sum() if 'Log2FC' in de_df.columns else 0
        }
    
    # Pathway Enrichment
    pathway_files = glob.glob("results/pathway_enrichment_fdr_*.tsv")
    if pathway_files:
        pathway_df = pd.read_csv(pathway_files[0], sep='\t')
        stats['Pathway'] = {
            'enriched_pathways': len(pathway_df),
            'top_pathway': pathway_df.iloc[0]['Term'] if len(pathway_df) > 0 and 'Term' in pathway_df.columns else 'N/A'
        }
    
    # ML Results
    ml_cv_file = "results/blood_ml_cv_results.tsv"
    if os.path.exists(ml_cv_file):
        ml_df = pd.read_csv(ml_cv_file, sep='\t')
        if 'Mean_AUC' in ml_df.columns:
            best_idx = ml_df['Mean_AUC'].idxmax()
            stats['ML'] = {
                'best_model': ml_df.loc[best_idx, 'Model'] if 'Model' in ml_df.columns else 'N/A',
                'best_auc': ml_df.loc[best_idx, 'Mean_AUC'],
                'best_std': ml_df.loc[best_idx, 'Std_AUC'] if 'Std_AUC' in ml_df.columns else 0
            }
    
    # EC Dysfunction
    ec_file = "results/ec_dysfunction_scores.tsv"
    if os.path.exists(ec_file):
        ec_df = pd.read_csv(ec_file, sep='\t')
        stats['EC'] = {
            'samples': len(ec_df),
            'mean_score': ec_df['EC_Dysfunction_Score'].mean() if 'EC_Dysfunction_Score' in ec_df.columns else 0
        }
    
    # Multi-Omics Integration
    risk_file = "results/risk_stratification.tsv"
    if os.path.exists(risk_file):
        risk_df = pd.read_csv(risk_file, sep='\t')
        stats['MultiOmics'] = {
            'samples': len(risk_df),
            'risk_categories': risk_df['Risk_Category'].value_counts().to_dict() if 'Risk_Category' in risk_df.columns else {}
        }
    
    return stats

def generate_comprehensive_report():
    """Generate comprehensive analysis report"""
    
    print("=" * 80)
    print("GENERATING COMPREHENSIVE ANALYSIS REPORT")
    print("=" * 80)
    
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("COMPREHENSIVE ANALYSIS REPORT")
    report_lines.append("Cerebral Aneurysm Genomics Project")
    report_lines.append("=" * 80)
    report_lines.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append("\n" + "=" * 80)
    
    # Executive Summary
    report_lines.append("\n## EXECUTIVE SUMMARY")
    report_lines.append("=" * 80)
    report_lines.append("\nThis report summarizes comprehensive multi-omics analyses performed on")
    report_lines.append("cerebral aneurysm datasets, including:")
    report_lines.append("  - Differential expression analysis (GSE122897)")
    report_lines.append("  - Pathway enrichment analysis")
    report_lines.append("  - Machine learning classification models")
    report_lines.append("  - Endothelial dysfunction analysis")
    report_lines.append("  - Multi-omics integration")
    report_lines.append("  - Gene ID conversion and cross-dataset integration")
    report_lines.append("  - Literature annotation")
    
    # Statistics Summary
    stats = load_summary_statistics()
    
    report_lines.append("\n\n## KEY STATISTICS")
    report_lines.append("=" * 80)
    
    if 'DE' in stats:
        de = stats['DE']
        report_lines.append(f"\n### Differential Expression Analysis")
        report_lines.append(f"  Total genes tested: {de['total_genes']:,}")
        report_lines.append(f"  Significant genes (P_adj < 0.05): {de['sig_genes_005']:,}")
        report_lines.append(f"  Highly significant (P_adj < 0.01): {de['sig_genes_001']:,}")
        report_lines.append(f"  Upregulated: {de['upregulated']:,}")
        report_lines.append(f"  Downregulated: {de['downregulated']:,}")
    
    if 'Pathway' in stats:
        pathway = stats['Pathway']
        report_lines.append(f"\n### Pathway Enrichment")
        report_lines.append(f"  Enriched pathways (FDR < 0.05): {pathway['enriched_pathways']}")
        report_lines.append(f"  Top pathway: {pathway['top_pathway']}")
    
    if 'ML' in stats:
        ml = stats['ML']
        report_lines.append(f"\n### Machine Learning Classification")
        report_lines.append(f"  Best model: {ml['best_model']}")
        report_lines.append(f"  Cross-validation AUC: {ml['best_auc']:.3f} ± {ml['best_std']:.3f}")
    
    if 'EC' in stats:
        ec = stats['EC']
        report_lines.append(f"\n### Endothelial Dysfunction")
        report_lines.append(f"  Samples analyzed: {ec['samples']}")
        report_lines.append(f"  Mean dysfunction score: {ec['mean_score']:.3f}")
    
    if 'MultiOmics' in stats:
        multi = stats['MultiOmics']
        report_lines.append(f"\n### Multi-Omics Integration")
        report_lines.append(f"  Samples stratified: {multi['samples']}")
        if multi['risk_categories']:
            report_lines.append(f"  Risk categories:")
            for category, count in multi['risk_categories'].items():
                report_lines.append(f"    {category}: {count}")
    
    # Output Files Inventory
    report_lines.append("\n\n## OUTPUT FILES INVENTORY")
    report_lines.append("=" * 80)
    
    file_categories = {
        "Differential Expression": [
            "results/gse122897_differential_expression.tsv",
            "results/gse122897_top_genes.tsv",
            "results/gse122897_sig_genes_fdr_0.05.tsv"
        ],
        "Pathway Enrichment": [
            "results/pathway_enrichment_fdr_0.05.tsv",
            "results/pathway_enrichment_summary.tsv"
        ],
        "Machine Learning": [
            "results/blood_ml_cv_results.tsv",
            "results/blood_ml_feature_importance.tsv",
            "results/blood_ml_summary.tsv"
        ],
        "Endothelial Dysfunction": [
            "results/ec_dysfunction_scores.tsv",
            "results/ec_dysfunction_markers.tsv"
        ],
        "Multi-Omics Integration": [
            "results/integrated_features.tsv",
            "results/risk_stratification.tsv",
            "results/feature_importance.tsv"
        ],
        "Gene ID Conversion": [
            "results/gene_id_mapping_unified.tsv",
            "results/integrated_cross_dataset_genes.tsv"
        ],
        "Literature Annotation": [
            "results/gse159610_top_genes_literature_annotated.tsv"
        ]
    }
    
    for category, files in file_categories.items():
        report_lines.append(f"\n### {category}")
        for file in files:
            exists, size = check_file_exists(file)
            if exists:
                report_lines.append(f"  [OK] {file} ({size:.1f} KB)")
            else:
                report_lines.append(f"  [MISSING] {file}")
    
    # Key Findings
    report_lines.append("\n\n## KEY FINDINGS")
    report_lines.append("=" * 80)
    
    # Try to load top genes
    top_genes_file = "results/gse122897_top_genes.tsv"
    if os.path.exists(top_genes_file):
        top_genes = pd.read_csv(top_genes_file, sep='\t')
        if len(top_genes) > 0:
            report_lines.append(f"\n### Top Differentially Expressed Genes")
            report_lines.append(f"  Number of top genes identified: {len(top_genes)}")
            if 'Gene' in top_genes.columns:
                top_10 = top_genes.head(10)['Gene'].tolist()
                report_lines.append(f"  Top 10 genes: {', '.join([str(g) for g in top_10[:10]])}")
    
    # Pathway findings
    pathway_file = "results/pathway_enrichment_fdr_0.05.tsv"
    if os.path.exists(pathway_file):
        pathway_df = pd.read_csv(pathway_file, sep='\t')
        if len(pathway_df) > 0:
            report_lines.append(f"\n### Enriched Pathways")
            for idx, row in pathway_df.head(5).iterrows():
                term = row.get('Term', 'N/A')
                pval = row.get('P_value', 'N/A')
                source = row.get('Source', 'N/A')
                report_lines.append(f"  - {term} ({source}, P={pval:.4e})")
    
    # ML findings
    if 'ML' in stats:
        ml = stats['ML']
        report_lines.append(f"\n### Machine Learning Performance")
        report_lines.append(f"  Best performing model: {ml['best_model']}")
        report_lines.append(f"  Cross-validation AUC: {ml['best_auc']:.3f} ± {ml['best_std']:.3f}")
        report_lines.append(f"  This indicates {'excellent' if ml['best_auc'] > 0.8 else 'good' if ml['best_auc'] > 0.7 else 'moderate'} predictive performance")
    
    # Next Steps
    report_lines.append("\n\n## RECOMMENDED NEXT STEPS")
    report_lines.append("=" * 80)
    report_lines.append("\n1. Validate top candidate genes in independent datasets")
    report_lines.append("2. Perform functional validation of enriched pathways")
    report_lines.append("3. Compare multi-omics model with clinical scores (PHASES)")
    report_lines.append("4. Generate publication-quality figures and visualizations")
    report_lines.append("5. Prepare manuscript with integrated findings")
    report_lines.append("6. Consider external validation cohort if available")
    
    # Write report
    report_text = "\n".join(report_lines)
    report_file = "results/FINAL_COMPREHENSIVE_ANALYSIS_REPORT.txt"
    
    os.makedirs("results", exist_ok=True)
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    print(f"\n[OK] Comprehensive report saved to {report_file}")
    print(f"\nReport summary:")
    print(f"  - {len(report_lines)} lines")
    print(f"  - {len([l for l in report_lines if l.startswith('###')])} major sections")
    
    return report_file

def generate_validation_comparison_report():
    """Generate validation and comparison analysis report"""
    
    print("\n" + "=" * 80)
    print("GENERATING VALIDATION AND COMPARISON REPORT")
    print("=" * 80)
    
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("VALIDATION AND COMPARISON ANALYSIS REPORT")
    report_lines.append("=" * 80)
    report_lines.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Cross-dataset validation
    report_lines.append("\n## CROSS-DATASET VALIDATION")
    report_lines.append("=" * 80)
    
    integrated_file = "results/integrated_cross_dataset_genes.tsv"
    if os.path.exists(integrated_file):
        integrated_df = pd.read_csv(integrated_file, sep='\t')
        report_lines.append(f"\nCross-dataset integration results:")
        report_lines.append(f"  Total unique genes: {len(integrated_df)}")
        report_lines.append(f"  In blood only: {(integrated_df['In_Blood'] & ~integrated_df['In_Tissue']).sum()}")
        report_lines.append(f"  In tissue only: {(integrated_df['In_Tissue'] & ~integrated_df['In_Blood']).sum()}")
        report_lines.append(f"  In both datasets: {(integrated_df['In_Blood'] & integrated_df['In_Tissue']).sum()}")
        
        # Cross-validated genes
        cross_validated = integrated_df[(integrated_df['In_Blood'] & integrated_df['In_Tissue'])].head(10)
        if len(cross_validated) > 0:
            report_lines.append(f"\n  Top cross-validated genes:")
            for idx, row in cross_validated.iterrows():
                symbol = row.get('Gene_Symbol', 'N/A')
                blood_fc = row.get('Blood_Log2FC', 0)
                tissue_fc = row.get('Tissue_Log2FC', 0)
                report_lines.append(f"    {symbol}: Blood FC={blood_fc:.2f}, Tissue FC={tissue_fc:.2f}")
    
    # Model comparison
    report_lines.append("\n## MODEL COMPARISON")
    report_lines.append("=" * 80)
    
    ml_file = "results/blood_ml_cv_results.tsv"
    if os.path.exists(ml_file):
        ml_df = pd.read_csv(ml_file, sep='\t')
        report_lines.append(f"\nMachine Learning Model Performance:")
        for idx, row in ml_df.iterrows():
            model = row.get('Model', 'N/A')
            auc = row.get('Mean_AUC', 0)
            std = row.get('Std_AUC', 0)
            report_lines.append(f"  {model}: AUC = {auc:.3f} ± {std:.3f}")
    
    # Write report
    report_text = "\n".join(report_lines)
    report_file = "results/VALIDATION_COMPARISON_REPORT.txt"
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    print(f"\n[OK] Validation report saved to {report_file}")
    
    return report_file

def main():
    """Main report generation"""
    
    # Generate comprehensive report
    comprehensive_report = generate_comprehensive_report()
    
    # Generate validation report
    validation_report = generate_validation_comparison_report()
    
    print("\n" + "=" * 80)
    print("REPORT GENERATION COMPLETE")
    print("=" * 80)
    print(f"\nGenerated reports:")
    print(f"  1. {comprehensive_report}")
    print(f"  2. {validation_report}")
    print("\nAll reports are available in the 'results/' directory")

if __name__ == "__main__":
    main()






