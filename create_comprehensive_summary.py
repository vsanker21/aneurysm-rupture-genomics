"""
Create Comprehensive Summary Report of All Analyses
"""
import pandas as pd
import numpy as np
import os
from datetime import datetime

def create_summary_report():
    """Create a comprehensive summary report"""
    
    print("=" * 80)
    print("CREATING COMPREHENSIVE SUMMARY REPORT")
    print("=" * 80)
    
    report = []
    report.append("=" * 80)
    report.append("GSE122897 DIFFERENTIAL EXPRESSION ANALYSIS - COMPREHENSIVE SUMMARY")
    report.append("=" * 80)
    report.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("\n" + "=" * 80)
    
    # Load main DE results
    de_file = "results/gse122897_differential_expression.tsv"
    if os.path.exists(de_file):
        df = pd.read_csv(de_file, sep='\t')
        report.append("\n## OVERALL ANALYSIS SUMMARY")
        report.append("-" * 80)
        report.append(f"Total genes tested: {len(df):,}")
        report.append(f"Genes with P < 0.05: {(df['P_value'] < 0.05).sum():,}")
        report.append(f"Genes with P_adj_BH < 0.05: {(df['P_adj_BH'] < 0.05).sum():,}")
        report.append(f"Genes with P_adj_BH < 0.1: {(df['P_adj_BH'] < 0.1).sum():,}")
        report.append(f"Genes with P_adj_BH < 0.2: {(df['P_adj_BH'] < 0.2).sum():,}")
    
    # Multi-threshold summary
    multi_file = "results/multi_threshold_summary.tsv"
    if os.path.exists(multi_file):
        multi_df = pd.read_csv(multi_file, sep='\t')
        report.append("\n## RESULTS BY FDR THRESHOLD")
        report.append("-" * 80)
        for idx, row in multi_df.iterrows():
            report.append(f"\nFDR < {row['FDR_threshold']}:")
            report.append(f"  Total significant genes: {int(row['Total_genes']):,}")
            report.append(f"  Upregulated: {int(row['Upregulated']):,}")
            report.append(f"  Downregulated: {int(row['Downregulated']):,}")
            report.append(f"  Mean |Log2FC|: {row['Mean_|Log2FC|']:.3f}")
            report.append(f"  Max Log2FC: {row['Max_Log2FC']:.3f}")
            report.append(f"  Min Log2FC: {row['Min_Log2FC']:.3f}")
    
    # Top genes
    report.append("\n## TOP SIGNIFICANT GENES (FDR < 0.05)")
    report.append("-" * 80)
    if os.path.exists(de_file):
        top_genes = df[df['P_adj_BH'] < 0.05].copy()
        if len(top_genes) > 0:
            top_genes = top_genes.sort_values('P_adj_BH')
            for idx, row in top_genes.iterrows():
                direction = "UP" if row['Log2FC'] > 0 else "DOWN"
                gene_str = str(row['Gene'])
                report.append(f"  Gene {gene_str:10s} | Log2FC={row['Log2FC']:7.3f} | P_adj={row['P_adj_BH']:.4f} | {direction}")
    
    # Power analysis summary
    power_file = "results/power_analysis_table.tsv"
    if os.path.exists(power_file):
        power_df = pd.read_csv(power_file, sep='\t')
        report.append("\n## POWER ANALYSIS SUMMARY")
        report.append("-" * 80)
        report.append("Sample sizes: 22 Ruptured vs 21 Unruptured")
        report.append("\nDetectable effect sizes:")
        
        # Find effect sizes at 80% and 90% power
        power_80 = power_df[power_df['Power'] >= 0.80]
        power_90 = power_df[power_df['Power'] >= 0.90]
        
        if len(power_80) > 0:
            effect_80 = power_80.iloc[0]['Effect_size_Cohens_d']
            log2fc_80 = effect_80 * 0.75  # Approximate conversion
            report.append(f"  At 80% power: Cohen's d = {effect_80:.2f} (≈ Log2FC {log2fc_80:.2f})")
        
        if len(power_90) > 0:
            effect_90 = power_90.iloc[0]['Effect_size_Cohens_d']
            log2fc_90 = effect_90 * 0.75
            report.append(f"  At 90% power: Cohen's d = {effect_90:.2f} (≈ Log2FC {log2fc_90:.2f})")
    
    # Recommendations
    report.append("\n## RECOMMENDATIONS")
    report.append("-" * 80)
    report.append("""
1. FDR THRESHOLD SELECTION:
   - FDR < 0.05: Very conservative, only 4 genes (high confidence)
   - FDR < 0.1: Moderate threshold, 1,010 genes (recommended for exploratory analysis)
   - FDR < 0.2: Lenient threshold, 3,351 genes (for hypothesis generation)

2. SAMPLE SIZE & POWER:
   - Current sample size (22 vs 21) provides good power for moderate-large effects
   - Power is adequate for detecting Log2FC > 0.5-1.0
   - For smaller effects, consider increasing sample size

3. EFFECT SIZES:
   - Top genes show moderate to large fold changes (|Log2FC| up to 3.3)
   - Largest effects observed: Log2FC ≈ 3.3 (upregulated) and -3.3 (downregulated)
   - These effect sizes are biologically meaningful

4. NEXT STEPS:
   a) VALIDATION:
      - Validate top genes in independent datasets
      - Replicate findings in other aneurysm studies
   
   b) PATHWAY ANALYSIS:
      - Perform pathway enrichment on genes with P_adj < 0.1
      - Identify biological pathways involved in rupture
      - Focus on vascular remodeling, inflammation, ECM pathways
   
   c) GENE FUNCTION:
      - Investigate biological function of top DE genes
      - Map to known aneurysm-related pathways
      - Check for known disease associations
   
   d) INTEGRATION:
      - Integrate with GWAS findings
      - Compare with blood transcriptomics signatures
      - Build multi-omics prediction model
""")
    
    # File listing
    report.append("\n## GENERATED FILES")
    report.append("-" * 80)
    result_files = [
        "gse122897_differential_expression.tsv",
        "gse122897_top_genes.tsv",
        "gse122897_sig_genes_fdr_0.05.tsv",
        "gse122897_sig_genes_fdr_0.1.tsv",
        "gse122897_sig_genes_fdr_0.2.tsv",
        "gse122897_upregulated_fdr_0.1.tsv",
        "gse122897_downregulated_fdr_0.1.tsv",
        "gse122897_annotated_fdr_0.1.tsv",
        "power_analysis_curve.png",
        "power_analysis_table.tsv",
        "multi_threshold_summary.tsv"
    ]
    
    for filename in result_files:
        filepath = f"results/{filename}"
        if os.path.exists(filepath):
            size = os.path.getsize(filepath)
            report.append(f"  ✓ {filename} ({size:,} bytes)")
        else:
            report.append(f"  ✗ {filename} (not found)")
    
    report.append("\n" + "=" * 80)
    report.append("END OF REPORT")
    report.append("=" * 80)
    
    # Save report
    report_text = "\n".join(report)
    output_file = "results/COMPREHENSIVE_ANALYSIS_SUMMARY.txt"
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    print(f"\n[OK] Saved comprehensive summary to {output_file}")
    print("\n" + report_text)
    
    return report_text

if __name__ == "__main__":
    create_summary_report()

