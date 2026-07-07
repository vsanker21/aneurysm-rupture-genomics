"""
Generate Final Comprehensive Summary Document
Compiles all findings into a manuscript-ready summary
"""
import pandas as pd
import numpy as np
import os
from datetime import datetime

def load_all_results():
    """Load all available result files"""
    results = {}
    
    # Multi-omics integrated
    files_to_check = {
        'multi_omics_integrated': 'results/multi_omics_comprehensive_prioritized_genes.tsv',
        'cross_validated': 'results/integrated_cross_dataset_genes.tsv',
        'literature_annotated': 'results/multi_omics_prioritized_genes_literature_annotated.tsv',
        'blood_de': 'results/gse159610_top_candidate_genes.tsv',
        'tissue_de': 'results/gse122897_differential_expression.tsv',
        'gwas_integration': 'results/gwas_integration_report.tsv',
        'pathway_enrichment': 'results/pathway_enrichment_fdr_0.05.tsv'
    }
    
    for key, filepath in files_to_check.items():
        if os.path.exists(filepath):
            try:
                results[key] = pd.read_csv(filepath, sep='\t')
                print(f"Loaded {key}: {len(results[key])} records")
            except Exception as e:
                print(f"Error loading {key}: {e}")
    
    return results

def generate_executive_summary(results):
    """Generate executive summary section"""
    summary = []
    summary.append("=" * 80)
    summary.append("EXECUTIVE SUMMARY")
    summary.append("=" * 80)
    summary.append("")
    summary.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    summary.append("")
    
    # Dataset overview
    summary.append("DATASET OVERVIEW")
    summary.append("-" * 80)
    
    if 'blood_de' in results:
        blood_df = results['blood_de']
        summary.append(f"Blood Transcriptomics (GSE159610):")
        summary.append(f"  - Total prioritized genes: {len(blood_df):,}")
        if 'Log2FC' in blood_df.columns:
            large_fc = (np.abs(blood_df['Log2FC']) > 2.0).sum()
            summary.append(f"  - Genes with |Log2FC| > 2.0: {large_fc:,}")
            summary.append(f"  - Upregulated: {(blood_df['Log2FC'] > 0).sum():,}")
            summary.append(f"  - Downregulated: {(blood_df['Log2FC'] < 0).sum():,}")
    
    if 'tissue_de' in results:
        tissue_df = results['tissue_de']
        summary.append(f"\nTissue Transcriptomics (GSE122897):")
        summary.append(f"  - Total genes analyzed: {len(tissue_df):,}")
        if 'P_adj_BH' in tissue_df.columns:
            sig_genes = (tissue_df['P_adj_BH'] < 0.05).sum()
            summary.append(f"  - Significant genes (FDR < 0.05): {sig_genes:,}")
        if 'Log2FC' in tissue_df.columns:
            summary.append(f"  - Upregulated: {(tissue_df['Log2FC'] > 0).sum():,}")
            summary.append(f"  - Downregulated: {(tissue_df['Log2FC'] < 0).sum():,}")
    
    # Multi-omics integration
    if 'multi_omics_integrated' in results:
        integrated_df = results['multi_omics_integrated']
        summary.append(f"\nMulti-Omics Integration:")
        summary.append(f"  - Total unique genes: {len(integrated_df):,}")
        if 'In_Blood' in integrated_df.columns:
            in_blood = integrated_df['In_Blood'].sum()
            in_tissue = integrated_df['In_Tissue'].sum()
            in_both = ((integrated_df['In_Blood'] == True) & (integrated_df['In_Tissue'] == True)).sum()
            summary.append(f"  - In blood dataset: {in_blood:,}")
            summary.append(f"  - In tissue dataset: {in_tissue:,}")
            summary.append(f"  - Cross-validated (in both): {in_both:,}")
    
    # Cross-validation
    if 'cross_validated' in results:
        cross_df = results['cross_validated']
        summary.append(f"\nCross-Validated Genes:")
        summary.append(f"  - Genes found in both datasets: {len(cross_df):,}")
        if 'Integrated_Priority_Score' in cross_df.columns:
            top_score = cross_df['Integrated_Priority_Score'].max()
            summary.append(f"  - Maximum priority score: {top_score:.1f}")
    
    # Literature support
    if 'literature_annotated' in results:
        lit_df = results['literature_annotated']
        summary.append(f"\nLiterature Support:")
        if 'Has_Literature_Support' in lit_df.columns:
            with_lit = lit_df['Has_Literature_Support'].sum()
            summary.append(f"  - Genes with literature support: {with_lit:,}")
        if 'PubMed_Count' in lit_df.columns:
            total_refs = lit_df['PubMed_Count'].sum()
            summary.append(f"  - Total PubMed references: {total_refs:,}")
    
    summary.append("")
    return "\n".join(summary)

def generate_top_genes_section(results):
    """Generate top prioritized genes section"""
    summary = []
    summary.append("=" * 80)
    summary.append("TOP PRIORITIZED GENES")
    summary.append("=" * 80)
    summary.append("")
    
    # Use literature annotated if available, otherwise use integrated
    df = None
    priority_col = None
    
    if 'literature_annotated' in results:
        df = results['literature_annotated']
        priority_col = 'Integrated_Priority_Score' if 'Integrated_Priority_Score' in df.columns else None
        if priority_col is None:
            for col in ['Multi_Omics_Priority_Score', 'Priority_Score']:
                if col in df.columns:
                    priority_col = col
                    break
    elif 'multi_omics_integrated' in results:
        df = results['multi_omics_integrated']
        priority_col = 'Multi_Omics_Priority_Score'
    elif 'cross_validated' in results:
        df = results['cross_validated']
        priority_col = 'Integrated_Priority_Score'
    
    if df is None or priority_col is None:
        summary.append("[No prioritized genes data available]")
        return "\n".join(summary)
    
    # Get top 20 genes
    top_20 = df.nlargest(20, priority_col).copy()
    
    summary.append(f"Top 20 Prioritized Genes (sorted by {priority_col}):")
    summary.append("-" * 80)
    summary.append("")
    
    # Determine available columns
    has_symbol = 'Gene_Symbol' in top_20.columns or 'Annotated_Gene_Symbol' in top_20.columns
    has_name = 'Gene_Name' in top_20.columns or 'Annotated_Gene_Name' in top_20.columns
    has_lit = 'Has_Literature_Support' in top_20.columns
    
    for idx, (row_idx, row) in enumerate(top_20.iterrows(), 1):
        # Get gene identifier
        symbol = row.get('Annotated_Gene_Symbol', row.get('Gene_Symbol', row.get('Gene_ID', 'Unknown')))
        if pd.isna(symbol) or symbol == '':
            symbol = str(row.get('Gene_ID', f'Gene_{idx}'))
        
        # Get priority score
        score = row.get(priority_col, 0)
        
        # Build info string
        info_parts = [f"Score: {score:.1f}"]
        
        # Blood data
        if 'Blood_Log2FC' in row and pd.notna(row['Blood_Log2FC']):
            info_parts.append(f"Blood FC: {row['Blood_Log2FC']:.2f}")
        
        # Tissue data
        if 'Tissue_Log2FC' in row and pd.notna(row['Tissue_Log2FC']):
            info_parts.append(f"Tissue FC: {row['Tissue_Log2FC']:.2f}")
        elif 'Log2FC' in row and pd.notna(row['Log2FC']):
            info_parts.append(f"FC: {row['Log2FC']:.2f}")
        
        # Literature
        if has_lit and row.get('Has_Literature_Support', False):
            pubmed_count = row.get('PubMed_Count', 0)
            info_parts.append(f"Lit: {pubmed_count} refs")
        
        # Gene name
        if has_name:
            name = row.get('Annotated_Gene_Name', row.get('Gene_Name', ''))
            if pd.notna(name) and name != '':
                symbol = f"{symbol} ({name})"
        
        info_str = " | ".join(info_parts)
        summary.append(f"{idx:2d}. {str(symbol):25s} | {info_str}")
    
    summary.append("")
    return "\n".join(summary)

def generate_methodology_section():
    """Generate methodology section"""
    summary = []
    summary.append("=" * 80)
    summary.append("METHODOLOGY")
    summary.append("=" * 80)
    summary.append("")
    
    summary.append("1. DATASET ACQUISITION")
    summary.append("-" * 80)
    summary.append("  - Blood Transcriptomics: GSE159610 (aneurysm vs control)")
    summary.append("  - Tissue Transcriptomics: GSE122897 (ruptured vs unruptured)")
    summary.append("  - GWAS Integration: Bakker et al. 2020 summary statistics")
    summary.append("")
    
    summary.append("2. DIFFERENTIAL EXPRESSION ANALYSIS")
    summary.append("-" * 80)
    summary.append("  - Statistical tests: Student's t-test with multiple testing correction")
    summary.append("  - Multiple testing correction: Benjamini-Hochberg FDR")
    summary.append("  - Effect size: Log2 fold change")
    summary.append("  - Significance threshold: FDR < 0.05")
    summary.append("")
    
    summary.append("3. GENE ID CONVERSION")
    summary.append("-" * 80)
    summary.append("  - Tool: MyGene.info API")
    summary.append("  - Conversion: Ensembl ↔ Entrez ID ↔ Gene Symbol")
    summary.append("  - Purpose: Enable cross-dataset integration")
    summary.append("")
    
    summary.append("4. MULTI-OMICS INTEGRATION")
    summary.append("-" * 80)
    summary.append("  - Priority scoring: Weighted combination of")
    summary.append("    * Effect size (Log2FC)")
    summary.append("    * Statistical significance (P-value, FDR)")
    summary.append("    * Cross-validation (presence in multiple datasets)")
    summary.append("    * Literature support")
    summary.append("    * Pathway enrichment")
    summary.append("")
    
    summary.append("5. LITERATURE ANNOTATION")
    summary.append("-" * 80)
    summary.append("  - Database: PubMed (via NCBI E-utilities)")
    summary.append("  - Search strategy: Gene symbol AND aneurysm-related terms")
    summary.append("  - Gene information: MyGene.info API")
    summary.append("")
    
    return "\n".join(summary)

def generate_key_findings_section(results):
    """Generate key findings section"""
    summary = []
    summary.append("=" * 80)
    summary.append("KEY FINDINGS")
    summary.append("=" * 80)
    summary.append("")
    
    # Find cross-validated genes
    if 'cross_validated' in results:
        cross_df = results['cross_validated']
        cross_df = cross_df[(cross_df.get('In_Blood', False)) & (cross_df.get('In_Tissue', False))]
        
        if len(cross_df) > 0:
            summary.append("1. CROSS-VALIDATED GENES")
            summary.append("-" * 80)
            summary.append(f"  Found {len(cross_df)} genes differentially expressed in both blood and tissue.")
            summary.append("  These genes represent the most robust biomarkers.")
            summary.append("")
            
            top_cross = cross_df.nlargest(5, 'Integrated_Priority_Score' if 'Integrated_Priority_Score' in cross_df.columns else 'Multi_Omics_Priority_Score')
            summary.append("  Top cross-validated genes:")
            for idx, (_, row) in enumerate(top_cross.iterrows(), 1):
                symbol = row.get('Gene_Symbol', row.get('Gene_ID', 'Unknown'))
                summary.append(f"    {idx}. {symbol}")
            summary.append("")
    
    # Large effect sizes
    if 'blood_de' in results:
        blood_df = results['blood_de']
        large_effect = blood_df[np.abs(blood_df['Log2FC']) > 3.0]
        if len(large_effect) > 0:
            summary.append("2. LARGE EFFECT SIZES")
            summary.append("-" * 80)
            summary.append(f"  {len(large_effect)} genes show extreme fold changes (|Log2FC| > 3.0) in blood.")
            summary.append("  These may represent key drivers of aneurysm pathophysiology.")
            summary.append("")
    
    # Literature support
    if 'literature_annotated' in results:
        lit_df = results['literature_annotated']
        if 'Has_Literature_Support' in lit_df.columns:
            with_lit = lit_df[lit_df['Has_Literature_Support'] == True]
            if len(with_lit) > 0:
                summary.append("3. LITERATURE SUPPORT")
                summary.append("-" * 80)
                summary.append(f"  {len(with_lit)} prioritized genes have existing literature support.")
                summary.append("  This validates our findings with prior aneurysm research.")
                summary.append("")
    
    return "\n".join(summary)

def generate_output_files_section():
    """Generate list of output files"""
    summary = []
    summary.append("=" * 80)
    summary.append("OUTPUT FILES")
    summary.append("=" * 80)
    summary.append("")
    
    output_files = {
        'Prioritized Genes': [
            'results/multi_omics_comprehensive_prioritized_genes.tsv',
            'results/integrated_cross_dataset_genes.tsv',
            'results/multi_omics_prioritized_top50.tsv'
        ],
        'Literature Annotated': [
            'results/multi_omics_prioritized_genes_literature_annotated.tsv',
            'results/top_genes_with_literature_support.tsv'
        ],
        'Visualizations': [
            'results/figures/volcano_plot_blood.png',
            'results/figures/volcano_plot_tissue.png',
            'results/figures/venn_blood_tissue.png',
            'results/figures/priority_score_distribution.png'
        ],
        'Reports': [
            'results/multi_omics_prioritization_report.txt',
            'results/comprehensive_gene_prioritization_report.txt'
        ]
    }
    
    for category, files in output_files.items():
        summary.append(f"{category}:")
        for file in files:
            if os.path.exists(file):
                summary.append(f"  [OK] {file}")
            else:
                summary.append(f"  [ ] {file}")
        summary.append("")
    
    return "\n".join(summary)

def main():
    """Generate final comprehensive summary document"""
    print("=" * 80)
    print("GENERATING FINAL COMPREHENSIVE SUMMARY DOCUMENT")
    print("=" * 80)
    
    # Load all results
    print("\nLoading all result files...")
    results = load_all_results()
    
    # Generate sections
    print("\nGenerating document sections...")
    
    document = []
    
    # Title
    document.append("=" * 80)
    document.append("COMPREHENSIVE MULTI-OMICS ANALYSIS SUMMARY")
    document.append("Cerebral Aneurysm Genomics Project")
    document.append("=" * 80)
    document.append("")
    
    # Executive summary
    document.append(generate_executive_summary(results))
    document.append("")
    
    # Key findings
    document.append(generate_key_findings_section(results))
    document.append("")
    
    # Top genes
    document.append(generate_top_genes_section(results))
    document.append("")
    
    # Methodology
    document.append(generate_methodology_section())
    document.append("")
    
    # Output files
    document.append(generate_output_files_section())
    document.append("")
    
    # Final summary
    document.append("=" * 80)
    document.append("CONCLUSIONS")
    document.append("=" * 80)
    document.append("")
    document.append("This comprehensive multi-omics integration has identified:")
    document.append("  - Prioritized gene candidates for aneurysm rupture prediction")
    document.append("  - Cross-validated biomarkers present in both blood and tissue")
    document.append("  - Genes with literature support validating prior research")
    document.append("  - Pathway-enriched gene sets for biological interpretation")
    document.append("")
    document.append("These findings provide a foundation for:")
    document.append("  - Development of blood-based biomarkers")
    document.append("  - Understanding molecular mechanisms of rupture")
    document.append("  - Identification of therapeutic targets")
    document.append("  - Clinical risk stratification models")
    document.append("")
    document.append("=" * 80)
    
    # Save document
    os.makedirs("results", exist_ok=True)
    output_file = "results/FINAL_COMPREHENSIVE_SUMMARY.txt"
    
    document_text = "\n".join(document)
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(document_text)
    
    print(f"\n[OK] Saved final summary document to {output_file}")
    
    # Print summary
    print("\n" + document_text[:2000] + "\n..." if len(document_text) > 2000 else document_text)
    
    print("\n" + "=" * 80)
    print("FINAL SUMMARY DOCUMENT GENERATION COMPLETE!")
    print("=" * 80)

if __name__ == "__main__":
    main()






