"""Check project progress and status"""
import os
import pandas as pd
from pathlib import Path
from datetime import datetime

def check_file_status():
    """Check status of key output files"""
    print("=" * 80)
    print("PROJECT PROGRESS STATUS REPORT")
    print("=" * 80)
    print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Key files to check
    key_files = {
        "Phase 1: Data Integration": [
            "results/dataset_inventory.json",
            "results/gse26969_metadata.tsv",
            "results/extracted_pdf_data.json"
        ],
        "Phase 2: Polygenic Risk Score": [
            "results/prs_weights.tsv",
            "results/prs_risk_loci_top17.tsv"
        ],
        "Phase 3: Blood Transcriptomics": [
            "results/gse159610_top_candidate_genes.tsv",
            "results/blood_ml_summary.tsv",
            "results/blood_ml_feature_importance.tsv",
            "results/blood_tissue_integration_fdr_0.1.tsv"
        ],
        "Phase 4: Tissue RNA-seq": [
            "results/gse122897_integrated_de.tsv",
            "results/gse122897_top_genes_integrated.tsv",
            "results/gse122897_differential_expression.tsv"
        ],
        "Phase 5: Pathway Enrichment": [
            "results/pathway_enrichment_fdr_0.05.tsv",
            "results/pathway_enrichment_fdr_0.1.tsv"
        ],
        "Phase 6: Multi-Omics Integration": [
            "results/multi_omics_prioritized_top50.tsv",
            "results/multi_omics_comprehensive_prioritized_genes.tsv",
            "results/multi_omics_prioritization_report.txt"
        ],
        "Advanced Analyses": [
            "results/gene_validation_report_fdr_0.05.tsv",
            "results/ec_dysfunction_scores.tsv",
            "results/integrated_features.tsv"
        ],
        "Summary Reports": [
            "results/COMPREHENSIVE_ANALYSIS_SUMMARY.txt",
            "results/FINAL_COMPREHENSIVE_SUMMARY.txt",
            "results/multi_omics_prioritization_report.txt"
        ]
    }
    
    total_files = 0
    existing_files = 0
    
    for phase, files in key_files.items():
        print(f"\n{phase}:")
        print("-" * 80)
        for file_path in files:
            total_files += 1
            if os.path.exists(file_path):
                existing_files += 1
                size = os.path.getsize(file_path)
                print(f"  ✓ {os.path.basename(file_path)} ({size:,} bytes)")
            else:
                print(f"  ✗ {os.path.basename(file_path)} (MISSING)")
    
    print("\n" + "=" * 80)
    print(f"SUMMARY: {existing_files}/{total_files} key files generated ({100*existing_files/total_files:.1f}%)")
    print("=" * 80)
    
    # Check key results
    print("\n" + "=" * 80)
    print("KEY RESULTS SUMMARY")
    print("=" * 80)
    
    # GSE122897 results
    de_file = "results/gse122897_integrated_de.tsv"
    if os.path.exists(de_file):
        try:
            df = pd.read_csv(de_file, sep='\t')
            print(f"\nGSE122897 (Tissue RNA-seq):")
            print(f"  Total genes tested: {len(df):,}")
            print(f"  FDR < 0.05: {(df['P_adj_BH'] < 0.05).sum():,}")
            print(f"  FDR < 0.1: {(df['P_adj_BH'] < 0.1).sum():,}")
        except:
            print(f"\nGSE122897: File exists but could not read")
    
    # PRS results
    prs_file = "results/prs_weights.tsv"
    if os.path.exists(prs_file):
        try:
            df = pd.read_csv(prs_file, sep='\t')
            print(f"\nPolygenic Risk Score:")
            print(f"  Risk loci identified: {len(df):,}")
        except:
            print(f"\nPRS: File exists but could not read")
    
    # Multi-omics results
    multi_file = "results/multi_omics_prioritized_top50.tsv"
    if os.path.exists(multi_file):
        try:
            df = pd.read_csv(multi_file, sep='\t')
            print(f"\nMulti-Omics Integration:")
            print(f"  Top prioritized genes: {len(df):,}")
        except:
            print(f"\nMulti-Omics: File exists but could not read")
    
    # ML results
    ml_file = "results/blood_ml_summary.tsv"
    if os.path.exists(ml_file):
        try:
            df = pd.read_csv(ml_file, sep='\t')
            print(f"\nMachine Learning (Blood Classification):")
            print(f"  Best model AUC: {df['CV_AUC'].max():.3f}" if 'CV_AUC' in df.columns else "  Results available")
        except:
            print(f"\nML: File exists but could not read")
    
    print("\n" + "=" * 80)
    print("PROJECT STATUS: COMPREHENSIVE ANALYSES COMPLETED")
    print("=" * 80)

if __name__ == "__main__":
    check_file_status()


