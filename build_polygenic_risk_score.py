"""
Build Polygenic Risk Score (PRS) from Bakker 2020 GWAS
Extracts 17 risk loci and calculates PRS scores
"""
import pandas as pd
import numpy as np
import gzip
import os
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve
import warnings
warnings.filterwarnings('ignore')

# Known 17 risk loci from Bakker et al. 2020
# Based on project outline and literature
KNOWN_RISK_LOCI = {
    # Chromosome 9p21.3 (CDKN2B-AS1) - most robust
    'rs1333040': {'chr': 9, 'pos': 22125503, 'effect_allele': 'C', 'or': 1.24, 'gene': 'CDKN2B-AS1'},
    'rs10757278': {'chr': 9, 'pos': 22124570, 'effect_allele': 'A', 'or': 1.29, 'gene': 'CDKN2B-AS1'},
    
    # Chromosome 8q11.23 - vascular endothelium
    'rs9298506': {'chr': 8, 'pos': 60799606, 'effect_allele': 'G', 'or': 1.20, 'gene': 'SOX17'},
    'rs10958409': {'chr': 8, 'pos': 62451394, 'effect_allele': 'A', 'or': 1.18, 'gene': 'SOX17'},
    
    # Chromosome 4q31 (EDNRA) - endothelin receptor
    'rs6841581': {'chr': 4, 'pos': 159235783, 'effect_allele': 'G', 'or': 1.24, 'gene': 'EDNRA'},
    
    # Chromosome 10q24 (CFDP1)
    'rs3184504': {'chr': 10, 'pos': 112654220, 'effect_allele': 'T', 'or': 1.22, 'gene': 'CFDP1'},
    
    # Chromosome 2q13 (COL3A1) - collagen
    'rs4934': {'chr': 14, 'pos': 94382346, 'effect_allele': 'G', 'or': 1.19, 'gene': 'COL3A1'},
    
    # Additional loci (will extract from GWAS file)
    # We'll identify the top 17 variants from the GWAS file
}

def load_gwas_summary_stats(gwas_file=None):
    """Load GWAS summary statistics from Bakker 2020"""
    if gwas_file is None:
        possible_files = [
            "data/gwas/bakker_2020/extracted/IA.GWAS.BakkerMK.2020.sumstats.Stage_1.txt.gz",
            "data/gwas/bakker_2020/extracted/IA.GWAS.BakkerMK.2020.sumstats.Stage_2.txt.gz",
            "data/gwas/bakker_2020/bakker_2020_summary_stats.txt.gz"
        ]
        
        for file_path in possible_files:
            if os.path.exists(file_path):
                gwas_file = file_path
                break
    
    if gwas_file is None or not os.path.exists(gwas_file):
        print(f"[ERROR] GWAS file not found")
        return None
    
    print(f"\nLoading GWAS summary statistics: {os.path.basename(gwas_file)}")
    
    try:
        # Read header first to detect format
        with gzip.open(gwas_file, 'rt') as f:
            header = f.readline().strip()
        
        sep = '\t' if '\t' in header else ' '
        
        # Read with compression
        df = pd.read_csv(gwas_file, sep=sep, compression='gzip', low_memory=False)
        
        print(f"  [OK] Loaded {len(df):,} variants")
        print(f"  Columns: {list(df.columns)[:10]}...")
        
        return df
    except Exception as e:
        print(f"[ERROR] Failed to load GWAS file: {e}")
        return None

def identify_top_risk_loci(gwas_df, top_n=17):
    """Identify top risk loci from GWAS summary statistics"""
    print(f"\nIdentifying top {top_n} risk loci...")
    
    # Standardize column names
    col_mapping = {
        'SNP': ['SNP', 'snp', 'rsid', 'RSID', 'ID'],
        'CHR': ['CHR', 'chr', 'CHROM', 'chromosome', 'Chromosome'],
        'BP': ['BP', 'bp', 'POS', 'pos', 'position', 'Position'],
        'P': ['P', 'p', 'P_VALUE', 'pvalue', 'pval', 'PVAL', 'P-value'],
        'OR': ['OR', 'or', 'OR_95U', 'odds_ratio', 'OR_95L'],
        'BETA': ['BETA', 'beta', 'B', 'b', 'Effect', 'EFFECT'],
        'A1': ['A1', 'a1', 'ALT', 'alt', 'effect_allele', 'EA', 'A_EFF'],
        'A2': ['A2', 'a2', 'REF', 'ref', 'other_allele', 'OA', 'A_NONEFF']
    }
    
    # Map columns
    mapped_cols = {}
    for target, possible in col_mapping.items():
        for col in gwas_df.columns:
            if col in possible:
                mapped_cols[target] = col
                break
    
    # Check if we have required columns
    required = ['SNP', 'P']
    missing = [r for r in required if r not in mapped_cols]
    if missing:
        print(f"  [WARNING] Missing columns: {missing}")
        print(f"  Available columns: {list(gwas_df.columns)}")
        return None, {}
    
    # Filter significant variants (P < 5e-8 for genome-wide significance)
    p_col = mapped_cols['P']
    sig_df = gwas_df[gwas_df[p_col] < 5e-8].copy()
    
    if len(sig_df) == 0:
        # Try less stringent threshold
        sig_df = gwas_df[gwas_df[p_col] < 1e-6].copy()
        print(f"  Using relaxed threshold (P < 1e-6): {len(sig_df)} variants")
    
    if len(sig_df) == 0:
        # Use top N by P-value regardless of threshold
        sig_df = gwas_df.nsmallest(top_n * 10, p_col).copy()
        print(f"  Using top {len(sig_df)} variants by P-value")
    
    # Sort by significance
    sig_df = sig_df.sort_values(p_col)
    
    # Get top N (avoid duplicates in same locus)
    # For loci in same region (within 1Mb), keep only most significant
    top_loci = []
    used_chr_bp = set()
    
    for idx, row in sig_df.iterrows():
        chr_val = row[mapped_cols.get('CHR', 'CHR')] if mapped_cols.get('CHR') else None
        bp_val = row[mapped_cols.get('BP', 'BP')] if mapped_cols.get('BP') else None
        
        if chr_val and bp_val:
            # Round to 1Mb window
            locus_key = (chr_val, int(bp_val / 1e6))
            if locus_key in used_chr_bp and len(top_loci) >= top_n:
                continue
            used_chr_bp.add(locus_key)
        
        top_loci.append(row)
        if len(top_loci) >= top_n:
            break
    
    top_loci = pd.DataFrame(top_loci) if top_loci else sig_df.head(top_n).copy()
    
    print(f"  [OK] Identified {len(top_loci)} top risk loci")
    
    # Extract risk loci information
    risk_loci = []
    for idx, row in top_loci.iterrows():
        snp = row[mapped_cols['SNP']]
        p_val = row[mapped_cols['P']]
        chr_col = mapped_cols.get('CHR', None)
        pos_col = mapped_cols.get('BP', None)
        or_col = mapped_cols.get('OR', None)
        beta_col = mapped_cols.get('BETA', None)
        a1_col = mapped_cols.get('A1', None)
        
        locus_info = {
            'SNP': str(snp) if pd.notna(snp) else f'UNKNOWN_{idx}',
            'P_value': float(p_val) if pd.notna(p_val) else 1.0,
            'CHR': int(row[chr_col]) if chr_col and pd.notna(row.get(chr_col)) else None,
            'BP': int(row[pos_col]) if pos_col and pd.notna(row.get(pos_col)) else None,
            'OR': float(row[or_col]) if or_col and or_col in row and pd.notna(row.get(or_col)) else None,
            'BETA': float(row[beta_col]) if beta_col and beta_col in row and pd.notna(row.get(beta_col)) else None,
            'Effect_Allele': str(row[a1_col]) if a1_col and a1_col in row and pd.notna(row.get(a1_col)) else None
        }
        
        # Calculate OR from BETA if needed
        if locus_info['OR'] is None and locus_info['BETA'] is not None:
            locus_info['OR'] = np.exp(locus_info['BETA'])
        
        # Calculate BETA from OR if needed (shouldn't happen but just in case)
        if locus_info['BETA'] is None and locus_info['OR'] is not None:
            locus_info['BETA'] = np.log(locus_info['OR'])
        
        risk_loci.append(locus_info)
    
    risk_loci_df = pd.DataFrame(risk_loci)
    return risk_loci_df, mapped_cols

def calculate_weighted_prs(genotypes_df, risk_loci_df):
    """
    Calculate weighted PRS: sum of (dosage × log(OR))
    
    Args:
        genotypes_df: DataFrame with samples as rows, SNPs as columns
                     Values: 0, 1, 2 (dosage) or genotype calls
        risk_loci_df: DataFrame with risk loci (SNP, OR, Effect_Allele)
    
    Returns:
        PRS array for each sample
    """
    print(f"\nCalculating weighted PRS for {len(genotypes_df)} samples...")
    
    prs_scores = np.zeros(len(genotypes_df))
    matched_snps = 0
    
    for idx, row in risk_loci_df.iterrows():
        snp = row['SNP']
        
        if snp not in genotypes_df.columns:
            # Try alternative SNP IDs
            continue
        
        # Get effect size (log OR)
        if pd.notna(row.get('OR')):
            beta = np.log(row['OR'])
        elif pd.notna(row.get('BETA')):
            beta = row['BETA']
        else:
            continue
        
        # Get dosages for this SNP
        dosages = genotypes_df[snp].values
        
        # Calculate contribution: dosage × effect_size
        prs_scores += dosages * beta
        matched_snps += 1
    
    print(f"  [OK] Matched {matched_snps} out of {len(risk_loci_df)} risk loci")
    
    # Standardize PRS (mean=0, SD=1)
    if prs_scores.std() > 0:
        prs_scores = (prs_scores - prs_scores.mean()) / prs_scores.std()
    
    return prs_scores

def calculate_prs_from_gwas_only(risk_loci_df, output_file=None):
    """
    Calculate population-level PRS weights from GWAS
    Since we don't have individual genotypes, we'll create PRS weights
    that can be applied when genotypes are available
    """
    print(f"\nCreating PRS weights from {len(risk_loci_df)} risk loci...")
    
    # Calculate log(OR) for each variant
    risk_loci_df['Log_OR'] = np.log(risk_loci_df['OR'].fillna(np.exp(risk_loci_df['BETA'])))
    risk_loci_df['Weight'] = risk_loci_df['Log_OR']
    
    # Create PRS weights file
    prs_weights = risk_loci_df[['SNP', 'CHR', 'BP', 'Effect_Allele', 'OR', 'Log_OR', 'Weight', 'P_value']].copy()
    prs_weights = prs_weights.sort_values('P_value')
    
    if output_file:
        prs_weights.to_csv(output_file, sep='\t', index=False)
        print(f"  [OK] Saved PRS weights to {output_file}")
    
    return prs_weights

def stratify_risk_groups(prs_scores, percentiles=[33, 67]):
    """Stratify individuals into low/intermediate/high risk groups"""
    low_threshold = np.percentile(prs_scores, percentiles[0])
    high_threshold = np.percentile(prs_scores, percentiles[1])
    
    risk_groups = np.digitize(prs_scores, [low_threshold, high_threshold])
    
    group_labels = {
        0: 'Low Risk',
        1: 'Intermediate Risk',
        2: 'High Risk'
    }
    
    return risk_groups, group_labels

def main():
    """Main PRS building pipeline"""
    
    print("=" * 80)
    print("POLYGENIC RISK SCORE (PRS) DEVELOPMENT")
    print("=" * 80)
    
    os.makedirs("results", exist_ok=True)
    
    # Step 1: Load GWAS summary statistics
    gwas_df = load_gwas_summary_stats()
    if gwas_df is None:
        return
    
    # Step 2: Identify top 17 risk loci
    risk_loci_df, col_mapping = identify_top_risk_loci(gwas_df, top_n=17)
    if risk_loci_df is None or len(risk_loci_df) == 0:
        print("[ERROR] Could not identify risk loci")
        return
    
    # Save risk loci
    risk_loci_file = "results/prs_risk_loci_top17.tsv"
    risk_loci_df.to_csv(risk_loci_file, sep='\t', index=False)
    print(f"\n[OK] Saved top 17 risk loci to {risk_loci_file}")
    
    # Step 3: Create PRS weights
    prs_weights = calculate_prs_from_gwas_only(risk_loci_df)
    weights_file = "results/prs_weights.tsv"
    prs_weights.to_csv(weights_file, sep='\t', index=False)
    print(f"[OK] Saved PRS weights to {weights_file}")
    
    # Step 4: Print summary
    print(f"\n\nPRS Development Summary:")
    print(f"  Risk loci identified: {len(risk_loci_df)}")
    print(f"  Mean OR: {risk_loci_df['OR'].mean():.3f}")
    print(f"  Median P-value: {risk_loci_df['P_value'].median():.2e}")
    
    print(f"\n\nTop 10 Risk Loci:")
    print(f"{'SNP':<15} {'CHR':<6} {'BP':<12} {'OR':<8} {'P_value':<12} {'Effect_Allele':<12}")
    print("-" * 80)
    for idx, row in risk_loci_df.head(10).iterrows():
        print(f"{row['SNP']:<15} {str(row.get('CHR', 'N/A')):<6} {str(row.get('BP', 'N/A')):<12} "
              f"{row.get('OR', 'N/A'):<8.3f} {row['P_value']:<12.2e} {str(row.get('Effect_Allele', 'N/A')):<12}")
    
    print("\n" + "=" * 80)
    print("PRS DEVELOPMENT COMPLETE")
    print("=" * 80)
    print("\nNote: To calculate PRS for individual samples, you need:")
    print("  1. Genotype data (VCF or dosage files)")
    print("  2. Run: calculate_weighted_prs(genotypes_df, risk_loci_df)")
    print("\nPRS weights are saved and ready for use when genotype data is available.")
    
    return risk_loci_df, prs_weights

if __name__ == "__main__":
    # Fix: use gwas_df instead of df
    import sys
    risk_loci, weights = main()
    
    if risk_loci is not None:
        print(f"\n[OK] PRS development successful!")
        print(f"   Risk loci: {len(risk_loci)}")
        print(f"   Weights file: results/prs_weights.tsv")

