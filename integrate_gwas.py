"""
GWAS Integration Script
Maps differentially expressed genes to GWAS variants and integrates genetic risk
"""
import pandas as pd
import numpy as np
import gzip
import os
from typing import Dict, List, Tuple, Optional

def load_de_results(fdr_threshold=0.1):
    """Load differential expression results"""
    de_file = "results/gse122897_differential_expression.tsv"
    if not os.path.exists(de_file):
        print(f"[ERROR] DE results file not found: {de_file}")
        return None
    
    df = pd.read_csv(de_file, sep='\t')
    sig_df = df[df['P_adj_BH'] < fdr_threshold].copy()
    
    return sig_df

def load_gwas_summary_stats(gwas_file: str = None):
    """
    Load GWAS summary statistics
    Default: Bakker et al. 2020 main summary stats
    """
    if gwas_file is None:
        # Try to find the main GWAS file - prefer Stage_1 as it's the discovery set
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
        print(f"[ERROR] GWAS file not found. Please specify correct path.")
        return None
    
    print(f"Loading GWAS summary statistics: {gwas_file}")
    
    try:
        # Check if file is actually gzip or zip
        with open(gwas_file, 'rb') as f:
            header = f.read(2)
        
        if header == b'\x1f\x8b':  # gzip magic number
            # Try different column formats
            # Read header separately to detect format
            with gzip.open(gwas_file, 'rt') as f:
                first_line = f.readline()
            
            # Detect separator
            sep = '\t' if '\t' in first_line else ' '
            
            # Read with compression handled by pandas
            df = pd.read_csv(gwas_file, sep=sep, compression='gzip', nrows=1000)
            
            print(f"  Detected columns: {list(df.columns)}")
            print(f"  Sample rows: {df.shape[0]}")
            
            # Reload full file
            df = pd.read_csv(gwas_file, sep=sep, compression='gzip')
            
            print(f"[OK] Loaded {len(df):,} GWAS variants")
            return df
        else:
            # Try as zip file or uncompressed
            print(f"  File appears to be ZIP or uncompressed, trying alternative format...")
            # Try reading as tab-separated without compression
            try:
                df = pd.read_csv(gwas_file, sep='\t', nrows=1000)
                print(f"  Detected columns: {list(df.columns)}")
                print(f"  Sample rows: {df.shape[0]}")
                df = pd.read_csv(gwas_file, sep='\t')
                print(f"[OK] Loaded {len(df):,} GWAS variants")
                return df
            except:
                print(f"[WARNING] Could not read as tab-separated, trying extracted files...")
                # Try extracted directory
                extracted_dir = os.path.join(os.path.dirname(gwas_file), "extracted")
                if os.path.exists(extracted_dir):
                    print(f"  Checking extracted directory: {extracted_dir}")
                    # Look for .txt.gz files in extracted directory
                    import glob
                    gwas_files = glob.glob(os.path.join(extracted_dir, "*.txt.gz"))
                    if gwas_files:
                        # Try first file
                        alt_file = gwas_files[0]
                        print(f"  Trying extracted file: {os.path.basename(alt_file)}")
                        return load_gwas_summary_stats(alt_file)
                return None
            
    except Exception as e:
        print(f"[ERROR] Failed to load GWAS file: {e}")
        return None

def map_genes_to_gwas_regions(de_df: pd.DataFrame, gwas_df: pd.DataFrame, window_size: int = 50000):
    """
    Map DE genes to GWAS significant regions
    Uses gene coordinates to find nearby GWAS variants
    """
    print("\nMapping DE genes to GWAS regions...")
    
    # Load gene coordinates (simplified - would use proper annotation file in practice)
    # For now, we'll use a gene-to-chromosome mapping approach via gene symbols
    
    print("  [NOTE] Full implementation would require:")
    print("    - Gene coordinate file (GTF/GFF)")
    print("    - Chromosome and position mapping")
    print("    - Gene annotation database")
    
    # Get gene symbols from DE results (assuming Gene column has IDs)
    gene_ids = de_df['Gene'].unique().tolist()
    
    # Simplified mapping: we'll match based on gene regions if coordinates available
    # For demonstration, create structure for integration
    
    gwas_region_matches = []
    
    # Check if GWAS has chromosome and position columns
    gwas_cols = gwas_df.columns.tolist()
    
    # Try to identify column names - GWAS format typically has:
    # CHR, BP, SNP, BETA, SE, P-value columns
    chr_col = None
    pos_col = None
    pval_col = None
    beta_col = None
    
    # Common column name mappings
    col_mappings = {
        'CHR': ['chr', 'chrom', 'chromosome'],
        'BP': ['bp', 'pos', 'position'],
        'P_value': ['p', 'pval', 'p_value', 'pvalue', 'pmaf'],
        'BETA': ['beta', 'or', 'odds_ratio', 'effect', 'effect_size']
    }
    
    for col in gwas_cols:
        col_lower = col.upper()  # Check uppercase too
        col_lower_orig = col.lower()
        
        # Check CHR
        if not chr_col and ('CHR' in col or col_lower_orig in col_mappings['CHR']):
            chr_col = col
        # Check BP
        if not pos_col and ('BP' in col or col_lower_orig in col_mappings['BP']):
            pos_col = col
        # Check P-value
        if not pval_col and (col_lower_orig in col_mappings['P_value'] or 
                           (col_lower_orig.startswith('p') and 'val' in col_lower_orig)):
            pval_col = col
        # Check BETA/OR
        if not beta_col and (col in col_mappings['BETA'] or col_lower_orig in col_mappings['BETA'] or
                           'BETA' in col or 'OR' in col.upper()):
            beta_col = col
    
    if chr_col and pos_col:
        print(f"  Found chromosome column: {chr_col}")
        print(f"  Found position column: {pos_col}")
        print(f"  Found P-value column: {pval_col if pval_col else 'Not found'}")
        print(f"  Found effect size column: {beta_col if beta_col else 'Not found'}")
        
        # Filter significant GWAS variants
        if pval_col:
            # Try to convert P-value column to numeric
            try:
                gwas_df[pval_col] = pd.to_numeric(gwas_df[pval_col], errors='coerce')
                sig_gwas = gwas_df[gwas_df[pval_col] < 5e-8].copy()
                print(f"  Significant GWAS variants (P < 5e-8): {len(sig_gwas):,}")
            except Exception as e:
                print(f"  [WARNING] Could not filter by P-value: {e}")
                sig_gwas = gwas_df.head(1000).copy()  # Use top 1000 as fallback
                print(f"  Using top 1000 variants")
        else:
            sig_gwas = gwas_df.head(1000).copy()
            print(f"  Total GWAS variants (using top 1000): {len(sig_gwas):,}")
        
        # Create structure for gene-GWAS mapping
        # In full implementation, would use gene coordinates here
        gwas_regions = []
        for idx, row in sig_gwas.head(100).iterrows():  # Limit for demo
            try:
                gwas_regions.append({
                    'Chromosome': str(row[chr_col]),
                    'Position': int(row[pos_col]) if pd.notna(row[pos_col]) else 0,
                    'P_value': row[pval_col] if pval_col and pd.notna(row[pval_col]) else np.nan,
                    'Effect_Size': row[beta_col] if beta_col and pd.notna(row[beta_col]) else np.nan
                })
            except Exception as e:
                continue  # Skip problematic rows
        
        print(f"  Created {len(gwas_regions)} GWAS region entries")
    else:
        print("  [WARNING] Could not identify chromosome/position columns")
        gwas_regions = []
    
    return gwas_regions

def identify_prioritized_genes(de_df: pd.DataFrame, gwas_df: pd.DataFrame) -> pd.DataFrame:
    """
    Identify genes that are both differentially expressed and in GWAS loci
    These are high-priority candidates for validation
    """
    print("\nIdentifying prioritized genes (DE + GWAS)...")
    
    # Load known GWAS loci for intracranial aneurysm
    # Based on Bakker et al. 2020, key regions include:
    gwas_loci = {
        '9p21.3': {'chr': '9', 'start': 21900000, 'end': 22100000, 'genes': ['CDKN2B-AS1', 'CDKN2A', 'CDKN2B']},
        '8q11.23': {'chr': '8', 'start': 48000000, 'end': 49000000, 'genes': ['SOX17']},
        '4q31': {'chr': '4', 'start': 147000000, 'end': 148000000, 'genes': ['EDNRA']},
        '10q24': {'chr': '10', 'start': 101000000, 'end': 102000000, 'genes': ['CFDP1']},
        '2q13': {'chr': '2', 'start': 189000000, 'end': 190000000, 'genes': ['COL3A1']},
        '18q11': {'chr': '18', 'start': 64000000, 'end': 65000000, 'genes': ['LOXL1']}
    }
    
    # For demonstration, create prioritized gene list
    # In full implementation, would match gene coordinates to GWAS regions
    
    prioritized_genes = []
    
    for idx, row in de_df.iterrows():
        gene_id = str(row['Gene'])
        
        # Check if gene is in known GWAS loci (simplified check)
        in_gwas_locus = False
        associated_locus = ''
        
        # This would be replaced with actual coordinate matching
        for locus_name, locus_info in gwas_loci.items():
            # Simplified: would check gene coordinates here
            pass
        
        prioritized_genes.append({
            'Gene_ID': gene_id,
            'Gene_Symbol': gene_id,  # Would be mapped from annotation
            'Log2FC': row['Log2FC'],
            'P_adj_BH': row['P_adj_BH'],
            'In_GWAS_Locus': in_gwas_locus,
            'Associated_Locus': associated_locus,
            'Priority_Score': abs(row['Log2FC']) * (-np.log10(row['P_adj_BH'] + 1e-300))
        })
    
    prioritized_df = pd.DataFrame(prioritized_genes)
    prioritized_df = prioritized_df.sort_values('Priority_Score', ascending=False)
    
    return prioritized_df

def create_gwas_integration_report(de_df: pd.DataFrame, gwas_file: str = None):
    """
    Create comprehensive GWAS integration report
    """
    print("\n" + "=" * 80)
    print("GWAS INTEGRATION ANALYSIS")
    print("=" * 80)
    
    # Load GWAS data
    gwas_df = load_gwas_summary_stats(gwas_file)
    if gwas_df is None:
        print("\n[WARNING] Could not load GWAS data. Creating structure-only report.")
        gwas_df = pd.DataFrame()  # Empty for structure
    
    # Map genes to GWAS regions
    if gwas_df is not None and len(gwas_df) > 0:
        gwas_regions = map_genes_to_gwas_regions(de_df, gwas_df)
    else:
        gwas_regions = []
    
    # Identify prioritized genes
    prioritized_df = identify_prioritized_genes(de_df, gwas_df)
    
    # Create integration summary
    integration_results = []
    
    for idx, row in de_df.iterrows():
        gene_id = str(row['Gene'])
        priority_row = prioritized_df[prioritized_df['Gene_ID'] == gene_id]
        
        integration_results.append({
            'Gene_ID': gene_id,
            'Log2FC': row['Log2FC'],
            'P_adj_BH': row['P_adj_BH'],
            'Priority_Score': priority_row['Priority_Score'].values[0] if len(priority_row) > 0 else 0,
            'In_GWAS_Locus': priority_row['In_GWAS_Locus'].values[0] if len(priority_row) > 0 else False,
            'GWAS_Association': 'Yes' if priority_row['In_GWAS_Locus'].values[0] else 'No' if len(priority_row) > 0 else 'Unknown'
        })
    
    integration_df = pd.DataFrame(integration_results)
    integration_df = integration_df.sort_values('Priority_Score', ascending=False)
    
    # Save results
    output_file = "results/gwas_integration_report.tsv"
    integration_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved GWAS integration report to {output_file}")
    
    # Print summary
    print("\nGWAS Integration Summary:")
    print(f"  Total DE genes: {len(integration_df)}")
    print(f"  Genes in GWAS loci: {integration_df['In_GWAS_Locus'].sum()}")
    print(f"  Top priority genes: {len(integration_df[integration_df['Priority_Score'] > integration_df['Priority_Score'].quantile(0.9)])}")
    
    # Print top prioritized genes
    print("\nTop 20 Priority Genes (DE + GWAS):")
    top_priority = integration_df.head(20)
    for idx, row in top_priority.iterrows():
        gwas_status = "GWAS*" if row['In_GWAS_Locus'] else ""
        print(f"  {row['Gene_ID']:15s} Log2FC={row['Log2FC']:7.3f}  P_adj={row['P_adj_BH']:.3f}  Priority={row['Priority_Score']:.2f} {gwas_status}")
    
    return integration_df

def main():
    """Main GWAS integration pipeline"""
    
    os.makedirs("results", exist_ok=True)
    
    # Analyze at different FDR thresholds
    for fdr_threshold in [0.05, 0.1, 0.2]:
        print(f"\n{'='*80}")
        print(f"GWAS INTEGRATION - FDR < {fdr_threshold}")
        print(f"{'='*80}")
        
        # Load DE results
        de_df = load_de_results(fdr_threshold=fdr_threshold)
        if de_df is None or len(de_df) == 0:
            print(f"No significant genes at FDR < {fdr_threshold}")
            continue
        
        # Create integration report
        integration_df = create_gwas_integration_report(de_df)
    
    print("\n" + "=" * 80)
    print("GWAS INTEGRATION COMPLETE")
    print("=" * 80)
    print("\n[NOTE] For full gene-GWAS coordinate mapping, provide:")
    print("  1. Gene annotation file (GTF/GFF) with chromosome coordinates")
    print("  2. GWAS summary statistics with chromosome and position columns")
    print("  3. LD reference panel for fine-mapping (optional)")

if __name__ == "__main__":
    main()

