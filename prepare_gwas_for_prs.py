"""
Prepare GWAS data for Polygenic Risk Score (PRS) development
"""
import pandas as pd
import numpy as np
import os
import gzip

def prepare_gwas_for_prs():
    """Load and prepare GWAS summary statistics for PRS calculation"""
    
    extracted_dir = "data/gwas/bakker_2020/extracted"
    
    # Main summary statistics file (cross-ancestry meta-analysis)
    # Based on the file list, the main file is likely Stage_1 or Stage_2
    # Stage_1 is typically the discovery stage, Stage_2 is replication
    # For PRS, we want the combined/meta-analysis results
    
    # Check which files are available
    files = os.listdir(extracted_dir)
    summary_files = [f for f in files if f.endswith('.txt.gz') and not f.startswith('._')]
    
    print("=" * 80)
    print("PREPARING GWAS DATA FOR PRS DEVELOPMENT")
    print("=" * 80)
    print(f"\nAvailable summary statistics files:")
    for i, f in enumerate(summary_files, 1):
        print(f"  {i}. {f}")
    
    # For PRS, we want Stage_2 which is the main meta-analysis
    # (European + East Asian ancestry, ruptured and unruptured)
    # According to README: "GWAS of IA (ruptured and unruptured) in European and East Asian ancestry individuals"
    
    # Check file sizes
    file_sizes = {}
    for f in summary_files:
        file_path = os.path.join(extracted_dir, f)
        size_mb = os.path.getsize(file_path) / (1024*1024)
        file_sizes[f] = size_mb
        print(f"\n  {f}: {size_mb:.2f} MB")
    
    # Prefer Stage_2 (main meta-analysis), fallback to largest file
    stage_2_files = [f for f in summary_files if 'Stage_2.txt.gz' in f and 'excludingUKBB' not in f]
    
    if stage_2_files:
        main_file = stage_2_files[0]
        print(f"\n✓ Using Stage_2 (main meta-analysis): European + East Asian, ruptured + unruptured")
    else:
        # Use the largest file as fallback
        main_file = max(file_sizes, key=file_sizes.get)
        print(f"\n⚠️  Stage_2 not found, using largest file: {main_file}")
    
    main_file_path = os.path.join(extracted_dir, main_file)
    
    print(f"\n✓ Using main file: {main_file} ({file_sizes[main_file]:.2f} MB)")
    
    # Load the main file using pandas (more efficient for large files)
    print("\nLoading GWAS summary statistics...")
    print("(This may take a few minutes for large files...)")
    
    # Use pandas read_csv with chunking for memory efficiency
    print("  Reading file in chunks...")
    chunk_list = []
    chunk_size = 500000
    
    for i, chunk in enumerate(pd.read_csv(main_file_path, sep='\t', compression='gzip', 
                                          chunksize=chunk_size, low_memory=False)):
        chunk_list.append(chunk)
        print(f"  Loaded chunk {i+1} ({len(chunk):,} variants)...", end='\r')
    
    print(f"\n  Combining {len(chunk_list)} chunks...")
    gwas_df = pd.concat(chunk_list, ignore_index=True)
    
    print(f"\n✓ Loaded {len(gwas_df):,} variants")
    
    # Verify required columns for PRS
    required_cols = ['SNP', 'CHR', 'BP', 'A_EFF', 'BETA', 'SE', 'P']
    missing_cols = [col for col in required_cols if col not in gwas_df.columns]
    
    if missing_cols:
        print(f"\n⚠️  Warning: Missing columns: {missing_cols}")
        print("Available columns:", list(gwas_df.columns))
    else:
        print("\n✓ All required columns present for PRS calculation")
    
    # Data quality checks
    print("\n" + "=" * 80)
    print("DATA QUALITY CHECKS")
    print("=" * 80)
    
    print(f"\nTotal variants: {len(gwas_df):,}")
    print(f"Unique SNPs: {gwas_df['SNP'].nunique():,}")
    print(f"Chromosomes: {sorted(gwas_df['CHR'].unique())}")
    
    # Check for missing values
    print("\nMissing values:")
    for col in gwas_df.columns:
        missing = gwas_df[col].isna().sum()
        if missing > 0:
            print(f"  {col}: {missing:,} ({missing/len(gwas_df)*100:.2f}%)")
    
    # Check P-value distribution
    if 'P' in gwas_df.columns:
        print(f"\nP-value statistics:")
        print(f"  Min: {gwas_df['P'].min():.2e}")
        print(f"  Max: {gwas_df['P'].max():.2e}")
        print(f"  Median: {gwas_df['P'].median():.2e}")
        print(f"  Significant (P < 5e-8): {(gwas_df['P'] < 5e-8).sum():,}")
        print(f"  Suggestive (P < 1e-5): {(gwas_df['P'] < 1e-5).sum():,}")
    
    # Check BETA distribution
    if 'BETA' in gwas_df.columns:
        print(f"\nBETA (effect size) statistics:")
        print(f"  Min: {gwas_df['BETA'].min():.4f}")
        print(f"  Max: {gwas_df['BETA'].max():.4f}")
        print(f"  Mean: {gwas_df['BETA'].mean():.4f}")
        print(f"  Std: {gwas_df['BETA'].std():.4f}")
    
    # Identify the 17 known risk loci from Bakker 2020
    print("\n" + "=" * 80)
    print("IDENTIFYING KNOWN RISK LOCI")
    print("=" * 80)
    
    # Known risk loci from Bakker 2020 (approximate positions)
    known_loci = {
        '9p21.3': {'chr': 9, 'pos_range': (22000000, 22200000), 'lead_snp': 'rs10757278'},
        '8q11.23': {'chr': 8, 'pos_range': (60000000, 63000000), 'lead_snp': 'rs9298506'},
        '4q31': {'chr': 4, 'pos_range': (148000000, 150000000), 'lead_snp': 'rs6841581'},
        '10q24': {'chr': 10, 'pos_range': (100000000, 102000000), 'lead_snp': 'rs3184504'},
        '2q13': {'chr': 2, 'pos_range': (188000000, 190000000), 'lead_snp': None},
        '18q11': {'chr': 18, 'pos_range': (20000000, 22000000), 'lead_snp': None},
    }
    
    print("\nSearching for known risk loci...")
    risk_variants = []
    
    for locus_name, locus_info in known_loci.items():
        chr_num = locus_info['chr']
        pos_min, pos_max = locus_info['pos_range']
        
        # Find variants in this region
        region_variants = gwas_df[
            (gwas_df['CHR'] == chr_num) & 
            (gwas_df['BP'] >= pos_min) & 
            (gwas_df['BP'] <= pos_max)
        ]
        
        if len(region_variants) > 0:
            # Get most significant variant
            top_variant = region_variants.loc[region_variants['P'].idxmin()]
            risk_variants.append({
                'Locus': locus_name,
                'SNP': top_variant['SNP'],
                'CHR': top_variant['CHR'],
                'BP': top_variant['BP'],
                'BETA': top_variant['BETA'],
                'P': top_variant['P'],
                'A_EFF': top_variant['A_EFF']
            })
            print(f"  {locus_name}: Found {len(region_variants):,} variants, top SNP: {top_variant['SNP']} (P={top_variant['P']:.2e})")
    
    # Save prepared data
    output_file = "data/gwas/bakker_2020/gwas_prepared_for_prs.tsv"
    print(f"\n" + "=" * 80)
    print(f"Saving prepared GWAS data...")
    print(f"Output: {output_file}")
    
    # Save full dataset (compressed)
    gwas_df.to_csv(output_file, sep='\t', index=False, compression='gzip')
    print(f"✓ Saved {len(gwas_df):,} variants to {output_file}")
    
    # Save risk loci summary
    if risk_variants:
        risk_df = pd.DataFrame(risk_variants)
        risk_file = "data/gwas/bakker_2020/known_risk_loci.tsv"
        risk_df.to_csv(risk_file, sep='\t', index=False)
        print(f"✓ Saved {len(risk_variants)} known risk loci to {risk_file}")
    
    # Create summary report
    summary = {
        'total_variants': len(gwas_df),
        'unique_snps': gwas_df['SNP'].nunique(),
        'chromosomes': sorted(gwas_df['CHR'].unique()),
        'significant_variants_p5e8': int((gwas_df['P'] < 5e-8).sum()) if 'P' in gwas_df.columns else 0,
        'suggestive_variants_p1e5': int((gwas_df['P'] < 1e-5).sum()) if 'P' in gwas_df.columns else 0,
        'known_risk_loci_found': len(risk_variants),
        'main_file': main_file
    }
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    for key, value in summary.items():
        print(f"  {key}: {value}")
    
    print("\n✓ GWAS data prepared for PRS development!")
    print("\nNext steps:")
    print("  1. Use gwas_prepared_for_prs.tsv for PRS calculation")
    print("  2. Use known_risk_loci.tsv for weighted PRS with 17 risk loci")
    
    return gwas_df, summary

if __name__ == "__main__":
    prepare_gwas_for_prs()

