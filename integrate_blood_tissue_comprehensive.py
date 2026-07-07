"""
Comprehensive integration of blood (GSE159610) and tissue (GSE122897) DE results
Compares substantial fold change genes from blood with tissue DE genes
"""
import pandas as pd
import numpy as np
import os

def load_blood_results():
    """Load blood DE results"""
    blood_file = "results/gse159610_substantial_fc_genes.tsv"
    
    if not os.path.exists(blood_file):
        print(f"[ERROR] Blood results not found: {blood_file}")
        return None
    
    df = pd.read_csv(blood_file, sep='\t')
    print(f"Loaded {len(df):,} blood genes with substantial FC")
    return df

def load_tissue_results():
    """Load tissue DE results"""
    tissue_file = "results/gse122897_differential_expression.tsv"
    
    if not os.path.exists(tissue_file):
        print(f"[ERROR] Tissue results not found: {tissue_file}")
        return None
    
    df = pd.read_csv(tissue_file, sep='\t')
    print(f"Loaded {len(df):,} tissue genes")
    return df

def map_gene_ids(blood_df, tissue_df):
    """
    Map gene IDs between Ensembl (blood) and Entrez (tissue) formats
    This is a simplified mapping - in practice would use proper ID conversion
    """
    print("\n" + "=" * 80)
    print("MAPPING GENE IDS")
    print("=" * 80)
    
    # Try to extract numeric part from Ensembl IDs for matching
    # Ensembl IDs like ENSG00000229873.1 -> try to match with Entrez IDs
    # This is a simplified approach - proper mapping would use biomaRt or similar
    
    # For now, we'll create a mapping based on common patterns
    # In practice, you'd use a proper gene ID conversion tool
    
    print("  Note: Using simplified gene ID matching")
    print("  For production, use biomaRt, mygene, or g:Profiler for proper ID conversion")
    
    # Extract numeric parts (this is a placeholder - real mapping needed)
    blood_df['Gene_ID_clean'] = blood_df['Gene_ID'].astype(str).str.split('.').str[0]
    
    # Try to match - this is simplified, real implementation would use proper mapping
    # For demonstration, we'll create a structure that can be filled with proper mapping
    
    return blood_df, tissue_df

def integrate_results(blood_df, tissue_df):
    """
    Integrate blood and tissue results
    Focus on genes with substantial changes in both or either dataset
    """
    print("\n" + "=" * 80)
    print("INTEGRATING BLOOD AND TISSUE RESULTS")
    print("=" * 80)
    
    # For now, create separate analyses since proper ID mapping requires external tools
    # We'll focus on comparing patterns and effect sizes
    
    print("\nBlood Results Summary:")
    print(f"  Total substantial FC genes: {len(blood_df):,}")
    print(f"  Upregulated: {(blood_df['Log2FC'] > 0).sum():,}")
    print(f"  Downregulated: {(blood_df['Log2FC'] < 0).sum():,}")
    print(f"  Large FC (|Log2FC| > 1.0): {(np.abs(blood_df['Log2FC']) > 1.0).sum():,}")
    print(f"  Very large FC (|Log2FC| > 1.585): {(np.abs(blood_df['Log2FC']) > 1.585).sum():,}")
    
    print("\nTissue Results Summary:")
    print(f"  Total genes: {len(tissue_df):,}")
    print(f"  FDR < 0.1: {(tissue_df['P_adj_BH'] < 0.1).sum():,}")
    print(f"  FDR < 0.2: {(tissue_df['P_adj_BH'] < 0.2).sum():,}")
    print(f"  Large FC (|Log2FC| > 1.0): {(np.abs(tissue_df['Log2FC']) > 1.0).sum():,}")
    print(f"  Very large FC (|Log2FC| > 1.585): {(np.abs(tissue_df['Log2FC']) > 1.585).sum():,}")
    
    # Compare effect size distributions
    print("\n" + "=" * 80)
    print("EFFECT SIZE COMPARISON")
    print("=" * 80)
    
    blood_fc = np.abs(blood_df['Log2FC'])
    tissue_fc = np.abs(tissue_df['Log2FC'])
    
    print(f"\nBlood FC statistics:")
    print(f"  Mean |Log2FC|: {blood_fc.mean():.3f}")
    print(f"  Median |Log2FC|: {blood_fc.median():.3f}")
    print(f"  Max |Log2FC|: {blood_fc.max():.3f}")
    
    print(f"\nTissue FC statistics:")
    print(f"  Mean |Log2FC|: {tissue_fc.mean():.3f}")
    print(f"  Median |Log2FC|: {tissue_fc.median():.3f}")
    print(f"  Max |Log2FC|: {tissue_fc.max():.3f}")
    
    # Identify top genes in each dataset
    print("\n" + "=" * 80)
    print("TOP GENES IN EACH DATASET")
    print("=" * 80)
    
    # Top blood genes by effect size
    blood_df['Effect_Score'] = np.abs(blood_df['Log2FC']) * (-np.log10(blood_df['P_value'] + 1e-300))
    top_blood = blood_df.nlargest(20, 'Effect_Score')
    
    print(f"\nTop 20 Blood Genes (by effect score):")
    print(top_blood[['Gene_ID', 'Log2FC', 'P_value', 'Effect_Score']].to_string(index=False))
    
    # Top tissue genes
    tissue_df['Effect_Score'] = np.abs(tissue_df['Log2FC']) * (-np.log10(tissue_df['P_value'] + 1e-300))
    top_tissue = tissue_df.nlargest(20, 'Effect_Score')
    
    print(f"\nTop 20 Tissue Genes (by effect score):")
    print(top_tissue[['Gene', 'Log2FC', 'P_value', 'P_adj_BH', 'Effect_Score']].to_string(index=False))
    
    # Create integration summary
    integration_summary = {
        'Blood_Total_Substantial_FC': len(blood_df),
        'Blood_Large_FC': (np.abs(blood_df['Log2FC']) > 1.0).sum(),
        'Blood_Very_Large_FC': (np.abs(blood_df['Log2FC']) > 1.585).sum(),
        'Tissue_Total_Genes': len(tissue_df),
        'Tissue_FDR_Significant_0.1': (tissue_df['P_adj_BH'] < 0.1).sum(),
        'Tissue_Large_FC': (np.abs(tissue_df['Log2FC']) > 1.0).sum(),
        'Tissue_Very_Large_FC': (np.abs(tissue_df['Log2FC']) > 1.585).sum(),
        'Blood_Mean_FC': blood_fc.mean(),
        'Tissue_Mean_FC': tissue_fc.mean(),
    }
    
    return integration_summary, top_blood, top_tissue

def save_integration_results(integration_summary, top_blood, top_tissue):
    """Save integration results"""
    os.makedirs("results", exist_ok=True)
    
    # Save summary
    summary_df = pd.DataFrame([integration_summary])
    summary_file = "results/blood_tissue_integration_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"\n[OK] Saved integration summary to {summary_file}")
    
    # Save top genes
    top_blood_file = "results/blood_top_20_genes.tsv"
    top_blood.to_csv(top_blood_file, sep='\t', index=False)
    print(f"[OK] Saved top 20 blood genes to {top_blood_file}")
    
    top_tissue_file = "results/tissue_top_20_genes.tsv"
    top_tissue.to_csv(top_tissue_file, sep='\t', index=False)
    print(f"[OK] Saved top 20 tissue genes to {top_tissue_file}")

def main():
    """Main integration pipeline"""
    print("=" * 80)
    print("BLOOD-TISSUE INTEGRATION ANALYSIS")
    print("=" * 80)
    
    # Load results
    blood_df = load_blood_results()
    tissue_df = load_tissue_results()
    
    if blood_df is None or tissue_df is None:
        print("[ERROR] Could not load required data files")
        return
    
    # Map gene IDs (simplified for now)
    blood_df, tissue_df = map_gene_ids(blood_df, tissue_df)
    
    # Integrate results
    integration_summary, top_blood, top_tissue = integrate_results(blood_df, tissue_df)
    
    # Save results
    save_integration_results(integration_summary, top_blood, top_tissue)
    
    print("\n" + "=" * 80)
    print("INTEGRATION COMPLETE!")
    print("=" * 80)
    print("\nNote: For proper gene ID mapping between Ensembl and Entrez IDs,")
    print("      use biomaRt (R), mygene (Python), or g:Profiler API")
    print("      to identify overlapping genes between datasets.")

if __name__ == "__main__":
    main()






