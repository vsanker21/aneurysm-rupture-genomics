"""
Analyze genes with substantial fold changes (|Log2FC| > 0.585)
Focus on effect sizes even if not FDR-significant
"""
import pandas as pd
import numpy as np
import os

def load_substantial_fc_genes():
    """Load genes with substantial fold changes"""
    fc_file = "results/gse159610_substantial_fc_genes.tsv"
    
    if not os.path.exists(fc_file):
        print(f"[ERROR] File not found: {fc_file}")
        return None
    
    df = pd.read_csv(fc_file, sep='\t')
    print(f"Loaded {len(df):,} genes with |Log2FC| > 0.585")
    return df

def categorize_genes(df):
    """Categorize genes by fold change magnitude and direction"""
    print("\n" + "=" * 80)
    print("CATEGORIZING GENES BY EFFECT SIZE")
    print("=" * 80)
    
    # By fold change magnitude
    df['FC_Category'] = 'Moderate'
    df.loc[np.abs(df['Log2FC']) > 1.0, 'FC_Category'] = 'Large'
    df.loc[np.abs(df['Log2FC']) > 1.585, 'FC_Category'] = 'Very Large'
    df.loc[np.abs(df['Log2FC']) > 2.0, 'FC_Category'] = 'Extreme'
    
    # By direction
    df['Direction'] = 'Upregulated'
    df.loc[df['Log2FC'] < 0, 'Direction'] = 'Downregulated'
    
    # By significance (nominal p-value)
    df['Significance'] = 'Not Significant'
    df.loc[df['P_value'] < 0.05, 'Significance'] = 'p < 0.05'
    df.loc[df['P_value'] < 0.01, 'Significance'] = 'p < 0.01'
    df.loc[df['P_value'] < 0.001, 'Significance'] = 'p < 0.001'
    
    # Summary
    print(f"\nBy Fold Change Magnitude:")
    print(df['FC_Category'].value_counts().to_string())
    
    print(f"\nBy Direction:")
    print(df['Direction'].value_counts().to_string())
    
    print(f"\nBy Significance (nominal p-value):")
    print(df['Significance'].value_counts().to_string())
    
    # Combined categories
    print(f"\nTop Priority Genes (Large FC + Low p-value):")
    top_priority = df[(np.abs(df['Log2FC']) > 1.0) & (df['P_value'] < 0.01)].copy()
    top_priority['Abs_Log2FC'] = np.abs(top_priority['Log2FC'])
    top_priority = top_priority.sort_values(['P_value', 'Abs_Log2FC'], ascending=[True, False])
    top_priority = top_priority.drop('Abs_Log2FC', axis=1)
    print(f"  Found {len(top_priority):,} genes")
    
    if len(top_priority) > 0:
        print(f"\n  Top 20:")
        print(top_priority.head(20)[['Gene_ID', 'Log2FC', 'P_value', 'FC_Category']].to_string(index=False))
    
    return df

def identify_top_candidates(df, top_n=100):
    """Identify top candidate genes for further analysis"""
    print("\n" + "=" * 80)
    print(f"TOP {top_n} CANDIDATE GENES")
    print("=" * 80)
    
    # Score genes by combining effect size and significance
    # Score = |Log2FC| * (-log10(p_value))
    df['Effect_Score'] = np.abs(df['Log2FC']) * (-np.log10(df['P_value'] + 1e-300))
    
    # Sort by score
    top_candidates = df.nlargest(top_n, 'Effect_Score').copy()
    
    print(f"\nTop {top_n} genes by combined effect size and significance:")
    print(top_candidates[['Gene_ID', 'Log2FC', 'P_value', 'Effect_Score', 'FC_Category']].to_string(index=False))
    
    # Separate upregulated and downregulated
    upregulated = top_candidates[top_candidates['Log2FC'] > 0]
    downregulated = top_candidates[top_candidates['Log2FC'] < 0]
    
    print(f"\n  Upregulated: {len(upregulated):,}")
    print(f"  Downregulated: {len(downregulated):,}")
    
    return top_candidates

def analyze_expression_levels(df):
    """Analyze expression levels of substantial FC genes"""
    print("\n" + "=" * 80)
    print("EXPRESSION LEVEL ANALYSIS")
    print("=" * 80)
    
    # Calculate mean expression across both groups
    df['Mean_Expression'] = (df['Mean_Aneurysm'] + df['Mean_Control']) / 2
    
    # Categorize by expression level
    df['Expression_Level'] = 'Low'
    df.loc[df['Mean_Expression'] > 10, 'Expression_Level'] = 'Medium'
    df.loc[df['Mean_Expression'] > 50, 'Expression_Level'] = 'High'
    df.loc[df['Mean_Expression'] > 100, 'Expression_Level'] = 'Very High'
    
    print(f"\nExpression level distribution:")
    print(df['Expression_Level'].value_counts().to_string())
    
    # High expression + large FC genes (most reliable)
    high_expr_large_fc = df[(df['Expression_Level'].isin(['High', 'Very High'])) & 
                            (np.abs(df['Log2FC']) > 1.0)]
    print(f"\nHigh expression + Large FC genes: {len(high_expr_large_fc):,}")
    
    if len(high_expr_large_fc) > 0:
        print(f"\n  Top 20:")
        top_high = high_expr_large_fc.nlargest(20, 'Effect_Score')
        print(top_high[['Gene_ID', 'Log2FC', 'P_value', 'Mean_Expression', 'Expression_Level']].to_string(index=False))
    
    return df

def save_prioritized_genes(df, top_candidates):
    """Save prioritized gene lists"""
    os.makedirs("results", exist_ok=True)
    
    # Save all substantial FC genes with categories
    output_file = "results/gse159610_substantial_fc_genes_categorized.tsv"
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved categorized genes to {output_file}")
    
    # Save top candidates
    output_file = "results/gse159610_top_candidate_genes.tsv"
    top_candidates.to_csv(output_file, sep='\t', index=False)
    print(f"[OK] Saved top candidates to {output_file}")
    
    # Save by category
    for category in ['Moderate', 'Large', 'Very Large', 'Extreme']:
        category_genes = df[df['FC_Category'] == category].copy()
        if len(category_genes) > 0:
            output_file = f"results/gse159610_fc_{category.lower().replace(' ', '_')}_genes.tsv"
            category_genes = category_genes.sort_values('P_value')
            category_genes.to_csv(output_file, sep='\t', index=False)
            print(f"[OK] Saved {len(category_genes):,} {category} FC genes to {output_file}")

def main():
    """Main analysis pipeline"""
    print("=" * 80)
    print("ANALYZING SUBSTANTIAL FOLD CHANGE GENES")
    print("=" * 80)
    
    # Load data
    df = load_substantial_fc_genes()
    if df is None:
        return
    
    # Categorize genes
    df = categorize_genes(df)
    
    # Identify top candidates
    top_candidates = identify_top_candidates(df, top_n=100)
    
    # Analyze expression levels
    df = analyze_expression_levels(df)
    
    # Save results
    save_prioritized_genes(df, top_candidates)
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)

if __name__ == "__main__":
    main()

