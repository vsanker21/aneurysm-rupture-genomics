"""
Review Blood Integration Results for Biomarker Candidates
Identifies genes that are dysregulated in both tissue and blood
"""
import pandas as pd
import numpy as np
import os

def review_biomarker_candidates():
    """Review blood-tissue integration results for biomarker candidates"""
    
    print("=" * 80)
    print("BIOMARKER CANDIDATES REVIEW")
    print("=" * 80)
    
    # Load blood-tissue integration results
    blood_file = "results/blood_tissue_integration_fdr_0.1.tsv"
    
    if not os.path.exists(blood_file):
        print(f"[ERROR] Blood integration file not found: {blood_file}")
        return None
    
    # Check file size
    file_size = os.path.getsize(blood_file)
    if file_size < 10:
        print(f"[WARNING] Blood integration file appears empty")
        print("  This may indicate that blood DE analysis was not performed")
        print("  due to metadata or sample grouping issues")
        return None
    
    try:
        blood_df = pd.read_csv(blood_file, sep='\t')
        
        if len(blood_df) == 0:
            print("[WARNING] Blood integration file is empty")
            return None
        
        print(f"\nLoaded blood-tissue integration results: {len(blood_df)} genes")
        print(f"Columns: {list(blood_df.columns)}")
        
        # Check for biomarker-related columns
        biomarker_cols = ['Biomarker_Score', 'Dysregulated_in_Both', 'Concordant_Direction',
                         'Log2FC_Tissue', 'Log2FC_Blood', 'P_adj_BH_Tissue', 'P_adj_BH_Blood']
        
        available_cols = [col for col in biomarker_cols if col in blood_df.columns]
        print(f"Available biomarker columns: {available_cols}")
        
        # Identify biomarker candidates
        candidates = []
        
        if 'Biomarker_Score' in blood_df.columns:
            # Sort by biomarker score
            blood_df = blood_df.sort_values('Biomarker_Score', ascending=False)
            
            # Top candidates based on biomarker score
            top_candidates = blood_df.head(50)
            
            print(f"\n\nTop 50 Biomarker Candidates (by Biomarker Score):")
            print(f"{'Rank':<6} {'Gene_ID':<12} {'Log2FC_Tissue':<15} {'Log2FC_Blood':<15} {'Score':<12} {'Both':<8} {'Concordant':<10}")
            print("-" * 100)
            
            for idx, (i, row) in enumerate(top_candidates.iterrows(), 1):
                gene_id = str(row.get('Gene_ID', 'N/A'))
                log2fc_t = row.get('Log2FC_Tissue', np.nan)
                log2fc_b = row.get('Log2FC_Blood', np.nan)
                score = row.get('Biomarker_Score', np.nan)
                both = "Yes" if row.get('Dysregulated_in_Both', False) else "No"
                concordant = "Yes" if row.get('Concordant_Direction', False) else "No"
                
                fc_t_str = f"{log2fc_t:.3f}" if pd.notna(log2fc_t) else "N/A"
                fc_b_str = f"{log2fc_b:.3f}" if pd.notna(log2fc_b) else "N/A"
                score_str = f"{score:.2f}" if pd.notna(score) else "N/A"
                
                print(f"{idx:<6} {gene_id:<12} {fc_t_str:>14}  {fc_b_str:>14}  {score_str:>11}  {both:<8} {concordant:<10}")
                
                candidates.append({
                    'Rank': idx,
                    'Gene_ID': gene_id,
                    'Log2FC_Tissue': log2fc_t,
                    'Log2FC_Blood': log2fc_b,
                    'Biomarker_Score': score,
                    'Dysregulated_in_Both': both == "Yes",
                    'Concordant_Direction': concordant == "Yes"
                })
        elif 'Log2FC_Tissue' in blood_df.columns:
            # Sort by absolute Log2FC in tissue
            blood_df['Abs_Log2FC_Tissue'] = blood_df['Log2FC_Tissue'].abs()
            blood_df = blood_df.sort_values('Abs_Log2FC_Tissue', ascending=False)
            
            top_candidates = blood_df.head(50)
            
            print(f"\n\nTop 50 Biomarker Candidates (by |Log2FC_Tissue|):")
            print(f"{'Rank':<6} {'Gene_ID':<12} {'Log2FC_Tissue':<15} {'P_adj_BH_Tissue':<15}")
            print("-" * 80)
            
            for idx, (i, row) in enumerate(top_candidates.iterrows(), 1):
                gene_id = str(row.get('Gene_ID', 'N/A'))
                log2fc_t = row.get('Log2FC_Tissue', np.nan)
                p_adj_t = row.get('P_adj_BH_Tissue', np.nan)
                
                fc_t_str = f"{log2fc_t:.3f}" if pd.notna(log2fc_t) else "N/A"
                p_adj_str = f"{p_adj_t:.4f}" if pd.notna(p_adj_t) else "N/A"
                
                print(f"{idx:<6} {gene_id:<12} {fc_t_str:>14}  {p_adj_str:>14}")
                
                candidates.append({
                    'Rank': idx,
                    'Gene_ID': gene_id,
                    'Log2FC_Tissue': log2fc_t,
                    'P_adj_BH_Tissue': p_adj_t
                })
        else:
            print("[WARNING] No biomarker-related columns found")
            print("  Displaying all genes:")
            print(blood_df.head(20).to_string())
            return blood_df
        
        # Save biomarker candidates
        if candidates:
            candidates_df = pd.DataFrame(candidates)
            output_file = "results/biomarker_candidates_review.tsv"
            candidates_df.to_csv(output_file, sep='\t', index=False)
            print(f"\n[OK] Saved biomarker candidates to {output_file}")
        
        # Summary statistics
        print(f"\n\nBiomarker Candidates Summary:")
        print(f"  Total genes analyzed: {len(blood_df)}")
        
        if 'Dysregulated_in_Both' in blood_df.columns:
            both_count = blood_df['Dysregulated_in_Both'].sum()
            print(f"  Dysregulated in both tissue and blood: {both_count}")
        
        if 'Concordant_Direction' in blood_df.columns:
            concordant_count = blood_df['Concordant_Direction'].sum()
            print(f"  Concordant direction of change: {concordant_count}")
        
        if 'In_Blood_Dataset' in blood_df.columns:
            in_blood = blood_df['In_Blood_Dataset'].sum()
            print(f"  Present in blood dataset: {in_blood}")
        
        return candidates_df if candidates else blood_df
        
    except Exception as e:
        print(f"[ERROR] Failed to review biomarker candidates: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    """Main function"""
    candidates = review_biomarker_candidates()
    
    print("\n" + "=" * 80)
    print("BIOMARKER REVIEW COMPLETE")
    print("=" * 80)
    
    if candidates is not None:
        print(f"\n[NOTE] Review results saved to results/biomarker_candidates_review.tsv")
        print("  Next steps:")
        print("  1. Validate top candidates in independent cohorts")
        print("  2. Check if genes are detectable in blood samples")
        print("  3. Assess biomarker potential for clinical use")

if __name__ == "__main__":
    main()






