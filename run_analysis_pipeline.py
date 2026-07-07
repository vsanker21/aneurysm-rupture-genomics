"""
Main analysis pipeline - runs GWAS preparation and rupture prediction analysis
"""
import os
import sys

def main():
    """Run the complete analysis pipeline"""
    
    print("=" * 80)
    print("CEREBRAL ANEURYSM GENOMICS PROJECT - ANALYSIS PIPELINE")
    print("=" * 80)
    
    # Create results directory
    os.makedirs("results", exist_ok=True)
    
    # Step 1: Prepare GWAS data for PRS
    print("\n" + "=" * 80)
    print("STEP 1: PREPARING GWAS DATA FOR PRS")
    print("=" * 80)
    
    try:
        from prepare_gwas_for_prs import prepare_gwas_for_prs
        gwas_df, summary = prepare_gwas_for_prs()
        print("\n✓ GWAS data preparation complete!")
    except Exception as e:
        print(f"\n⚠️  Error in GWAS preparation: {e}")
        print("Continuing with RNA-seq analysis...")
        gwas_df = None
    
    # Step 2: Rupture prediction analysis
    print("\n" + "=" * 80)
    print("STEP 2: RUPTURE PREDICTION ANALYSIS (GSE122897)")
    print("=" * 80)
    
    try:
        from rupture_prediction_analysis import (
            load_gse122897, quality_control, normalize_data, 
            differential_expression, identify_top_rupture_genes
        )
        
        # Load data
        expr_df, metadata = load_gse122897()
        
        if expr_df is not None:
            # Quality control
            expr_df_qc, qc_metrics = quality_control(expr_df, metadata)
            
            # Normalization
            expr_norm, norm_method = normalize_data(expr_df_qc)
            
            # For GSE122897, we need to manually annotate samples
            # Based on the paper: "RNA Sequencing Analysis of Intracranial Aneurysm Walls 
            # Reveals Involvement of Lysosomes and Immunoglobulins in Rupture"
            # The dataset compares ruptured vs unruptured aneurysms
            
            # Try to infer or create sample groups
            # Since we don't have metadata, we'll need to check the paper or 
            # use a different approach - for now, let's create a placeholder
            print("\n⚠️  Sample group information not available in files.")
            print("Creating sample groups based on typical GSE122897 structure...")
            print("(In practice, you would need to check the publication or GEO page)")
            
            # For demonstration, split samples in half
            # In reality, you need the actual group labels from the publication
            samples = expr_df_qc.columns.tolist()
            n_samples = len(samples)
            
            # Create a simple metadata structure
            # NOTE: This is a placeholder - you need actual group labels!
            metadata_simple = pd.DataFrame({
                'Sample': samples,
                'Group': ['Ruptured' if i < n_samples//2 else 'Unruptured' 
                         for i in range(n_samples)]
            })
            
            print(f"\n⚠️  Using placeholder groups (check publication for actual labels):")
            print(f"  Ruptured: {sum(metadata_simple['Group'] == 'Ruptured')} samples")
            print(f"  Unruptured: {sum(metadata_simple['Group'] == 'Unruptured')} samples")
            
            # Update qc_metrics with groups
            qc_metrics = qc_metrics.merge(metadata_simple, on='Sample', how='left')
            
            # Differential expression
            de_results = differential_expression(expr_norm, qc_metrics)
            
            if de_results is not None:
                # Identify top genes
                top_genes = identify_top_rupture_genes(de_results, top_n=50)
                
                # Save results
                output_file = "results/gse122897_differential_expression.tsv"
                de_results.to_csv(output_file, sep='\t', index=False)
                print(f"\n✓ Saved full DE results to {output_file}")
                
                print("\n" + "=" * 80)
                print("ANALYSIS COMPLETE!")
                print("=" * 80)
                print("\n⚠️  IMPORTANT: Sample groups are placeholders!")
                print("   Please verify actual ruptured/unruptured labels from:")
                print("   - GSE122897 publication")
                print("   - GEO database metadata")
                print("   - Contact dataset authors if needed")
            else:
                print("\n⚠️  Could not complete differential expression analysis")
        else:
            print("\n⚠️  Could not load GSE122897 data")
            
    except Exception as e:
        print(f"\n⚠️  Error in rupture prediction analysis: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "=" * 80)
    print("PIPELINE COMPLETE")
    print("=" * 80)
    print("\nNext steps:")
    print("  1. Verify sample group labels for GSE122897")
    print("  2. Re-run DE analysis with correct group labels")
    print("  3. Proceed with PRS development using prepared GWAS data")
    print("  4. Integrate results for multi-omics model")

if __name__ == "__main__":
    import pandas as pd
    main()






