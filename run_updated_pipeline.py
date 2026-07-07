"""
Main Pipeline Runner for Updated Analysis
Runs preprocessing and analysis for all datasets with updated pipeline
"""
import os
import sys
from pathlib import Path

def run_preprocessing():
    """Run preprocessing for all datasets"""
    print("=" * 80)
    print("STEP 1: PREPROCESSING")
    print("=" * 80)
    
    # Preprocess GSE26969 microarray data
    print("\n--- Preprocessing GSE26969 (Microarray) ---")
    try:
        from preprocess_gse26969 import main as preprocess_gse26969
        preprocess_gse26969()
    except Exception as e:
        print(f"[WARNING] GSE26969 preprocessing failed: {e}")
        print("   This is expected if R/affy is not installed")
        print("   You can still proceed with GSE122897 analysis")
    
    print("\n[OK] Preprocessing step complete")

def run_integrated_analysis():
    """Run integrated analysis pipeline"""
    print("\n" + "=" * 80)
    print("STEP 2: INTEGRATED ANALYSIS")
    print("=" * 80)
    
    try:
        from analyze_integrated_datasets import main as integrated_analysis
        integrated_analysis()
    except Exception as e:
        print(f"[ERROR] Integrated analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True

def setup_single_cell_framework():
    """Set up single-cell integration framework"""
    print("\n" + "=" * 80)
    print("STEP 3: SINGLE-CELL INTEGRATION FRAMEWORK")
    print("=" * 80)
    
    try:
        from prepare_single_cell_integration import main as setup_sc
        setup_sc()
    except Exception as e:
        print(f"[WARNING] Single-cell framework setup failed: {e}")

def main():
    """Main pipeline runner"""
    
    print("=" * 80)
    print("UPDATED ANALYSIS PIPELINE")
    print("=" * 80)
    print("\nThis pipeline will:")
    print("1. Preprocess GSE26969 microarray data (if R is available)")
    print("2. Run integrated analysis (GSE122897 + GSE26969)")
    print("3. Set up single-cell atlas integration framework")
    print("\n" + "=" * 80)
    
    # Create results directory
    os.makedirs("results", exist_ok=True)
    
    # Step 1: Preprocessing
    run_preprocessing()
    
    # Step 2: Integrated analysis
    success = run_integrated_analysis()
    
    # Step 3: Single-cell framework
    setup_single_cell_framework()
    
    # Summary
    print("\n" + "=" * 80)
    print("PIPELINE EXECUTION COMPLETE")
    print("=" * 80)
    
    if success:
        print("\n[OK] Analysis completed successfully!")
        print("\nOutput files:")
        print("  - results/*_integrated_de.tsv (differential expression results)")
        print("  - results/*_top_genes_integrated.tsv (top significant genes)")
    else:
        print("\n[WARNING] Some steps may have failed")
        print("   Check error messages above for details")
    
    print("\nNext steps:")
    print("1. Review differential expression results")
    print("2. Update GSE26969 metadata if needed")
    print("3. When single-cell atlas is available, run deconvolution")
    print("4. Integrate findings across all datasets")

if __name__ == "__main__":
    main()






