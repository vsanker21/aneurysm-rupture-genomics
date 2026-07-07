"""
Master script to run all advanced analyses:
1. Convert Ensembl IDs to gene symbols
2. Pathway enrichment with symbols
3. Literature annotation
4. Machine learning classification
5. DESeq2 setup (if available)
"""
import subprocess
import sys
import os

def run_script(script_name, description):
    """Run a Python script and report results"""
    print("\n" + "=" * 80)
    print(f"RUNNING: {description}")
    print("=" * 80)
    
    if not os.path.exists(script_name):
        print(f"[ERROR] Script not found: {script_name}")
        return False
    
    try:
        result = subprocess.run([sys.executable, script_name],
                              capture_output=True, text=True, timeout=1800)
        
        if result.returncode == 0:
            print("[OK] Script completed successfully")
            if result.stdout:
                print(result.stdout[-1000:])  # Last 1000 chars
            return True
        else:
            print("[ERROR] Script failed")
            if result.stderr:
                print("STDERR:", result.stderr[-500:])
            return False
    except subprocess.TimeoutExpired:
        print("[ERROR] Script timed out (>30 minutes)")
        return False
    except Exception as e:
        print(f"[ERROR] Failed to run script: {e}")
        return False

def main():
    """Run all advanced analyses"""
    print("=" * 80)
    print("RUNNING ALL ADVANCED ANALYSES")
    print("=" * 80)
    
    analyses = [
        ("convert_ensembl_to_symbols.py", "Gene ID Conversion (Ensembl → Symbols)"),
        ("pathway_enrichment_blood_genes_symbols.py", "Pathway Enrichment with Symbols"),
        ("annotate_top_genes_literature.py", "Literature Annotation"),
        ("ml_classification_blood_genes.py", "Machine Learning Classification"),
        ("setup_deseq2_analysis.py", "DESeq2 Setup and Analysis")
    ]
    
    results = {}
    
    for script, description in analyses:
        success = run_script(script, description)
        results[description] = success
        
        if not success and "DESeq2" not in description:
            print(f"\n[WARNING] {description} failed, but continuing...")
    
    # Summary
    print("\n" + "=" * 80)
    print("ANALYSIS SUMMARY")
    print("=" * 80)
    
    for description, success in results.items():
        status = "✓ COMPLETE" if success else "✗ FAILED"
        print(f"{status}: {description}")
    
    print("\n" + "=" * 80)
    print("ALL ANALYSES COMPLETE!")
    print("=" * 80)
    
    print("\nGenerated files:")
    print("1. results/gse159610_*_annotated.tsv - Genes with symbols and annotations")
    print("2. results/blood_pathway_enrichment_symbols_*.tsv - Pathway enrichment results")
    print("3. results/gse159610_top_genes_literature_annotated.tsv - Literature annotations")
    print("4. results/blood_ml_*.tsv - Machine learning results")
    print("5. results/gse159610_deseq2_results.tsv - DESeq2 results (if available)")

if __name__ == "__main__":
    main()






