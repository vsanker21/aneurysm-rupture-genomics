"""
Master script to run complete validation and integration pipeline:
1. Pathway enrichment
2. Gene validation
3. Functional annotation
4. GWAS integration
5. Blood transcriptomics integration
"""
import os
import sys
import subprocess
import time

def run_script(script_name, description):
    """Run a Python script and report results"""
    print("\n" + "=" * 80)
    print(f"RUNNING: {description}")
    print("=" * 80)
    
    try:
        result = subprocess.run(
            [sys.executable, script_name],
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        
        if result.returncode == 0:
            print(f"[OK] {description} completed successfully")
            if result.stdout:
                print(result.stdout[-2000:])  # Last 2000 chars
            return True
        else:
            print(f"[ERROR] {description} failed with return code {result.returncode}")
            if result.stderr:
                print("Error output:")
                print(result.stderr[-1000:])  # Last 1000 chars
            return False
            
    except subprocess.TimeoutExpired:
        print(f"[ERROR] {description} timed out after 10 minutes")
        return False
    except Exception as e:
        print(f"[ERROR] {description} failed: {e}")
        return False

def main():
    """Run complete validation pipeline"""
    
    print("=" * 80)
    print("COMPLETE VALIDATION AND INTEGRATION PIPELINE")
    print("=" * 80)
    print("\nThis pipeline will run:")
    print("  1. Pathway enrichment analysis")
    print("  2. Gene validation")
    print("  3. Functional annotation (enhanced)")
    print("  4. GWAS integration")
    print("  5. Blood transcriptomics integration")
    print("\n[NOTE] This may take 30-60 minutes depending on API rates")
    
    # Ensure results directory exists
    os.makedirs("results", exist_ok=True)
    
    # Check if DE results exist
    de_file = "results/gse122897_differential_expression.tsv"
    if not os.path.exists(de_file):
        print(f"\n[ERROR] DE results file not found: {de_file}")
        print("   Please run analyze_gse122897.py first")
        return
    
    # Run each analysis
    scripts = [
        ("pathway_enrichment_analysis.py", "Pathway Enrichment Analysis"),
        ("gene_function_annotation.py", "Functional Annotation"),
        ("validate_de_genes.py", "Gene Validation"),
        ("integrate_gwas.py", "GWAS Integration"),
        ("integrate_blood_transcriptomics.py", "Blood Transcriptomics Integration")
    ]
    
    results = {}
    
    for script, description in scripts:
        if os.path.exists(script):
            success = run_script(script, description)
            results[description] = success
            time.sleep(2)  # Brief pause between scripts
        else:
            print(f"\n[WARNING] Script not found: {script}")
            results[description] = False
    
    # Print summary
    print("\n" + "=" * 80)
    print("PIPELINE SUMMARY")
    print("=" * 80)
    
    for description, success in results.items():
        status = "[OK]" if success else "[FAILED]"
        print(f"  {status} {description}")
    
    successful = sum(1 for s in results.values() if s)
    total = len(results)
    
    print(f"\nCompleted: {successful}/{total} analyses")
    
    if successful == total:
        print("\n[SUCCESS] All analyses completed successfully!")
        print("\nOutput files are in the 'results/' directory:")
        print("  - pathway_enrichment_*.tsv")
        print("  - gene_validation_report_*.tsv")
        print("  - gse122897_annotated_*.tsv")
        print("  - gwas_integration_report.tsv")
        print("  - blood_tissue_integration_*.tsv")
        print("  - biomarker_candidates_*.tsv")
    else:
        print(f"\n[WARNING] {total - successful} analysis(es) failed. Check error messages above.")

if __name__ == "__main__":
    main()






