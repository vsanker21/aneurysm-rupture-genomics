"""Check prerequisites and run the updated pipeline"""
import subprocess
import sys
import os
from pathlib import Path

def check_r_installation():
    """Check if R is installed"""
    print("=" * 80)
    print("CHECKING PREREQUISITES")
    print("=" * 80)
    
    r_available = False
    rscript_available = False
    
    # Check for R
    try:
        result = subprocess.run(['R', '--version'], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            r_available = True
            print("\n[OK] R is installed")
            print(f"   Version info: {result.stdout.split(chr(10))[0] if result.stdout else 'N/A'}")
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    
    # Check for Rscript
    try:
        result = subprocess.run(['Rscript', '--version'], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            rscript_available = True
            print("[OK] Rscript is available")
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    
    if not r_available and not rscript_available:
        print("\n[WARNING] R is not installed or not in PATH")
        print("   Microarray normalization will be skipped")
        print("   To install R: https://cran.r-project.org/")
        print("   After installing R, run: Rscript install_r_packages.R")
        return False
    
    # Check if packages are installed
    if rscript_available:
        print("\nChecking for required R packages...")
        check_script = """
        if (requireNamespace("affy", quietly = TRUE) && 
            requireNamespace("Biobase", quietly = TRUE)) {
            cat("SUCCESS: Packages installed\\n")
        } else {
            cat("WARNING: Packages not installed\\n")
            cat("Run: Rscript install_r_packages.R\\n")
        }
        """
        try:
            result = subprocess.run(['Rscript', '-e', check_script],
                                  capture_output=True, text=True, timeout=10)
            print(result.stdout)
            if "SUCCESS" in result.stdout:
                return True
        except Exception as e:
            print(f"   Could not check packages: {e}")
    
    return False

def run_preprocessing():
    """Run preprocessing step"""
    print("\n" + "=" * 80)
    print("STEP 1: PREPROCESSING")
    print("=" * 80)
    
    try:
        from preprocess_gse26969 import main as preprocess_main
        preprocess_main()
        return True
    except Exception as e:
        print(f"[WARNING] Preprocessing encountered issues: {e}")
        print("   This is expected if R is not available")
        print("   Continuing with RNA-seq analysis only...")
        return False

def run_integrated_analysis():
    """Run integrated analysis"""
    print("\n" + "=" * 80)
    print("STEP 2: INTEGRATED ANALYSIS")
    print("=" * 80)
    
    try:
        from analyze_integrated_datasets import main as integrated_main
        integrated_main()
        return True
    except Exception as e:
        print(f"[ERROR] Integrated analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Main function"""
    print("=" * 80)
    print("UPDATED PIPELINE EXECUTION")
    print("=" * 80)
    
    # Check prerequisites
    r_available = check_r_installation()
    
    # Run preprocessing (will handle R not being available)
    preprocessing_success = run_preprocessing()
    
    # Run integrated analysis
    analysis_success = run_integrated_analysis()
    
    # Summary
    print("\n" + "=" * 80)
    print("PIPELINE EXECUTION SUMMARY")
    print("=" * 80)
    
    print(f"\nPrerequisites:")
    print(f"  R available: {'Yes' if r_available else 'No (microarray normalization skipped)'}")
    
    print(f"\nSteps completed:")
    print(f"  Preprocessing: {'Success' if preprocessing_success else 'Skipped/Partial'}")
    print(f"  Integrated analysis: {'Success' if analysis_success else 'Failed'}")
    
    if analysis_success:
        print("\n[OK] Analysis complete! Check results in results/ directory:")
        print("  - results/gse122897_integrated_de.tsv")
        print("  - results/gse26969_integrated_de.tsv (if microarray data processed)")
        print("  - results/*_top_genes_integrated.tsv")
    else:
        print("\n[WARNING] Some steps failed. Check error messages above.")

if __name__ == "__main__":
    main()






