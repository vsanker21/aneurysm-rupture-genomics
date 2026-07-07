"""
Run Full Multi-Omics Pipeline
Executes all scripts in sequence and reports results
"""
import subprocess
import sys
import os
import pandas as pd

def run_script(script_name, description):
    """Run a Python script and capture output"""
    print("\n" + "=" * 80)
    print(f"RUNNING: {description}")
    print("=" * 80)
    print(f"Script: {script_name}")
    
    try:
        result = subprocess.run(
            [sys.executable, script_name],
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        
        if result.returncode == 0:
            print(f"\n[OK] {description} completed successfully")
            return True
        else:
            print(f"\n[ERROR] {description} failed with return code {result.returncode}")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"\n[ERROR] {description} timed out")
        return False
    except Exception as e:
        print(f"\n[ERROR] {description} failed: {e}")
        return False

def check_results():
    """Check what results were generated"""
    print("\n" + "=" * 80)
    print("CHECKING RESULTS")
    print("=" * 80)
    
    results_files = {
        "results/integrated_features_full.tsv": "Integrated features",
        "results/full_multi_omics_results.tsv": "Model performance results",
        "results/feature_importance_shap.tsv": "Feature importance",
        "results/risk_stratification_full.tsv": "Risk stratification",
        "results/risk_category_statistics.tsv": "Risk category statistics"
    }
    
    for filepath, description in results_files.items():
        if os.path.exists(filepath):
            print(f"\n[OK] {description}: {filepath}")
            try:
                df = pd.read_csv(filepath, sep='\t', nrows=5)
                print(f"  Shape: {df.shape if hasattr(df, 'shape') else 'N/A'}")
                if 'Mean_AUC' in df.columns:
                    print(f"  Best AUC: {df['Mean_AUC'].max():.3f}")
            except:
                pass
        else:
            print(f"\n[SKIP] {description}: {filepath} not found")

def main():
    """Run full pipeline"""
    print("=" * 80)
    print("FULL MULTI-OMICS PIPELINE EXECUTION")
    print("=" * 80)
    
    scripts = [
        ("build_full_multi_omics_model.py", "Multi-Omics Integration"),
        ("shap_interpretation.py", "SHAP Interpretation"),
        ("risk_stratification.py", "Risk Stratification")
    ]
    
    results = {}
    for script, description in scripts:
        success = run_script(script, description)
        results[description] = success
    
    # Check results
    check_results()
    
    # Summary
    print("\n" + "=" * 80)
    print("PIPELINE EXECUTION SUMMARY")
    print("=" * 80)
    
    for description, success in results.items():
        status = "[OK]" if success else "[FAILED]"
        print(f"{status} {description}")
    
    all_success = all(results.values())
    if all_success:
        print("\n[SUCCESS] All scripts completed successfully!")
    else:
        print("\n[WARNING] Some scripts failed. Check output above.")

if __name__ == "__main__":
    main()


