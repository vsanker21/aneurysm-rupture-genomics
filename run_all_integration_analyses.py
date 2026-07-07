"""
Run All Integration Analyses
Master script to run all integration analyses in sequence
"""
import subprocess
import sys
import os

def run_script(script_name, description):
    """Run a Python script and report results"""
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
        
        if result.returncode == 0:
            print(f"\n✅ {description} - SUCCESS")
            if result.stdout:
                # Print last 20 lines of output
                lines = result.stdout.strip().split('\n')
                print("\nLast output lines:")
                for line in lines[-20:]:
                    print(f"  {line}")
            return True
        else:
            print(f"\n❌ {description} - FAILED")
            if result.stderr:
                print("\nError output:")
                print(result.stderr[:500])
            return False
    except subprocess.TimeoutExpired:
        print(f"\n⏱️ {description} - TIMEOUT (exceeded 10 minutes)")
        return False
    except Exception as e:
        print(f"\n❌ {description} - ERROR: {e}")
        return False

def main():
    """Run all integration analyses"""
    
    print("=" * 80)
    print("MULTI-OMICS INTEGRATION - COMPLETE PIPELINE")
    print("=" * 80)
    
    scripts = [
        ("build_polygenic_risk_score.py", "Phase 2: Polygenic Risk Score Development"),
        ("complete_blood_integration.py", "Phase 3: Blood Integration - Signature Finalization"),
        ("endothelial_dysfunction_analysis.py", "Phase 4: Endothelial Dysfunction Analysis"),
        ("multi_omics_integration.py", "Phase 5: Multi-Omics Integration")
    ]
    
    results = {}
    
    for script, description in scripts:
        if os.path.exists(script):
            success = run_script(script, description)
            results[script] = success
        else:
            print(f"\n⚠️ Script not found: {script}")
            results[script] = False
    
    # Summary
    print("\n" + "=" * 80)
    print("EXECUTION SUMMARY")
    print("=" * 80)
    
    for script, description in scripts:
        status = "✅ SUCCESS" if results.get(script, False) else "❌ FAILED/NOT RUN"
        print(f"{description}: {status}")
    
    successful = sum(1 for v in results.values() if v)
    total = len(results)
    
    print(f"\nOverall: {successful}/{total} scripts completed successfully")
    
    # Check for output files
    print("\n" + "=" * 80)
    print("OUTPUT FILES CHECK")
    print("=" * 80)
    
    expected_files = [
        "results/prs_risk_loci_top17.tsv",
        "results/prs_weights.tsv",
        "results/blood_signature_25genes.tsv",
        "results/blood_signature_30genes.tsv",
        "results/ec_dysfunction_scores.tsv",
        "results/integrated_features.tsv",
        "results/risk_stratification.tsv"
    ]
    
    for file_path in expected_files:
        if os.path.exists(file_path):
            size = os.path.getsize(file_path)
            print(f"✅ {file_path} ({size:,} bytes)")
        else:
            print(f"❌ {file_path} - NOT FOUND")
    
    print("\n" + "=" * 80)
    print("PIPELINE EXECUTION COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()






