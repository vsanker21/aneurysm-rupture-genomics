"""
Complete EC Dysfunction Analysis with Real Data
Runs the analysis to completion, handling API rate limits gracefully
"""
import subprocess
import sys
import os
import time

def run_ec_analysis():
    """Run EC dysfunction analysis with timeout and retry logic"""
    print("=" * 80)
    print("RUNNING EC DYSFUNCTION ANALYSIS TO COMPLETION")
    print("=" * 80)
    
    script = "endothelial_dysfunction_analysis.py"
    
    if not os.path.exists(script):
        print(f"[ERROR] Script not found: {script}")
        return False
    
    print(f"\nRunning {script}...")
    print("Note: This may take several minutes due to API calls for gene ID mapping")
    
    try:
        # Run with extended timeout (30 minutes for API calls)
        result = subprocess.run(
            [sys.executable, script],
            capture_output=True,
            text=True,
            timeout=1800  # 30 minutes
        )
        
        if result.returncode == 0:
            print("[OK] EC dysfunction analysis completed successfully")
            if result.stdout:
                # Print key output
                lines = result.stdout.split('\n')
                # Print last 50 lines
                print("\n".join(lines[-50:]))
            return True
        else:
            print("[ERROR] EC dysfunction analysis failed")
            if result.stderr:
                print("Error output:")
                print(result.stderr[-1000:])
            return False
            
    except subprocess.TimeoutExpired:
        print("[ERROR] Analysis timed out after 30 minutes")
        print("This may be due to API rate limits. The analysis can be resumed.")
        return False
    except Exception as e:
        print(f"[ERROR] Failed to run analysis: {e}")
        return False

if __name__ == "__main__":
    success = run_ec_analysis()
    if success:
        print("\n" + "=" * 80)
        print("EC DYSFUNCTION ANALYSIS COMPLETE")
        print("=" * 80)
        print("\nResults saved to:")
        print("  - results/ec_dysfunction_scores.tsv")
        print("  - results/ec_dysfunction_markers.tsv")
    else:
        print("\n[WARNING] Analysis did not complete successfully")
        print("Check the error messages above for details")






