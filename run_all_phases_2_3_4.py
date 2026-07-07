"""
Run Phases 2, 3, and 4 in sequence
All phases using REAL data only
"""
import subprocess
import sys
import os

print("=" * 80)
print("RUNNING PHASES 2, 3, AND 4")
print("=" * 80)

scripts = [
    ("Phase 2: Polygenic Risk Score", "phase2_population_prs.py"),
    ("Phase 4: EC Dysfunction", "phase4_ec_dysfunction_bulk.py"),
    ("Phase 3: Blood Signature Validation", "phase3_validate_blood_signature.py")
]

for phase_name, script in scripts:
    print(f"\n{'='*80}")
    print(f"Running {phase_name}...")
    print(f"{'='*80}")
    
    if not os.path.exists(script):
        print(f"✗ Error: Script not found: {script}")
        continue
    
    try:
        result = subprocess.run(
            [sys.executable, script],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='replace',
            timeout=600  # 10 minutes per script
        )
        
        if result.returncode == 0:
            print(f"✅ {phase_name} completed successfully")
            # Print last 20 lines of output
            if result.stdout:
                lines = result.stdout.strip().split('\n')
                for line in lines[-20:]:
                    print(line)
        else:
            print(f"✗ {phase_name} failed with return code {result.returncode}")
            if result.stderr:
                print("Error output:")
                print(result.stderr[-500:])
    except subprocess.TimeoutExpired:
        print(f"✗ {phase_name} timed out after 10 minutes")
    except Exception as e:
        print(f"✗ Error running {phase_name}: {e}")

print("\n" + "=" * 80)
print("PHASES 2, 3, AND 4 EXECUTION COMPLETE")
print("=" * 80)

