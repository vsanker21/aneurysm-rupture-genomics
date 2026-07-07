"""
Simplified script to run all analyses 3-8 quickly
Runs each analysis independently with immediate file saving
"""
import sys
import os
sys.stdout.reconfigure(encoding='utf-8')  # Ensure UTF-8 output

# Run each analysis as a separate import/execution
print("=" * 80)
print("RUNNING ALL ANALYSES 3-8 (SIMPLIFIED FAST VERSIONS)")
print("=" * 80)

# Import the main script
exec(open('run_analyses_3_to_8_fast.py', encoding='utf-8').read())

print("\n" + "=" * 80)
print("SCRIPT EXECUTION COMPLETE")
print("=" * 80)

