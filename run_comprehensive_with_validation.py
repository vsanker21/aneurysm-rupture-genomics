"""
Run comprehensive integration with validation and error handling
Ensures all results are successfully generated
"""
import sys
import os
import traceback

# Force UTF-8
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

print("=" * 80)
print("COMPREHENSIVE MULTI-OMICS INTEGRATION - EXECUTION WITH VALIDATION")
print("=" * 80, flush=True)

try:
    # Execute the main script with proper encoding
    with open('execute_comprehensive_integration.py', 'r', encoding='utf-8', errors='replace') as f:
        exec(f.read())
    
    print("\n" + "=" * 80, flush=True)
    print("VALIDATING RESULTS...", flush=True)
    print("=" * 80, flush=True)
    
    # Check for result files
    expected_files = [
        'results/comprehensive_multi_omics_model_summary.tsv',
        'results/ensemble_tissue_predictions.tsv',
        'results/ensemble_blood_predictions.tsv',
        'results/common_enriched_pathways.tsv',
        'results/entrez_to_ensembl_mapping.tsv'
    ]
    
    print("\nChecking for result files:", flush=True)
    for file_path in expected_files:
        if os.path.exists(file_path):
            size = os.path.getsize(file_path)
            print(f"  ✓ {file_path} ({size} bytes)", flush=True)
        else:
            print(f"  ✗ {file_path} - MISSING", flush=True)
    
    print("\n" + "=" * 80, flush=True)
    print("EXECUTION COMPLETE", flush=True)
    print("=" * 80, flush=True)
    
except Exception as e:
    print(f"\nERROR: {e}", file=sys.stderr, flush=True)
    traceback.print_exc()
    sys.exit(1)

