"""
Test pipeline execution - minimal version to verify it works
"""
import sys
import os

# Force UTF-8
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

print("=" * 80, flush=True)
print("TESTING PIPELINE EXECUTION", flush=True)
print("=" * 80, flush=True)

# Test imports
try:
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.preprocessing import StandardScaler
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import SVC
    from sklearn.model_selection import cross_val_predict, StratifiedKFold
    from sklearn.metrics import roc_curve, auc, roc_auc_score
    from imblearn.over_sampling import SMOTE
    print("✓ All imports successful", flush=True)
except Exception as e:
    print(f"✗ Import error: {e}", flush=True)
    sys.exit(1)

# Test data loading
try:
    tissue_expr_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    if os.path.exists(tissue_expr_file):
        tissue_expr = pd.read_csv(tissue_expr_file, sep='\t', compression='gzip', index_col=0, nrows=10)
        print(f"✓ Tissue data accessible: {tissue_expr.shape}", flush=True)
    else:
        print("✗ Tissue data file not found", flush=True)
        sys.exit(1)
except Exception as e:
    print(f"✗ Data loading error: {e}", flush=True)
    sys.exit(1)

# Test directory creation
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/gene_lists", exist_ok=True)
os.makedirs("results/risk_scores", exist_ok=True)
print("✓ Directories created", flush=True)

# Test figure creation
try:
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot([0, 1], [0, 1], 'b-', label='Test')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Test Figure')
    ax.legend()
    plt.savefig('results/figures/test_figure.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("✓ Figure creation successful", flush=True)
except Exception as e:
    print(f"✗ Figure creation error: {e}", flush=True)
    sys.exit(1)

print("\n" + "=" * 80, flush=True)
print("✅ ALL TESTS PASSED - READY TO RUN FULL PIPELINE", flush=True)
print("=" * 80, flush=True)

# Now try to import and run the full pipeline
print("\nAttempting to run full pipeline...", flush=True)
try:
    # Read the full pipeline script
    with open('comprehensive_analysis_pipeline.py', 'r', encoding='utf-8', errors='replace') as f:
        code = f.read()
    
    # Execute it
    exec(code, {'__name__': '__main__'})
    print("\n✅ Full pipeline executed successfully!", flush=True)
except Exception as e:
    import traceback
    print(f"\n✗ Pipeline execution error: {e}", flush=True)
    traceback.print_exc()
    sys.exit(1)

