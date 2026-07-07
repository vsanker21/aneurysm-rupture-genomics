"""
Direct execution of comprehensive analysis pipeline
Ensures all outputs are generated successfully
"""
import sys
import os

# Force UTF-8 encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

print("=" * 80)
print("COMPREHENSIVE ANALYSIS PIPELINE - DIRECT EXECUTION")
print("=" * 80)
sys.stdout.flush()

try:
    # Import required packages
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.preprocessing import StandardScaler
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.svm import SVC
    from sklearn.model_selection import cross_val_predict, StratifiedKFold
    from sklearn.metrics import roc_curve, auc, roc_auc_score
    from imblearn.over_sampling import SMOTE
    print("✓ All imports successful")
    sys.stdout.flush()
except Exception as e:
    print(f"✗ Import error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Create output directories
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/gene_lists", exist_ok=True)
os.makedirs("results/risk_scores", exist_ok=True)
print("✓ Directories created")
sys.stdout.flush()

# Now execute the main pipeline script
print("\nExecuting main pipeline script...")
sys.stdout.flush()

try:
    # Read and execute the comprehensive pipeline
    with open('comprehensive_analysis_pipeline.py', 'r', encoding='utf-8', errors='replace') as f:
        pipeline_code = f.read()
    
    # Execute it
    exec(pipeline_code, {'__name__': '__main__', '__builtins__': __builtins__})
    
    print("\n" + "=" * 80)
    print("✅ PIPELINE EXECUTION COMPLETE")
    print("=" * 80)
    sys.stdout.flush()
    
except Exception as e:
    import traceback
    print(f"\n✗ Execution error: {e}")
    traceback.print_exc()
    sys.exit(1)

