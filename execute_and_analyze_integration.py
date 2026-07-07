"""
Execute integration and analyze results
"""
import sys
import os
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve, classification_report

# Force UTF-8 encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

print("=" * 80)
print("EXECUTING MULTI-OMICS INTEGRATION")
print("=" * 80)

# Import and run integration
try:
    from integrate_real_data_multi_omics import main
    print("\nRunning integration pipeline...")
    main()
    print("\n[OK] Integration script completed")
except Exception as e:
    print(f"\n[ERROR] Integration failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Analyze results
print("\n" + "=" * 80)
print("ANALYZING INTEGRATION RESULTS")
print("=" * 80)

# Check for output files
output_files = {
    'predictions': 'results/multi_omics_predictions_real_data.tsv',
    'features': 'results/multi_omics_features_real_data.tsv',
    'signature': 'results/multi_omics_signature_genes.tsv',
    'importance': 'results/multi_omics_feature_importance_real_data.tsv'
}

for name, path in output_files.items():
    if os.path.exists(path):
        print(f"\n[OK] {name.capitalize()} file exists: {path}")
        if name == 'predictions':
            df = pd.read_csv(path, sep='\t')
            print(f"  Samples: {len(df)}")
            if 'Label' in df.columns:
                print(f"  Positive: {df['Label'].sum()}")
                print(f"  Negative: {(df['Label'] == 0).sum()}")
            if 'GB_Predicted_Proba' in df.columns:
                if len(df['Label'].unique()) > 1:
                    auc = roc_auc_score(df['Label'], df['GB_Predicted_Proba'])
                    print(f"  GB AUC: {auc:.3f}")
        elif name == 'importance':
            df = pd.read_csv(path, sep='\t')
            print(f"  Features: {len(df)}")
            print("\n  Top 10 most important features:")
            for idx, row in df.head(10).iterrows():
                print(f"    {row['Feature']:<40s} GB: {row['GB_Importance']:.4f}")
        elif name == 'signature':
            df = pd.read_csv(path, sep='\t')
            print(f"  Signature genes: {len(df)}")
            if 'Source' in df.columns:
                print(df['Source'].value_counts().to_string())
    else:
        print(f"\n[WARNING] {name.capitalize()} file not found: {path}")

# Check existing results
print("\n" + "=" * 80)
print("CHECKING EXISTING RESULTS")
print("=" * 80)

existing_results = {
    'multi_omics_prioritized': 'results/multi_omics_prioritized_top50.tsv',
    'blood_tissue_integration': 'results/blood_tissue_integration_fdr_0.05.tsv',
    'full_multi_omics': 'results/full_multi_omics_results.tsv'
}

for name, path in existing_results.items():
    if os.path.exists(path):
        print(f"\n[OK] Existing {name} results found")
        df = pd.read_csv(path, sep='\t')
        print(f"  Rows: {len(df)}")
        print(f"  Columns: {list(df.columns[:10])}")

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)


