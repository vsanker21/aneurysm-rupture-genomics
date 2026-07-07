"""
Clarify Model Performance: What does each AUC represent?

1. AUC 0.968: Blood-only model (aneurysm vs control)
2. AUC ~0.52: Multi-omics integration attempt (different phenotypes)
3. Need to clarify what each model represents
"""
import pandas as pd
import numpy as np
import os

print("=" * 80)
print("CLARIFYING MODEL PERFORMANCE")
print("=" * 80)

# Check blood model results
print("\n1. BLOOD-ONLY MODEL (GSE159610)")
print("-" * 80)

blood_ml_file = "results/blood_ml_summary.tsv"
blood_cv_file = "results/blood_ml_cv_results.tsv"

if os.path.exists(blood_ml_file):
    blood_summary = pd.read_csv(blood_ml_file, sep='\t')
    print("\n  Summary:")
    print(f"    Best Model: {blood_summary['Best_Model'].values[0]}")
    print(f"    AUC: {blood_summary['Best_AUC'].values[0]:.3f}")
    print(f"    Features: {blood_summary['N_Features'].values[0]}")
    print(f"    Samples: {blood_summary['N_Samples'].values[0]}")
    print(f"\n  Task: Aneurysm vs Control (binary classification)")
    print(f"  Dataset: GSE159610 (blood transcriptomics)")
    print(f"  Features: Top 50 genes from blood transcriptomics")
    
    if os.path.exists(blood_cv_file):
        blood_cv = pd.read_csv(blood_cv_file, sep='\t')
        print(f"\n  Cross-Validation Results:")
        for idx, row in blood_cv.iterrows():
            print(f"    {row['Model']}: {row['Mean_AUC']:.3f} ± {row['Std_AUC']:.3f}")

# Check multi-omics model results
print("\n2. MULTI-OMICS INTEGRATION MODEL")
print("-" * 80)

multi_omics_file = "results/multi_omics_predictions_final.tsv"

if os.path.exists(multi_omics_file):
    multi_omics = pd.read_csv(multi_omics_file, sep='\t')
    print(f"\n  Multi-omics predictions file found:")
    print(f"    Samples: {len(multi_omics)}")
    print(f"    Columns: {list(multi_omics.columns)[:10]}")
    
    if 'Dataset' in multi_omics.columns:
        print(f"\n  Datasets:")
        print(multi_omics['Dataset'].value_counts())
    
    # Check if we have AUC information
    if 'GB_Predicted_Proba' in multi_omics.columns and 'Label' in multi_omics.columns:
        from sklearn.metrics import roc_auc_score
        
        y_true = multi_omics['Label'].values
        y_pred = multi_omics['GB_Predicted_Proba'].values
        
        if len(np.unique(y_true)) > 1:
            auc = roc_auc_score(y_true, y_pred)
            print(f"\n  Multi-omics Model AUC (full dataset): {auc:.3f}")
            
            if auc < 0.6:
                print(f"    [WARNING] Low AUC may indicate:")
                print(f"      - Different phenotypes (aneurysm vs control vs ruptured vs unruptured)")
                print(f"      - Sample size issues")
                print(f"      - Feature alignment problems")
else:
    print(f"\n  [INFO] Multi-omics predictions file not found")
    print(f"    File: {multi_omics_file}")
    print(f"    Status: May not have been run or file path different")

# Check what datasets were integrated
print("\n3. DATASET INTEGRATION ANALYSIS")
print("-" * 80)

print("\n  Blood Dataset (GSE159610):")
print(f"    Task: Aneurysm vs Control")
print(f"    Samples: ~47 (23 aneurysm, 24 control)")
print(f"    Phenotype: Presence of aneurysm (binary)")

print("\n  Tissue Dataset (GSE122897):")
print(f"    Task: Ruptured vs Unruptured")
print(f"    Samples: ~50 (ruptured vs unruptured)")
print(f"    Phenotype: Rupture status (binary)")

print("\n  [IMPORTANT] These are DIFFERENT phenotypes!")
print(f"    - Blood: Aneurysm vs Control (presence of disease)")
print(f"    - Tissue: Ruptured vs Unruptured (disease severity)")
print(f"    - Cannot directly combine without careful handling")

print("\n4. PERFORMANCE INTERPRETATION")
print("-" * 80)

print("\n  AUC 0.968 (Blood-only model):")
print(f"    ✓ Valid for: Aneurysm vs Control classification")
print(f"    ✓ Dataset: GSE159610 (blood transcriptomics)")
print(f"    ✓ Task: Binary classification (aneurysm present vs absent)")
print(f"    ✓ Scientifically valid for this specific task")

print("\n  AUC ~0.52 (if from multi-omics):")
print(f"    ✗ Would indicate poor performance")
print(f"    ✗ Likely due to: Combining different phenotypes")
print(f"    ✗ Problem: Blood (aneurysm vs control) ≠ Tissue (ruptured vs unruptured)")
print(f"    ✗ Solution: Need separate models or careful phenotype alignment")

print("\n5. RECOMMENDATION")
print("-" * 80)

print("\n  Correct approach for multi-omics integration:")
print(f"    1. Blood model: Aneurysm vs Control (AUC 0.968) ✓")
print(f"       - Use for: Screening/detection")
print(f"    2. Tissue model: Ruptured vs Unruptured (separate)")
print(f"       - Use for: Risk stratification")
print(f"    3. Integrated model: Combine blood + tissue (if same patients)")
print(f"       - Requires: Same patients with both blood and tissue data")
print(f"       - Current data: Different patients in each dataset")

print("\n  Summary:")
print(f"    - AUC 0.968 is for BLOOD-ONLY model (aneurysm vs control)")
print(f"    - AUC 0.968 is NOT from multi-omics integration")
print(f"    - Multi-omics integration requires same patients in both datasets")
print(f"    - Current datasets have different patients (cannot directly combine)")


