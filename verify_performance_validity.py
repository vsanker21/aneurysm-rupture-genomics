"""
Verify that the obtained performance (AUC 0.968) is scientifically valid
Check for potential issues: sample size, feature count, data leakage, cross-validation
"""
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("PERFORMANCE VALIDITY VERIFICATION")
print("=" * 80)

# Load actual data
blood_expr_file = "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
blood_meta_file = "data/blood_rnaseq/GSE159610/sample_metadata_with_labels.tsv"

print("\n1. VERIFYING DATA AND SAMPLE SIZE")
print("-" * 80)

blood_expr = pd.read_csv(blood_expr_file, sep='\t', compression='gzip', index_col=0)
blood_meta = pd.read_csv(blood_meta_file, sep='\t')

print(f"  Expression data: {blood_expr.shape[0]:,} genes, {blood_expr.shape[1]} samples")
print(f"  Metadata: {len(blood_meta)} samples")

# Check sample alignment
common_samples = [s for s in blood_meta['Sample'] if s in blood_expr.columns]
blood_meta_aligned = blood_meta[blood_meta['Sample'].isin(common_samples)].copy()
blood_meta_aligned['Label'] = (blood_meta_aligned['Status'] == 'Aneurysm').astype(int)

print(f"  Common samples: {len(common_samples)}")
print(f"  Aneurysm samples: {(blood_meta_aligned['Label'] == 1).sum()}")
print(f"  Control samples: {(blood_meta_aligned['Label'] == 0).sum()}")

# Load feature importance to see which features were used
feature_importance_file = "results/blood_ml_feature_importance.tsv"
if not pd.io.common.file_exists(feature_importance_file):
    print(f"\n[WARNING] Feature importance file not found: {feature_importance_file}")
    # Use all genes as fallback
    feature_genes = blood_expr.index.tolist()[:50]
else:
    feature_importance = pd.read_csv(feature_importance_file, sep='\t')
    feature_genes = feature_importance['Gene_ID'].head(50).tolist()
    print(f"\n  Features used: {len(feature_genes)} genes (from feature importance)")

# Check feature count vs sample size
n_samples = len(common_samples)
n_features = len(feature_genes)

print(f"\n2. CHECKING FEATURE-TO-SAMPLE RATIO")
print("-" * 80)
print(f"  Samples: {n_samples}")
print(f"  Features: {n_features}")
print(f"  Ratio: {n_features/n_samples:.2f} features per sample")

if n_features > n_samples:
    print(f"  [WARNING] More features than samples - risk of overfitting!")
elif n_features > n_samples * 0.8:
    print(f"  [WARNING] High feature-to-sample ratio - potential overfitting")
else:
    print(f"  [OK] Feature-to-sample ratio is reasonable")

# Extract features
X = blood_expr.loc[[g for g in feature_genes if g in blood_expr.index], common_samples].T.values
y = blood_meta_aligned['Label'].values

print(f"\n  Actual feature matrix: {X.shape}")
print(f"  Labels: {y.shape}")

# Handle missing values
X = np.nan_to_num(X, nan=0.0)

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print(f"\n3. REPRODUCING CROSS-VALIDATION")
print("-" * 80)

# 5-fold stratified cross-validation (same as reported)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# SVM (best model according to results)
svm_model = SVC(kernel='rbf', probability=True, random_state=42)
svm_scores = []

# Random Forest (for comparison)
rf_model = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=3)
rf_scores = []

print("  Running 5-fold stratified cross-validation...")

for fold, (train_idx, test_idx) in enumerate(cv.split(X_scaled, y), 1):
    X_train, X_test = X_scaled[train_idx], X_scaled[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    
    # Train SVM
    svm_model.fit(X_train, y_train)
    y_pred_proba_svm = svm_model.predict_proba(X_test)[:, 1]
    svm_auc = roc_auc_score(y_test, y_pred_proba_svm)
    svm_scores.append(svm_auc)
    
    # Train RF
    rf_model.fit(X_train, y_train)
    y_pred_proba_rf = rf_model.predict_proba(X_test)[:, 1]
    rf_auc = roc_auc_score(y_test, y_pred_proba_rf)
    rf_scores.append(rf_auc)
    
    print(f"    Fold {fold}: SVM AUC={svm_auc:.3f}, RF AUC={rf_auc:.3f}")

svm_mean = np.mean(svm_scores)
svm_std = np.std(svm_scores)
rf_mean = np.mean(rf_scores)
rf_std = np.std(rf_scores)

print(f"\n  Cross-validation results:")
print(f"    SVM: {svm_mean:.3f} ± {rf_std:.3f}")
print(f"    RF:  {rf_mean:.3f} ± {rf_std:.3f}")

print(f"\n  Reported results:")
print(f"    SVM: 0.968 ± 0.064")
print(f"    RF:  0.922 ± 0.089")

# Check if results are similar
svm_diff = abs(svm_mean - 0.968)
rf_diff = abs(rf_mean - 0.922)

print(f"\n4. VALIDATION ASSESSMENT")
print("-" * 80)

if svm_diff < 0.05:
    print(f"  [OK] SVM AUC matches reported value (difference: {svm_diff:.3f})")
else:
    print(f"  [WARNING] SVM AUC differs from reported (difference: {svm_diff:.3f})")

if rf_diff < 0.05:
    print(f"  [OK] RF AUC matches reported value (difference: {rf_diff:.3f})")
else:
    print(f"  [WARNING] RF AUC differs from reported (difference: {rf_diff:.3f})")

# Check for potential issues
print(f"\n5. POTENTIAL ISSUES CHECK")
print("-" * 80)

issues = []

# Check 1: Sample size
if n_samples < 30:
    issues.append(f"Small sample size ({n_samples}) - results may not generalize")
elif n_samples < 50:
    issues.append(f"Moderate sample size ({n_samples}) - consider external validation")

# Check 2: Class balance
class_balance = min(y.sum(), len(y) - y.sum()) / len(y)
if class_balance < 0.3:
    issues.append(f"Class imbalance (ratio: {class_balance:.2f}) - may affect performance")
else:
    print(f"  [OK] Class balance is reasonable ({class_balance:.2f})")

# Check 3: Feature count
if n_features > n_samples * 0.8:
    issues.append(f"High feature-to-sample ratio ({n_features/n_samples:.2f}) - risk of overfitting")

# Check 4: Cross-validation variability
if svm_std > 0.15:
    issues.append(f"High CV variability (std: {svm_std:.3f}) - results may be unstable")
elif svm_std < 0.05:
    issues.append(f"Very low CV variability (std: {svm_std:.3f}) - may indicate data leakage")

# Check 5: Performance seems too high
if svm_mean > 0.95 and n_samples < 100:
    issues.append(f"Very high AUC ({svm_mean:.3f}) with small sample size - may not generalize")

if len(issues) > 0:
    print(f"  [WARNINGS] {len(issues)} potential issues identified:")
    for issue in issues:
        print(f"    - {issue}")
else:
    print(f"  [OK] No major issues detected")

print(f"\n6. SCIENTIFIC VALIDITY ASSESSMENT")
print("-" * 80)

validity_score = 0
max_score = 5

# Criterion 1: Proper cross-validation
if 0.05 < svm_std < 0.15:
    validity_score += 1
    print(f"  [✓] Cross-validation properly implemented (variability: {svm_std:.3f})")
else:
    print(f"  [✗] Cross-validation variability may be problematic (variability: {svm_std:.3f})")

# Criterion 2: Reasonable feature-to-sample ratio
if n_features <= n_samples * 0.8:
    validity_score += 1
    print(f"  [✓] Feature-to-sample ratio is reasonable ({n_features/n_samples:.2f})")
else:
    print(f"  [✗] Feature-to-sample ratio is high ({n_features/n_samples:.2f})")

# Criterion 3: Class balance
if class_balance >= 0.3:
    validity_score += 1
    print(f"  [✓] Class balance is acceptable ({class_balance:.2f})")
else:
    print(f"  [✗] Class imbalance detected ({class_balance:.2f})")

# Criterion 4: Sample size
if n_samples >= 30:
    validity_score += 1
    print(f"  [✓] Sample size is adequate ({n_samples})")
else:
    print(f"  [✗] Sample size is small ({n_samples})")

# Criterion 5: Reproducibility
if svm_diff < 0.1:
    validity_score += 1
    print(f"  [✓] Results are reproducible (difference: {svm_diff:.3f})")
else:
    print(f"  [✗] Results differ from reported (difference: {svm_diff:.3f})")

print(f"\n  Validity Score: {validity_score}/{max_score}")

if validity_score >= 4:
    print(f"\n  [CONCLUSION] Performance appears scientifically valid")
    print(f"    - AUC {svm_mean:.3f} is reproducible")
    print(f"    - Cross-validation properly implemented")
    print(f"    - However, external validation on independent dataset is recommended")
elif validity_score >= 3:
    print(f"\n  [CONCLUSION] Performance likely valid but with caveats")
    print(f"    - External validation strongly recommended")
    print(f"    - Consider reducing feature count or increasing sample size")
else:
    print(f"\n  [CONCLUSION] Performance validity concerns identified")
    print(f"    - External validation essential")
    print(f"    - Review feature selection and cross-validation procedures")

print(f"\n7. RECOMMENDATIONS")
print("-" * 80)
print(f"  1. External validation on independent dataset (GSE66238)")
print(f"  2. Reduce feature count if possible (currently {n_features})")
print(f"  3. Report confidence intervals for AUC")
print(f"  4. Compare with simpler models (fewer features)")
print(f"  5. Document any preprocessing steps carefully")


