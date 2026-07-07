"""
Phase 7: Validation and Comparison (REAL DATA ONLY)
Comprehensive evaluation using only real data - NO SIMULATIONS
"""
import pandas as pd
import numpy as np
import pickle
import json
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score, roc_curve, precision_recall_curve, auc,
    confusion_matrix, classification_report, 
    accuracy_score, precision_score, recall_score, f1_score
)

print("=" * 80)
print("PHASE 7: VALIDATION AND COMPARISON (REAL DATA ONLY)")
print("NO SIMULATED DATA - USING ONLY ACTUAL DATA FROM DATASETS")
print("=" * 80)

# Load data
features_df = pd.read_csv('results/integrated_features.tsv', sep='\t')
features_df = features_df[features_df['Group'].isin(['Ruptured', 'Unruptured'])].copy()
labels = (features_df['Group'] == 'Ruptured').astype(int)

# Load optimized model
with open('results/phase5_optimized_model.pkl', 'rb') as f:
    optimized_model = pickle.load(f)

with open('results/phase5_optimized_scaler.pkl', 'rb') as f:
    scaler = pickle.load(f)

with open('results/phase5_best_model_summary.json', 'r') as f:
    model_summary = json.load(f)

optimized_features = model_summary['features_used']
X_optimized = features_df[optimized_features].fillna(0).values
X_optimized_scaled = scaler.transform(X_optimized)

# Get optimized model predictions
optimized_predictions = optimized_model.predict_proba(X_optimized_scaled)[:, 1]

print("\n" + "=" * 80)
print("1. OPTIMIZED MODEL EVALUATION (REAL DATA)")
print("=" * 80)

# Evaluate optimized model on real data
optimized_auc = roc_auc_score(labels, optimized_predictions)
fpr_opt, tpr_opt, _ = roc_curve(labels, optimized_predictions)
precision_opt, recall_opt, _ = precision_recall_curve(labels, optimized_predictions)
pr_auc_opt = auc(recall_opt, precision_opt)

# Optimal threshold
def find_optimal_threshold(y_true, y_scores):
    fpr, tpr, thresholds = roc_curve(y_true, y_scores)
    youden = tpr - fpr
    optimal_idx = np.argmax(youden)
    return thresholds[optimal_idx]

opt_threshold = find_optimal_threshold(labels, optimized_predictions)
opt_pred_binary = (optimized_predictions >= opt_threshold).astype(int)

print(f"\nOptimized Model Performance:")
print(f"  ROC AUC: {optimized_auc:.3f}")
print(f"  PR AUC: {pr_auc_opt:.3f}")
print(f"  Accuracy: {accuracy_score(labels, opt_pred_binary):.3f}")
print(f"  Precision: {precision_score(labels, opt_pred_binary, zero_division=0):.3f}")
print(f"  Recall: {recall_score(labels, opt_pred_binary, zero_division=0):.3f}")
print(f"  F1-Score: {f1_score(labels, opt_pred_binary, zero_division=0):.3f}")
print(f"  Optimal Threshold: {opt_threshold:.3f}")

cm_opt = confusion_matrix(labels, opt_pred_binary)
print(f"\nConfusion Matrix:")
print(f"  True Negative:  {cm_opt[0, 0]:3d} | False Positive: {cm_opt[0, 1]:3d}")
print(f"  False Negative: {cm_opt[1, 0]:3d} | True Positive:  {cm_opt[1, 1]:3d}")

print("\n" + "=" * 80)
print("2. BASELINE COMPARISON (USING AVAILABLE REAL FEATURES)")
print("=" * 80)

# Instead of simulating PHASES, compare with baseline models using real features
# Use Clinical_Score as a simple baseline since it's available in the data
if 'Clinical_Score' in features_df.columns:
    clinical_scores = features_df['Clinical_Score'].fillna(0).values
    clinical_scores_normalized = (clinical_scores - clinical_scores.min()) / (clinical_scores.max() - clinical_scores.min() + 1e-10)
    clinical_auc = roc_auc_score(labels, clinical_scores_normalized)
    
    print(f"\nClinical Score Baseline (Real Data):")
    print(f"  ROC AUC: {clinical_auc:.3f}")
else:
    # If Clinical_Score not available, use DE_Signature_Score as single-feature baseline
    de_scores = features_df['DE_Signature_Score'].fillna(0).values
    de_scores_normalized = (de_scores - de_scores.min()) / (de_scores.max() - de_scores.min() + 1e-10)
    clinical_auc = roc_auc_score(labels, de_scores_normalized)
    
    print(f"\nDE_Signature_Score Single-Feature Baseline (Real Data):")
    print(f"  ROC AUC: {clinical_auc:.3f}")

improvement = optimized_auc - clinical_auc
relative_improvement = (improvement / clinical_auc) * 100 if clinical_auc > 0 else 0

print(f"\nPerformance Comparison (Real Data):")
print(f"  Optimized Model AUC: {optimized_auc:.3f}")
print(f"  Baseline AUC:        {clinical_auc:.3f}")
print(f"  Improvement:         +{improvement:.3f} AUC points")
print(f"  Relative Improvement: +{relative_improvement:.1f}%")

print("\n" + "=" * 80)
print("3. DETAILED METRICS (REAL DATA)")
print("=" * 80)

# Detailed metrics for optimized model only (we don't have real PHASES data)
metrics_df = pd.DataFrame({
    'Metric': ['ROC_AUC', 'PR_AUC', 'Accuracy', 'Precision', 'Recall', 'F1_Score', 
               'True_Negative', 'False_Positive', 'False_Negative', 'True_Positive'],
    'Optimized_Model': [
        optimized_auc,
        pr_auc_opt,
        accuracy_score(labels, opt_pred_binary),
        precision_score(labels, opt_pred_binary, zero_division=0),
        recall_score(labels, opt_pred_binary, zero_division=0),
        f1_score(labels, opt_pred_binary, zero_division=0),
        float(cm_opt[0, 0]),
        float(cm_opt[0, 1]),
        float(cm_opt[1, 0]),
        float(cm_opt[1, 1])
    ]
})

print("\nOptimized Model Metrics:")
print(metrics_df.to_string(index=False))

metrics_df.to_csv('results/phase7_optimized_model_metrics.tsv', sep='\t', index=False)

print("\n" + "=" * 80)
print("4. ROC AND PR CURVES (REAL DATA)")
print("=" * 80)

# Save ROC and PR curves
roc_data = pd.DataFrame({
    'FPR': fpr_opt,
    'TPR': tpr_opt
})
roc_data.to_csv('results/phase7_roc_curve_optimized.tsv', sep='\t', index=False)

pr_data = pd.DataFrame({
    'Recall': recall_opt,
    'Precision': precision_opt
})
pr_data.to_csv('results/phase7_pr_curve_optimized.tsv', sep='\t', index=False)

print(f"\nROC AUC: {optimized_auc:.3f}")
print(f"PR AUC: {pr_auc_opt:.3f}")
print(f"[OK] ROC and PR curve data saved")

print("\n" + "=" * 80)
print("5. RISK STRATIFICATION (REAL DATA)")
print("=" * 80)

# Risk stratification
opt_risk_categories = pd.cut(optimized_predictions, bins=[0, 0.33, 0.67, 1.0], 
                             labels=['Low', 'Medium', 'High'])

risk_strat = pd.DataFrame({
    'Sample': features_df.index,
    'True_Label': labels,
    'Optimized_Score': optimized_predictions,
    'Optimized_Risk': opt_risk_categories
})

risk_strat.to_csv('results/phase7_risk_stratification.tsv', sep='\t', index=False)

print("\nOptimized Model Risk Distribution:")
print(opt_risk_categories.value_counts().sort_index())

# Calculate rupture rates by risk category
print("\nRupture Rates by Risk Category:")
for risk_level in ['Low', 'Medium', 'High']:
    mask = opt_risk_categories == risk_level
    if mask.sum() > 0:
        rupture_rate = labels[mask].mean()
        n_samples = mask.sum()
        print(f"  {risk_level:8s} Risk: {rupture_rate:.1%} ({labels[mask].sum()}/{n_samples})")

print("\n" + "=" * 80)
print("6. COMPREHENSIVE SUMMARY (REAL DATA ONLY)")
print("=" * 80)

summary = {
    'dataset': {
        'n_samples': len(features_df),
        'n_ruptured': int(labels.sum()),
        'n_unruptured': int((labels == 0).sum()),
        'rupture_rate': float(labels.mean())
    },
    'optimized_model': {
        'model_type': model_summary['best_model'],
        'features': optimized_features,
        'auc': float(optimized_auc),
        'pr_auc': float(pr_auc_opt),
        'accuracy': float(accuracy_score(labels, opt_pred_binary)),
        'precision': float(precision_score(labels, opt_pred_binary, zero_division=0)),
        'recall': float(recall_score(labels, opt_pred_binary, zero_division=0)),
        'f1': float(f1_score(labels, opt_pred_binary, zero_division=0)),
        'optimal_threshold': float(opt_threshold)
    },
    'baseline_comparison': {
        'baseline_type': 'Clinical_Score' if 'Clinical_Score' in features_df.columns else 'DE_Signature_Score',
        'baseline_auc': float(clinical_auc),
        'optimized_auc': float(optimized_auc),
        'improvement': float(improvement),
        'relative_improvement_percent': float(relative_improvement)
    },
    'note': 'All results based on REAL DATA from datasets. No simulated or synthetic data used.'
}

with open('results/phase7_comprehensive_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print("\n[OK] All results saved (REAL DATA ONLY):")
print("  - phase7_optimized_model_metrics.tsv: Detailed metrics")
print("  - phase7_roc_curve_optimized.tsv: ROC curve data")
print("  - phase7_pr_curve_optimized.tsv: PR curve data")
print("  - phase7_risk_stratification.tsv: Risk stratification")
print("  - phase7_comprehensive_summary.json: Complete summary")

print("\n" + "=" * 80)
print("KEY FINDINGS (REAL DATA ONLY)")
print("=" * 80)
print(f"\nOptimized Model AUC: {optimized_auc:.3f}")
print(f"Baseline AUC:        {clinical_auc:.3f}")
print(f"Improvement:         +{improvement:.3f} AUC points (+{relative_improvement:.1f}%)")
print(f"\nNOTE: All results are based on REAL DATA from the actual datasets.")
print(f"      No simulated, synthetic, or randomly generated data was used.")

print("\n" + "=" * 80)
print("PHASE 7 COMPLETE (REAL DATA ONLY)")
print("=" * 80)
