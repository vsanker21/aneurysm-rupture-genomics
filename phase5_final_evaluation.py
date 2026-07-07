"""Final evaluation of optimized model"""
import pandas as pd
import numpy as np
import pickle
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, auc, confusion_matrix, classification_report
import json

# Load data
features = pd.read_csv('results/integrated_features.tsv', sep='\t')
feature_cols = ['Blood_RNA_Score', 'EC_Dysfunction_Score', 'Inflammation_Score', 
                'ECM_Remodeling_Score', 'DE_Signature_Score', 'Clinical_Score']

features = features[features['Group'].isin(['Ruptured', 'Unruptured'])].copy()
labels = (features['Group'] == 'Ruptured').astype(int)

# Load best model summary
with open('results/phase5_best_model_summary.json', 'r') as f:
    best_summary = json.load(f)

best_features = best_summary['features_used']
feature_indices = [feature_cols.index(f) for f in best_features]

# Prepare data with selected features
X = features[best_features].fillna(0).values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Load optimized model
with open('results/phase5_optimized_model.pkl', 'rb') as f:
    model = pickle.load(f)

# Generate predictions
predictions = model.predict_proba(X_scaled)[:, 1]
predicted_labels = (predictions >= 0.5).astype(int)

# Calculate metrics
roc_auc = roc_auc_score(labels, predictions)
fpr, tpr, _ = roc_curve(labels, predictions)
precision, recall, _ = precision_recall_curve(labels, predictions)
pr_auc = auc(recall, precision)
cm = confusion_matrix(labels, predicted_labels)

# Save results
roc_data = pd.DataFrame({'FPR': fpr, 'TPR': tpr})
roc_data.to_csv('results/phase5_optimized_roc_curve.tsv', sep='\t', index=False)

pr_data = pd.DataFrame({'Recall': recall, 'Precision': precision})
pr_data.to_csv('results/phase5_optimized_pr_curve.tsv', sep='\t', index=False)

predictions_df = pd.DataFrame({
    'Sample': features.index,
    'True_Label': labels,
    'Predicted_Probability': predictions,
    'Predicted_Label': predicted_labels
})
predictions_df.to_csv('results/phase5_optimized_predictions.tsv', sep='\t', index=False)

metrics_df = pd.DataFrame({
    'Metric': ['ROC_AUC', 'PR_AUC', 'True_Negative', 'False_Positive', 
               'False_Negative', 'True_Positive', 'Accuracy', 'Precision', 'Recall', 'F1_Score'],
    'Value': [
        roc_auc,
        pr_auc,
        float(cm[0, 0]),
        float(cm[0, 1]),
        float(cm[1, 0]),
        float(cm[1, 1]),
        float((cm[0, 0] + cm[1, 1]) / cm.sum()),
        float(cm[1, 1] / (cm[1, 1] + cm[0, 1])) if (cm[1, 1] + cm[0, 1]) > 0 else 0,
        float(cm[1, 1] / (cm[1, 1] + cm[1, 0])) if (cm[1, 1] + cm[1, 0]) > 0 else 0,
        float(2 * cm[1, 1] / (2 * cm[1, 1] + cm[0, 1] + cm[1, 0])) if (2 * cm[1, 1] + cm[0, 1] + cm[1, 0]) > 0 else 0
    ]
})
metrics_df.to_csv('results/phase5_optimized_metrics.tsv', sep='\t', index=False)

print("=" * 80)
print("OPTIMIZED MODEL FINAL EVALUATION")
print("=" * 80)
print(f"\nModel: {best_summary['best_model']}")
print(f"Features: {', '.join(best_features)}")
print(f"\nTraining Metrics:")
print(f"  ROC AUC: {best_summary['train_auc']:.3f}")
print(f"  CV AUC: {best_summary['cv_auc_mean']:.3f} ± {best_summary['cv_auc_std']:.3f}")
print(f"  Overfitting Gap: {best_summary['overfitting_gap']:.3f}")
print(f"\nFull Dataset Metrics:")
print(f"  ROC AUC: {roc_auc:.3f}")
print(f"  PR AUC: {pr_auc:.3f}")
print(f"  Accuracy: {metrics_df[metrics_df['Metric']=='Accuracy']['Value'].values[0]:.3f}")
print(f"\nConfusion Matrix:")
print(f"  True Negative:  {int(cm[0, 0])}")
print(f"  False Positive: {int(cm[0, 1])}")
print(f"  False Negative: {int(cm[1, 0])}")
print(f"  True Positive:  {int(cm[1, 1])}")
print(f"\nImprovement over original:")
print(f"  Original CV AUC: 0.425")
print(f"  Optimized CV AUC: {best_summary['cv_auc_mean']:.3f}")
print(f"  Improvement: +{best_summary['improvement']:.3f} ({best_summary['improvement']/0.425*100:.1f}%)")
print("\n[OK] Final evaluation complete. Files saved.")
print("=" * 80)
