"""Phase 5: SHAP interpretation and model analysis"""
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.inspection import permutation_importance
import pickle
import json

# Load data
features = pd.read_csv('results/integrated_features.tsv', sep='\t')
feature_cols = ['Blood_RNA_Score', 'EC_Dysfunction_Score', 'Inflammation_Score', 
                'ECM_Remodeling_Score', 'DE_Signature_Score', 'Clinical_Score']
features = features[features['Group'].isin(['Ruptured', 'Unruptured'])].copy()
labels = (features['Group'] == 'Ruptured').astype(int)
X = features[feature_cols].fillna(0).values

# Scale features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Load or train model
model = GradientBoostingClassifier(n_estimators=200, max_depth=5, 
                                   learning_rate=0.05, random_state=42)
model.fit(X_scaled, labels)

import sys
sys.stdout.flush()
sys.stderr.flush()

print("=" * 80, flush=True)
print("PHASE 5: MODEL INTERPRETATION AND SHAP ANALYSIS", flush=True)
print("=" * 80, flush=True)

# Try SHAP if available, otherwise use permutation importance
try:
    import shap
    print("\n[OK] SHAP available - calculating SHAP values...")
    
    # Create SHAP explainer
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_scaled)
    
    # Handle binary classification (SHAP returns array for each class)
    if isinstance(shap_values, list):
        shap_values = shap_values[1]  # Use positive class
    
    # Calculate mean absolute SHAP values per feature
    mean_abs_shap = np.abs(shap_values).mean(axis=0)
    
    shap_importance = pd.DataFrame({
        'Feature': feature_cols,
        'Mean_ABS_SHAP': mean_abs_shap,
        'Percent_Contribution': mean_abs_shap / mean_abs_shap.sum() * 100
    }).sort_values('Mean_ABS_SHAP', ascending=False)
    
    shap_importance.to_csv('results/phase5_shap_importance.tsv', sep='\t', index=False)
    
    print("\nFeature Contribution (SHAP values):")
    for _, row in shap_importance.iterrows():
        print(f"  {row['Feature']:25s}: {row['Mean_ABS_SHAP']:.4f} ({row['Percent_Contribution']:.1f}%)")
    
    # Save SHAP values for all samples
    shap_df = pd.DataFrame(shap_values, columns=feature_cols, index=features.index)
    shap_df['Sample'] = features.index
    shap_df['True_Label'] = labels.values
    shap_df = shap_df[['Sample', 'True_Label'] + feature_cols]
    shap_df.to_csv('results/phase5_shap_values.tsv', sep='\t', index=False)
    
    print("\n[OK] SHAP values saved to results/phase5_shap_values.tsv")
    
except ImportError as e:
    print(f"\n[WARNING] SHAP not available - using permutation importance instead...")
    print(f"ImportError: {e}", flush=True)
    
    # Use permutation importance as alternative
    try:
        perm_importance = permutation_importance(model, X_scaled, labels, 
                                                n_repeats=10, random_state=42, 
                                                scoring='roc_auc')
    except Exception as e2:
        print(f"Error in permutation importance: {e2}", flush=True)
        raise
    
    perm_df = pd.DataFrame({
        'Feature': feature_cols,
        'Importance_Mean': perm_importance.importances_mean,
        'Importance_Std': perm_importance.importances_std,
        'Percent_Contribution': perm_importance.importances_mean / 
                               perm_importance.importances_mean.sum() * 100
    }).sort_values('Importance_Mean', ascending=False)
    
    perm_df.to_csv('results/phase5_permutation_importance.tsv', sep='\t', index=False)
    
    print("\nFeature Contribution (Permutation Importance):")
    for _, row in perm_df.iterrows():
        print(f"  {row['Feature']:25s}: {row['Importance_Mean']:.4f} ± {row['Importance_Std']:.4f} ({row['Percent_Contribution']:.1f}%)")

# Save model and scaler for later use
with open('results/phase5_model.pkl', 'wb') as f:
    pickle.dump(model, f)
    
with open('results/phase5_scaler.pkl', 'wb') as f:
    pickle.dump(scaler, f)

# Create summary
summary = {
    'n_samples': len(features),
    'n_ruptured': labels.sum(),
    'n_unruptured': (labels == 0).sum(),
    'feature_cols': feature_cols,
    'model_type': 'GradientBoostingClassifier',
    'model_params': {
        'n_estimators': 200,
        'max_depth': 5,
        'learning_rate': 0.05
    }
}

with open('results/phase5_model_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print("\n[OK] Model and scaler saved to results/")
print("[OK] Model summary saved to results/phase5_model_summary.json")
print("\n" + "=" * 80)
print("PHASE 5 ANALYSIS COMPLETE")
print("=" * 80)
