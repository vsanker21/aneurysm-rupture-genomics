"""
Phase 5: Model Optimization to Fix Overfitting and Improve Performance
Addressing the low CV AUC (0.42) and severe overfitting (train AUC=1.0)
"""
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, classification_report, confusion_matrix
import pickle
import json

print("=" * 80)
print("PHASE 5: MODEL OPTIMIZATION")
print("Fixing overfitting and improving CV performance")
print("=" * 80)

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

# Cross-validation setup
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

print(f"\nDataset: {len(features)} samples ({labels.sum()} ruptured, {(labels==0).sum()} unruptured)")
print(f"Features: {len(feature_cols)}")
print("\n" + "=" * 80)
print("MODEL COMPARISON")
print("=" * 80)

results = []

# 1. Regularized Gradient Boosting (reduced complexity)
print("\n1. Regularized Gradient Boosting (max_depth=2, subsample=0.8)...")
model1 = GradientBoostingClassifier(
    n_estimators=100,  # Reduced from 200
    max_depth=2,       # Reduced from 5
    learning_rate=0.1, # Increased from 0.05
    subsample=0.8,     # Add subsampling for regularization
    max_features='sqrt',  # Feature subsampling
    random_state=42
)
scores1 = cross_val_score(model1, X_scaled, labels, cv=cv, scoring='roc_auc')
train_model1 = model1.fit(X_scaled, labels)
train_pred1 = train_model1.predict_proba(X_scaled)[:, 1]
train_auc1 = roc_auc_score(labels, train_pred1)

results.append({
    'Model': 'GB_Regularized',
    'Mean_CV_AUC': scores1.mean(),
    'Std_CV_AUC': scores1.std(),
    'Train_AUC': train_auc1,
    'Overfit_Gap': train_auc1 - scores1.mean(),
    'Model_Obj': train_model1
})

print(f"  CV AUC: {scores1.mean():.3f} ± {scores1.std():.3f}")
print(f"  Train AUC: {train_auc1:.3f}")
print(f"  Overfitting gap: {train_auc1 - scores1.mean():.3f}")

# 2. Regularized Gradient Boosting (even more conservative)
print("\n2. Conservative Gradient Boosting (max_depth=2, subsample=0.7)...")
model2 = GradientBoostingClassifier(
    n_estimators=50,
    max_depth=2,
    learning_rate=0.1,
    subsample=0.7,
    max_features='sqrt',
    random_state=42
)
scores2 = cross_val_score(model2, X_scaled, labels, cv=cv, scoring='roc_auc')
train_model2 = model2.fit(X_scaled, labels)
train_pred2 = train_model2.predict_proba(X_scaled)[:, 1]
train_auc2 = roc_auc_score(labels, train_pred2)

results.append({
    'Model': 'GB_Conservative',
    'Mean_CV_AUC': scores2.mean(),
    'Std_CV_AUC': scores2.std(),
    'Train_AUC': train_auc2,
    'Overfit_Gap': train_auc2 - scores2.mean(),
    'Model_Obj': train_model2
})

print(f"  CV AUC: {scores2.mean():.3f} ± {scores2.std():.3f}")
print(f"  Train AUC: {train_auc2:.3f}")
print(f"  Overfitting gap: {train_auc2 - scores2.mean():.3f}")

# 3. Feature selection: Remove zero-importance features
print("\n3. Feature Selection: Using only top features (DE, EC, ECM)...")
top_features = ['DE_Signature_Score', 'EC_Dysfunction_Score', 'ECM_Remodeling_Score']
feature_indices = [feature_cols.index(f) for f in top_features]
X_selected = X_scaled[:, feature_indices]

model3 = GradientBoostingClassifier(
    n_estimators=100,
    max_depth=2,
    learning_rate=0.1,
    subsample=0.8,
    max_features='sqrt',
    random_state=42
)
scores3 = cross_val_score(model3, X_selected, labels, cv=cv, scoring='roc_auc')
train_model3 = model3.fit(X_selected, labels)
train_pred3 = train_model3.predict_proba(X_selected)[:, 1]
train_auc3 = roc_auc_score(labels, train_pred3)

results.append({
    'Model': 'GB_FeatureSelected',
    'Mean_CV_AUC': scores3.mean(),
    'Std_CV_AUC': scores3.std(),
    'Train_AUC': train_auc3,
    'Overfit_Gap': train_auc3 - scores3.mean(),
    'Model_Obj': train_model3,
    'Features_Used': top_features
})

print(f"  CV AUC: {scores3.mean():.3f} ± {scores3.std():.3f}")
print(f"  Train AUC: {train_auc3:.3f}")
print(f"  Overfitting gap: {train_auc3 - scores3.mean():.3f}")
print(f"  Features: {top_features}")

# 4. Regularized Random Forest
print("\n4. Regularized Random Forest...")
model4 = RandomForestClassifier(
    n_estimators=100,
    max_depth=3,  # Limited depth
    min_samples_split=5,
    min_samples_leaf=2,
    max_features='sqrt',
    random_state=42,
    class_weight='balanced'
)
scores4 = cross_val_score(model4, X_scaled, labels, cv=cv, scoring='roc_auc')
train_model4 = model4.fit(X_scaled, labels)
train_pred4 = train_model4.predict_proba(X_scaled)[:, 1]
train_auc4 = roc_auc_score(labels, train_pred4)

results.append({
    'Model': 'RF_Regularized',
    'Mean_CV_AUC': scores4.mean(),
    'Std_CV_AUC': scores4.std(),
    'Train_AUC': train_auc4,
    'Overfit_Gap': train_auc4 - scores4.mean(),
    'Model_Obj': train_model4
})

print(f"  CV AUC: {scores4.mean():.3f} ± {scores4.std():.3f}")
print(f"  Train AUC: {train_auc4:.3f}")
print(f"  Overfitting gap: {train_auc4 - scores4.mean():.3f}")

# 5. Logistic Regression with L2 regularization
print("\n5. Logistic Regression (L2 regularization)...")
model5 = LogisticRegression(
    C=1.0,  # Inverse of regularization strength
    penalty='l2',
    solver='liblinear',
    max_iter=1000,
    random_state=42,
    class_weight='balanced'
)
scores5 = cross_val_score(model5, X_scaled, labels, cv=cv, scoring='roc_auc')
train_model5 = model5.fit(X_scaled, labels)
train_pred5 = train_model5.predict_proba(X_scaled)[:, 1]
train_auc5 = roc_auc_score(labels, train_pred5)

results.append({
    'Model': 'Logistic_L2',
    'Mean_CV_AUC': scores5.mean(),
    'Std_CV_AUC': scores5.std(),
    'Train_AUC': train_auc5,
    'Overfit_Gap': train_auc5 - scores5.mean(),
    'Model_Obj': train_model5
})

print(f"  CV AUC: {scores5.mean():.3f} ± {scores5.std():.3f}")
print(f"  Train AUC: {train_auc5:.3f}")
print(f"  Overfitting gap: {train_auc5 - scores5.mean():.3f}")

# 6. Logistic Regression with feature selection
print("\n6. Logistic Regression with Feature Selection...")
model6 = LogisticRegression(
    C=1.0,
    penalty='l2',
    solver='liblinear',
    max_iter=1000,
    random_state=42,
    class_weight='balanced'
)
scores6 = cross_val_score(model6, X_selected, labels, cv=cv, scoring='roc_auc')
train_model6 = model6.fit(X_selected, labels)
train_pred6 = train_model6.predict_proba(X_selected)[:, 1]
train_auc6 = roc_auc_score(labels, train_pred6)

results.append({
    'Model': 'Logistic_FeatureSelected',
    'Mean_CV_AUC': scores6.mean(),
    'Std_CV_AUC': scores6.std(),
    'Train_AUC': train_auc6,
    'Overfit_Gap': train_auc6 - scores6.mean(),
    'Model_Obj': train_model6,
    'Features_Used': top_features
})

print(f"  CV AUC: {scores6.mean():.3f} ± {scores6.std():.3f}")
print(f"  Train AUC: {train_auc6:.3f}")
print(f"  Overfitting gap: {train_auc6 - scores6.mean():.3f}")

# 7. Ensemble: Average predictions from top models
print("\n7. Ensemble (Average of top 3 models)...")
# Select top 3 models by CV AUC
sorted_results = sorted(results, key=lambda x: x['Mean_CV_AUC'], reverse=True)
top3_models = sorted_results[:3]

ensemble_preds_train = np.zeros(len(labels))
for res in top3_models:
    if 'Features_Used' in res:
        X_used = X_selected
    else:
        X_used = X_scaled
    ensemble_preds_train += res['Model_Obj'].predict_proba(X_used)[:, 1]
ensemble_preds_train /= len(top3_models)
ensemble_train_auc = roc_auc_score(labels, ensemble_preds_train)

# Cross-validation for ensemble
ensemble_cv_scores = []
for train_idx, val_idx in cv.split(X_scaled, labels):
    X_train_cv, X_val_cv = X_scaled[train_idx], X_scaled[val_idx]
    X_train_cv_sel, X_val_cv_sel = X_selected[train_idx], X_selected[val_idx]
    y_train_cv, y_val_cv = labels.values[train_idx], labels.values[val_idx]
    
    ensemble_preds_val = np.zeros(len(y_val_cv))
    for res in top3_models:
        if res['Model'] == 'GB_FeatureSelected' or res['Model'] == 'Logistic_FeatureSelected':
            X_train_use = X_train_cv_sel
            X_val_use = X_val_cv_sel
        else:
            X_train_use = X_train_cv
            X_val_use = X_val_cv
        
        # Retrain model for this fold
        if 'GB' in res['Model']:
            if res['Model'] == 'GB_FeatureSelected':
                m = GradientBoostingClassifier(n_estimators=100, max_depth=2, learning_rate=0.1, 
                                             subsample=0.8, max_features='sqrt', random_state=42)
                m.fit(X_train_use, y_train_cv)
                ensemble_preds_val += m.predict_proba(X_val_use)[:, 1]
            elif res['Model'] == 'GB_Regularized':
                m = GradientBoostingClassifier(n_estimators=100, max_depth=2, learning_rate=0.1, 
                                             subsample=0.8, max_features='sqrt', random_state=42)
                m.fit(X_train_use, y_train_cv)
                ensemble_preds_val += m.predict_proba(X_val_use)[:, 1]
            else:
                m = GradientBoostingClassifier(n_estimators=50, max_depth=2, learning_rate=0.1, 
                                             subsample=0.7, max_features='sqrt', random_state=42)
                m.fit(X_train_use, y_train_cv)
                ensemble_preds_val += m.predict_proba(X_val_use)[:, 1]
        elif 'Logistic' in res['Model']:
            m = LogisticRegression(C=1.0, penalty='l2', solver='liblinear', max_iter=1000, 
                                  random_state=42, class_weight='balanced')
            m.fit(X_train_use, y_train_cv)
            ensemble_preds_val += m.predict_proba(X_val_use)[:, 1]
        else:  # Random Forest
            m = RandomForestClassifier(n_estimators=100, max_depth=3, min_samples_split=5, 
                                     min_samples_leaf=2, max_features='sqrt', random_state=42, 
                                     class_weight='balanced')
            m.fit(X_train_use, y_train_cv)
            ensemble_preds_val += m.predict_proba(X_val_use)[:, 1]
    
    ensemble_preds_val /= len(top3_models)
    ensemble_cv_scores.append(roc_auc_score(y_val_cv, ensemble_preds_val))

ensemble_cv_auc = np.mean(ensemble_cv_scores)
ensemble_cv_std = np.std(ensemble_cv_scores)

results.append({
    'Model': 'Ensemble_Top3',
    'Mean_CV_AUC': ensemble_cv_auc,
    'Std_CV_AUC': ensemble_cv_std,
    'Train_AUC': ensemble_train_auc,
    'Overfit_Gap': ensemble_train_auc - ensemble_cv_auc,
    'Model_Obj': top3_models  # List of models
})

print(f"  CV AUC: {ensemble_cv_auc:.3f} ± {ensemble_cv_std:.3f}")
print(f"  Train AUC: {ensemble_train_auc:.3f}")
print(f"  Overfitting gap: {ensemble_train_auc - ensemble_cv_auc:.3f}")
print(f"  Ensemble members: {[r['Model'] for r in top3_models]}")

# Summary
print("\n" + "=" * 80)
print("MODEL COMPARISON SUMMARY")
print("=" * 80)

summary_df = pd.DataFrame([{
    'Model': r['Model'],
    'CV_AUC_Mean': r['Mean_CV_AUC'],
    'CV_AUC_Std': r['Std_CV_AUC'],
    'Train_AUC': r['Train_AUC'],
    'Overfitting_Gap': r['Overfit_Gap']
} for r in results])

summary_df = summary_df.sort_values('CV_AUC_Mean', ascending=False)
print(summary_df.to_string(index=False))

# Select best model (highest CV AUC with reasonable overfitting gap)
best_idx = summary_df.index[0]
best_model_name = summary_df.iloc[0]['Model']
best_result = next(r for r in results if r['Model'] == best_model_name)

print(f"\n{'='*80}")
print(f"BEST MODEL: {best_model_name}")
print(f"{'='*80}")
print(f"CV AUC: {best_result['Mean_CV_AUC']:.3f} ± {best_result['Std_CV_AUC']:.3f}")
print(f"Train AUC: {best_result['Train_AUC']:.3f}")
print(f"Overfitting gap: {best_result['Overfit_Gap']:.3f}")

# Save best model and results
if 'Features_Used' in best_result:
    best_features = best_result['Features_Used']
    best_X = X_selected
    best_scaler = StandardScaler()
    best_scaler.fit(X[:, [feature_cols.index(f) for f in best_features]])
else:
    best_features = feature_cols
    best_X = X_scaled
    best_scaler = scaler

# Save best model
with open('results/phase5_optimized_model.pkl', 'wb') as f:
    pickle.dump(best_result['Model_Obj'], f)

with open('results/phase5_optimized_scaler.pkl', 'wb') as f:
    pickle.dump(best_scaler, f)

# Save comparison results
summary_df.to_csv('results/phase5_model_comparison.tsv', sep='\t', index=False)

# Save detailed results
detailed_results = []
for r in results:
    if 'Features_Used' in r:
        detailed_results.append({
            'Model': r['Model'],
            'CV_AUC_Mean': r['Mean_CV_AUC'],
            'CV_AUC_Std': r['Std_CV_AUC'],
            'Train_AUC': r['Train_AUC'],
            'Overfitting_Gap': r['Overfit_Gap'],
            'Features_Used': ', '.join(r['Features_Used'])
        })
    else:
        detailed_results.append({
            'Model': r['Model'],
            'CV_AUC_Mean': r['Mean_CV_AUC'],
            'CV_AUC_Std': r['Std_CV_AUC'],
            'Train_AUC': r['Train_AUC'],
            'Overfitting_Gap': r['Overfit_Gap'],
            'Features_Used': 'All'
        })

pd.DataFrame(detailed_results).to_csv('results/phase5_optimization_results.tsv', sep='\t', index=False)

# Save best model summary
best_summary = {
    'best_model': best_model_name,
    'cv_auc_mean': float(best_result['Mean_CV_AUC']),
    'cv_auc_std': float(best_result['Std_CV_AUC']),
    'train_auc': float(best_result['Train_AUC']),
    'overfitting_gap': float(best_result['Overfit_Gap']),
    'features_used': best_features if 'Features_Used' in best_result else feature_cols,
    'improvement': float(best_result['Mean_CV_AUC'] - 0.425)  # Improvement over original 0.425
}

with open('results/phase5_best_model_summary.json', 'w') as f:
    json.dump(best_summary, f, indent=2)

print(f"\n[OK] Best model saved: results/phase5_optimized_model.pkl")
print(f"[OK] Comparison results saved: results/phase5_model_comparison.tsv")
print(f"[OK] Detailed results saved: results/phase5_optimization_results.tsv")
print(f"[OK] Summary saved: results/phase5_best_model_summary.json")

improvement = best_result['Mean_CV_AUC'] - 0.425
print(f"\nIMPROVEMENT: {improvement:+.3f} AUC points over original model (0.425)")
print(f"Relative improvement: {(improvement/0.425)*100:+.1f}%")

print("\n" + "=" * 80)
print("OPTIMIZATION COMPLETE")
print("=" * 80)
