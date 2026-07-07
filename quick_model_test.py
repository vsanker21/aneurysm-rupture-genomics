"""Quick model test to generate results"""
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
import sys

print("Loading data...")
data = pd.read_csv('results/integrated_features_full.tsv', sep='\t', index_col=0)
features = data.drop('Ruptured', axis=1)
labels = data['Ruptured']

print(f"Samples: {len(features)}, Features: {list(features.columns)}")
print(f"Ruptured: {labels.sum()}, Unruptured: {(labels==0).sum()}")

X = features.values
y = labels.values
X = np.nan_to_num(X, nan=0.0)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

models = {
    'Gradient Boosting': GradientBoostingClassifier(n_estimators=200, max_depth=5, learning_rate=0.05, random_state=42),
    'Random Forest': RandomForestClassifier(n_estimators=200, max_depth=10, random_state=42),
    'XGBoost': XGBClassifier(n_estimators=200, max_depth=5, learning_rate=0.05, random_state=42, eval_metric='logloss', use_label_encoder=False)
}

print("\nTraining models...")
results = {}
for name, model in models.items():
    print(f"  {name}...")
    scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc', n_jobs=-1)
    results[name] = {
        'mean_auc': scores.mean(),
        'std_auc': scores.std(),
        'min_auc': scores.min(),
        'max_auc': scores.max()
    }

print("\nResults:")
for name, r in results.items():
    print(f"  {name}: AUC = {r['mean_auc']:.3f} ± {r['std_auc']:.3f} (range: {r['min_auc']:.3f} - {r['max_auc']:.3f})")

best_name = max(results.items(), key=lambda x: x[1]['mean_auc'])[0]
best_auc = results[best_name]['mean_auc']
print(f"\nBest model: {best_name} (AUC = {best_auc:.3f})")

# Save results
results_df = pd.DataFrame({
    'Model': list(results.keys()),
    'Mean_AUC': [r['mean_auc'] for r in results.values()],
    'Std_AUC': [r['std_auc'] for r in results.values()],
    'Min_AUC': [r['min_auc'] for r in results.values()],
    'Max_AUC': [r['max_auc'] for r in results.values()]
})

results_df.to_csv('results/full_multi_omics_results.tsv', sep='\t', index=False)
print("\nSaved results to results/full_multi_omics_results.tsv")

sys.stdout.flush()


