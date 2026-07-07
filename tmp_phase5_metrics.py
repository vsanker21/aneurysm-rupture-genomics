
import pandas as pd

import multi_omics_integration as moi

metadata_dir = 'data/sc_rnaseq/GSE122897'

metadata_file = metadata_dir + '/sample_metadata_with_clinical.tsv'

try:

    metadata = pd.read_csv(metadata_file, sep='\t')

except FileNotFoundError:

    metadata = pd.read_csv(metadata_dir + '/sample_metadata.tsv', sep='\t')

expr_df = moi.load_tissue_expression()

blood_sig = moi.load_blood_signature()

ec_scores = moi.load_ec_dysfunction_scores()

features, feature_cols = moi.prepare_integrated_features(metadata, expr_df, blood_sig, ec_scores)

labels = (features['Group'] == 'Ruptured').astype(int)

model, scaler, cv_scores, feature_importance = moi.build_integrated_model(features, feature_cols, labels)

import json
import os

cv_mean = float(cv_scores.mean())
cv_std = float(cv_scores.std())

print('CV_mean', cv_mean)
print('CV_std', cv_std)
print('feature_cols', feature_cols)
print('feature_importance', feature_importance.to_dict('records'))

# Save metrics to JSON file
os.makedirs('results', exist_ok=True)
metrics = {
    'CV_mean': cv_mean,
    'CV_std': cv_std,
    'CV_scores': cv_scores.tolist(),
    'feature_cols': feature_cols,
    'feature_importance': feature_importance.to_dict('records'),
    'n_samples': len(features),
    'n_ruptured': int(labels.sum()),
    'n_unruptured': int((labels == 0).sum())
}

with open('results/phase5_metrics.json', 'w') as f:
    json.dump(metrics, f, indent=2)

print('\n[OK] Metrics saved to results/phase5_metrics.json')

