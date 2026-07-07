"""
Build final multi-omics prediction model using real data
Integrates blood transcriptomics, tissue transcriptomics, and creates prediction model
"""
import sys
import os
import pandas as pd
import numpy as np

# Force UTF-8 encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve, classification_report

print("=" * 80)
print("BUILDING FINAL MULTI-OMICS PREDICTION MODEL")
print("=" * 80)

# Load existing integration results
print("\n[Step 1] Loading existing integration data...")

blood_tissue_file = "results/blood_tissue_integration_fdr_0.05.tsv"
if os.path.exists(blood_tissue_file):
    print(f"  Loading: {blood_tissue_file}")
    integration_df = pd.read_csv(blood_tissue_file, sep='\t')
    print(f"  [OK] Loaded {len(integration_df)} genes")
else:
    print(f"  [ERROR] File not found: {blood_tissue_file}")
    sys.exit(1)

# Load expression data
print("\n[Step 2] Loading expression data...")

# Blood expression
blood_expr_file = "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
blood_meta_file = "data/blood_rnaseq/GSE159610/sample_metadata.tsv"

if os.path.exists(blood_expr_file) and os.path.exists(blood_meta_file):
    print(f"  Loading blood expression: {blood_expr_file}")
    blood_expr = pd.read_csv(blood_expr_file, sep='\t', compression='gzip', index_col=0)
    blood_meta = pd.read_csv(blood_meta_file, sep='\t')
    
    # Create status labels
    blood_meta['Status'] = blood_meta['Sample'].apply(
        lambda x: 1 if int(x.replace('GSM4835', '')) <= 5782 else 0
    )
    print(f"  [OK] Blood: {blood_expr.shape[0]:,} genes, {blood_expr.shape[1]} samples")
else:
    print("  [WARNING] Blood data not found")
    blood_expr = None
    blood_meta = None

# Tissue expression  
tissue_expr_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
tissue_meta_file = "data/sc_rnaseq/GSE122897/sample_metadata.tsv"

if os.path.exists(tissue_expr_file) and os.path.exists(tissue_meta_file):
    print(f"  Loading tissue expression: {tissue_expr_file}")
    tissue_expr = pd.read_csv(tissue_expr_file, sep='\t', compression='gzip', index_col=0)
    tissue_meta = pd.read_csv(tissue_meta_file, sep='\t')
    
    # Filter for ruptured vs unruptured
    if 'Group' in tissue_meta.columns:
        tissue_meta_aneurysm = tissue_meta[tissue_meta['Group'].isin(['Ruptured', 'Unruptured'])].copy()
        tissue_meta_aneurysm['Status'] = (tissue_meta_aneurysm['Group'] == 'Ruptured').astype(int)
        print(f"  [OK] Tissue: {tissue_expr.shape[0]:,} genes, {len(tissue_meta_aneurysm)} aneurysm samples")
    else:
        tissue_meta_aneurysm = None
else:
    print("  [WARNING] Tissue data not found")
    tissue_expr = None
    tissue_meta_aneurysm = None

# Get top signature genes
print("\n[Step 3] Selecting signature genes...")

# Prioritize genes found in both blood and tissue, or highly significant
if 'In_Blood' in integration_df.columns and 'In_Tissue' in integration_df.columns:
    # Genes in both datasets
    both_genes = integration_df[
        (integration_df['In_Blood'] == True) & (integration_df['In_Tissue'] == True)
    ].copy()
    
    # Highly significant genes
    sig_genes = integration_df[
        (integration_df.get('Blood_P_adj', 1) < 0.05) | 
        (integration_df.get('Tissue_P_adj', 1) < 0.05)
    ].copy()
    
    # Combine and sort by priority score
    top_genes = pd.concat([both_genes, sig_genes]).drop_duplicates()
    
    if 'Multi_Omics_Priority_Score' in top_genes.columns:
        top_genes = top_genes.sort_values('Multi_Omics_Priority_Score', ascending=False)
    
    selected_genes = top_genes.head(50)  # Top 50 genes
    
    print(f"  Selected {len(selected_genes)} signature genes")
    print(f"    In both datasets: {(selected_genes['In_Blood'] & selected_genes['In_Tissue']).sum()}")
    print(f"    Blood only: {(selected_genes['In_Blood'] & ~selected_genes['In_Tissue']).sum()}")
    print(f"    Tissue only: {(selected_genes['In_Tissue'] & ~selected_genes['In_Blood']).sum()}")
    
    # Get gene identifiers
    gene_ids = []
    for idx, row in selected_genes.iterrows():
        if pd.notna(row.get('Gene_Symbol')):
            gene_ids.append(str(row['Gene_Symbol']))
        elif pd.notna(row.get('Blood_Gene_ID')):
            gene_ids.append(str(row['Blood_Gene_ID']))
        elif pd.notna(row.get('Tissue_Gene_ID')):
            gene_ids.append(str(row['Tissue_Gene_ID']))
        elif pd.notna(row.get('Entrez_ID')):
            gene_ids.append(str(int(row['Entrez_ID'])))
else:
    # Fallback: use top genes by priority
    selected_genes = integration_df.head(50)
    print(f"  Selected top {len(selected_genes)} genes by priority")
    gene_ids = selected_genes[selected_genes.columns[0]].astype(str).tolist()

# Create feature matrix
print("\n[Step 4] Creating feature matrix...")

features_list = []

# Add blood samples
if blood_expr is not None and blood_meta is not None:
    for sample in blood_expr.columns:
        if sample in blood_meta['Sample'].values:
            sample_meta = blood_meta[blood_meta['Sample'] == sample].iloc[0]
            features = {
                'Sample_ID': sample,
                'Dataset': 'Blood_GSE159610',
                'Label': int(sample_meta['Status'])
            }
            
            # Extract expression for signature genes
            for gene_id in gene_ids:
                # Try different ID formats
                if gene_id in blood_expr.index:
                    features[f'Gene_{gene_id}'] = blood_expr.loc[gene_id, sample]
                else:
                    # Try matching by gene symbol or other identifiers
                    features[f'Gene_{gene_id}'] = 0
            
            features_list.append(features)
    
    print(f"  [OK] Added {len([f for f in features_list if f['Dataset'] == 'Blood_GSE159610'])} blood samples")

# Add tissue samples
if tissue_expr is not None and tissue_meta_aneurysm is not None:
    for sample in tissue_expr.columns:
        if sample in tissue_meta_aneurysm['Sample'].values:
            sample_meta = tissue_meta_aneurysm[tissue_meta_aneurysm['Sample'] == sample].iloc[0]
            features = {
                'Sample_ID': sample,
                'Dataset': 'Tissue_GSE122897',
                'Label': int(sample_meta['Status'])
            }
            
            # Extract expression for signature genes
            for gene_id in gene_ids:
                # Try matching by numeric ID or symbol
                if str(gene_id) in tissue_expr.index or str(int(float(gene_id))) in tissue_expr.index:
                    key = str(gene_id) if str(gene_id) in tissue_expr.index else str(int(float(gene_id)))
                    features[f'Gene_{gene_id}'] = tissue_expr.loc[key, sample]
                else:
                    features[f'Gene_{gene_id}'] = 0
            
            features_list.append(features)
    
    print(f"  [OK] Added {len([f for f in features_list if f['Dataset'] == 'Tissue_GSE122897'])} tissue samples")

if len(features_list) == 0:
    print("  [ERROR] No features created")
    sys.exit(1)

# Create DataFrame
features_df = pd.DataFrame(features_list)

# Ensure all gene columns exist
for gene_id in gene_ids:
    col = f'Gene_{gene_id}'
    if col not in features_df.columns:
        features_df[col] = 0

# Fill missing values
features_df = features_df.fillna(0)

print(f"\n[OK] Feature matrix created: {len(features_df)} samples, {len(gene_ids)} features")

# Build model
print("\n[Step 5] Building prediction models...")

feature_cols = [col for col in features_df.columns if col.startswith('Gene_')]
X = features_df[feature_cols].values
y = features_df['Label'].values

# Handle missing/infinite values
X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print(f"  Features: {len(feature_cols)}")
print(f"  Samples: {len(X_scaled)}")
print(f"  Positive: {y.sum()}")
print(f"  Negative: {(y == 0).sum()}")

# Train models
print("\n  Training Gradient Boosting...")
gb_model = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=42)
gb_model.fit(X_scaled, y)

print("  Training Random Forest...")
rf_model = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)
rf_model.fit(X_scaled, y)

# Cross-validation
print("\n  Performing 5-fold cross-validation...")
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

if len(np.unique(y)) > 1:
    gb_cv_scores = cross_val_score(gb_model, X_scaled, y, cv=cv, scoring='roc_auc')
    rf_cv_scores = cross_val_score(rf_model, X_scaled, y, cv=cv, scoring='roc_auc')
    
    print(f"\n  Cross-validation AUC:")
    print(f"    Gradient Boosting: {gb_cv_scores.mean():.3f} ± {gb_cv_scores.std():.3f}")
    print(f"    Random Forest: {rf_cv_scores.mean():.3f} ± {rf_cv_scores.std():.3f}")
    
    # Predictions
    gb_pred_proba = gb_model.predict_proba(X_scaled)[:, 1]
    rf_pred_proba = rf_model.predict_proba(X_scaled)[:, 1]
    
    gb_auc = roc_auc_score(y, gb_pred_proba)
    rf_auc = roc_auc_score(y, rf_pred_proba)
    
    print(f"\n  Full dataset AUC:")
    print(f"    Gradient Boosting: {gb_auc:.3f}")
    print(f"    Random Forest: {rf_auc:.3f}")
else:
    print("  [WARNING] Only one class present, skipping AUC calculation")
    gb_cv_scores = [0]
    rf_cv_scores = [0]

# Feature importance
feature_importance = pd.DataFrame({
    'Feature': feature_cols,
    'Gene_ID': [col.replace('Gene_', '') for col in feature_cols],
    'GB_Importance': gb_model.feature_importances_,
    'RF_Importance': rf_model.feature_importances_
}).sort_values('GB_Importance', ascending=False)

print(f"\n  Top 10 most important features:")
for idx, row in feature_importance.head(10).iterrows():
    print(f"    {row['Gene_ID']:<20s} GB: {row['GB_Importance']:.4f}  RF: {row['RF_Importance']:.4f}")

# Save results
print("\n[Step 6] Saving results...")

results_df = features_df.copy()
results_df['GB_Predicted_Proba'] = gb_model.predict_proba(X_scaled)[:, 1]
results_df['RF_Predicted_Proba'] = rf_model.predict_proba(X_scaled)[:, 1]
results_df['GB_Predicted'] = gb_model.predict(X_scaled)
results_df['RF_Predicted'] = rf_model.predict(X_scaled)

os.makedirs("results", exist_ok=True)

results_file = "results/multi_omics_predictions_final.tsv"
results_df.to_csv(results_file, sep='\t', index=False)
print(f"  [OK] Saved predictions: {results_file}")

importance_file = "results/multi_omics_feature_importance_final.tsv"
feature_importance.to_csv(importance_file, sep='\t', index=False)
print(f"  [OK] Saved feature importance: {importance_file}")

# Summary
print("\n" + "=" * 80)
print("FINAL MULTI-OMICS MODEL COMPLETE")
print("=" * 80)
print(f"\nSummary:")
print(f"  Signature genes: {len(selected_genes)}")
print(f"  Samples: {len(features_df)}")
print(f"    Blood: {(features_df['Dataset'] == 'Blood_GSE159610').sum()}")
print(f"    Tissue: {(features_df['Dataset'] == 'Tissue_GSE122897').sum()}")
print(f"  Features: {len(feature_cols)}")
if len(np.unique(y)) > 1:
    print(f"  GB Cross-validation AUC: {gb_cv_scores.mean():.3f} ± {gb_cv_scores.std():.3f}")
    print(f"  RF Cross-validation AUC: {rf_cv_scores.mean():.3f} ± {rf_cv_scores.std():.3f}")
print(f"\nTarget AUC: 0.85-0.92 (integrated model)")
print(f"\nOutput files:")
print(f"  - {results_file}")
print(f"  - {importance_file}")


