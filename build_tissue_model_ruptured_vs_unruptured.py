"""
Build Tissue Model for Ruptured vs Unruptured Classification
Using GSE122897 tissue transcriptomics data
"""
import sys
import os

# Force UTF-8 encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve, classification_report, confusion_matrix
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("TISSUE MODEL: RUPTURED VS UNRUPTURED CLASSIFICATION")
print("=" * 80)

# ============================================================================
# LOAD TISSUE DATA
# ============================================================================

def load_tissue_data():
    """Load tissue expression data and metadata"""
    print("\n[Step 1] Loading tissue transcriptomics data...")
    
    # Load expression data
    tissue_expr_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    tissue_meta_file = "data/sc_rnaseq/GSE122897/sample_metadata.tsv"
    
    if not os.path.exists(tissue_expr_file):
        print(f"  [ERROR] Expression file not found: {tissue_expr_file}")
        return None, None
    
    if not os.path.exists(tissue_meta_file):
        print(f"  [ERROR] Metadata file not found: {tissue_meta_file}")
        return None, None
    
    print(f"  Loading expression data...")
    tissue_expr = pd.read_csv(tissue_expr_file, sep='\t', compression='gzip', index_col=0)
    print(f"    Expression: {tissue_expr.shape[0]:,} genes, {tissue_expr.shape[1]} samples")
    
    print(f"  Loading metadata...")
    tissue_meta = pd.read_csv(tissue_meta_file, sep='\t')
    print(f"    Metadata: {len(tissue_meta)} samples")
    
    # Filter for ruptured vs unruptured only
    tissue_meta_aneurysm = tissue_meta[
        tissue_meta['Group'].isin(['Ruptured', 'Unruptured'])
    ].copy()
    
    print(f"\n  Filtered for ruptured vs unruptured:")
    print(f"    Ruptured: {(tissue_meta_aneurysm['Group'] == 'Ruptured').sum()}")
    print(f"    Unruptured: {(tissue_meta_aneurysm['Group'] == 'Unruptured').sum()}")
    print(f"    Total: {len(tissue_meta_aneurysm)} samples")
    
    # Create labels
    tissue_meta_aneurysm['Label'] = (tissue_meta_aneurysm['Group'] == 'Ruptured').astype(int)
    
    # Align samples
    common_samples = [s for s in tissue_meta_aneurysm['Sample'] if s in tissue_expr.columns]
    tissue_meta_aligned = tissue_meta_aneurysm[tissue_meta_aneurysm['Sample'].isin(common_samples)].copy()
    
    print(f"\n  Common samples: {len(common_samples)}")
    
    return tissue_expr, tissue_meta_aligned

# ============================================================================
# FEATURE SELECTION
# ============================================================================

def select_tissue_features(tissue_expr, tissue_meta):
    """Select top differentially expressed genes for features"""
    print("\n[Step 2] Selecting tissue-specific features...")
    
    # Load pre-computed differential expression if available
    de_file = "results/gse122897_differential_expression.tsv"
    
    if os.path.exists(de_file):
        print(f"  Loading pre-computed differential expression...")
        de_results = pd.read_csv(de_file, sep='\t')
        
        # Filter for significant genes
        if 'P_adj_BH' in de_results.columns:
            sig_genes = de_results[de_results['P_adj_BH'] < 0.05].copy()
        elif 'P_value' in de_results.columns:
            sig_genes = de_results[de_results['P_value'] < 0.01].copy()
        else:
            sig_genes = de_results.copy()
        
        # Sort by absolute log2FC
        if 'Log2FC' in sig_genes.columns:
            sig_genes['Abs_Log2FC'] = np.abs(sig_genes['Log2FC'])
            sig_genes = sig_genes.sort_values('Abs_Log2FC', ascending=False)
        
        # Get top genes
        top_genes = sig_genes.head(50)
        
        # Extract gene IDs
        if 'Gene' in top_genes.columns:
            feature_genes = top_genes['Gene'].astype(str).tolist()
        elif 'Gene_ID' in top_genes.columns:
            feature_genes = top_genes['Gene_ID'].astype(str).tolist()
        else:
            feature_genes = top_genes.index.astype(str).tolist()
        
        print(f"  Selected top {len(feature_genes)} differentially expressed genes")
        print(f"    Upregulated in ruptured: {(top_genes['Log2FC'] > 0).sum() if 'Log2FC' in top_genes.columns else 'N/A'}")
        print(f"    Downregulated in ruptured: {(top_genes['Log2FC'] < 0).sum() if 'Log2FC' in top_genes.columns else 'N/A'}")
        
    else:
        print(f"  [INFO] DE results not found, using top variable genes")
        # Use top variable genes
        common_samples = [s for s in tissue_meta['Sample'] if s in tissue_expr.columns]
        gene_vars = tissue_expr.loc[:, common_samples].var(axis=1)
        feature_genes = gene_vars.nlargest(50).index.astype(str).tolist()
        print(f"  Selected top {len(feature_genes)} variable genes")
    
    # Filter for genes present in expression data
    feature_genes = [g for g in feature_genes if g in tissue_expr.index]
    print(f"  Genes present in expression data: {len(feature_genes)}")
    
    return feature_genes

# ============================================================================
# BUILD MODEL
# ============================================================================

def build_tissue_model(tissue_expr, tissue_meta, feature_genes):
    """Build tissue model for ruptured vs unruptured classification"""
    print("\n[Step 3] Building tissue model...")
    
    # Extract features
    common_samples = [s for s in tissue_meta['Sample'] if s in tissue_expr.columns]
    X = tissue_expr.loc[feature_genes, common_samples].T.values
    y = tissue_meta[tissue_meta['Sample'].isin(common_samples)]['Label'].values
    
    print(f"  Feature matrix: {X.shape}")
    print(f"  Labels: Ruptured={y.sum()}, Unruptured={(y == 0).sum()}")
    
    # Handle missing values
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Cross-validation
    print("\n  Performing 5-fold cross-validation...")
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    # Models
    models = {
        'Random Forest': RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42),
        'Gradient Boosting': GradientBoostingClassifier(n_estimators=100, max_depth=3, learning_rate=0.1, random_state=42),
        'SVM': SVC(kernel='rbf', probability=True, random_state=42, C=1.0, gamma='scale')
    }
    
    results = {}
    
    for name, model in models.items():
        print(f"\n  {name}:")
        try:
            cv_scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')
            results[name] = {
                'mean_auc': cv_scores.mean(),
                'std_auc': cv_scores.std(),
                'scores': cv_scores,
                'min_auc': cv_scores.min(),
                'max_auc': cv_scores.max()
            }
            print(f"    CV AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
            print(f"    Range: {cv_scores.min():.3f} - {cv_scores.max():.3f}")
            
            # Train on full data for feature importance
            model.fit(X_scaled, y)
            results[name]['model'] = model
            
        except Exception as e:
            print(f"    [ERROR] {e}")
            results[name] = None
    
    # Find best model
    valid_results = {k: v for k, v in results.items() if v is not None}
    if valid_results:
        best_model_name = max(valid_results, key=lambda x: valid_results[x]['mean_auc'])
        best_auc = valid_results[best_model_name]['mean_auc']
        print(f"\n  Best model: {best_model_name} (AUC = {best_auc:.3f})")
        
        # Feature importance
        best_model = valid_results[best_model_name]['model']
        if hasattr(best_model, 'feature_importances_'):
            feature_importance = pd.DataFrame({
                'Gene_ID': feature_genes,
                'Importance': best_model.feature_importances_
            }).sort_values('Importance', ascending=False)
            
            print(f"\n  Top 10 most important features:")
            for idx, row in feature_importance.head(10).iterrows():
                print(f"    {row['Gene_ID']:<20s} {row['Importance']:.4f}")
            
            # Save feature importance
            importance_file = "results/tissue_model_feature_importance.tsv"
            feature_importance.to_csv(importance_file, sep='\t', index=False)
            print(f"\n  [OK] Saved feature importance to {importance_file}")
    
    return results, scaler, X_scaled, y, feature_genes

# ============================================================================
# SAVE RESULTS
# ============================================================================

def save_tissue_model_results(results, tissue_meta, scaler, X_scaled, y, feature_genes):
    """Save tissue model results"""
    print("\n[Step 4] Saving results...")
    
    os.makedirs("results", exist_ok=True)
    
    # Model summary
    summary_data = []
    for name, result in results.items():
        if result is not None:
            summary_data.append({
                'Model': name,
                'Mean_AUC': result['mean_auc'],
                'Std_AUC': result['std_auc'],
                'Min_AUC': result['min_auc'],
                'Max_AUC': result['max_auc'],
                'N_Features': len(feature_genes),
                'N_Samples': len(y),
                'N_Ruptured': y.sum(),
                'N_Unruptured': (y == 0).sum()
            })
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = "results/tissue_model_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"  [OK] Saved model summary to {summary_file}")
    
    # Best model predictions
    valid_results = {k: v for k, v in results.items() if v is not None}
    if valid_results:
        best_model_name = max(valid_results, key=lambda x: valid_results[x]['mean_auc'])
        best_model = valid_results[best_model_name]['model']
        
        predictions = best_model.predict_proba(X_scaled)[:, 1]
        
        predictions_df = tissue_meta.copy()
        predictions_df['Rupture_Probability'] = predictions
        predictions_df['Predicted'] = best_model.predict(X_scaled)
        
        predictions_file = "results/tissue_model_predictions.tsv"
        predictions_df.to_csv(predictions_file, sep='\t', index=False)
        print(f"  [OK] Saved predictions to {predictions_file}")
    
    return summary_df

# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main pipeline"""
    
    # Load data
    tissue_expr, tissue_meta = load_tissue_data()
    if tissue_expr is None or tissue_meta is None:
        print("\n[ERROR] Failed to load tissue data")
        return
    
    if len(tissue_meta) < 10:
        print(f"\n[WARNING] Small sample size ({len(tissue_meta)}) may affect model performance")
    
    # Select features
    feature_genes = select_tissue_features(tissue_expr, tissue_meta)
    if len(feature_genes) == 0:
        print("\n[ERROR] No features selected")
        return
    
    # Build model
    results, scaler, X_scaled, y, feature_genes = build_tissue_model(
        tissue_expr, tissue_meta, feature_genes
    )
    
    # Save results
    summary_df = save_tissue_model_results(
        results, tissue_meta, scaler, X_scaled, y, feature_genes
    )
    
    # Summary
    print("\n" + "=" * 80)
    print("TISSUE MODEL COMPLETE")
    print("=" * 80)
    print(f"\nTask: Ruptured vs Unruptured classification")
    print(f"Dataset: GSE122897 (tissue transcriptomics)")
    print(f"Samples: {len(y)} (Ruptured={y.sum()}, Unruptured={(y == 0).sum()})")
    print(f"Features: {len(feature_genes)} genes")
    
    if len(summary_df) > 0:
        print(f"\nModel Performance:")
        for idx, row in summary_df.iterrows():
            print(f"  {row['Model']:<20s} AUC: {row['Mean_AUC']:.3f} ± {row['Std_AUC']:.3f}")
        
        best_auc = summary_df['Mean_AUC'].max()
        print(f"\nBest AUC: {best_auc:.3f}")
        print(f"Expected performance: Better than random (0.50)")
        print(f"Target: AUC > 0.70 for clinical utility")
        
        if best_auc > 0.70:
            print(f"  [OK] Meets target for clinical utility")
        elif best_auc > 0.60:
            print(f"  [INFO] Moderate performance, may benefit from more samples")
        else:
            print(f"  [WARNING] Performance may be limited by sample size or features")
    
    print(f"\nOutput files:")
    print(f"  - results/tissue_model_summary.tsv")
    print(f"  - results/tissue_model_feature_importance.tsv")
    print(f"  - results/tissue_model_predictions.tsv")

if __name__ == "__main__":
    main()


