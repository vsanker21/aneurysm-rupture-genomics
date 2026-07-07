"""
Enhanced Tissue Model for Ruptured vs Unruptured Classification
Integrates GSE122897 and GSE66238 for improved prediction
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
from sklearn.model_selection import cross_val_score, StratifiedKFold, train_test_split
from sklearn.metrics import roc_auc_score, roc_curve, classification_report, confusion_matrix, accuracy_score
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("ENHANCED TISSUE MODEL: RUPTURED VS UNRUPTURED CLASSIFICATION")
print("Integrating GSE122897 and GSE66238")
print("=" * 80)

# ============================================================================
# LOAD TISSUE DATA - GSE122897 (Primary Dataset)
# ============================================================================

def load_gse122897_data():
    """Load GSE122897 tissue expression data and metadata"""
    print("\n[Step 1] Loading GSE122897 tissue transcriptomics data...")
    
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
    tissue_meta_aneurysm['Dataset'] = 'GSE122897'
    
    # Align samples
    common_samples = [s for s in tissue_meta_aneurysm['Sample'] if s in tissue_expr.columns]
    tissue_meta_aligned = tissue_meta_aneurysm[tissue_meta_aneurysm['Sample'].isin(common_samples)].copy()
    
    print(f"\n  Common samples: {len(common_samples)}")
    
    return tissue_expr, tissue_meta_aligned

# ============================================================================
# LOAD TISSUE DATA - GSE66238 (Secondary Dataset)
# ============================================================================

def load_gse66238_data():
    """Load GSE66238 tissue expression data if available"""
    print("\n[Step 2] Attempting to load GSE66238 tissue transcriptomics data...")
    
    tissue_expr_file = "data/sc_rnaseq/GSE66238/GSE66238_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    
    if not os.path.exists(tissue_expr_file):
        print(f"  [INFO] GSE66238 expression file not found: {tissue_expr_file}")
        print(f"  [INFO] Continuing with GSE122897 only")
        return None, None
    
    try:
        print(f"  Loading expression data...")
        tissue_expr = pd.read_csv(tissue_expr_file, sep='\t', compression='gzip', index_col=0)
        print(f"    Expression: {tissue_expr.shape[0]:,} genes, {tissue_expr.shape[1]} samples")
        
        # Try to extract metadata from sample names or check for metadata file
        metadata_file = "data/sc_rnaseq/GSE66238/sample_metadata.tsv"
        
        if os.path.exists(metadata_file):
            tissue_meta = pd.read_csv(metadata_file, sep='\t')
            print(f"    Metadata: {len(tissue_meta)} samples")
        else:
            print(f"  [INFO] Metadata file not found, will need manual annotation")
            # Create placeholder metadata
            tissue_meta = pd.DataFrame({
                'Sample': tissue_expr.columns,
                'Group': ['Unknown'] * len(tissue_expr.columns),
                'Label': [0] * len(tissue_expr.columns),
                'Dataset': 'GSE66238'
            })
        
        # Filter for ruptured vs unruptured if available
        if 'Group' in tissue_meta.columns:
            tissue_meta_aneurysm = tissue_meta[
                tissue_meta['Group'].isin(['Ruptured', 'Unruptured'])
            ].copy()
            
            if len(tissue_meta_aneurysm) > 0:
                tissue_meta_aneurysm['Label'] = (tissue_meta_aneurysm['Group'] == 'Ruptured').astype(int)
                print(f"    Ruptured vs Unruptured: {len(tissue_meta_aneurysm)} samples")
                return tissue_expr, tissue_meta_aneurysm
        
        print(f"  [INFO] No ruptured/unruptured labels available for GSE66238")
        return None, None
        
    except Exception as e:
        print(f"  [WARNING] Error loading GSE66238: {e}")
        print(f"  [INFO] Continuing with GSE122897 only")
        return None, None

# ============================================================================
# FEATURE SELECTION
# ============================================================================

def select_tissue_features(tissue_expr, tissue_meta):
    """Select top differentially expressed genes for features"""
    print("\n[Step 3] Selecting tissue-specific features...")
    
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
        top_genes = sig_genes.head(100)
        
        # Extract gene IDs - try multiple column names and index
        feature_genes = []
        if 'Gene' in top_genes.columns:
            feature_genes = top_genes['Gene'].astype(str).tolist()
        elif 'Gene_ID' in top_genes.columns:
            feature_genes = top_genes['Gene_ID'].astype(str).tolist()
        elif 'gene' in top_genes.columns:
            feature_genes = top_genes['gene'].astype(str).tolist()
        else:
            # Try using index
            feature_genes = top_genes.index.astype(str).tolist()
        
        print(f"  Selected top {len(feature_genes)} differentially expressed genes")
        if 'Log2FC' in top_genes.columns:
            print(f"    Upregulated in ruptured: {(top_genes['Log2FC'] > 0).sum()}")
            print(f"    Downregulated in ruptured: {(top_genes['Log2FC'] < 0).sum()}")
        
        # Check what gene IDs look like in expression data
        expr_gene_ids = list(tissue_expr.index[:10])
        de_gene_ids = feature_genes[:10] if len(feature_genes) > 10 else feature_genes
        print(f"  Sample expression gene IDs: {expr_gene_ids}")
        print(f"  Sample DE gene IDs: {de_gene_ids}")
        
    else:
        print(f"  [INFO] DE results not found, using top variable genes")
        # Use top variable genes
        common_samples = [s for s in tissue_meta['Sample'] if s in tissue_expr.columns]
        if len(common_samples) > 0:
            gene_vars = tissue_expr.loc[:, common_samples].var(axis=1)
            feature_genes = gene_vars.nlargest(100).index.astype(str).tolist()
            print(f"  Selected top {len(feature_genes)} variable genes")
        else:
            print(f"  [ERROR] No common samples found!")
            return []
    
    # Filter for genes present in expression data - try exact match first
    feature_genes_in_expr = [g for g in feature_genes if g in tissue_expr.index]
    
    # If no matches, try converting to int or removing prefix
    if len(feature_genes_in_expr) == 0 and len(feature_genes) > 0:
        print(f"  [INFO] No exact matches, trying alternative matching...")
        # Try converting to int
        try:
            expr_index_int = [str(int(x)) for x in tissue_expr.index if str(x).isdigit()]
            feature_genes_int = [str(int(x)) for x in feature_genes if str(x).isdigit()]
            feature_genes_in_expr = [g for g in feature_genes_int if g in expr_index_int]
            if len(feature_genes_in_expr) > 0:
                print(f"  [INFO] Found {len(feature_genes_in_expr)} matches after int conversion")
                # Map back to original expression index
                feature_genes_in_expr = [tissue_expr.index[i] for i, g in enumerate([str(int(x)) if str(x).isdigit() else x for x in tissue_expr.index]) if g in feature_genes_int]
        except:
            pass
    
    # If still no matches, use top variable genes instead
    if len(feature_genes_in_expr) == 0:
        print(f"  [WARNING] No feature genes matched, using top variable genes instead")
        common_samples = [s for s in tissue_meta['Sample'] if s in tissue_expr.columns]
        if len(common_samples) > 0:
            gene_vars = tissue_expr.loc[:, common_samples].var(axis=1)
            feature_genes_in_expr = gene_vars.nlargest(100).index.astype(str).tolist()
            print(f"  Selected top {len(feature_genes_in_expr)} variable genes")
        else:
            print(f"  [ERROR] No common samples found!")
            return []
    
    feature_genes = feature_genes_in_expr
    print(f"  Genes present in expression data: {len(feature_genes)}")
    
    return feature_genes

# ============================================================================
# BUILD MODEL
# ============================================================================

def build_tissue_model(tissue_expr, tissue_meta, feature_genes):
    """Build tissue model for ruptured vs unruptured classification"""
    print("\n[Step 4] Building enhanced tissue model...")
    
    # Extract features
    common_samples = [s for s in tissue_meta['Sample'] if s in tissue_expr.columns]
    
    if len(common_samples) == 0:
        print(f"  [ERROR] No common samples found between expression and metadata!")
        return None, None, None, None, None
    
    if len(feature_genes) == 0:
        print(f"  [ERROR] No features selected!")
        return None, None, None, None, None
    
    X = tissue_expr.loc[feature_genes, common_samples].T.values
    y = tissue_meta[tissue_meta['Sample'].isin(common_samples)]['Label'].values
    
    print(f"  Feature matrix: {X.shape}")
    print(f"  Labels: Ruptured={y.sum()}, Unruptured={(y == 0).sum()}")
    
    # Check for class imbalance
    if y.sum() < 3 or (y == 0).sum() < 3:
        print(f"  [WARNING] Very small sample size may affect model performance")
    
    # Handle missing values
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Cross-validation
    print("\n  Performing 5-fold cross-validation...")
    
    # Use fewer folds if small sample size
    n_splits = min(5, min(y.sum(), (y == 0).sum()))
    if n_splits < 3:
        print(f"  [WARNING] Sample size too small for cross-validation, using leave-one-out")
        from sklearn.model_selection import LeaveOneOut
        cv = LeaveOneOut()
    else:
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    
    # Models
    models = {
        'Random Forest': RandomForestClassifier(
            n_estimators=200, 
            max_depth=5, 
            min_samples_split=2,
            min_samples_leaf=1,
            random_state=42,
            class_weight='balanced'
        ),
        'Gradient Boosting': GradientBoostingClassifier(
            n_estimators=200, 
            max_depth=3, 
            learning_rate=0.05,
            random_state=42,
            subsample=0.8
        ),
        'SVM': SVC(
            kernel='rbf', 
            probability=True, 
            random_state=42, 
            C=1.0, 
            gamma='scale',
            class_weight='balanced'
        )
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
            import traceback
            traceback.print_exc()
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
            
            print(f"\n  Top 15 most important features:")
            for idx, row in feature_importance.head(15).iterrows():
                print(f"    {row['Gene_ID']:<20s} {row['Importance']:.4f}")
        
    return results, scaler, X_scaled, y, feature_genes

# ============================================================================
# SAVE RESULTS
# ============================================================================

def save_tissue_model_results(results, tissue_meta, scaler, X_scaled, y, feature_genes):
    """Save tissue model results"""
    print("\n[Step 5] Saving results...")
    
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
    summary_file = "results/enhanced_tissue_model_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"  [OK] Saved model summary to {summary_file}")
    
    # Best model predictions
    valid_results = {k: v for k, v in results.items() if v is not None}
    if valid_results:
        best_model_name = max(valid_results, key=lambda x: valid_results[x]['mean_auc'])
        best_model = valid_results[best_model_name]['model']
        
        predictions = best_model.predict_proba(X_scaled)[:, 1]
        predicted_labels = best_model.predict(X_scaled)
        
        predictions_df = tissue_meta.copy()
        predictions_df['Rupture_Probability'] = predictions
        predictions_df['Predicted'] = predicted_labels
        predictions_df['Correct'] = (predictions_df['Predicted'] == predictions_df['Label']).astype(int)
        
        predictions_file = "results/enhanced_tissue_model_predictions.tsv"
        predictions_df.to_csv(predictions_file, sep='\t', index=False)
        print(f"  [OK] Saved predictions to {predictions_file}")
        
        # Feature importance
        if hasattr(best_model, 'feature_importances_'):
            feature_importance = pd.DataFrame({
                'Gene_ID': feature_genes,
                'Importance': best_model.feature_importances_
            }).sort_values('Importance', ascending=False)
            
            importance_file = "results/enhanced_tissue_model_feature_importance.tsv"
            feature_importance.to_csv(importance_file, sep='\t', index=False)
            print(f"  [OK] Saved feature importance to {importance_file}")
            
            # Calculate accuracy
            accuracy = accuracy_score(y, predicted_labels)
            print(f"\n  Model Performance:")
            print(f"    Accuracy: {accuracy:.3f}")
            print(f"    AUC: {valid_results[best_model_name]['mean_auc']:.3f}")
            
            # Classification report
            print(f"\n  Classification Report:")
            print(classification_report(y, predicted_labels, 
                                      target_names=['Unruptured', 'Ruptured']))
    
    return summary_df

# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main pipeline"""
    
    # Load primary dataset
    tissue_expr, tissue_meta = load_gse122897_data()
    if tissue_expr is None or tissue_meta is None:
        print("\n[ERROR] Failed to load GSE122897 data")
        return
    
    if len(tissue_meta) < 10:
        print(f"\n[WARNING] Small sample size ({len(tissue_meta)}) may affect model performance")
    
    # Try to load secondary dataset
    tissue_expr2, tissue_meta2 = load_gse66238_data()
    
    # For now, use GSE122897 only (can extend later to integrate both)
    if tissue_expr2 is not None and tissue_meta2 is not None and len(tissue_meta2) > 0:
        print("\n[INFO] GSE66238 data available but not integrated yet")
        print("  [INFO] Using GSE122897 for now, can integrate later")
    
    # Select features
    feature_genes = select_tissue_features(tissue_expr, tissue_meta)
    if len(feature_genes) == 0:
        print("\n[ERROR] No features selected")
        return
    
    # Build model
    results, scaler, X_scaled, y, feature_genes = build_tissue_model(
        tissue_expr, tissue_meta, feature_genes
    )
    
    if results is None:
        print("\n[ERROR] Failed to build model")
        return
    
    # Save results
    summary_df = save_tissue_model_results(
        results, tissue_meta, scaler, X_scaled, y, feature_genes
    )
    
    # Summary
    print("\n" + "=" * 80)
    print("ENHANCED TISSUE MODEL COMPLETE")
    print("=" * 80)
    print(f"\nTask: Ruptured vs Unruptured classification")
    print(f"Primary Dataset: GSE122897 (tissue transcriptomics)")
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
    print(f"  - results/enhanced_tissue_model_summary.tsv")
    print(f"  - results/enhanced_tissue_model_feature_importance.tsv")
    print(f"  - results/enhanced_tissue_model_predictions.tsv")

if __name__ == "__main__":
    main()

