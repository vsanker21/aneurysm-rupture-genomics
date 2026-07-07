"""
Multi-Omics Integration with Real Data
Integrates PRS, Blood Transcriptomics (GSE159610), and Tissue Transcriptomics (GSE122897)
for aneurysm rupture prediction using real patient data
"""
import sys
import os

# Force UTF-8 encoding for Windows
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve, classification_report, confusion_matrix
import warnings
warnings.filterwarnings('ignore')

def load_blood_expression_data():
    """Load GSE159610 blood RNA-seq expression data"""
    print("\n" + "=" * 80)
    print("LOADING BLOOD TRANSCRIPTOMICS DATA (GSE159610)")
    print("=" * 80)
    
    # Load raw counts
    raw_counts_file = "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
    
    if not os.path.exists(raw_counts_file):
        print(f"[ERROR] Blood expression file not found: {raw_counts_file}")
        return None, None
    
    print(f"Loading: {raw_counts_file}")
    expr_df = pd.read_csv(raw_counts_file, sep='\t', compression='gzip', index_col=0)
    
    print(f"[OK] Loaded blood expression data:")
    print(f"  Genes: {expr_df.shape[0]:,}")
    print(f"  Samples: {expr_df.shape[1]}")
    
    # Load metadata
    metadata_file = "data/blood_rnaseq/GSE159610/sample_metadata.tsv"
    if os.path.exists(metadata_file):
        metadata = pd.read_csv(metadata_file, sep='\t')
        # Extract status from title (aneurysm vs control based on angiography)
        # Samples DSA1-DSA23 are aneurysm positive, DSA24-DSA47 are control
        metadata['Status'] = metadata['Sample'].apply(
            lambda x: 'Aneurysm' if int(x.replace('GSM4835', '')) <= 5782 else 'Control'
        )
        print(f"\n[OK] Loaded metadata:")
        print(f"  Aneurysm: {(metadata['Status'] == 'Aneurysm').sum()}")
        print(f"  Control: {(metadata['Status'] == 'Control').sum()}")
    else:
        print("[WARNING] Metadata file not found, creating placeholder")
        metadata = pd.DataFrame({
            'Sample': expr_df.columns,
            'Status': ['Aneurysm'] * 23 + ['Control'] * 24
        })
    
    return expr_df, metadata

def load_tissue_expression_data():
    """Load GSE122897 tissue RNA-seq expression data"""
    print("\n" + "=" * 80)
    print("LOADING TISSUE TRANSCRIPTOMICS DATA (GSE122897)")
    print("=" * 80)
    
    # Load TPM-normalized data
    tpm_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    
    if not os.path.exists(tpm_file):
        print(f"[ERROR] Tissue expression file not found: {tpm_file}")
        return None, None
    
    print(f"Loading: {tpm_file}")
    expr_df = pd.read_csv(tpm_file, sep='\t', compression='gzip', index_col=0)
    
    print(f"[OK] Loaded tissue expression data:")
    print(f"  Genes: {expr_df.shape[0]:,}")
    print(f"  Samples: {expr_df.shape[1]}")
    
    # Load metadata
    metadata_file = "data/sc_rnaseq/GSE122897/sample_metadata.tsv"
    if os.path.exists(metadata_file):
        metadata = pd.read_csv(metadata_file, sep='\t')
        if 'Group' in metadata.columns:
            print(f"\n[OK] Loaded metadata:")
            print(metadata['Group'].value_counts().to_string())
        else:
            print("[WARNING] No 'Group' column in metadata")
            metadata['Group'] = 'Unknown'
    else:
        print("[WARNING] Metadata file not found")
        metadata = pd.DataFrame({
            'Sample': expr_df.columns,
            'Group': ['Unknown'] * len(expr_df.columns)
        })
    
    return expr_df, metadata

def load_differential_expression_results():
    """Load pre-computed differential expression results"""
    print("\n" + "=" * 80)
    print("LOADING DIFFERENTIAL EXPRESSION RESULTS")
    print("=" * 80)
    
    blood_de_file = "results/gse159610_full_de_results.tsv"
    tissue_de_file = "results/gse122897_differential_expression.tsv"
    
    blood_de = None
    tissue_de = None
    
    if os.path.exists(blood_de_file):
        print(f"Loading blood DE results: {blood_de_file}")
        blood_de = pd.read_csv(blood_de_file, sep='\t')
        print(f"  [OK] Loaded {len(blood_de):,} genes")
        if 'P_adj_BH' in blood_de.columns:
            sig_blood = (blood_de['P_adj_BH'] < 0.05).sum()
            print(f"  Significant (FDR < 0.05): {sig_blood:,} genes")
    
    if os.path.exists(tissue_de_file):
        print(f"Loading tissue DE results: {tissue_de_file}")
        tissue_de = pd.read_csv(tissue_de_file, sep='\t')
        print(f"  [OK] Loaded {len(tissue_de):,} genes")
        if 'P_adj_BH' in tissue_de.columns:
            sig_tissue = (tissue_de['P_adj_BH'] < 0.05).sum()
            print(f"  Significant (FDR < 0.05): {sig_tissue:,} genes")
    
    return blood_de, tissue_de

def identify_rupture_signature_genes(blood_de, tissue_de, fdr_threshold=0.05, top_n=50):
    """Identify rupture-associated signature genes from both blood and tissue"""
    print(f"\n" + "=" * 80)
    print(f"IDENTIFYING RUPTURE SIGNATURE GENES")
    print("=" * 80)
    
    signature_genes = []
    
    # Blood signature (aneurysm vs control)
    if blood_de is not None:
        # Check which columns exist
        log2fc_col = None
        for col in ['Log2FC', 'log2FC', 'Log2_FC']:
            if col in blood_de.columns:
                log2fc_col = col
                break
        
        p_adj_col = 'P_adj_BH' if 'P_adj_BH' in blood_de.columns else 'P_value'
        gene_col = 'Gene_ID' if 'Gene_ID' in blood_de.columns else blood_de.columns[0]
        
        if p_adj_col in blood_de.columns:
            # Use P-value threshold instead of FDR if FDR is too strict
            if fdr_threshold == 0.05:
                p_threshold = 0.01  # More lenient for small sample sizes
            else:
                p_threshold = fdr_threshold
            
            blood_sig = blood_de[
                (blood_de[p_adj_col] < p_threshold) &
                (np.abs(blood_de[log2fc_col] if log2fc_col else 0) > 0.5 if log2fc_col else True)
            ].copy()
            
            if len(blood_sig) > 0:
                blood_sig['Source'] = 'Blood'
                blood_sig['Comparison'] = 'Aneurysm_vs_Control'
                signature_genes.append(blood_sig)
                
                log2fc_vals = blood_sig[log2fc_col] if log2fc_col else 0
                print(f"\nBlood signature (P < {p_threshold}):")
                print(f"  Total significant: {len(blood_sig):,} genes")
                if log2fc_col:
                    print(f"  Upregulated in aneurysm: {(log2fc_vals > 0).sum():,}")
                    print(f"  Downregulated in aneurysm: {(log2fc_vals < 0).sum():,}")
    
    # Tissue signature (ruptured vs unruptured)
    if tissue_de is not None:
        log2fc_col = 'Log2FC' if 'Log2FC' in tissue_de.columns else None
        p_adj_col = 'P_adj_BH' if 'P_adj_BH' in tissue_de.columns else 'P_value'
        gene_col = 'Gene' if 'Gene' in tissue_de.columns else tissue_de.columns[0]
        
        if p_adj_col in tissue_de.columns:
            tissue_sig = tissue_de[
                (tissue_de[p_adj_col] < fdr_threshold) &
                (np.abs(tissue_de[log2fc_col] if log2fc_col else 0) > 0.5 if log2fc_col else True)
            ].copy()
            
            if len(tissue_sig) > 0:
                tissue_sig['Source'] = 'Tissue'
                tissue_sig['Comparison'] = 'Ruptured_vs_Unruptured'
                signature_genes.append(tissue_sig)
                
                log2fc_vals = tissue_sig[log2fc_col] if log2fc_col else 0
                print(f"\nTissue signature (FDR < {fdr_threshold}):")
                print(f"  Total significant: {len(tissue_sig):,} genes")
                if log2fc_col:
                    print(f"  Upregulated in ruptured: {(log2fc_vals > 0).sum():,}")
                    print(f"  Downregulated in ruptured: {(log2fc_vals < 0).sum():,}")
    
    if signature_genes:
        combined = pd.concat(signature_genes, ignore_index=True)
        return combined
    else:
        return pd.DataFrame()

def create_multi_omics_features(blood_expr, blood_meta, tissue_expr, tissue_meta, signature_genes):
    """Create integrated multi-omics feature matrix"""
    print(f"\n" + "=" * 80)
    print("CREATING MULTI-OMICS FEATURE MATRIX")
    print("=" * 80)
    
    # Get top signature genes
    if len(signature_genes) == 0:
        print("[WARNING] No signature genes found, using top DE genes")
        return pd.DataFrame()
    
    # Prioritize genes found in both datasets or highly significant
    top_genes = signature_genes.nlargest(100, 'P_adj_BH' if 'P_adj_BH' in signature_genes.columns else 'P_value')
    
    # Extract gene symbols/IDs
    gene_col = 'Gene' if 'Gene' in signature_genes.columns else signature_genes.columns[0]
    selected_genes = top_genes[gene_col].unique()[:50]  # Top 50 genes
    
    print(f"\nSelected {len(selected_genes)} signature genes for feature extraction")
    
    # Create feature matrix for blood samples
    blood_features = []
    blood_labels = []
    
    if blood_expr is not None and blood_meta is not None:
        print("\nExtracting blood features...")
        blood_expr_filtered = blood_expr[
            blood_expr.index.isin(selected_genes)
        ]
        
        for sample in blood_expr.columns:
            if sample in blood_meta['Sample'].values:
                sample_meta = blood_meta[blood_meta['Sample'] == sample].iloc[0]
                features = {
                    'Sample_ID': sample,
                    'Dataset': 'Blood_GSE159610',
                    'Label': 1 if sample_meta['Status'] == 'Aneurysm' else 0
                }
                
                # Add gene expression features
                for gene in selected_genes:
                    if gene in blood_expr_filtered.index:
                        features[f'Blood_{gene}'] = blood_expr_filtered.loc[gene, sample]
                    else:
                        features[f'Blood_{gene}'] = 0
                
                blood_features.append(features)
                blood_labels.append(features['Label'])
        
        print(f"  [OK] Created features for {len(blood_features)} blood samples")
    
    # Create feature matrix for tissue samples
    tissue_features = []
    
    if tissue_expr is not None and tissue_meta is not None:
        print("\nExtracting tissue features...")
        tissue_expr_filtered = tissue_expr[
            tissue_expr.index.isin(selected_genes)
        ]
        
        # Filter for ruptured vs unruptured comparison
        tissue_meta_aneurysm = tissue_meta[
            tissue_meta['Group'].isin(['Ruptured', 'Unruptured'])
        ].copy()
        
        for sample in tissue_expr.columns:
            if sample in tissue_meta_aneurysm['Sample'].values:
                sample_meta = tissue_meta_aneurysm[tissue_meta_aneurysm['Sample'] == sample].iloc[0]
                features = {
                    'Sample_ID': sample,
                    'Dataset': 'Tissue_GSE122897',
                    'Label': 1 if sample_meta['Group'] == 'Ruptured' else 0
                }
                
                # Add gene expression features
                for gene in selected_genes:
                    if gene in tissue_expr_filtered.index:
                        features[f'Tissue_{gene}'] = tissue_expr_filtered.loc[gene, sample]
                    else:
                        features[f'Tissue_{gene}'] = 0
                
                tissue_features.append(features)
        
        print(f"  [OK] Created features for {len(tissue_features)} tissue samples")
    
    # Combine features
    if blood_features and tissue_features:
        # Align features
        all_gene_cols = [f'Blood_{g}' for g in selected_genes] + [f'Tissue_{g}' for g in selected_genes]
        
        # Add missing columns
        for feat in blood_features:
            for col in all_gene_cols:
                if col not in feat:
                    feat[col] = 0
        
        for feat in tissue_features:
            for col in all_gene_cols:
                if col not in feat:
                    feat[col] = 0
        
        combined_features = pd.DataFrame(blood_features + tissue_features)
    elif blood_features:
        combined_features = pd.DataFrame(blood_features)
    elif tissue_features:
        combined_features = pd.DataFrame(tissue_features)
    else:
        return pd.DataFrame()
    
    print(f"\n[OK] Combined feature matrix:")
    print(f"  Total samples: {len(combined_features)}")
    print(f"  Total features: {len(combined_features.columns)}")
    print(f"  Positive (Aneurysm/Ruptured): {combined_features['Label'].sum()}")
    print(f"  Negative (Control/Unruptured): {(combined_features['Label'] == 0).sum()}")
    
    return combined_features

def build_prediction_model(features_df):
    """Build and validate multi-omics prediction model"""
    print(f"\n" + "=" * 80)
    print("BUILDING MULTI-OMICS PREDICTION MODEL")
    print("=" * 80)
    
    if len(features_df) == 0:
        print("[ERROR] No features available for model building")
        return None, None
    
    # Prepare data
    feature_cols = [col for col in features_df.columns 
                   if col not in ['Sample_ID', 'Dataset', 'Label']]
    X = features_df[feature_cols].values
    y = features_df['Label'].values
    
    # Handle missing values
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    print(f"\nModel features: {len(feature_cols)}")
    print(f"  Samples: {len(X_scaled)}")
    print(f"  Positive class: {y.sum()}")
    print(f"  Negative class: {(y == 0).sum()}")
    
    # Split data
    if len(features_df) < 20:
        print("[WARNING] Small sample size, using all data for training")
        X_train, X_test = X_scaled, X_scaled
        y_train, y_test = y, y
    else:
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled, y, test_size=0.2, random_state=42, stratify=y
        )
    
    # Train Gradient Boosting Classifier
    print("\nTraining Gradient Boosting Classifier...")
    gb_model = GradientBoostingClassifier(
        n_estimators=100,
        learning_rate=0.1,
        max_depth=3,
        random_state=42
    )
    gb_model.fit(X_train, y_train)
    
    # Train Random Forest Classifier
    print("Training Random Forest Classifier...")
    rf_model = RandomForestClassifier(
        n_estimators=100,
        max_depth=5,
        random_state=42
    )
    rf_model.fit(X_train, y_train)
    
    # Cross-validation
    print("\nPerforming 5-fold cross-validation...")
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    gb_cv_scores = cross_val_score(gb_model, X_scaled, y, cv=cv, scoring='roc_auc')
    rf_cv_scores = cross_val_score(rf_model, X_scaled, y, cv=cv, scoring='roc_auc')
    
    print(f"\nCross-validation AUC:")
    print(f"  Gradient Boosting: {gb_cv_scores.mean():.3f} ± {gb_cv_scores.std():.3f}")
    print(f"  Random Forest: {rf_cv_scores.mean():.3f} ± {rf_cv_scores.std():.3f}")
    
    # Test set performance
    gb_pred_proba = gb_model.predict_proba(X_test)[:, 1]
    rf_pred_proba = rf_model.predict_proba(X_test)[:, 1]
    
    if len(np.unique(y_test)) > 1:
        gb_auc = roc_auc_score(y_test, gb_pred_proba)
        rf_auc = roc_auc_score(y_test, rf_pred_proba)
        
        print(f"\nTest set AUC:")
        print(f"  Gradient Boosting: {gb_auc:.3f}")
        print(f"  Random Forest: {rf_auc:.3f}")
    
    # Feature importance
    feature_importance = pd.DataFrame({
        'Feature': feature_cols,
        'GB_Importance': gb_model.feature_importances_,
        'RF_Importance': rf_model.feature_importances_
    }).sort_values('GB_Importance', ascending=False)
    
    print(f"\nTop 10 most important features:")
    for idx, row in feature_importance.head(10).iterrows():
        print(f"  {row['Feature']:<30s} GB: {row['GB_Importance']:.4f}  RF: {row['RF_Importance']:.4f}")
    
    # Save results
    results_df = features_df.copy()
    results_df['GB_Predicted_Proba'] = gb_model.predict_proba(X_scaled)[:, 1]
    results_df['RF_Predicted_Proba'] = rf_model.predict_proba(X_scaled)[:, 1]
    results_df['GB_Predicted'] = gb_model.predict(X_scaled)
    results_df['RF_Predicted'] = rf_model.predict(X_scaled)
    
    return results_df, {
        'gb_model': gb_model,
        'rf_model': rf_model,
        'scaler': scaler,
        'feature_cols': feature_cols,
        'gb_cv_auc': gb_cv_scores.mean(),
        'rf_cv_auc': rf_cv_scores.mean()
    }

def main():
    """Main integration pipeline"""
    
    print("=" * 80)
    print("MULTI-OMICS INTEGRATION WITH REAL DATA")
    print("=" * 80)
    print("\nThis script integrates:")
    print("  1. Blood transcriptomics (GSE159610: aneurysm vs control)")
    print("  2. Tissue transcriptomics (GSE122897: ruptured vs unruptured)")
    print("  3. Signature genes from differential expression analysis")
    print("  4. Machine learning prediction models")
    
    os.makedirs("results", exist_ok=True)
    
    # Step 1: Load expression data
    blood_expr, blood_meta = load_blood_expression_data()
    tissue_expr, tissue_meta = load_tissue_expression_data()
    
    # Step 2: Load DE results
    blood_de, tissue_de = load_differential_expression_results()
    
    # Step 3: Identify signature genes
    signature_genes = identify_rupture_signature_genes(blood_de, tissue_de, fdr_threshold=0.05)
    
    if len(signature_genes) > 0:
        signature_file = "results/multi_omics_signature_genes.tsv"
        signature_genes.to_csv(signature_file, sep='\t', index=False)
        print(f"\n[OK] Saved signature genes to {signature_file}")
    
    # Step 4: Create multi-omics features
    features_df = create_multi_omics_features(
        blood_expr, blood_meta, tissue_expr, tissue_meta, signature_genes
    )
    
    if len(features_df) == 0:
        print("\n[ERROR] Could not create feature matrix")
        return
    
    features_file = "results/multi_omics_features_real_data.tsv"
    features_df.to_csv(features_file, sep='\t', index=False)
    print(f"\n[OK] Saved features to {features_file}")
    
    # Step 5: Build prediction models
    results_df, model_info = build_prediction_model(features_df)
    
    if results_df is not None:
        results_file = "results/multi_omics_predictions_real_data.tsv"
        results_df.to_csv(results_file, sep='\t', index=False)
        print(f"\n[OK] Saved predictions to {results_file}")
        
        # Save feature importance
        feature_importance = pd.DataFrame({
            'Feature': model_info['feature_cols'],
            'GB_Importance': model_info['gb_model'].feature_importances_,
            'RF_Importance': model_info['rf_model'].feature_importances_
        }).sort_values('GB_Importance', ascending=False)
        
        importance_file = "results/multi_omics_feature_importance_real_data.tsv"
        feature_importance.to_csv(importance_file, sep='\t', index=False)
        print(f"[OK] Saved feature importance to {importance_file}")
        
        # Summary
        print("\n" + "=" * 80)
        print("MULTI-OMICS INTEGRATION COMPLETE")
        print("=" * 80)
        print(f"\nSummary:")
        print(f"  Signature genes: {len(signature_genes)}")
        print(f"  Samples analyzed: {len(features_df)}")
        print(f"  Features per sample: {len(model_info['feature_cols'])}")
        print(f"  GB Cross-validation AUC: {model_info['gb_cv_auc']:.3f}")
        print(f"  RF Cross-validation AUC: {model_info['rf_cv_auc']:.3f}")
        print(f"\nTarget AUC: 0.85-0.92 (integrated model)")

if __name__ == "__main__":
    try:
        main()
        print("\n[SUCCESS] Integration complete!")
    except Exception as e:
        print(f"\n[ERROR] Integration failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

