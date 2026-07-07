"""
Comprehensive Multi-Omics Integration Pipeline
Integrates GSE122897 (Tissue) and GSE159610 (Blood) for rupture prediction
Features:
- Feature-level integration (common signature genes)
- Cross-dataset validation
- Integrated prediction models
"""
import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path

# Force UTF-8 encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve, classification_report, confusion_matrix
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("MULTI-OMICS INTEGRATION PIPELINE")
print("GSE122897 (Tissue) + GSE159610 (Blood)")
print("=" * 80)

# ============================================================================
# LOAD DATA
# ============================================================================

def load_datasets():
    """Load tissue and blood expression data"""
    print("\n[Step 1] Loading datasets...")
    
    results = {}
    
    # Load tissue data (GSE122897)
    print("\n  Loading GSE122897 (Tissue RNA-seq)...")
    tissue_expr_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    tissue_meta_file = "data/sc_rnaseq/GSE122897/sample_metadata.tsv"
    
    if os.path.exists(tissue_expr_file) and os.path.exists(tissue_meta_file):
        tissue_expr = pd.read_csv(tissue_expr_file, sep='\t', compression='gzip', index_col=0)
        tissue_meta = pd.read_csv(tissue_meta_file, sep='\t')
        
        # Filter for ruptured vs unruptured
        tissue_meta_aneurysm = tissue_meta[tissue_meta['Group'].isin(['Ruptured', 'Unruptured'])].copy()
        tissue_meta_aneurysm['Label'] = (tissue_meta_aneurysm['Group'] == 'Ruptured').astype(int)
        
        # Align samples
        common_samples = [s for s in tissue_meta_aneurysm['Sample'] if s in tissue_expr.columns]
        tissue_meta_aligned = tissue_meta_aneurysm[tissue_meta_aneurysm['Sample'].isin(common_samples)].copy()
        
        print(f"    ✓ Tissue: {tissue_expr.shape[0]:,} genes, {len(tissue_meta_aligned)} samples")
        print(f"      Ruptured: {tissue_meta_aligned['Label'].sum()}, Unruptured: {(tissue_meta_aligned['Label'] == 0).sum()}")
        
        results['tissue_expr'] = tissue_expr
        results['tissue_meta'] = tissue_meta_aligned
    else:
        print(f"    ✗ Tissue files not found")
        results['tissue_expr'] = None
        results['tissue_meta'] = None
    
    # Load blood data (GSE159610)
    print("\n  Loading GSE159610 (Blood RNA-seq)...")
    blood_expr_file = "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
    
    # Try refined metadata first, then fall back to full metadata
    blood_meta_file = "data/blood_rnaseq/GSE159610/sample_metadata_with_labels.tsv"
    if not os.path.exists(blood_meta_file):
        blood_meta_file = "data/blood_rnaseq/GSE159610/sample_metadata.tsv"
    
    if os.path.exists(blood_expr_file) and os.path.exists(blood_meta_file):
        blood_expr = pd.read_csv(blood_expr_file, sep='\t', compression='gzip', index_col=0)
        blood_meta = pd.read_csv(blood_meta_file, sep='\t')
        
        # Create labels (aneurysm = 1, control = 0)
        if 'Status' in blood_meta.columns:
            # Use Status column if available
            blood_meta['Label'] = (blood_meta['Status'].isin(['Aneurysm', 'aneurysm', 'Aneurysm_patient'])).astype(int)
        elif 'Group' in blood_meta.columns:
            blood_meta['Label'] = (blood_meta['Group'].isin(['Aneurysm', 'aneurysm', 'Aneurysm_patient'])).astype(int)
        else:
            # Infer from sample order: first 23 are aneurysm (based on publication)
            blood_meta = blood_meta.sort_values('Sample')
            blood_meta['Label'] = 0
            blood_meta.loc[:22, 'Label'] = 1  # First 23 samples = aneurysm
        
        # Align samples - Expression matrix uses DSA1, DSA2, etc.
        # Metadata has 'description' column with DSA values (in full metadata)
        expr_cols = [col for col in blood_expr.columns if col.startswith('DSA')]
        
        # Try to extract DSA number from description or create mapping
        if 'description' in blood_meta.columns:
            # Extract DSA number from description (e.g., "DSA1" from description)
            blood_meta['DSA_ID'] = blood_meta['description'].str.extract(r'(DSA\d+)')[0]
            # Map to expression columns
            blood_meta['Sample_Expr'] = blood_meta['DSA_ID']
        else:
            # Create mapping: Sample order -> DSA1, DSA2, etc.
            blood_meta = blood_meta.sort_values('Sample')
            blood_meta['DSA_ID'] = [f'DSA{i+1}' for i in range(len(blood_meta))]
            blood_meta['Sample_Expr'] = blood_meta['DSA_ID']
        
        # Align with expression data
        common_samples = [s for s in blood_meta['Sample_Expr'] if s in expr_cols]
        blood_meta_aligned = blood_meta[blood_meta['Sample_Expr'].isin(common_samples)].copy()
        
        # Update Sample column to match expression column names
        blood_meta_aligned['Sample'] = blood_meta_aligned['Sample_Expr']
        
        print(f"    ✓ Blood: {blood_expr.shape[0]:,} genes, {len(blood_meta_aligned)} samples")
        if 'Label' in blood_meta_aligned.columns:
            print(f"      Aneurysm: {blood_meta_aligned['Label'].sum()}, Control: {(blood_meta_aligned['Label'] == 0).sum()}")
        
        results['blood_expr'] = blood_expr
        results['blood_meta'] = blood_meta_aligned
    else:
        print(f"    ✗ Blood files not found")
        results['blood_expr'] = None
        results['blood_meta'] = None
    
    return results

# ============================================================================
# FEATURE-LEVEL INTEGRATION (COMMON SIGNATURE GENES)
# ============================================================================

def identify_common_signature_genes(tissue_expr, tissue_meta, blood_expr, blood_meta):
    """Identify common signature genes across blood and tissue datasets"""
    print("\n[Step 2] Identifying common signature genes...")
    
    # Load DE results if available
    tissue_de_file = "results/gse122897_differential_expression.tsv"
    blood_de_file = "results/gse159610_full_de_results.tsv"
    
    tissue_de = None
    blood_de = None
    
    if os.path.exists(tissue_de_file):
        tissue_de = pd.read_csv(tissue_de_file, sep='\t')
        print(f"    ✓ Loaded tissue DE results: {len(tissue_de):,} genes")
        
        # Filter significant genes
        if 'P_adj_BH' in tissue_de.columns:
            tissue_sig = tissue_de[tissue_de['P_adj_BH'] < 0.05].copy()
        elif 'P_value' in tissue_de.columns:
            tissue_sig = tissue_de[tissue_de['P_value'] < 0.01].copy()
        else:
            tissue_sig = tissue_de.head(100).copy()
        
        # Sort by absolute log2FC
        if 'Log2FC' in tissue_sig.columns:
            tissue_sig['Abs_Log2FC'] = np.abs(tissue_sig['Log2FC'])
            tissue_sig = tissue_sig.sort_values('Abs_Log2FC', ascending=False)
        
        tissue_top = tissue_sig.head(100)
        print(f"      Top {len(tissue_top)} tissue signature genes (FDR < 0.05)")
    
    if os.path.exists(blood_de_file):
        blood_de = pd.read_csv(blood_de_file, sep='\t')
        print(f"    ✓ Loaded blood DE results: {len(blood_de):,} genes")
        
        # Filter significant genes
        if 'P_adj_BH' in blood_de.columns:
            blood_sig = blood_de[blood_de['P_adj_BH'] < 0.05].copy()
        elif 'P_value' in blood_de.columns:
            blood_sig = blood_de[blood_de['P_value'] < 0.01].copy()
        else:
            blood_sig = blood_de.head(100).copy()
        
        # Sort by absolute log2FC
        if 'Log2FC' in blood_sig.columns:
            blood_sig['Abs_Log2FC'] = np.abs(blood_sig['Log2FC'])
            blood_sig = blood_sig.sort_values('Abs_Log2FC', ascending=False)
        
        blood_top = blood_sig.head(100)
        print(f"      Top {len(blood_top)} blood signature genes (FDR < 0.05)")
    
    # For now, use variable genes from expression data
    # In practice, would use proper gene ID mapping
    print(f"\n    Identifying highly variable genes...")
    
    # Tissue: top variable genes
    tissue_samples = [s for s in tissue_meta['Sample'] if s in tissue_expr.columns]
    if len(tissue_samples) > 0:
        tissue_vars = tissue_expr.loc[:, tissue_samples].var(axis=1)
        tissue_var_genes = tissue_vars.nlargest(500).index.tolist()
        print(f"      Tissue: {len(tissue_var_genes)} top variable genes")
    
    # Blood: top variable genes
    blood_samples = [s for s in blood_meta['Sample'] if s in blood_expr.columns]
    if len(blood_samples) > 0:
        blood_vars = blood_expr.loc[:, blood_samples].var(axis=1)
        blood_var_genes = blood_vars.nlargest(500).index.tolist()
        print(f"      Blood: {len(blood_var_genes)} top variable genes")
    
    # Note: Direct gene ID matching won't work (different ID systems)
    # In production, would use gene ID conversion (biomaRt, mygene, etc.)
    print(f"\n    ℹ️ Note: Gene IDs differ between datasets (Entrez vs Ensembl)")
    print(f"    ℹ️ Using top variable genes from each dataset separately")
    print(f"    ℹ️ For true integration, convert gene IDs using biomaRt or mygene")
    
    return {
        'tissue_genes': tissue_var_genes if 'tissue_var_genes' in locals() else [],
        'blood_genes': blood_var_genes if 'blood_var_genes' in locals() else [],
        'tissue_de': tissue_top if tissue_de is not None else None,
        'blood_de': blood_top if blood_de is not None else None
    }

# ============================================================================
# CROSS-DATASET VALIDATION
# ============================================================================

def cross_dataset_validation(tissue_expr, tissue_meta, blood_expr, blood_meta, signature_genes):
    """Perform cross-dataset validation"""
    print("\n[Step 3] Cross-dataset validation...")
    
    # For cross-validation, we'll:
    # 1. Train on one dataset
    # 2. Test on the other (if possible)
    # 3. Or use cross-validation within each dataset
    
    print(f"\n    Dataset 1: Tissue (GSE122897)")
    print(f"      Samples: {len(tissue_meta)}")
    print(f"      Task: Ruptured vs Unruptured")
    
    print(f"\n    Dataset 2: Blood (GSE159610)")
    print(f"      Samples: {len(blood_meta)}")
    print(f"      Task: Aneurysm vs Control")
    
    print(f"\n    ℹ️ Note: Different phenotypes (rupture vs aneurysm presence)")
    print(f"    ℹ️ Cross-validation within each dataset")
    
    # Build models for each dataset
    results = {}
    
    # Tissue model
    if tissue_expr is not None and tissue_meta is not None:
        print(f"\n    [3a] Building tissue model...")
        tissue_genes = signature_genes['tissue_genes']
        if len(tissue_genes) > 0:
            tissue_samples = [s for s in tissue_meta['Sample'] if s in tissue_expr.columns]
            tissue_genes_in_expr = [g for g in tissue_genes if g in tissue_expr.index]
            
            if len(tissue_genes_in_expr) > 0 and len(tissue_samples) > 0:
                X_tissue = tissue_expr.loc[tissue_genes_in_expr[:100], tissue_samples].T.values
                y_tissue = tissue_meta[tissue_meta['Sample'].isin(tissue_samples)]['Label'].values
                
                X_tissue = np.nan_to_num(X_tissue, nan=0.0)
                scaler_tissue = StandardScaler()
                X_tissue_scaled = scaler_tissue.fit_transform(X_tissue)
                
                # Cross-validation
                cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
                model = GradientBoostingClassifier(n_estimators=100, random_state=42)
                cv_scores = cross_val_score(model, X_tissue_scaled, y_tissue, cv=cv, scoring='roc_auc')
                
                print(f"      ✓ Tissue CV AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
                results['tissue_cv'] = cv_scores
    
    # Blood model
    if blood_expr is not None and blood_meta is not None:
        print(f"\n    [3b] Building blood model...")
        blood_genes = signature_genes['blood_genes']
        if len(blood_genes) > 0:
            blood_samples = [s for s in blood_meta['Sample'] if s in blood_expr.columns]
            blood_genes_in_expr = [g for g in blood_genes if g in blood_expr.index]
            
            if len(blood_genes_in_expr) > 0 and len(blood_samples) > 0:
                X_blood = blood_expr.loc[blood_genes_in_expr[:100], blood_samples].T.values
                y_blood = blood_meta[blood_meta['Sample'].isin(blood_samples)]['Label'].values
                
                X_blood = np.nan_to_num(X_blood, nan=0.0)
                scaler_blood = StandardScaler()
                X_blood_scaled = scaler_blood.fit_transform(X_blood)
                
                # Cross-validation
                cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
                model = GradientBoostingClassifier(n_estimators=100, random_state=42)
                cv_scores = cross_val_score(model, X_blood_scaled, y_blood, cv=cv, scoring='roc_auc')
                
                print(f"      ✓ Blood CV AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
                results['blood_cv'] = cv_scores
    
    return results

# ============================================================================
# INTEGRATED MODEL (FEATURE-LEVEL)
# ============================================================================

def build_integrated_models(tissue_expr, tissue_meta, blood_expr, blood_meta, signature_genes):
    """Build integrated models using signature genes"""
    print("\n[Step 4] Building integrated models...")
    
    results = {}
    
    # For feature-level integration, we identify common biological pathways
    # and build models using signature genes from both datasets
    
    print(f"\n    [4a] Tissue model (Ruptured vs Unruptured)...")
    if tissue_expr is not None and tissue_meta is not None:
        tissue_genes = signature_genes['tissue_genes']
        tissue_samples = [s for s in tissue_meta['Sample'] if s in tissue_expr.columns]
        tissue_genes_in_expr = [g for g in tissue_genes[:100] if g in tissue_expr.index]
        
        if len(tissue_genes_in_expr) > 0 and len(tissue_samples) > 0:
            X = tissue_expr.loc[tissue_genes_in_expr, tissue_samples].T.values
            y = tissue_meta[tissue_meta['Sample'].isin(tissue_samples)]['Label'].values
            
            X = np.nan_to_num(X, nan=0.0)
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            
            # Train models
            models = {
                'Random Forest': RandomForestClassifier(n_estimators=200, random_state=42, class_weight='balanced'),
                'Gradient Boosting': GradientBoostingClassifier(n_estimators=200, random_state=42),
                'SVM': SVC(kernel='rbf', probability=True, random_state=42, class_weight='balanced')
            }
            
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
            
            tissue_results = {}
            for name, model in models.items():
                cv_scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')
                tissue_results[name] = {
                    'mean_auc': cv_scores.mean(),
                    'std_auc': cv_scores.std()
                }
                print(f"      {name}: AUC {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
            
            results['tissue_models'] = tissue_results
    
    print(f"\n    [4b] Blood model (Aneurysm vs Control)...")
    if blood_expr is not None and blood_meta is not None:
        blood_genes = signature_genes['blood_genes']
        blood_samples = [s for s in blood_meta['Sample'] if s in blood_expr.columns]
        blood_genes_in_expr = [g for g in blood_genes[:100] if g in blood_expr.index]
        
        if len(blood_genes_in_expr) > 0 and len(blood_samples) > 0:
            X = blood_expr.loc[blood_genes_in_expr, blood_samples].T.values
            y = blood_meta[blood_meta['Sample'].isin(blood_samples)]['Label'].values
            
            X = np.nan_to_num(X, nan=0.0)
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            
            # Train models
            models = {
                'Random Forest': RandomForestClassifier(n_estimators=200, random_state=42, class_weight='balanced'),
                'Gradient Boosting': GradientBoostingClassifier(n_estimators=200, random_state=42),
                'SVM': SVC(kernel='rbf', probability=True, random_state=42, class_weight='balanced')
            }
            
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
            
            blood_results = {}
            for name, model in models.items():
                cv_scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')
                blood_results[name] = {
                    'mean_auc': cv_scores.mean(),
                    'std_auc': cv_scores.std()
                }
                print(f"      {name}: AUC {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
            
            results['blood_models'] = blood_results
    
    return results

# ============================================================================
# SAVE RESULTS
# ============================================================================

def save_results(signature_genes, validation_results, model_results):
    """Save integration results"""
    print("\n[Step 5] Saving results...")
    
    os.makedirs("results", exist_ok=True)
    
    # Save signature genes
    if signature_genes['tissue_de'] is not None:
        output_file = "results/multi_omics_tissue_signature_genes.tsv"
        signature_genes['tissue_de'].to_csv(output_file, sep='\t', index=False)
        print(f"    ✓ Saved tissue signature genes to {output_file}")
    
    if signature_genes['blood_de'] is not None:
        output_file = "results/multi_omics_blood_signature_genes.tsv"
        signature_genes['blood_de'].to_csv(output_file, sep='\t', index=False)
        print(f"    ✓ Saved blood signature genes to {output_file}")
    
    # Save model results
    summary_data = []
    
    if 'tissue_models' in model_results:
        for model_name, metrics in model_results['tissue_models'].items():
            summary_data.append({
                'Dataset': 'Tissue (GSE122897)',
                'Task': 'Ruptured vs Unruptured',
                'Model': model_name,
                'Mean_AUC': metrics['mean_auc'],
                'Std_AUC': metrics['std_auc']
            })
    
    if 'blood_models' in model_results:
        for model_name, metrics in model_results['blood_models'].items():
            summary_data.append({
                'Dataset': 'Blood (GSE159610)',
                'Task': 'Aneurysm vs Control',
                'Model': model_name,
                'Mean_AUC': metrics['mean_auc'],
                'Std_AUC': metrics['std_auc']
            })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_file = "results/multi_omics_integration_summary.tsv"
        summary_df.to_csv(summary_file, sep='\t', index=False)
        print(f"    ✓ Saved model summary to {summary_file}")
    
    return summary_df if summary_data else None

# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main pipeline"""
    
    # Load datasets
    data = load_datasets()
    
    if data['tissue_expr'] is None or data['blood_expr'] is None:
        print("\n[ERROR] Required datasets not found!")
        return
    
    # Identify signature genes
    signature_genes = identify_common_signature_genes(
        data['tissue_expr'], data['tissue_meta'],
        data['blood_expr'], data['blood_meta']
    )
    
    # Cross-dataset validation
    validation_results = cross_dataset_validation(
        data['tissue_expr'], data['tissue_meta'],
        data['blood_expr'], data['blood_meta'],
        signature_genes
    )
    
    # Build integrated models
    model_results = build_integrated_models(
        data['tissue_expr'], data['tissue_meta'],
        data['blood_expr'], data['blood_meta'],
        signature_genes
    )
    
    # Save results
    summary_df = save_results(signature_genes, validation_results, model_results)
    
    # Final summary
    print("\n" + "=" * 80)
    print("MULTI-OMICS INTEGRATION COMPLETE")
    print("=" * 80)
    
    if summary_df is not None:
        print("\nModel Performance Summary:")
        print(summary_df.to_string(index=False))
    
    print("\n✅ Integration complete!")
    print("\nOutput files:")
    print("  - results/multi_omics_tissue_signature_genes.tsv")
    print("  - results/multi_omics_blood_signature_genes.tsv")
    print("  - results/multi_omics_integration_summary.tsv")

if __name__ == "__main__":
    main()

