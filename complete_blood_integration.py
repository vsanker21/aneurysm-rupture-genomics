"""
Complete Blood Integration Analysis
Finalize 20-30 gene signature and cross-validate performance
"""
import pandas as pd
import numpy as np
import os
from scipy import stats
from scipy.stats import ttest_ind
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve, classification_report
import warnings
warnings.filterwarnings('ignore')

def load_blood_data():
    """Load blood RNA-seq expression data"""
    blood_file = "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
    
    if not os.path.exists(blood_file):
        print(f"[ERROR] Blood transcriptomics file not found: {blood_file}")
        return None, None
    
    print(f"\nLoading blood transcriptomics data...")
    expr_df = pd.read_csv(blood_file, sep='\t', compression='gzip', index_col=0)
    print(f"  [OK] Expression matrix: {expr_df.shape[0]} genes, {expr_df.shape[1]} samples")
    
    # Load refined metadata
    metadata_file = "data/blood_rnaseq/GSE159610/sample_metadata_with_labels.tsv"
    if not os.path.exists(metadata_file):
        print(f"[ERROR] Refined metadata not found: {metadata_file}")
        return None, None
    
    metadata = pd.read_csv(metadata_file, sep='\t')
    print(f"  [OK] Metadata: {metadata.shape}")
    
    # Check sample groups
    if 'Status' in metadata.columns:
        status_counts = metadata['Status'].value_counts()
        print(f"  Sample groups: {status_counts.to_dict()}")
    
    return expr_df, metadata

def identify_rupture_signature(expr_df, metadata, min_fc=1.5, fdr_threshold=0.1):
    """
    Identify rupture-associated gene signature from blood RNA-seq
    """
    print(f"\n" + "=" * 80)
    print("IDENTIFYING RUPTURE-ASSOCIATED GENE SIGNATURE")
    print("=" * 80)
    
    # Get sample groups
    if 'Status' not in metadata.columns:
        print("[ERROR] No 'Status' column in metadata")
        return None
    
    aneurysm_samples_meta = metadata[metadata['Status'] == 'Aneurysm']['Sample'].tolist()
    control_samples_meta = metadata[metadata['Status'] == 'Control']['Sample'].tolist()
    
    if not aneurysm_samples_meta or not control_samples_meta:
        print("[ERROR] Could not identify aneurysm and control samples")
        return None
    
    print(f"\nSample groups (from metadata):")
    print(f"  Aneurysm: {len(aneurysm_samples_meta)} samples")
    print(f"  Control: {len(control_samples_meta)} samples")
    
    # Filter expression matrix to only data columns (remove annotation columns)
    expr_cols = list(expr_df.columns)
    data_cols = [c for c in expr_cols if c.startswith('DSA') or c.startswith('GSM') or 
                 (c not in ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'gene_name', 'gene_type', 'Geneid'])]
    
    if len(data_cols) == 0:
        print("[ERROR] No data columns found in expression matrix")
        return None
    
    print(f"  Expression matrix has {len(data_cols)} sample columns")
    print(f"  Sample column names: {data_cols[:5]}... (showing first 5)")
    
    # Try direct matching first
    aneurysm_samples = [s for s in aneurysm_samples_meta if s in expr_df.columns]
    control_samples = [s for s in control_samples_meta if s in expr_df.columns]
    
    # If direct matching fails, map by position
    if len(aneurysm_samples) == 0 or len(control_samples) == 0:
        print(f"  [INFO] Direct name matching failed. Attempting position-based matching...")
        
        # Sort data columns to ensure consistent order (DSA1, DSA2, ...)
        data_cols_sorted = sorted([c for c in data_cols if c.startswith('DSA')], 
                                  key=lambda x: int(x.replace('DSA', '')) if x.replace('DSA', '').isdigit() else 999)
        
        # Also include GSM columns if any
        gsm_cols = [c for c in data_cols if c.startswith('GSM')]
        data_cols_sorted = gsm_cols + data_cols_sorted
        
        # Match by position: metadata row i corresponds to expression column i
        # Ensure metadata is in the same order as DSA columns
        if len(data_cols_sorted) == len(metadata):
            print(f"  [INFO] Sample counts match. Mapping by position...")
            aneurysm_samples = [data_cols_sorted[i] for i in range(len(metadata)) 
                               if metadata.iloc[i]['Status'] == 'Aneurysm']
            control_samples = [data_cols_sorted[i] for i in range(len(metadata)) 
                             if metadata.iloc[i]['Status'] == 'Control']
            print(f"  [OK] Position-based mapping successful!")
        elif len(data_cols) >= len(metadata):
            # Try to match first N data columns to metadata
            print(f"  [WARNING] Sample count mismatch: {len(data_cols)} expression cols vs {len(metadata)} metadata rows")
            print(f"  Attempting to match first {len(metadata)} expression columns...")
            data_cols_to_use = data_cols_sorted[:len(metadata)]
            aneurysm_samples = [data_cols_to_use[i] for i in range(len(metadata)) 
                               if metadata.iloc[i]['Status'] == 'Aneurysm']
            control_samples = [data_cols_to_use[i] for i in range(len(metadata)) 
                             if metadata.iloc[i]['Status'] == 'Control']
        else:
            print(f"  [ERROR] Not enough expression columns ({len(data_cols)}) for metadata samples ({len(metadata)})")
            return None
    
    print(f"  After alignment: Aneurysm={len(aneurysm_samples)}, Control={len(control_samples)}")
    
    # Verify we have samples
    if len(aneurysm_samples) == 0 or len(control_samples) == 0:
        print("[ERROR] Could not align samples between metadata and expression matrix")
        return None
    
    # Subset expression data to aligned samples only
    all_aligned_samples = aneurysm_samples + control_samples
    expr_data = expr_df[all_aligned_samples].copy()
    
    # Perform differential expression analysis
    print(f"\nPerforming differential expression analysis...")
    print(f"  Testing {len(expr_data)} genes")
    print(f"  Using {len(aneurysm_samples) + len(control_samples)} samples")
    
    de_results = []
    
    # Limit to genes with actual expression data (skip annotation rows)
    gene_ids = [idx for idx in expr_data.index if str(idx).isdigit() or len(str(idx)) > 3]
    
    for gene_id in gene_ids[:10000]:  # Limit to first 10000 for speed
        # Get expression values
        try:
            aneurysm_expr = expr_data.loc[gene_id, aneurysm_samples].values.astype(float)
            control_expr = expr_data.loc[gene_id, control_samples].values.astype(float)
        except:
            continue
        
        # Filter out genes with too many zeros
        if np.sum(aneurysm_expr > 0) < 3 or np.sum(control_expr > 0) < 3:
            continue
        
        # Log transform
        aneurysm_expr_log = np.log2(aneurysm_expr + 1)
        control_expr_log = np.log2(control_expr + 1)
        
        # T-test
        try:
            t_stat, p_value = ttest_ind(aneurysm_expr_log, control_expr_log)
            
            # Calculate fold change
            mean_aneurysm = np.mean(aneurysm_expr)
            mean_control = np.mean(control_expr)
            log2fc = np.log2((mean_aneurysm + 1e-6) / (mean_control + 1e-6))
            
            de_results.append({
                'Gene_ID': gene_id,
                'Mean_Aneurysm': mean_aneurysm,
                'Mean_Control': mean_control,
                'Log2FC': log2fc,
                'P_value': p_value,
                'T_statistic': t_stat
            })
        except:
            continue
    
    de_df = pd.DataFrame(de_results)
    
    if len(de_df) == 0:
        print("[ERROR] No DE results generated")
        return None
    
    # Multiple testing correction
    from statsmodels.stats.multitest import multipletests
    _, p_adj, _, _ = multipletests(de_df['P_value'], method='fdr_bh')
    de_df['P_adj_BH'] = p_adj
    
    # Filter significant genes
    sig_genes = de_df[
        (de_df['P_adj_BH'] < fdr_threshold) & 
        (np.abs(de_df['Log2FC']) > np.log2(min_fc))
    ].copy()
    
    sig_genes = sig_genes.sort_values('P_adj_BH')
    
    print(f"\n[OK] Identified {len(sig_genes)} significant genes")
    print(f"  Upregulated: {(sig_genes['Log2FC'] > 0).sum()}")
    print(f"  Downregulated: {(sig_genes['Log2FC'] < 0).sum()}")
    
    return sig_genes

def select_signature_genes(de_df, n_genes=25):
    """Select top N genes for signature"""
    # Sort by significance and fold change
    de_df['Abs_Log2FC'] = np.abs(de_df['Log2FC'])
    de_df = de_df.sort_values(['P_adj_BH', 'Abs_Log2FC'], ascending=[True, False])
    
    # Select top N
    signature = de_df.head(n_genes).copy()
    
    return signature

def cross_validate_signature(expr_df, metadata, signature_genes):
    """
    Cross-validate rupture signature using machine learning
    """
    print(f"\n" + "=" * 80)
    print("CROSS-VALIDATING SIGNATURE PERFORMANCE")
    print("=" * 80)
    
    # Get sample groups from metadata
    aneurysm_samples_meta = metadata[metadata['Status'] == 'Aneurysm']['Sample'].tolist()
    control_samples_meta = metadata[metadata['Status'] == 'Control']['Sample'].tolist()
    
    # Get data columns from expression matrix
    expr_cols = list(expr_df.columns)
    data_cols = [c for c in expr_cols if c.startswith('DSA') or c.startswith('GSM') or 
                 (c not in ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'gene_name', 'gene_type'])]
    
    # Try direct matching first
    aneurysm_samples = [s for s in aneurysm_samples_meta if s in expr_df.columns]
    control_samples = [s for s in control_samples_meta if s in expr_df.columns]
    
    # If direct matching fails, map by position
    if len(aneurysm_samples) == 0 or len(control_samples) == 0:
        # Sort data columns to ensure consistent order (DSA1, DSA2, ...)
        data_cols_sorted = sorted([c for c in data_cols if c.startswith('DSA')], 
                                  key=lambda x: int(x.replace('DSA', '')) if x.replace('DSA', '').isdigit() else 999)
        gsm_cols = [c for c in data_cols if c.startswith('GSM')]
        data_cols_sorted = gsm_cols + data_cols_sorted
        
        if len(data_cols_sorted) >= len(metadata):
            data_cols_to_use = data_cols_sorted[:len(metadata)]
            aneurysm_samples = [data_cols_to_use[i] for i in range(len(metadata)) 
                               if metadata.iloc[i]['Status'] == 'Aneurysm']
            control_samples = [data_cols_to_use[i] for i in range(len(metadata)) 
                             if metadata.iloc[i]['Status'] == 'Control']
    
    # Combine all samples
    all_samples = aneurysm_samples + control_samples
    
    # Create labels (1 = aneurysm, 0 = control)
    labels = []
    for sample in all_samples:
        if sample in aneurysm_samples:
            labels.append(1)
        else:
            labels.append(0)
    
    labels = np.array(labels)
    
    # Get signature gene expression
    signature_ids = signature_genes['Gene_ID'].tolist()
    available_genes = [g for g in signature_ids if g in expr_df.index]
    
    if len(available_genes) < 5:
        print(f"[WARNING] Only {len(available_genes)} signature genes available")
        return None
    
    print(f"\nUsing {len(available_genes)} signature genes for prediction")
    
    # Extract expression matrix
    X = expr_df.loc[available_genes, all_samples].T.values
    y = labels
    
    # Log transform
    X = np.log2(X + 1)
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Random forest classifier
    rf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
    
    # 5-fold cross-validation
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = cross_val_score(rf, X_scaled, y, cv=cv, scoring='roc_auc')
    
    print(f"\nCross-Validation Results:")
    print(f"  AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
    print(f"  Range: {cv_scores.min():.3f} - {cv_scores.max():.3f}")
    
    # Train on full data for feature importance
    rf.fit(X_scaled, y)
    
    # Feature importance
    feature_importance = pd.DataFrame({
        'Gene_ID': available_genes,
        'Importance': rf.feature_importances_
    }).sort_values('Importance', ascending=False)
    
    print(f"\nTop 10 Most Important Genes:")
    for idx, row in feature_importance.head(10).iterrows():
        print(f"  {row['Gene_ID']}: {row['Importance']:.4f}")
    
    return {
        'cv_auc_mean': cv_scores.mean(),
        'cv_auc_std': cv_scores.std(),
        'cv_scores': cv_scores,
        'feature_importance': feature_importance,
        'model': rf,
        'scaler': scaler
    }

def main():
    """Main blood integration pipeline"""
    
    print("=" * 80)
    print("COMPLETE BLOOD INTEGRATION - SIGNATURE FINALIZATION")
    print("=" * 80)
    
    os.makedirs("results", exist_ok=True)
    
    # Step 1: Load blood data
    expr_df, metadata = load_blood_data()
    if expr_df is None or metadata is None:
        return
    
    # Step 2: Identify rupture signature
    de_results = identify_rupture_signature(expr_df, metadata, min_fc=1.5, fdr_threshold=0.1)
    if de_results is None or len(de_results) == 0:
        print("[ERROR] Could not identify signature genes")
        return
    
    # Save DE results
    de_file = "results/blood_de_results.tsv"
    de_results.to_csv(de_file, sep='\t', index=False)
    print(f"\n[OK] Saved DE results to {de_file}")
    
    # Step 3: Select signature genes (20-30)
    signature_25 = select_signature_genes(de_results, n_genes=25)
    signature_30 = select_signature_genes(de_results, n_genes=30)
    
    # Save signatures
    sig_25_file = "results/blood_signature_25genes.tsv"
    sig_30_file = "results/blood_signature_30genes.tsv"
    signature_25.to_csv(sig_25_file, sep='\t', index=False)
    signature_30.to_csv(sig_30_file, sep='\t', index=False)
    print(f"[OK] Saved 25-gene signature to {sig_25_file}")
    print(f"[OK] Saved 30-gene signature to {sig_30_file}")
    
    # Step 4: Cross-validate signature
    cv_results_25 = cross_validate_signature(expr_df, metadata, signature_25)
    cv_results_30 = cross_validate_signature(expr_df, metadata, signature_30)
    
    # Print summary
    print(f"\n\n" + "=" * 80)
    print("BLOOD INTEGRATION SUMMARY")
    print("=" * 80)
    print(f"\n25-gene signature:")
    if cv_results_25:
        print(f"  Cross-validation AUC: {cv_results_25['cv_auc_mean']:.3f} ± {cv_results_25['cv_auc_std']:.3f}")
    
    print(f"\n30-gene signature:")
    if cv_results_30:
        print(f"  Cross-validation AUC: {cv_results_30['cv_auc_mean']:.3f} ± {cv_results_30['cv_auc_std']:.3f}")
    
    # Select best signature
    if cv_results_25 and cv_results_30:
        if cv_results_25['cv_auc_mean'] >= cv_results_30['cv_auc_mean']:
            best_signature = signature_25
            best_cv = cv_results_25
            n_genes = 25
        else:
            best_signature = signature_30
            best_cv = cv_results_30
            n_genes = 30
        
        best_file = f"results/blood_signature_final_{n_genes}genes.tsv"
        best_signature.to_csv(best_file, sep='\t', index=False)
        print(f"\n[OK] Selected {n_genes}-gene signature (AUC: {best_cv['cv_auc_mean']:.3f})")
        print(f"  Saved to {best_file}")
    
    print("\n" + "=" * 80)
    print("BLOOD INTEGRATION COMPLETE")
    print("=" * 80)
    
    return signature_25, signature_30, cv_results_25, cv_results_30

if __name__ == "__main__":
    main()

