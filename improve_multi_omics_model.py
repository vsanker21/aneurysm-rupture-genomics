"""
Immediate Improvements to Multi-Omics Model
Implements quick wins to improve AUC from 0.528
"""
import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score
from xgboost import XGBClassifier
import warnings
warnings.filterwarnings('ignore')

# Extended EC markers for better dysfunction scoring
EXTENDED_EC_MARKERS = {
    # Existing markers
    'ANGPT2': {'weight': 1.2, 'direction': 'up'},
    'ICAM1': {'weight': 1.1, 'direction': 'up'},
    'VCAM1': {'weight': 1.1, 'direction': 'up'},
    'THBS1': {'weight': 0.5, 'direction': 'up'},
    'PERP': {'weight': 0.5, 'direction': 'up'},
    'CXCL8': {'weight': 1.0, 'direction': 'up'},
    'CDH5': {'weight': -1.0, 'direction': 'down'},
    'ENG': {'weight': -1.0, 'direction': 'down'},
    'KDR': {'weight': -0.5, 'direction': 'down'},
    'TEK': {'weight': -0.5, 'direction': 'down'},
    'VWF': {'weight': -0.8, 'direction': 'down'},
    'PECAM1': {'weight': 0.3, 'direction': 'up'},
    'LYVE1': {'weight': 0.3, 'direction': 'up'},
    # Additional markers
    'SELE': {'weight': 0.8, 'direction': 'up'},  # E-selectin
    'SELP': {'weight': 0.8, 'direction': 'up'},  # P-selectin
    'NOS3': {'weight': -0.6, 'direction': 'down'},  # eNOS
    'NOS2': {'weight': 0.7, 'direction': 'up'},  # iNOS
}

# Inflammation genes
INFLAMMATION_GENES = ['IL6', 'TNF', 'IL1B', 'CRP', 'CXCL8', 'IL8', 'CXCL1', 'CXCL2']

# ECM remodeling genes
ECM_GENES = ['MMP9', 'MMP2', 'TIMP1', 'COL3A1', 'COL1A1', 'ELN', 'FBN1', 'LAMA3']

# Top DE genes from our analysis
TOP_DE_GENES = {
    'SMS': 6611,  # Spermine synthase
    'Ccd': 124220,  # Strongly upregulated
    'CPXM2': 119587,  # Downregulated
    'MARCO': 8685,  # Macrophage receptor
    'LAMA3': 3909,  # Laminin
    'TGFBI': 7045,  # Transforming growth factor
}

def load_tissue_expression():
    """Load tissue expression data"""
    expr_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    if not os.path.exists(expr_file):
        print(f"[ERROR] Expression file not found: {expr_file}")
        return None, None
    
    expr_df = pd.read_csv(expr_file, sep='\t', compression='gzip', index_col=0)
    
    # Load metadata
    metadata_file = "data/sc_rnaseq/GSE122897/sample_metadata.tsv"
    metadata = None
    if os.path.exists(metadata_file):
        metadata = pd.read_csv(metadata_file, sep='\t')
    
    return expr_df, metadata

def load_gene_mapping():
    """Load gene symbol to Entrez ID mapping"""
    mapping_file = "results/gene_symbol_to_entrez_mapping.tsv"
    if os.path.exists(mapping_file):
        mapping_df = pd.read_csv(mapping_file, sep='\t')
        return dict(zip(mapping_df['Gene_Symbol'], mapping_df['Entrez_ID'].astype(int)))
    return {}

def calculate_enhanced_ec_dysfunction(expr_df):
    """Calculate enhanced EC dysfunction score with more markers"""
    print("\nCalculating enhanced EC dysfunction score...")
    
    # Load gene mapping
    gene_mapping = load_gene_mapping()
    expr_ids = set(expr_df.index.astype(str))
    expr_ids_int = set(expr_df.index)
    
    # Find available markers
    available_markers = {}
    for marker, info in EXTENDED_EC_MARKERS.items():
        # Try to find by symbol first
        if marker in expr_df.index:
            available_markers[marker] = info
        elif marker in gene_mapping:
            # Try to find by Entrez ID from mapping
            entrez_id = gene_mapping[marker]
            if entrez_id in expr_ids_int:
                available_markers[entrez_id] = info
    
    print(f"  Found {len(available_markers)}/{len(EXTENDED_EC_MARKERS)} EC markers")
    
    if len(available_markers) < 5:
        print("  [WARNING] Limited markers found, using existing EC dysfunction scores")
        ec_file = "results/ec_dysfunction_scores.tsv"
        if os.path.exists(ec_file):
            ec_df = pd.read_csv(ec_file, sep='\t')
            return ec_df['EC_Dysfunction_Score'].values if 'EC_Dysfunction_Score' in ec_df.columns else None
        return None
    
    # Calculate enhanced score
    dysfunction_scores = np.zeros(expr_df.shape[1])
    
    for marker_id, info in available_markers.items():
        # marker_id could be symbol or Entrez ID (as int)
        # Try to find in index (which is numeric)
        try:
            if marker_id in expr_df.index:
                expr_values = np.log2(expr_df.loc[marker_id].values + 1)
            else:
                continue
                
            weight = info['weight']
            
            if info['direction'] == 'up':
                dysfunction_scores += weight * expr_values
            else:
                dysfunction_scores -= abs(weight) * expr_values
        except (KeyError, ValueError) as e:
            continue
    
    # Normalize
    dysfunction_scores = (dysfunction_scores - dysfunction_scores.mean()) / (dysfunction_scores.std() + 1e-6)
    
    return dysfunction_scores

def calculate_inflammation_score(expr_df):
    """Calculate inflammation score from tissue expression"""
    print("\nCalculating inflammation score...")
    
    # Load gene mapping
    gene_mapping = load_gene_mapping()
    expr_ids = set(expr_df.index.astype(str))
    expr_ids_int = set(expr_df.index)
    
    # Map gene symbols to Entrez IDs
    available_genes = []
    for gene_symbol in INFLAMMATION_GENES:
        if gene_symbol in expr_df.index:
            available_genes.append(gene_symbol)
        elif gene_symbol in gene_mapping:
            entrez_id = gene_mapping[gene_symbol]
            if entrez_id in expr_ids_int:
                available_genes.append(entrez_id)
    
    if len(available_genes) < 3:
        print(f"  [WARNING] Only {len(available_genes)} inflammation genes found")
        return None
    
    print(f"  Found {len(available_genes)}/{len(INFLAMMATION_GENES)} inflammation genes")
    
    # Calculate mean expression
    inflammation_expr = expr_df.loc[available_genes]
    inflammation_score = np.log2(inflammation_expr + 1).mean(axis=0).values
    
    # Normalize
    inflammation_score = (inflammation_score - inflammation_score.mean()) / (inflammation_score.std() + 1e-6)
    
    return inflammation_score

def calculate_ecm_score(expr_df):
    """Calculate ECM remodeling score"""
    print("\nCalculating ECM remodeling score...")
    
    # Load gene mapping
    gene_mapping = load_gene_mapping()
    expr_ids = set(expr_df.index.astype(str))
    expr_ids_int = set(expr_df.index)
    
    # Map gene symbols to Entrez IDs
    available_genes = []
    for gene_symbol in ECM_GENES:
        if gene_symbol in expr_df.index:
            available_genes.append(gene_symbol)
        elif gene_symbol in gene_mapping:
            entrez_id = gene_mapping[gene_symbol]
            if entrez_id in expr_ids_int:
                available_genes.append(entrez_id)
    
    if len(available_genes) < 3:
        print(f"  [WARNING] Only {len(available_genes)} ECM genes found")
        return None
    
    print(f"  Found {len(available_genes)}/{len(ECM_GENES)} ECM genes")
    
    ecm_expr = expr_df.loc[available_genes]
    ecm_score = np.log2(ecm_expr + 1).mean(axis=0).values
    
    # Normalize
    ecm_score = (ecm_score - ecm_score.mean()) / (ecm_score.std() + 1e-6)
    
    return ecm_score

def calculate_de_signature_score(expr_df):
    """Calculate signature score from top DE genes"""
    print("\nCalculating top DE gene signature score...")
    
    # Try to find top DE genes by ID
    available_genes = []
    for symbol, gene_id in TOP_DE_GENES.items():
        gene_id_str = str(gene_id)
        if gene_id_str in expr_df.index.astype(str):
            idx = expr_df.index[expr_df.index.astype(str) == gene_id_str][0]
            available_genes.append(idx)
    
    if len(available_genes) < 2:
        print(f"  [WARNING] Only {len(available_genes)} top DE genes found")
        return None
    
    print(f"  Found {len(available_genes)} top DE genes")
    
    de_expr = expr_df.loc[available_genes]
    de_score = np.log2(de_expr + 1).mean(axis=0).values
    
    # Normalize
    de_score = (de_score - de_score.mean()) / (de_score.std() + 1e-6)
    
    return de_score

def prepare_enhanced_features(expr_df, metadata):
    """Prepare enhanced feature set"""
    print("\n" + "=" * 80)
    print("PREPARING ENHANCED FEATURES")
    print("=" * 80)
    
    features = {}
    
    # 1. Enhanced EC dysfunction score
    ec_score = calculate_enhanced_ec_dysfunction(expr_df)
    if ec_score is not None:
        features['EC_Dysfunction_Score'] = ec_score
        print("  [OK] EC dysfunction score calculated")
    
    # 2. Inflammation score
    inflam_score = calculate_inflammation_score(expr_df)
    if inflam_score is not None:
        features['Inflammation_Score'] = inflam_score
        print("  [OK] Inflammation score calculated")
    
    # 3. ECM remodeling score
    ecm_score = calculate_ecm_score(expr_df)
    if ecm_score is not None:
        features['ECM_Remodeling_Score'] = ecm_score
        print("  [OK] ECM remodeling score calculated")
    
    # 4. Top DE gene signature
    de_score = calculate_de_signature_score(expr_df)
    if de_score is not None:
        features['DE_Signature_Score'] = de_score
        print("  [OK] DE signature score calculated")
    
    # Create DataFrame
    if len(features) == 0:
        print("  [ERROR] No features calculated")
        return None
    
    feature_df = pd.DataFrame(features, index=expr_df.columns)
    
    print(f"\n  [OK] Prepared {len(features)} features for {len(feature_df)} samples")
    
    return feature_df

def build_improved_model(features, labels):
    """Build improved model with multiple algorithms"""
    print("\n" + "=" * 80)
    print("BUILDING IMPROVED MODELS")
    print("=" * 80)
    
    X = features.values
    y = labels.values
    
    # Handle missing values
    X = np.nan_to_num(X, nan=0.0)
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Cross-validation
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    models = {
        'Gradient Boosting': GradientBoostingClassifier(
            n_estimators=200, max_depth=5, learning_rate=0.05, random_state=42
        ),
        'Random Forest': RandomForestClassifier(
            n_estimators=200, max_depth=10, random_state=42
        ),
        'XGBoost': XGBClassifier(
            n_estimators=200, max_depth=5, learning_rate=0.05, random_state=42
        )
    }
    
    results = {}
    
    for name, model in models.items():
        print(f"\n{name}:")
        try:
            scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')
            results[name] = {
                'mean_auc': scores.mean(),
                'std_auc': scores.std(),
                'scores': scores
            }
            print(f"  AUC: {scores.mean():.3f} ± {scores.std():.3f}")
            print(f"  Range: {scores.min():.3f} - {scores.max():.3f}")
        except Exception as e:
            print(f"  [ERROR] {e}")
            results[name] = None
    
    # Find best model
    valid_results = {k: v for k, v in results.items() if v is not None}
    if valid_results:
        best_model_name = max(valid_results, key=lambda x: valid_results[x]['mean_auc'])
        best_auc = valid_results[best_model_name]['mean_auc']
        print(f"\n  Best model: {best_model_name} (AUC = {best_auc:.3f})")
    
    return results, scaler

def main():
    """Main improvement pipeline"""
    print("=" * 80)
    print("IMPROVING MULTI-OMICS MODEL")
    print("=" * 80)
    
    os.makedirs("results", exist_ok=True)
    
    # Load data
    print("\nLoading tissue expression data...")
    expr_df, metadata = load_tissue_expression()
    if expr_df is None:
        return
    
    print(f"  Expression matrix: {expr_df.shape[0]} genes, {expr_df.shape[1]} samples")
    
    # Prepare enhanced features
    features = prepare_enhanced_features(expr_df, metadata)
    if features is None:
        return
    
    # Get labels
    if metadata is not None and 'Group' in metadata.columns:
        # Merge with metadata
        features = features.merge(
            metadata[['Sample', 'Group']].set_index('Sample'),
            left_index=True,
            right_index=True,
            how='left'
        )
        
        labels = (features['Group'] == 'Ruptured').astype(int)
        features = features.drop('Group', axis=1)
        
        print(f"\n  Labels: Ruptured={labels.sum()}, Unruptured={(labels==0).sum()}")
    else:
        print("\n  [ERROR] No group labels available")
        return
    
    # Build improved models
    results, scaler = build_improved_model(features, labels)
    
    # Save results
    results_df = pd.DataFrame({
        'Model': list(results.keys()),
        'Mean_AUC': [r['mean_auc'] if r else np.nan for r in results.values()],
        'Std_AUC': [r['std_auc'] if r else np.nan for r in results.values()]
    })
    
    results_file = "results/improved_multi_omics_results.tsv"
    results_df.to_csv(results_file, sep='\t', index=False)
    print(f"\n[OK] Saved results to {results_file}")
    
    # Save enhanced features
    features_file = "results/enhanced_features.tsv"
    features.to_csv(features_file, sep='\t')
    print(f"[OK] Saved enhanced features to {features_file}")
    
    # Summary
    print("\n" + "=" * 80)
    print("IMPROVEMENT SUMMARY")
    print("=" * 80)
    print(f"\nOriginal AUC: 0.528")
    if results:
        best_auc = max([r['mean_auc'] for r in results.values() if r is not None])
        improvement = best_auc - 0.528
        print(f"Improved AUC: {best_auc:.3f}")
        print(f"Improvement: +{improvement:.3f} ({improvement/0.528*100:.1f}% relative)")
    
    print("\n" + "=" * 80)
    print("IMPROVEMENT COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()





