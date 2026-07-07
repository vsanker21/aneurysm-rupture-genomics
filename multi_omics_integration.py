"""
Multi-Omics Integration
Combines PRS + Blood RNA + EC Dysfunction + Clinical factors
Trains integrated model with SHAP interpretation
"""
import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve, classification_report, confusion_matrix
import warnings
from typing import List, Optional
warnings.filterwarnings('ignore')

try:
    import shap
    SHAP_AVAILABLE = True
except (ImportError, KeyboardInterrupt, Exception) as e:
    SHAP_AVAILABLE = False
    print(f"[WARNING] SHAP not available or failed to load: {e}")
    print("  Continuing without SHAP interpretation")

def load_prs_weights():
    """Load PRS weights"""
    prs_file = "results/prs_weights.tsv"
    if os.path.exists(prs_file):
        return pd.read_csv(prs_file, sep='\t')
    return None

def load_blood_signature():
    """Load blood signature genes"""
    sig_file = "results/blood_signature_final_25genes.tsv"
    if not os.path.exists(sig_file):
        sig_file = "results/blood_signature_25genes.tsv"
    if not os.path.exists(sig_file):
        sig_file = "results/blood_signature_30genes.tsv"
    
    if os.path.exists(sig_file):
        return pd.read_csv(sig_file, sep='\t')
    return None

def load_ec_dysfunction_scores():
    """Load EC dysfunction scores"""
    ec_file = "results/ec_dysfunction_scores.tsv"
    if os.path.exists(ec_file):
        return pd.read_csv(ec_file, sep='\t')
    return None

def load_tissue_expression():
    """Load tissue expression for blood signature genes"""
    expr_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    if os.path.exists(expr_file):
        return pd.read_csv(expr_file, sep='\t', compression='gzip', index_col=0)
    return None

# -----------------------------
# New: Tissue-derived feature builders
# -----------------------------

INFLAMMATION_GENES = ['IL6', 'TNF', 'IL1B', 'CRP', 'CXCL8', 'IL8', 'CXCL1', 'CXCL2']
ECM_GENES = ['MMP9', 'MMP2', 'TIMP1', 'COL3A1', 'COL1A1', 'ELN', 'FBN1', 'LAMA3']
# Top DE gene IDs observed in results
TOP_DE_GENE_IDS = [6611, 124220, 119587, 8685, 3909, 7045]

GENE_MAPPING_FILE = os.path.join("results", "gene_symbol_to_entrez_mapping.tsv")
_GENE_SYMBOL_TO_ENTREZ = None


def _load_local_gene_mapping():
    """Load cached gene symbol → Entrez ID mappings from disk if available."""
    global _GENE_SYMBOL_TO_ENTREZ
    if _GENE_SYMBOL_TO_ENTREZ is not None:
        return _GENE_SYMBOL_TO_ENTREZ

    mapping = {}
    if os.path.exists(GENE_MAPPING_FILE):
        try:
            mapping_df = pd.read_csv(GENE_MAPPING_FILE, sep='\t')
            if 'Gene_Symbol' in mapping_df.columns and 'Entrez_ID' in mapping_df.columns:
                mapping_df = mapping_df.dropna(subset=['Gene_Symbol', 'Entrez_ID'])
                mapping = {
                    str(sym).strip(): str(int(entrez))
                    for sym, entrez in zip(mapping_df['Gene_Symbol'], mapping_df['Entrez_ID'])
                    if str(sym).strip() != ''
                }
                print(f"  [OK] Loaded {len(mapping)} gene symbol mappings from {GENE_MAPPING_FILE}")
        except Exception as exc:
            print(f"  [WARNING] Failed to load local gene mapping ({exc})")

    _GENE_SYMBOL_TO_ENTREZ = mapping
    return _GENE_SYMBOL_TO_ENTREZ


ENHANCED_EC_MARKERS = {
    'ANGPT2': 1.2, 'ICAM1': 1.1, 'VCAM1': 1.1, 'THBS1': 0.5, 'PERP': 0.5, 'CXCL8': 1.0,
    'CDH5': -1.0, 'ENG': -1.0, 'KDR': -0.5, 'TEK': -0.5, 'VWF': -0.8, 'PECAM1': 0.3, 'LYVE1': 0.3,
    'SELE': 0.8, 'SELP': 0.8, 'NOS3': -0.6, 'NOS2': 0.7
}

def _map_symbols_to_expr_ids(expr_df: pd.DataFrame, symbols: List[str]) -> List[str]:
    """
    Map gene symbols to IDs present in expr_df.index using available resources.
    Priority order:
      1. Direct symbol match (if the expression matrix already uses symbols)
      2. Local curated mapping file (results/gene_symbol_to_entrez_mapping.tsv)
      3. Online MyGene.info lookup (if network is available)
    Returns a list of index entries (as strings) that exist in expr_df.
    """
    if expr_df is None or len(expr_df.index) == 0:
        return []

    idx = expr_df.index
    idx_str = idx.astype(str).tolist()
    hits = []

    # 1) Direct matches
    for symbol in symbols:
        if symbol in idx:
            hits.append(symbol)

    # Prepare set for quick lookup
    idx_str_set = set(idx_str)

    # 2) Local mapping fallback (offline-friendly)
    local_mapping = _load_local_gene_mapping()
    for symbol in symbols:
        if symbol in local_mapping:
            entrez = local_mapping[symbol]
            if entrez in idx_str_set:
                # Capture the exact index token to preserve dtype
                match_positions = [i for i, token in enumerate(idx_str) if token == entrez]
                if match_positions:
                    hits.append(idx[match_positions[0]])

    # Determine which symbols still need mapping
    mapped_symbols = set(symbols) & set([str(h) for h in hits])
    remaining = [s for s in symbols if s not in mapped_symbols]

    # 3) Online lookup (best effort, skipped if unavailable)
    if remaining:
        try:
            import requests
            import time
            url = "https://mygene.info/v3/query"
            for symbol in remaining:
                params = {
                    'q': symbol,
                    'scopes': 'symbol',
                    'fields': 'entrezgene',
                    'species': 'human'
                }
                resp = requests.get(url, params=params, timeout=10)
                if resp.status_code == 200:
                    data = resp.json()
                    if isinstance(data, dict):
                        data = [data]
                    for item in data:
                        eg = item.get('entrezgene')
                        if eg is None:
                            continue
                        eg_str = str(eg)
                        if eg_str in idx_str_set:
                            match_positions = [i for i, token in enumerate(idx_str) if token == eg_str]
                            if match_positions:
                                hits.append(idx[match_positions[0]])
                                break
                time.sleep(0.1)
        except Exception as exc:
            print(f"  [WARNING] Online gene mapping lookup failed ({exc})")

    # Ensure unique and preserve order encountered
    unique_hits = []
    seen = set()
    for h in hits:
        h_str = str(h)
        if h_str in idx_str_set and h_str not in seen:
            unique_hits.append(h)
            seen.add(h_str)

    return unique_hits

def _mean_score(expr_df: pd.DataFrame, gene_list: List[str]) -> Optional[pd.Series]:
    # Map provided symbols to actual expr index entries (Entrez IDs or symbols)
    available = _map_symbols_to_expr_ids(expr_df, gene_list)
    if len(available) < 2:
        return None
    score = np.log2(expr_df.loc[available] + 1).mean(axis=0)
    return (score - score.mean()) / (score.std() + 1e-6)

def build_inflammation_score(expr_df: pd.DataFrame) -> Optional[pd.Series]:
    return _mean_score(expr_df, INFLAMMATION_GENES)

def build_ecm_score(expr_df: pd.DataFrame) -> Optional[pd.Series]:
    return _mean_score(expr_df, ECM_GENES)

def build_de_signature_score(expr_df: pd.DataFrame) -> Optional[pd.Series]:
    # indices may be Entrez IDs; attempt by string matching
    idx_str_set = set(expr_df.index.astype(str))
    hit_ids = []
    for gid in TOP_DE_GENE_IDS:
        gid_str = str(gid)
        if gid_str in idx_str_set:
            # Find the actual index value that matches (could be string or numeric)
            matches = [idx for idx in expr_df.index if str(idx) == gid_str]
            if matches:
                hit_ids.append(matches[0])
    if len(hit_ids) < 2:
        return None
    score = np.log2(expr_df.loc[hit_ids] + 1).mean(axis=0)
    return (score - score.mean()) / (score.std() + 1e-6)

def build_enhanced_ec_dysfunction(expr_df: pd.DataFrame) -> Optional[pd.Series]:
    accum = None
    found_any = False
    # Map all marker symbols to actual expr indices first
    mapped = {}
    for gene in ENHANCED_EC_MARKERS.keys():
        hits = _map_symbols_to_expr_ids(expr_df, [gene])
        if hits:
            mapped[gene] = hits[0]  # use the first valid mapping
    for gene, weight in ENHANCED_EC_MARKERS.items():
        idx_token = mapped.get(gene, None)
        if idx_token is not None:
            vals = np.log2(expr_df.loc[idx_token] + 1)
            term = weight * vals if weight > 0 else -abs(weight) * vals
            accum = term if accum is None else accum + term
            found_any = True
    if not found_any:
        return None
    return (accum - accum.mean()) / (accum.std() + 1e-6)

def prepare_integrated_features(metadata, expr_df=None, blood_sig=None, ec_scores=None):
    """
    Prepare integrated feature matrix combining all omics layers
    """
    print(f"\n" + "=" * 80)
    print("PREPARING INTEGRATED FEATURES")
    print("=" * 80)
    
    # Start with metadata
    features = metadata.copy()
    debug_rows = []  # collect mapping diagnostics
    
    # Layer 1: PRS - Only include if we have real PRS data
    # PRS requires genotype data which we don't have for these samples
    # We'll exclude PRS from the model since we don't have real genotype data
    if 'PRS' not in features.columns:
        print(f"  [NOTE] PRS not available - requires genotype data (excluded from model)")
        # Don't add placeholder - we'll work without PRS
    
    # Layer 2: Blood RNA Signature Score
    if expr_df is not None and blood_sig is not None:
        print(f"\n  Calculating Blood RNA Signature Score...")
        blood_genes = blood_sig['Gene_ID'].tolist()
        available_genes = [g for g in blood_genes if str(g) in expr_df.index.astype(str)]
        
        if len(available_genes) > 0:
            print(f"    Found {len(available_genes)} out of {len(blood_genes)} signature genes")
            
            # Get expression for signature genes
            blood_expr = expr_df.loc[available_genes].T
            blood_expr = np.log2(blood_expr + 1)
            
            # Calculate signature score as mean of signature genes
            features['Blood_RNA_Score'] = blood_expr.mean(axis=1).values
            print(f"    [OK] Blood RNA Score calculated")
        else:
            features['Blood_RNA_Score'] = 0
            print(f"    [WARNING] No signature genes found in expression data")
    else:
        features['Blood_RNA_Score'] = 0
        print(f"  [WARNING] Blood signature not available")
    
    # Layer 3: EC Dysfunction Score
    if ec_scores is not None:
        print(f"\n  Adding EC Dysfunction Score...")
        # Merge EC scores
        features = features.merge(
            ec_scores[['Sample', 'EC_Dysfunction_Score']],
            on='Sample',
            how='left'
        )
        features['EC_Dysfunction_Score'] = features['EC_Dysfunction_Score'].fillna(0)
        print(f"    [OK] EC Dysfunction Score added")
    else:
        # Try to compute enhanced EC dysfunction directly from expression
        if expr_df is not None:
            print(f"\n  Computing Enhanced EC Dysfunction Score from expression...")
            ec_enhanced = build_enhanced_ec_dysfunction(expr_df)
            if ec_enhanced is not None:
                ec_df = ec_enhanced.reset_index()
                ec_df.columns = ['Sample', 'EC_Dysfunction_Score']
                features = features.merge(ec_df, on='Sample', how='left')
                features['EC_Dysfunction_Score'] = features['EC_Dysfunction_Score'].fillna(0)
                print(f"    [OK] Enhanced EC Dysfunction Score added")
            else:
                features['EC_Dysfunction_Score'] = 0
                print(f"  [WARNING] Could not compute EC dysfunction from expression")
        else:
            features['EC_Dysfunction_Score'] = 0
            print(f"  [WARNING] EC dysfunction scores not available")

    # New: add tissue-derived features if expression is available
    if expr_df is not None:
        print(f"\n  Calculating additional tissue-derived scores...")
        # Mapping diagnostics
        infl_hits = _map_symbols_to_expr_ids(expr_df, INFLAMMATION_GENES)
        ecm_hits = _map_symbols_to_expr_ids(expr_df, ECM_GENES)
        # DE hits via Entrez IDs present
        idx_str_set = set(expr_df.index.astype(str))
        de_hits = [str(gid) for gid in TOP_DE_GENE_IDS if str(gid) in idx_str_set]
        # EC enhanced mapped
        ec_enh_mapped = []
        for g in ENHANCED_EC_MARKERS.keys():
            h = _map_symbols_to_expr_ids(expr_df, [g])
            if h:
                ec_enh_mapped.append(h[0])
        # Append to debug_rows
        debug_rows.append({'Feature': 'Inflammation_Score', 'Requested': len(INFLAMMATION_GENES), 'Mapped': len(infl_hits), 'Mapped_Names': ",".join(map(str, infl_hits))})
        debug_rows.append({'Feature': 'ECM_Remodeling_Score', 'Requested': len(ECM_GENES), 'Mapped': len(ecm_hits), 'Mapped_Names': ",".join(map(str, ecm_hits))})
        debug_rows.append({'Feature': 'DE_Signature_Score', 'Requested': len(TOP_DE_GENE_IDS), 'Mapped': len(de_hits), 'Mapped_Names': ",".join(de_hits)})
        debug_rows.append({'Feature': 'EC_Dysfunction_Score_Enhanced', 'Requested': len(ENHANCED_EC_MARKERS), 'Mapped': len(ec_enh_mapped), 'Mapped_Names': ",".join(map(str, ec_enh_mapped))})

        # Ensure columns align: expr_df columns are samples
        infl = build_inflammation_score(expr_df)
        if infl is not None:
            infl_df = infl.reset_index()
            infl_df.columns = ['Sample', 'Inflammation_Score']
            features = features.merge(infl_df, on='Sample', how='left')
            print(f"    [OK] Inflammation Score added")
        ecm = build_ecm_score(expr_df)
        if ecm is not None:
            ecm_df = ecm.reset_index()
            ecm_df.columns = ['Sample', 'ECM_Remodeling_Score']
            features = features.merge(ecm_df, on='Sample', how='left')
            print(f"    [OK] ECM Remodeling Score added")
        de_sig = build_de_signature_score(expr_df)
        if de_sig is not None:
            de_df = de_sig.reset_index()
            de_df.columns = ['Sample', 'DE_Signature_Score']
            features = features.merge(de_df, on='Sample', how='left')
            print(f"    [OK] DE Signature Score added")
        # Save mapping diagnostics
        try:
            os.makedirs("results", exist_ok=True)
            debug_df = pd.DataFrame(debug_rows)
            # Write TSV
            debug_df.to_csv("results/integration_feature_debug.tsv", sep='\t', index=False)
            # Also append a human-readable section to integration_analysis_report.txt
            with open("results/integration_analysis_report.txt", "a", encoding="utf-8") as f:
                f.write("\n=== Feature Mapping Diagnostics ===\n")
                for _, row in debug_df.iterrows():
                    f.write(f"{row['Feature']}: requested={row['Requested']}, mapped={row['Mapped']}\n")
                f.write("Mapped names:\n")
                for _, row in debug_df.iterrows():
                    f.write(f"  {row['Feature']}: {row['Mapped_Names']}\n")
            print(f"    [OK] Saved feature mapping diagnostics to results/integration_feature_debug.tsv and appended to integration_analysis_report.txt")
        except Exception as e:
            print(f"    [WARNING] Could not save mapping diagnostics: {e}")
    
    # Layer 4: Clinical Factors (if available)
    clinical_cols = []
    if 'Age' in features.columns:
        features['Clinical_Age'] = features['Age']
        clinical_cols.append('Clinical_Age')
    if 'Sex' in features.columns or 'Gender' in features.columns:
        sex_col = 'Sex' if 'Sex' in features.columns else 'Gender'
        features['Clinical_Female'] = (features[sex_col] == 'Female').astype(int)
        clinical_cols.append('Clinical_Female')
    if 'Hypertension' in features.columns:
        # Try to map textual HTN to binary
        if features['Hypertension'].dtype == object:
            features['Clinical_HTN'] = features['Hypertension'].str.contains('1|yes|true|hypertension', case=False, na=False).astype(int)
        else:
            features['Clinical_HTN'] = features['Hypertension'].fillna(0).astype(int)
        clinical_cols.append('Clinical_HTN')
    if 'Size' in features.columns:
        # Attempt to parse numeric size if string like '5 mm'
        size_numeric = pd.to_numeric(features['Size'].astype(str).str.extract(r'([0-9]+\\.?[0-9]*)')[0], errors='coerce')
        features['Clinical_SizeMM'] = size_numeric
        clinical_cols.append('Clinical_SizeMM')
    
    if len(clinical_cols) > 0:
        features['Clinical_Score'] = features[clinical_cols].mean(axis=1)
    else:
        features['Clinical_Score'] = 0
    
    # Select feature columns (exclude PRS if not available - we don't use placeholders)
    feature_cols = ['Blood_RNA_Score', 'EC_Dysfunction_Score', 'Inflammation_Score', 'ECM_Remodeling_Score', 'DE_Signature_Score', 'Clinical_Score']
    # Only include PRS if it's real data (not placeholder)
    if 'PRS' in features.columns and not features['PRS'].isna().all():
        # Check if PRS looks like real data (not all zeros or random)
        if features['PRS'].std() > 0.1:  # Real PRS should have variation
            feature_cols.insert(0, 'PRS')
    
    available_features = [f for f in feature_cols if f in features.columns and features[f].notna().any()]
    
    print(f"\n  [OK] Integrated features prepared:")
    print(f"    Available features: {available_features}")
    
    return features, available_features

def build_integrated_model(features, feature_cols, labels):
    """
    Train integrated gradient boosting model
    """
    print(f"\n" + "=" * 80)
    print("TRAINING INTEGRATED MODEL")
    print("=" * 80)
    
    # Prepare data
    X = features[feature_cols].values
    y = labels.values if hasattr(labels, 'values') else np.array(labels)
    
    # Handle missing values
    X = np.nan_to_num(X, nan=0.0)
    
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    print(f"\n  Training on {len(X)} samples with {len(feature_cols)} features")
    print(f"    Features: {', '.join(feature_cols)}")
    print(f"    Classes: {np.unique(y, return_counts=True)}")
    
    # Gradient boosting classifier
    # Note: GradientBoostingClassifier doesn't support class_weight
    # Use sample_weight in fit() if needed for imbalanced classes
    model = GradientBoostingClassifier(
        n_estimators=200,
        max_depth=5,
        learning_rate=0.05,
        random_state=42
    )
    
    # Cross-validation
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')
    
    print(f"\n  Cross-Validation Results:")
    print(f"    AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
    print(f"    Range: {cv_scores.min():.3f} - {cv_scores.max():.3f}")
    
    # Train on full data
    model.fit(X_scaled, y)
    
    # Feature importance
    feature_importance = pd.DataFrame({
        'Feature': feature_cols,
        'Importance': model.feature_importances_
    }).sort_values('Importance', ascending=False)
    
    print(f"\n  Feature Importance:")
    for idx, row in feature_importance.iterrows():
        print(f"    {row['Feature']}: {row['Importance']:.4f}")
    
    return model, scaler, cv_scores, feature_importance

def shap_interpretation(model, scaler, X, feature_cols):
    """SHAP interpretation of model predictions"""
    if not SHAP_AVAILABLE:
        print(f"\n[WARNING] SHAP not available for interpretation")
        return None
    
    print(f"\n" + "=" * 80)
    print("SHAP INTERPRETATION")
    print("=" * 80)
    
    try:
        # Scale data
        X_scaled = scaler.transform(X)
        
        # Create SHAP explainer
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X_scaled)
        
        # Calculate mean absolute SHAP values
        if isinstance(shap_values, list):
            shap_values = shap_values[1]  # For binary classification
        
        mean_shap = np.abs(shap_values).mean(axis=0)
        
        shap_df = pd.DataFrame({
            'Feature': feature_cols,
            'Mean_|SHAP|': mean_shap
        }).sort_values('Mean_|SHAP|', ascending=False)
        
        print(f"\n  Mean Absolute SHAP Values:")
        for idx, row in shap_df.iterrows():
            pct = row['Mean_|SHAP|'] / mean_shap.sum() * 100
            print(f"    {row['Feature']}: {row['Mean_|SHAP|']:.4f} ({pct:.1f}%)")
        
        return shap_values, shap_df
    except Exception as e:
        print(f"[ERROR] SHAP interpretation failed: {e}")
        return None

def risk_stratification(model, scaler, features, feature_cols):
    """Stratify patients into risk categories"""
    print(f"\n" + "=" * 80)
    print("RISK STRATIFICATION")
    print("=" * 80)
    
    X = features[feature_cols].values
    X = np.nan_to_num(X, nan=0.0)
    X_scaled = scaler.transform(X)
    
    # Get predictions
    predictions = model.predict_proba(X_scaled)[:, 1]
    
    # Risk categories
    low_threshold = 0.33
    high_threshold = 0.67
    
    risk_categories = pd.cut(
        predictions,
        bins=[0, low_threshold, high_threshold, 1.0],
        labels=['Low Risk', 'Intermediate Risk', 'High Risk']
    )
    
    stratification = pd.DataFrame({
        'Sample': features['Sample'].values if 'Sample' in features.columns else range(len(predictions)),
        'Rupture_Probability': predictions,
        'Risk_Category': risk_categories
    })
    
    print(f"\n  Risk Stratification:")
    print(stratification['Risk_Category'].value_counts().to_string())
    
    return stratification

def main():
    """Main multi-omics integration pipeline"""
    
    print("=" * 80)
    print("MULTI-OMICS INTEGRATION")
    print("=" * 80)
    
    os.makedirs("results", exist_ok=True)
    
    # Load data
    print(f"\nLoading data components...")
    
    # Load metadata (from GSE122897) - prefer clinical-augmented file if available
    metadata_dir = "data/sc_rnaseq/GSE122897"
    metadata_file_clin = os.path.join(metadata_dir, "sample_metadata_with_clinical.tsv")
    metadata_file_base = os.path.join(metadata_dir, "sample_metadata.tsv")
    if os.path.exists(metadata_file_clin):
        metadata = pd.read_csv(metadata_file_clin, sep='\t')
        print(f"  [OK] Clinical metadata: {metadata.shape}")
    elif os.path.exists(metadata_file_base):
        metadata = pd.read_csv(metadata_file_base, sep='\t')
        print(f"  [OK] Metadata: {metadata.shape}")
    else:
        print(f"  [ERROR] Metadata not found")
        return
    
    # Load expression
    expr_df = load_tissue_expression()
    if expr_df is not None:
        print(f"  [OK] Expression data: {expr_df.shape}")
    
    # Load blood signature
    blood_sig = load_blood_signature()
    if blood_sig is not None:
        print(f"  [OK] Blood signature: {len(blood_sig)} genes")
    
    # Load EC dysfunction scores
    ec_scores = load_ec_dysfunction_scores()
    if ec_scores is not None:
        print(f"  [OK] EC dysfunction scores: {ec_scores.shape}")
    
    # Prepare integrated features
    features, feature_cols = prepare_integrated_features(
        metadata, expr_df, blood_sig, ec_scores
    )
    
    # Create labels (1 = Ruptured, 0 = Unruptured)
    if 'Group' in features.columns:
        labels = (features['Group'] == 'Ruptured').astype(int)
    else:
        print(f"[ERROR] No 'Group' column for labels")
        return
    
    # Build model
    model, scaler, cv_scores, feature_importance = build_integrated_model(
        features, feature_cols, labels
    )
    
    # SHAP interpretation
    X = features[feature_cols].values
    X = np.nan_to_num(X, nan=0.0)
    shap_results = shap_interpretation(model, scaler, X, feature_cols)
    
    # Risk stratification
    stratification = risk_stratification(model, scaler, features, feature_cols)
    
    # Save results
    features_file = "results/integrated_features.tsv"
    features.to_csv(features_file, sep='\t', index=False)
    print(f"\n[OK] Saved integrated features to {features_file}")
    
    stratification_file = "results/risk_stratification.tsv"
    stratification.to_csv(stratification_file, sep='\t', index=False)
    print(f"[OK] Saved risk stratification to {stratification_file}")
    
    importance_file = "results/feature_importance.tsv"
    feature_importance.to_csv(importance_file, sep='\t', index=False)
    print(f"[OK] Saved feature importance to {importance_file}")
    
    # Summary
    print(f"\n\n" + "=" * 80)
    print("INTEGRATION SUMMARY")
    print("=" * 80)
    print(f"\n  Cross-validation AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
    print(f"\n  Feature Importance:")
    for idx, row in feature_importance.iterrows():
        print(f"    {row['Feature']}: {row['Importance']:.4f}")
    
    print(f"\n  Risk Categories:")
    print(stratification['Risk_Category'].value_counts().to_string())
    
    print("\n" + "=" * 80)
    print("MULTI-OMICS INTEGRATION COMPLETE")
    print("=" * 80)
    
    return model, scaler, features, stratification

if __name__ == "__main__":
    main()

