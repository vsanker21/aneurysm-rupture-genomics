"""
Machine Learning Classification using Top Blood DE Genes
Builds classification models to distinguish aneurysm vs control samples
"""
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve, classification_report, confusion_matrix
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
import warnings
import os
warnings.filterwarnings('ignore')

def load_expression_data():
    """Load blood RNA-seq expression data"""
    expr_file = "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
    metadata_file = "data/blood_rnaseq/GSE159610/sample_metadata_with_labels.tsv"
    
    if not os.path.exists(expr_file) or not os.path.exists(metadata_file):
        print("[ERROR] Required data files not found")
        return None, None
    
    print("Loading expression data...")
    expr_df = pd.read_csv(expr_file, sep='\t', compression='gzip', index_col=0)
    print(f"  Expression matrix: {expr_df.shape[0]} genes, {expr_df.shape[1]} samples")
    
    metadata = pd.read_csv(metadata_file, sep='\t')
    print(f"  Metadata: {metadata.shape}")
    
    return expr_df, metadata

def select_feature_genes(expr_df, metadata, n_genes=50):
    """Select top DE genes as features"""
    print(f"\nSelecting top {n_genes} genes as features...")
    
    # Get sample columns
    data_cols = sorted([c for c in expr_df.columns if c.startswith('DSA')], 
                       key=lambda x: int(x.replace('DSA', '')))
    
    # Identify samples
    aneurysm_samples = [data_cols[i] for i in range(len(metadata)) 
                        if metadata.iloc[i]['Status'] == 'Aneurysm']
    control_samples = [data_cols[i] for i in range(len(metadata)) 
                       if metadata.iloc[i]['Status'] == 'Control']
    
    # Load DE results to get top genes
    de_file = "results/gse159610_top_candidate_genes.tsv"
    if os.path.exists(de_file):
        de_df = pd.read_csv(de_file, sep='\t')
        # Take top N genes
        top_genes = de_df.head(n_genes)['Gene_ID'].tolist()
        print(f"  Selected {len(top_genes)} genes from DE results")
    else:
        # Fallback: use genes with highest variance
        print("  [WARNING] DE results not found, using high-variance genes")
        gene_variance = expr_df.var(axis=1)
        top_genes = gene_variance.nlargest(n_genes).index.tolist()
    
    return top_genes, aneurysm_samples, control_samples

def prepare_features(expr_df, feature_genes, samples):
    """Prepare feature matrix"""
    # Subset to feature genes and samples
    feature_df = expr_df.loc[feature_genes, samples].T
    
    # Log2 transform
    feature_df = np.log2(feature_df + 1)
    
    return feature_df

def build_classification_models(X, y, feature_names):
    """Build and evaluate multiple classification models"""
    print("\n" + "=" * 80)
    print("BUILDING CLASSIFICATION MODELS")
    print("=" * 80)
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    models = {
        'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42, max_depth=10),
        'SVM': SVC(kernel='rbf', probability=True, random_state=42),
        'Logistic Regression': LogisticRegression(max_iter=1000, random_state=42)
    }
    
    results = {}
    
    for name, model in models.items():
        print(f"\n{name}:")
        
        # Train
        if name == 'SVM' or name == 'Logistic Regression':
            model.fit(X_train_scaled, y_train)
            y_pred = model.predict(X_test_scaled)
            y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        else:
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            y_pred_proba = model.predict_proba(X_test)[:, 1]
        
        # Evaluate
        auc = roc_auc_score(y_test, y_pred_proba)
        results[name] = {
            'AUC': auc,
            'Model': model,
            'Predictions': y_pred,
            'Probabilities': y_pred_proba
        }
        
        print(f"  Test AUC: {auc:.3f}")
        
        # Classification report
        print(f"  Classification Report:")
        report = classification_report(y_test, y_pred, target_names=['Control', 'Aneurysm'])
        print(report)
        
        # Feature importance (for Random Forest)
        if name == 'Random Forest':
            feature_importance = pd.DataFrame({
                'Gene_ID': feature_names,
                'Importance': model.feature_importances_
            }).sort_values('Importance', ascending=False)
            
            print(f"\n  Top 10 Most Important Features:")
            print(feature_importance.head(10).to_string(index=False))
            results[name]['Feature_Importance'] = feature_importance
    
    return results, X_test, y_test

def cross_validate_models(X, y, feature_names):
    """Perform cross-validation"""
    print("\n" + "=" * 80)
    print("CROSS-VALIDATION")
    print("=" * 80)
    
    # Scale features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    models = {
        'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42, max_depth=10),
        'SVM': SVC(kernel='rbf', probability=True, random_state=42),
        'Logistic Regression': LogisticRegression(max_iter=1000, random_state=42)
    }
    
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_results = {}
    
    for name, model in models.items():
        print(f"\n{name}:")
        
        if name == 'SVM' or name == 'Logistic Regression':
            scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')
        else:
            scores = cross_val_score(model, X, y, cv=cv, scoring='roc_auc')
        
        cv_results[name] = {
            'Mean_AUC': scores.mean(),
            'Std_AUC': scores.std(),
            'Scores': scores
        }
        
        print(f"  CV AUC: {scores.mean():.3f} (+/- {scores.std() * 2:.3f})")
        print(f"  Individual fold scores: {scores}")
    
    return cv_results

def save_results(results, cv_results, feature_names):
    """Save classification results"""
    os.makedirs("results", exist_ok=True)
    
    # Save feature importance
    if 'Random Forest' in results and 'Feature_Importance' in results['Random Forest']:
        importance_file = "results/blood_ml_feature_importance.tsv"
        results['Random Forest']['Feature_Importance'].to_csv(importance_file, sep='\t', index=False)
        print(f"\n[OK] Saved feature importance to {importance_file}")
    
    # Save CV results
    cv_df = pd.DataFrame({
        'Model': list(cv_results.keys()),
        'Mean_AUC': [cv_results[m]['Mean_AUC'] for m in cv_results.keys()],
        'Std_AUC': [cv_results[m]['Std_AUC'] for m in cv_results.keys()]
    })
    cv_file = "results/blood_ml_cv_results.tsv"
    cv_df.to_csv(cv_file, sep='\t', index=False)
    print(f"[OK] Saved CV results to {cv_file}")
    
    # Save summary
    summary = {
        'Best_Model': max(cv_results, key=lambda x: cv_results[x]['Mean_AUC']),
        'Best_AUC': max([cv_results[m]['Mean_AUC'] for m in cv_results.keys()]),
        'N_Features': len(feature_names),
        'N_Samples': len(feature_names)  # This will be updated with actual sample count
    }
    
    summary_file = "results/blood_ml_summary.tsv"
    pd.DataFrame([summary]).to_csv(summary_file, sep='\t', index=False)
    print(f"[OK] Saved summary to {summary_file}")

def main():
    """Main ML classification pipeline"""
    print("=" * 80)
    print("MACHINE LEARNING CLASSIFICATION")
    print("=" * 80)
    
    # Load data
    expr_df, metadata = load_expression_data()
    if expr_df is None:
        return
    
    # Select feature genes
    feature_genes, aneurysm_samples, control_samples = select_feature_genes(
        expr_df, metadata, n_genes=50
    )
    
    # Prepare features
    all_samples = aneurysm_samples + control_samples
    X = prepare_features(expr_df, feature_genes, all_samples)
    
    # Prepare labels (0 = Control, 1 = Aneurysm)
    y = np.array([0] * len(control_samples) + [1] * len(aneurysm_samples))
    
    print(f"\nFeature matrix shape: {X.shape}")
    print(f"  Samples: {X.shape[0]}")
    print(f"  Features: {X.shape[1]}")
    print(f"  Classes: Control={np.sum(y==0)}, Aneurysm={np.sum(y==1)}")
    
    # Build models
    results, X_test, y_test = build_classification_models(X.values, y, feature_genes)
    
    # Cross-validation
    cv_results = cross_validate_models(X.values, y, feature_genes)
    
    # Save results
    save_results(results, cv_results, feature_genes)
    
    print("\n" + "=" * 80)
    print("MACHINE LEARNING ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"\nBest model: {max(cv_results, key=lambda x: cv_results[x]['Mean_AUC'])}")
    print(f"Best CV AUC: {max([cv_results[m]['Mean_AUC'] for m in cv_results.keys()]):.3f}")

if __name__ == "__main__":
    main()






