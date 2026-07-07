"""
Risk Stratification for Multi-Omics Model
Creates clinically actionable risk categories
"""
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import seaborn as sns

def load_data():
    """Load features and labels"""
    print("\n" + "=" * 80)
    print("LOADING DATA")
    print("=" * 80)
    
    features_file = "results/integrated_features_full.tsv"
    if not os.path.exists(features_file):
        print(f"[ERROR] Features file not found: {features_file}")
        return None, None
    
    data = pd.read_csv(features_file, sep='\t', index_col=0)
    
    if 'Ruptured' not in data.columns:
        print("[ERROR] No 'Ruptured' column")
        return None, None
    
    features = data.drop('Ruptured', axis=1)
    labels = data['Ruptured']
    
    print(f"\n[OK] Loaded {len(features)} samples")
    print(f"  Features: {list(features.columns)}")
    print(f"  Ruptured: {labels.sum()}, Unruptured: {(labels==0).sum()}")
    
    return features, labels

def train_model_and_predict(features, labels):
    """Train model and get predictions"""
    print("\n" + "=" * 80)
    print("TRAINING MODEL AND GENERATING PREDICTIONS")
    print("=" * 80)
    
    X = features.values
    y = labels.values
    X = np.nan_to_num(X, nan=0.0)
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Use best model (XGBoost)
    model = XGBClassifier(
        n_estimators=200, max_depth=5, learning_rate=0.05,
        random_state=42, eval_metric='logloss', use_label_encoder=False
    )
    
    # Get cross-validated predictions
    print("\nGenerating cross-validated predictions...")
    from sklearn.model_selection import StratifiedKFold
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    y_pred_proba = cross_val_predict(model, X_scaled, y, cv=cv, method='predict_proba')[:, 1]
    
    # Also train on full data for final model
    model.fit(X_scaled, y)
    y_pred_proba_full = model.predict_proba(X_scaled)[:, 1]
    
    # Calculate AUC
    auc_cv = roc_auc_score(y, y_pred_proba)
    auc_full = roc_auc_score(y, y_pred_proba_full)
    
    print(f"\n[OK] Model performance:")
    print(f"  CV AUC: {auc_cv:.3f}")
    print(f"  Full model AUC: {auc_full:.3f}")
    
    return model, scaler, y_pred_proba, y_pred_proba_full

def stratify_risk(predictions, thresholds=[0.33, 0.67]):
    """
    Stratify patients into risk categories
    
    Parameters:
    - thresholds: [low_threshold, high_threshold]
      Low risk: < thresholds[0]
      Intermediate risk: thresholds[0] <= x < thresholds[1]
      High risk: >= thresholds[1]
    """
    print("\n" + "=" * 80)
    print("RISK STRATIFICATION")
    print("=" * 80)
    
    low_threshold = thresholds[0]
    high_threshold = thresholds[1]
    
    risk_categories = pd.cut(
        predictions,
        bins=[0, low_threshold, high_threshold, 1.0],
        labels=['Low Risk', 'Intermediate Risk', 'High Risk'],
        include_lowest=True
    )
    
    stratification = pd.DataFrame({
        'Sample': range(len(predictions)),
        'Rupture_Probability': predictions,
        'Risk_Category': risk_categories
    })
    
    print(f"\nRisk Stratification:")
    print(stratification['Risk_Category'].value_counts().sort_index())
    
    print(f"\nRisk Thresholds:")
    print(f"  Low Risk: < {low_threshold:.2f}")
    print(f"  Intermediate Risk: {low_threshold:.2f} - {high_threshold:.2f}")
    print(f"  High Risk: >= {high_threshold:.2f}")
    
    return stratification

def calculate_risk_statistics(stratification, labels):
    """Calculate statistics for each risk category"""
    print("\n" + "=" * 80)
    print("RISK CATEGORY STATISTICS")
    print("=" * 80)
    
    stats = []
    for category in ['Low Risk', 'Intermediate Risk', 'High Risk']:
        mask = stratification['Risk_Category'] == category
        n_samples = mask.sum()
        n_ruptured = labels[mask].sum()
        rupture_rate = n_ruptured / n_samples if n_samples > 0 else 0
        mean_prob = stratification.loc[mask, 'Rupture_Probability'].mean()
        
        stats.append({
            'Risk_Category': category,
            'N_Samples': n_samples,
            'N_Ruptured': n_ruptured,
            'Rupture_Rate': rupture_rate,
            'Mean_Probability': mean_prob
        })
        
        print(f"\n{category}:")
        print(f"  Samples: {n_samples}")
        print(f"  Ruptured: {n_ruptured} ({rupture_rate*100:.1f}%)")
        print(f"  Mean probability: {mean_prob:.3f}")
    
    stats_df = pd.DataFrame(stats)
    return stats_df

def plot_risk_distribution(stratification, labels):
    """Plot risk distribution"""
    print("\n" + "=" * 80)
    print("GENERATING RISK DISTRIBUTION PLOTS")
    print("=" * 80)
    
    os.makedirs("results/figures", exist_ok=True)
    
    # Add labels to stratification
    stratification['Ruptured'] = labels.values
    
    # 1. Histogram of risk probabilities
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # By rupture status
    axes[0].hist(
        stratification[stratification['Ruptured']==0]['Rupture_Probability'],
        bins=20, alpha=0.7, label='Unruptured', color='blue'
    )
    axes[0].hist(
        stratification[stratification['Ruptured']==1]['Rupture_Probability'],
        bins=20, alpha=0.7, label='Ruptured', color='red'
    )
    axes[0].axvline(0.33, color='green', linestyle='--', label='Low/Int Threshold')
    axes[0].axvline(0.67, color='orange', linestyle='--', label='Int/High Threshold')
    axes[0].set_xlabel('Rupture Probability', fontsize=12)
    axes[0].set_ylabel('Frequency', fontsize=12)
    axes[0].set_title('Distribution of Rupture Probabilities', fontsize=14, fontweight='bold')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # By risk category
    risk_order = ['Low Risk', 'Intermediate Risk', 'High Risk']
    category_counts = stratification['Risk_Category'].value_counts()
    colors = ['green', 'orange', 'red']
    axes[1].bar(range(len(risk_order)), [category_counts.get(cat, 0) for cat in risk_order], color=colors)
    axes[1].set_xticks(range(len(risk_order)))
    axes[1].set_xticklabels(risk_order, rotation=45, ha='right')
    axes[1].set_ylabel('Number of Samples', fontsize=12)
    axes[1].set_title('Sample Distribution by Risk Category', fontsize=14, fontweight='bold')
    axes[1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig("results/figures/risk_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("\n[OK] Saved: results/figures/risk_distribution.png")
    
    # 2. Box plot by risk category
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=stratification, x='Risk_Category', y='Rupture_Probability', 
                order=risk_order, palette=['green', 'orange', 'red'])
    plt.xlabel('Risk Category', fontsize=12)
    plt.ylabel('Rupture Probability', fontsize=12)
    plt.title('Rupture Probability by Risk Category', fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig("results/figures/risk_category_boxplot.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("[OK] Saved: results/figures/risk_category_boxplot.png")

def clinical_interpretation(stratification, labels):
    """Provide clinical interpretation"""
    print("\n" + "=" * 80)
    print("CLINICAL INTERPRETATION")
    print("=" * 80)
    
    print("\nClinical Recommendations:")
    print("\n  LOW RISK (< 0.33):")
    print("    - Conservative management")
    print("    - Annual imaging follow-up")
    print("    - Monitor for changes")
    
    print("\n  INTERMEDIATE RISK (0.33 - 0.67):")
    print("    - Close monitoring")
    print("    - 6-month imaging follow-up")
    print("    - Consider intervention if other risk factors present")
    
    print("\n  HIGH RISK (>= 0.67):")
    print("    - Strongly consider intervention")
    print("    - 3-month imaging follow-up if not treating")
    print("    - Discuss treatment options with patient")
    
    # Identify high-risk patients
    high_risk = stratification[stratification['Risk_Category'] == 'High Risk']
    print(f"\n  High-risk patients identified: {len(high_risk)}")
    if len(high_risk) > 0:
        high_risk_ruptured = labels[high_risk.index].sum()
        print(f"    - {high_risk_ruptured} actually ruptured ({high_risk_ruptured/len(high_risk)*100:.1f}%)")

def main():
    """Main risk stratification pipeline"""
    print("=" * 80)
    print("RISK STRATIFICATION FOR MULTI-OMICS MODEL")
    print("=" * 80)
    
    os.makedirs("results/figures", exist_ok=True)
    
    # Load data
    features, labels = load_data()
    if features is None or labels is None:
        print("\n[ERROR] Could not load data. Exiting.")
        return
    
    # Train model and get predictions
    model, scaler, y_pred_cv, y_pred_full = train_model_and_predict(features, labels)
    
    # Stratify risk
    stratification = stratify_risk(y_pred_full, thresholds=[0.33, 0.67])
    
    # Calculate statistics
    stats_df = calculate_risk_statistics(stratification, labels)
    
    # Add sample IDs if available
    if hasattr(features, 'index'):
        stratification['Sample_ID'] = features.index.values
    else:
        stratification['Sample_ID'] = range(len(stratification))
    
    # Add actual rupture status
    stratification['Actual_Ruptured'] = labels.values
    
    # Save results
    stratification_file = "results/risk_stratification_full.tsv"
    stratification.to_csv(stratification_file, sep='\t', index=False)
    print(f"\n[OK] Saved risk stratification to {stratification_file}")
    
    stats_file = "results/risk_category_statistics.tsv"
    stats_df.to_csv(stats_file, sep='\t', index=False)
    print(f"[OK] Saved risk statistics to {stats_file}")
    
    # Generate plots
    plot_risk_distribution(stratification, labels)
    
    # Clinical interpretation
    clinical_interpretation(stratification, labels)
    
    print("\n" + "=" * 80)
    print("RISK STRATIFICATION COMPLETE")
    print("=" * 80)
    
    return stratification, stats_df

if __name__ == "__main__":
    main()


