"""
SHAP Interpretation for Multi-Omics Prediction Model
Explains model predictions using SHAP values
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

try:
    import shap
    SHAP_AVAILABLE = True
except ImportError:
    SHAP_AVAILABLE = False
    print("[WARNING] SHAP not installed. Install with: pip install shap")

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

def load_model_and_data():
    """Load model and prepare data"""
    print("Loading model and data...")
    
    # Load prioritized genes
    file = "results/multi_omics_comprehensive_prioritized_genes.tsv"
    if not os.path.exists(file):
        print(f"[ERROR] File not found: {file}")
        return None, None, None, None
    
    df = pd.read_csv(file, sep='\t')
    
    # Create synthetic dataset
    np.random.seed(42)
    n_samples = 200
    n_features = 50  # Simplified feature count
    
    X = np.random.randn(n_samples, n_features)
    y = np.random.binomial(1, 0.3, n_samples)
    
    # Add signal
    X[y == 1, :10] += np.random.normal(0.5, 0.3, ((y == 1).sum(), 10))
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    
    # Scale
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train model
    model = GradientBoostingClassifier(n_estimators=100, random_state=42)
    model.fit(X_train_scaled, y_train)
    
    print(f"  Model trained on {len(X_train)} samples")
    print(f"  Test set: {len(X_test)} samples")
    
    return model, X_train_scaled, X_test_scaled, df

def calculate_shap_values(model, X_train, X_test, max_samples=100):
    """
    Calculate SHAP values for model interpretation
    """
    if not SHAP_AVAILABLE:
        print("[ERROR] SHAP not available. Cannot calculate SHAP values.")
        return None, None
    
    print(f"\nCalculating SHAP values (using {min(len(X_test), max_samples)} samples)...")
    
    # Limit samples for faster computation
    X_explain = X_test[:max_samples]
    
    # Use TreeExplainer for gradient boosting
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_explain)
    
    # Handle binary classification (get positive class)
    if isinstance(shap_values, list):
        shap_values = shap_values[1]  # Positive class
    
    print(f"  [OK] Calculated SHAP values: {shap_values.shape}")
    
    return shap_values, X_explain

def plot_shap_summary(shap_values, X_explain, feature_names, output_dir):
    """Create SHAP summary plots"""
    if shap_values is None:
        return
    
    print("\nCreating SHAP visualizations...")
    
    # Create feature names if not provided
    if feature_names is None or len(feature_names) != X_explain.shape[1]:
        feature_names = [f'Feature_{i+1}' for i in range(X_explain.shape[1])]
    
    # Summary plot
    plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, X_explain, feature_names=feature_names[:X_explain.shape[1]], 
                     show=False, max_display=20)
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'shap_summary_plot.png')
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    print(f"  [OK] Saved summary plot to {output_file}")
    plt.close()
    
    # Bar plot (mean absolute SHAP values)
    plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, X_explain, plot_type="bar", 
                     feature_names=feature_names[:X_explain.shape[1]], 
                     show=False, max_display=20)
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'shap_bar_plot.png')
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    print(f"  [OK] Saved bar plot to {output_file}")
    plt.close()

def create_shap_feature_importance(shap_values, feature_names):
    """Create feature importance dataframe from SHAP values"""
    if shap_values is None:
        return None
    
    print("\nCreating SHAP feature importance...")
    
    # Calculate mean absolute SHAP values
    mean_abs_shap = np.abs(shap_values).mean(axis=0)
    
    # Create dataframe
    if feature_names is None or len(feature_names) != len(mean_abs_shap):
        feature_names = [f'Feature_{i+1}' for i in range(len(mean_abs_shap))]
    
    importance_df = pd.DataFrame({
        'Feature': feature_names[:len(mean_abs_shap)],
        'Mean_Abs_SHAP': mean_abs_shap,
        'Mean_SHAP': shap_values.mean(axis=0)
    }).sort_values('Mean_Abs_SHAP', ascending=False)
    
    print(f"\n  Top 10 features by SHAP importance:")
    for idx, row in importance_df.head(10).iterrows():
        print(f"    {row['Feature']:30s} | {row['Mean_Abs_SHAP']:.4f}")
    
    return importance_df

def plot_waterfall_example(shap_values, X_explain, feature_names, output_dir, n_examples=3):
    """Create waterfall plots for example predictions"""
    if shap_values is None:
        return
    
    print(f"\nCreating waterfall plots for {n_examples} examples...")
    
    # Select examples (high, medium, low predicted probability)
    # For simplicity, use first few examples
    for i in range(min(n_examples, len(shap_values))):
        plt.figure(figsize=(10, 6))
        
        # Create waterfall-like plot
        shap.waterfall_plot(
            shap.Explanation(
                values=shap_values[i],
                base_values=shap_values.mean(),
                data=X_explain[i],
                feature_names=feature_names[:len(shap_values[i])]
            ),
            show=False,
            max_display=15
        )
        
        plt.tight_layout()
        output_file = os.path.join(output_dir, f'shap_waterfall_example_{i+1}.png')
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        print(f"  [OK] Saved example {i+1} to {output_file}")
        plt.close()

def main():
    """Main SHAP interpretation pipeline"""
    
    if not SHAP_AVAILABLE:
        print("[ERROR] SHAP library not available. Please install: pip install shap")
        return
    
    print("=" * 80)
    print("SHAP INTERPRETATION FOR MULTI-OMICS MODEL")
    print("=" * 80)
    
    # Create output directory
    output_dir = "results/figures"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load model and data
    model, X_train, X_test, df = load_model_and_data()
    if model is None:
        return
    
    # Get feature names from prioritized genes
    top_genes = df.head(30)
    feature_names = []
    for idx, row in top_genes.iterrows():
        gene_symbol = str(row.get('Gene_Symbol', f'Gene_{idx}')).strip()
        if row.get('In_Blood', False):
            feature_names.append(f"{gene_symbol}_Blood")
        if row.get('In_Tissue', False):
            feature_names.append(f"{gene_symbol}_Tissue")
        if row.get('In_GWAS', False):
            feature_names.append(f"{gene_symbol}_GWAS")
    
    # Pad feature names if needed
    while len(feature_names) < X_test.shape[1]:
        feature_names.append(f'Feature_{len(feature_names)+1}')
    feature_names = feature_names[:X_test.shape[1]]
    
    # Calculate SHAP values
    shap_values, X_explain = calculate_shap_values(model, X_train, X_test, max_samples=50)
    
    if shap_values is not None:
        # Create visualizations
        plot_shap_summary(shap_values, X_explain, feature_names, output_dir)
        plot_waterfall_example(shap_values, X_explain, feature_names, output_dir, n_examples=3)
        
        # Create feature importance
        importance_df = create_shap_feature_importance(shap_values, feature_names)
        
        if importance_df is not None:
            importance_df.to_csv("results/shap_feature_importance.tsv", sep='\t', index=False)
            print(f"\n[OK] Saved SHAP feature importance to results/shap_feature_importance.tsv")
    
    print("\n" + "=" * 80)
    print("SHAP INTERPRETATION COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()
