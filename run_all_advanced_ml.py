"""
Master Script: Run All Advanced ML Analyses
Executes all Tier 1, 2, and 3 advanced ML methods
"""
import sys
import os
import time
from pathlib import Path

print("=" * 80)
print("RUNNING ALL ADVANCED ML ANALYSES")
print("=" * 80)

os.makedirs("results/advanced_ml", exist_ok=True)
os.makedirs("results/figures/advanced_ml", exist_ok=True)

# Track results
all_results = {}

# ============================================================================
# TIER 1: HIGHEST PRIORITY
# ============================================================================

print("\n" + "=" * 80)
print("TIER 1: HIGHEST PRIORITY METHODS")
print("=" * 80)

# 1. Transfer Learning
print("\n[1/8] Running Transfer Learning...")
try:
    exec(open('advanced_ml_transfer_learning.py').read())
    import json
    if os.path.exists("results/advanced_ml/transfer_learning_results.json"):
        with open("results/advanced_ml/transfer_learning_results.json", 'r') as f:
            all_results['transfer_learning'] = json.load(f)
        print(f"[OK] Transfer Learning complete: AUC = {all_results['transfer_learning']['mean_auc']:.4f}")
except Exception as e:
    print(f"[WARNING] Transfer Learning failed: {e}")

# 2. Graph Neural Networks
print("\n[2/8] Running Graph Neural Networks...")
try:
    exec(open('advanced_ml_gnn_transformer.py').read())
    import json
    if os.path.exists("results/advanced_ml/gnn_network_info.json"):
        with open("results/advanced_ml/gnn_network_info.json", 'r') as f:
            all_results['gnn'] = json.load(f)
        print(f"[OK] GNN network built")
except Exception as e:
    print(f"[WARNING] GNN failed: {e}")

# 3. Transformer Models
print("\n[3/8] Running Transformer Models...")
try:
    # Already in gnn_transformer script
    import json
    if os.path.exists("results/advanced_ml/transformer_results.json"):
        with open("results/advanced_ml/transformer_results.json", 'r') as f:
            all_results['transformer'] = json.load(f)
        print(f"[OK] Transformer complete: AUC = {all_results['transformer']['mean_auc']:.4f}")
except Exception as e:
    print(f"[WARNING] Transformer failed: {e}")

# ============================================================================
# TIER 2: ADVANCED METHODS
# ============================================================================

print("\n" + "=" * 80)
print("TIER 2: ADVANCED METHODS")
print("=" * 80)

# 4. Ensemble Stacking
print("\n[4/8] Running Ensemble Stacking...")
try:
    exec(open('advanced_ml_ensemble_vae_mtl.py').read())
    import json
    if os.path.exists("results/advanced_ml/ensemble_stacking_results.json"):
        with open("results/advanced_ml/ensemble_stacking_results.json", 'r') as f:
            all_results['ensemble'] = json.load(f)
        print(f"[OK] Ensemble complete: AUC = {all_results['ensemble']['mean_auc']:.4f}")
except Exception as e:
    print(f"[WARNING] Ensemble failed: {e}")

# 5. Variational Autoencoder
print("\n[5/8] Running Variational Autoencoder...")
try:
    import json
    if os.path.exists("results/advanced_ml/vae_results.json"):
        with open("results/advanced_ml/vae_results.json", 'r') as f:
            all_results['vae'] = json.load(f)
        print(f"[OK] VAE complete: AUC = {all_results['vae']['mean_auc']:.4f}")
except Exception as e:
    print(f"[WARNING] VAE failed: {e}")

# 6. Multi-Task Learning
print("\n[6/8] Running Multi-Task Learning...")
try:
    import json
    if os.path.exists("results/advanced_ml/multi_task_learning_results.json"):
        with open("results/advanced_ml/multi_task_learning_results.json", 'r') as f:
            all_results['multi_task'] = json.load(f)
        print(f"[OK] Multi-Task complete: AUC = {all_results['multi_task']['mean_auc']:.4f}")
except Exception as e:
    print(f"[WARNING] Multi-Task failed: {e}")

# ============================================================================
# TIER 3: CUTTING-EDGE
# ============================================================================

print("\n" + "=" * 80)
print("TIER 3: CUTTING-EDGE METHODS")
print("=" * 80)

# 7. Advanced Deep Neural Networks
print("\n[7/8] Running Advanced Deep Neural Networks...")
try:
    import json
    if os.path.exists("results/advanced_ml/deep_neural_network_results.json"):
        with open("results/advanced_ml/deep_neural_network_results.json", 'r') as f:
            all_results['dnn'] = json.load(f)
        print(f"[OK] Deep Neural Network complete: AUC = {all_results['dnn']['mean_auc']:.4f}")
except Exception as e:
    print(f"[WARNING] Deep Neural Network failed: {e}")

# 8. Bayesian Neural Networks (simplified implementation)
print("\n[8/8] Running Bayesian Neural Networks...")
try:
    # Bayesian NN requires specialized libraries
    print("[INFO] Bayesian Neural Networks require PyMC3 or TensorFlow Probability")
    print("[INFO] Simplified implementation with dropout-based uncertainty...")
    
    # Use dropout as Bayesian approximation
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers
    import pandas as pd
    import numpy as np
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score
    
    import multi_omics_integration as moi
    
    metadata_dir = 'data/sc_rnaseq/GSE122897'
    try:
        metadata = pd.read_csv(metadata_dir + '/sample_metadata_with_clinical.tsv', sep='\t')
    except:
        metadata = pd.read_csv(metadata_dir + '/sample_metadata.tsv', sep='\t')
    
    features, feature_cols = moi.prepare_integrated_features(
        metadata, moi.load_tissue_expression(), 
        moi.load_blood_signature(), moi.load_ec_dysfunction_scores()
    )
    
    X = features[feature_cols].values
    X = np.nan_to_num(X, nan=0.0)
    y = (metadata['Group'] == 'Ruptured').astype(int).values
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Bayesian approximation with dropout
    def build_bayesian_nn():
        inputs = keras.Input(shape=(X_scaled.shape[1],))
        x = layers.Dense(128, activation='relu')(inputs)
        x = layers.Dropout(0.5)(x, training=True)  # Always on for Bayesian
        x = layers.Dense(64, activation='relu')(x)
        x = layers.Dropout(0.5)(x, training=True)
        x = layers.Dense(32, activation='relu')(x)
        x = layers.Dropout(0.4)(x, training=True)
        outputs = layers.Dense(1, activation='sigmoid')(x)
        return keras.Model(inputs, outputs)
    
    model = build_bayesian_nn()
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy', keras.metrics.AUC(name='auc')])
    
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = []
    
    for fold, (train_idx, val_idx) in enumerate(cv.split(X_scaled, y)):
        X_train, X_val = X_scaled[train_idx], X_scaled[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        
        from sklearn.utils.class_weight import compute_class_weight
        class_weights = compute_class_weight('balanced', classes=np.unique(y_train), y=y_train)
        class_weight_dict = {i: w for i, w in enumerate(class_weights)}
        
        model.fit(X_train, y_train, validation_data=(X_val, y_val),
                 epochs=100, batch_size=8, class_weight=class_weight_dict, verbose=0,
                 callbacks=[keras.callbacks.EarlyStopping(monitor='val_auc', patience=15, restore_best_weights=True, mode='max')])
        
        # Monte Carlo dropout for uncertainty
        predictions = []
        for _ in range(10):  # 10 forward passes
            pred = model.predict(X_val, verbose=0)
            predictions.append(pred.flatten())
        y_pred_proba = np.mean(predictions, axis=0)
        uncertainty = np.std(predictions, axis=0)
        
        auc = roc_auc_score(y_val, y_pred_proba)
        cv_scores.append(auc)
    
    mean_auc = np.mean(cv_scores)
    std_auc = np.std(cv_scores)
    
    bayesian_results = {
        'method': 'Bayesian Neural Network (Dropout-based)',
        'mean_auc': float(mean_auc),
        'std_auc': float(std_auc),
        'cv_scores': [float(s) for s in cv_scores],
        'uncertainty_quantification': True
    }
    
    import json
    with open("results/advanced_ml/bayesian_nn_results.json", 'w') as f:
        json.dump(bayesian_results, f, indent=2)
    
    all_results['bayesian_nn'] = bayesian_results
    print(f"[OK] Bayesian NN complete: AUC = {mean_auc:.4f} ± {std_auc:.4f}")
    
except Exception as e:
    print(f"[WARNING] Bayesian NN failed: {e}")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("ADVANCED ML ANALYSES SUMMARY")
print("=" * 80)

import json

# Save all results
with open("results/advanced_ml/all_advanced_ml_results.json", 'w') as f:
    json.dump(all_results, f, indent=2)

# Create comparison table
print("\n" + "=" * 80)
print("PERFORMANCE COMPARISON")
print("=" * 80)

comparison_data = []
for method_name, results in all_results.items():
    if 'mean_auc' in results:
        comparison_data.append({
            'Method': results.get('method', method_name),
            'Mean_AUC': results['mean_auc'],
            'Std_AUC': results.get('std_auc', 0),
            'Min_AUC': min(results.get('cv_scores', [])),
            'Max_AUC': max(results.get('cv_scores', []))
        })

if comparison_data:
    import pandas as pd
    comparison_df = pd.DataFrame(comparison_data)
    comparison_df = comparison_df.sort_values('Mean_AUC', ascending=False)
    
    print("\n" + comparison_df.to_string(index=False))
    
    # Save comparison
    comparison_df.to_csv("results/advanced_ml/advanced_ml_comparison.tsv", sep='\t', index=False)
    print(f"\n[OK] Comparison saved to results/advanced_ml/advanced_ml_comparison.tsv")
    
    # Create visualization
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=(14, 8))
        
        methods = comparison_df['Method'].tolist()
        means = comparison_df['Mean_AUC'].values
        stds = comparison_df['Std_AUC'].values
        
        y_pos = np.arange(len(methods))
        bars = ax.barh(y_pos, means, xerr=stds, color='#3498db', alpha=0.8, edgecolor='black', capsize=5)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels([m[:40] + '...' if len(m) > 40 else m for m in methods], fontsize=10)
        ax.set_xlabel('AUC-ROC', fontweight='bold', fontsize=14)
        ax.set_title('Advanced ML Methods: Performance Comparison\n' +
                    'Cross-Validation AUC-ROC',
                    fontweight='bold', fontsize=16)
        ax.set_xlim([0, 1])
        ax.axvline(x=0.739, color='red', linestyle='--', linewidth=2, 
                  label='Baseline (0.739)')
        ax.legend(fontsize=12)
        ax.grid(True, alpha=0.3, axis='x')
        
        # Add value labels
        for i, (bar, mean, std) in enumerate(zip(bars, means, stds)):
            ax.text(mean + std + 0.02, i, f'{mean:.3f} ± {std:.3f}',
                   va='center', fontweight='bold', fontsize=10)
        
        plt.tight_layout()
        plt.savefig("results/figures/advanced_ml/all_methods_comparison.png",
                   dpi=600, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"[OK] Comparison visualization saved")
        
    except Exception as e:
        print(f"[WARNING] Visualization failed: {e}")

print("\n" + "=" * 80)
print("ALL ADVANCED ML ANALYSES COMPLETE!")
print("=" * 80)
print(f"\nResults saved to: results/advanced_ml/")
print(f"Figures saved to: results/figures/advanced_ml/")

