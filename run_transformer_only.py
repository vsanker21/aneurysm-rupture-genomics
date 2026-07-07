"""
Run Transformer Model Only
Focused script to ensure transformer works correctly with real data
"""
import pandas as pd
import numpy as np
import os
import json
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("TRANSFORMER MODEL FOR GENE EXPRESSION")
print("=" * 80)

# Import libraries
try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers, callbacks
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score
    from sklearn.utils.class_weight import compute_class_weight
    print("[OK] All required libraries imported")
except ImportError as e:
    print(f"[ERROR] Missing library: {e}")
    print("Please install: pip install tensorflow scikit-learn")
    exit(1)

os.makedirs("results/advanced_ml", exist_ok=True)
os.makedirs("results/figures/advanced_ml", exist_ok=True)

# ============================================================================
# LOAD DATA (REAL DATA ONLY)
# ============================================================================

print("\n" + "=" * 80)
print("LOADING REAL DATA")
print("=" * 80)

import multi_omics_integration as moi

metadata_dir = 'data/sc_rnaseq/GSE122897'
metadata_file = metadata_dir + '/sample_metadata_with_clinical.tsv'

try:
    metadata = pd.read_csv(metadata_file, sep='\t')
    print(f"[OK] Metadata loaded: {metadata.shape}")
except FileNotFoundError:
    metadata = pd.read_csv(metadata_dir + '/sample_metadata.tsv', sep='\t')
    print(f"[OK] Metadata loaded: {metadata.shape}")

# Load REAL expression data
expr_df = moi.load_tissue_expression()
if expr_df is None:
    print("[ERROR] Expression data not available")
    exit(1)

labels = (metadata['Group'] == 'Ruptured').astype(int).values
print(f"[OK] Labels: {labels.sum()} ruptured, {(labels==0).sum()} unruptured")

# Prepare REAL expression data
expr_df = expr_df.T  # Samples x Genes
expr_df = expr_df.reindex(metadata['Sample'])
expr_df = expr_df.dropna()

aligned_labels = labels[metadata['Sample'].isin(expr_df.index)]

# Select top variable genes (real selection from real data)
gene_vars = expr_df.var(axis=0)
top_genes = gene_vars.nlargest(5000).index.tolist()
expr_subset = expr_df[top_genes].values
expr_subset = np.log2(expr_subset + 1)  # Log2 transform real data

print(f"[OK] Expression data prepared: {expr_subset.shape}")
print(f"  Samples: {expr_subset.shape[0]}")
print(f"  Top {len(top_genes)} variable genes selected")

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(expr_subset)
print(f"[OK] Data standardized")

# ============================================================================
# BUILD TRANSFORMER MODEL
# ============================================================================

print("\n" + "=" * 80)
print("BUILDING TRANSFORMER MODEL")
print("=" * 80)

n_genes = X_scaled.shape[1]
d_model = 128  # Reduced for smaller dataset

print(f"Configuration:")
print(f"  Sequence length (genes): {n_genes}")
print(f"  Model dimension: {d_model}")
print(f"  Attention heads: 4")
print(f"  Transformer blocks: 2")

# Transformer encoder function
def transformer_encoder(inputs, num_heads=4, ff_dim=256, dropout=0.1):
    """Transformer encoder block"""
    # Multi-head self-attention
    attention_output = layers.MultiHeadAttention(
        num_heads=num_heads, key_dim=d_model, dropout=dropout
    )(inputs, inputs)
    attention_output = layers.Dropout(dropout)(attention_output)
    attention_output = layers.LayerNormalization(epsilon=1e-6)(inputs + attention_output)
    
    # Feed-forward network
    ffn_output = layers.Dense(ff_dim, activation='relu')(attention_output)
    ffn_output = layers.Dense(d_model)(ffn_output)
    ffn_output = layers.Dropout(dropout)(ffn_output)
    output = layers.LayerNormalization(epsilon=1e-6)(attention_output + ffn_output)
    
    return output

# Build model
inputs = keras.Input(shape=(n_genes,), name='gene_expression')

# Reshape and project to d_model
x = layers.Reshape((n_genes, 1))(inputs)
x = layers.Dense(d_model, activation=None)(x)
x = layers.Dropout(0.1)(x)

# Positional encoding (learnable)
pos_emb = layers.Embedding(n_genes, d_model)(tf.range(n_genes))
pos_emb = tf.expand_dims(pos_emb, 0)  # (1, n_genes, d_model)
x = x + pos_emb

# Transformer encoder blocks
x = transformer_encoder(x, num_heads=4, ff_dim=256, dropout=0.2)
x = transformer_encoder(x, num_heads=4, ff_dim=256, dropout=0.2)

# Global pooling
x = layers.GlobalAveragePooling1D()(x)

# Classification head
x = layers.Dense(128, activation='relu')(x)
x = layers.Dropout(0.4)(x)
x = layers.Dense(64, activation='relu')(x)
x = layers.Dropout(0.3)(x)
outputs = layers.Dense(1, activation='sigmoid', name='classification')(x)

transformer_model = keras.Model(inputs, outputs, name='transformer_classifier')

transformer_model.compile(
    optimizer=keras.optimizers.Adam(learning_rate=0.0001),
    loss='binary_crossentropy',
    metrics=['accuracy', keras.metrics.AUC(name='auc')]
)

print(f"\n[OK] Transformer model created")
print(f"  Total parameters: {transformer_model.count_params():,}")
print(f"  Input shape: {transformer_model.input_shape}")
print(f"  Output shape: {transformer_model.output_shape}")

# Test forward pass
test_input = X_scaled[:2]
test_output = transformer_model.predict(test_input, verbose=0)
print(f"[OK] Forward pass test successful")
print(f"  Test output shape: {test_output.shape}")
print(f"  Test output range: [{test_output.min():.4f}, {test_output.max():.4f}]")

# ============================================================================
# CROSS-VALIDATION
# ============================================================================

print("\n" + "=" * 80)
print("CROSS-VALIDATION WITH TRANSFORMER")
print("=" * 80)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_scores = []
cv_predictions = []

for fold, (train_idx, val_idx) in enumerate(cv.split(X_scaled, aligned_labels)):
    print(f"\n  Fold {fold + 1}/5...")
    
    X_train, X_val = X_scaled[train_idx], X_scaled[val_idx]
    y_train, y_val = aligned_labels[train_idx], aligned_labels[val_idx]
    
    print(f"    Train: {len(X_train)} samples ({y_train.sum()} ruptured)")
    print(f"    Val: {len(X_val)} samples ({y_val.sum()} ruptured)")
    
    # Class weights
    class_weights = compute_class_weight('balanced', classes=np.unique(y_train), y=y_train)
    class_weight_dict = {i: float(w) for i, w in enumerate(class_weights)}
    
    # Early stopping
    early_stopping = callbacks.EarlyStopping(
        monitor='val_auc',
        patience=15,
        restore_best_weights=True,
        mode='max',
        verbose=0
    )
    
    # Reduce learning rate
    lr_scheduler = callbacks.ReduceLROnPlateau(
        monitor='val_loss',
        factor=0.5,
        patience=5,
        min_lr=1e-7,
        verbose=0
    )
    
    # Adaptive batch size
    batch_size = min(8, max(2, len(X_train) // 5))
    
    print(f"    Training (batch_size={batch_size})...")
    
    # Train
    history = transformer_model.fit(
        X_train, y_train,
        validation_data=(X_val, y_val),
        epochs=100,
        batch_size=batch_size,
        class_weight=class_weight_dict,
        callbacks=[early_stopping, lr_scheduler],
        verbose=0
    )
    
    # Evaluate
    y_pred_proba = transformer_model.predict(X_val, verbose=0).flatten()
    auc = roc_auc_score(y_val, y_pred_proba)
    cv_scores.append(auc)
    cv_predictions.append((y_val, y_pred_proba))
    
    print(f"    AUC: {auc:.4f}")
    print(f"    Best val AUC: {max(history.history['val_auc']):.4f}")

mean_auc = np.mean(cv_scores)
std_auc = np.std(cv_scores)

print(f"\n[OK] Transformer Cross-Validation Complete")
print(f"  Mean AUC: {mean_auc:.4f} ± {std_auc:.4f}")
print(f"  Range: {min(cv_scores):.4f} - {max(cv_scores):.4f}")

# Save results
results = {
    'method': 'Transformer Model for Gene Expression',
    'mean_auc': float(mean_auc),
    'std_auc': float(std_auc),
    'cv_scores': [float(s) for s in cv_scores],
    'n_samples': int(len(aligned_labels)),
    'n_genes': int(n_genes),
    'architecture': {
        'd_model': int(d_model),
        'num_heads': 4,
        'num_blocks': 2,
        'seq_length': int(n_genes),
        'total_parameters': int(transformer_model.count_params())
    },
    'data_source': 'GSE122897 (real RNA-seq data)',
    'note': 'No simulated or placeholder data used'
}

with open("results/advanced_ml/transformer_results.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n[OK] Results saved to results/advanced_ml/transformer_results.json")

# Visualization
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    fold_labels = [f'Fold {i+1}' for i in range(len(cv_scores))]
    bars = ax.bar(fold_labels, cv_scores, color='#e74c3c', alpha=0.8, edgecolor='black')
    ax.axhline(y=mean_auc, color='blue', linestyle='--', linewidth=2,
              label=f'Mean: {mean_auc:.4f} ± {std_auc:.4f}')
    ax.fill_between(range(len(cv_scores)+1), 
                    mean_auc - std_auc, mean_auc + std_auc, 
                    alpha=0.2, color='blue', label=f'±1 SD')
    
    for bar, score in zip(bars, cv_scores):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{score:.3f}', ha='center', va='bottom',
               fontweight='bold', fontsize=11)
    
    ax.set_ylabel('AUC-ROC', fontweight='bold', fontsize=14)
    ax.set_title('Transformer Model: Cross-Validation Performance\n' +
                f'Mean AUC: {mean_auc:.4f} ± {std_auc:.4f}',
                fontweight='bold', fontsize=16)
    ax.set_ylim([0, 1])
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig("results/figures/advanced_ml/transformer_cv.png",
               dpi=600, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"[OK] Visualization saved to results/figures/advanced_ml/transformer_cv.png")
    
except Exception as e:
    print(f"[WARNING] Visualization failed: {e}")

print("\n" + "=" * 80)
print("TRANSFORMER MODEL COMPLETE!")
print("=" * 80)
print(f"\nResults:")
print(f"  Mean AUC: {mean_auc:.4f} ± {std_auc:.4f}")
print(f"  Improvement vs baseline (0.739): {((mean_auc - 0.739) / 0.739 * 100):+.1f}%")
print(f"\nNote: All analysis used REAL data - no simulated or placeholder data")

