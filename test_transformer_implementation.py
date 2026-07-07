"""
Test Transformer Model Implementation
Verify it works correctly with real data
"""
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("TESTING TRANSFORMER MODEL IMPLEMENTATION")
print("=" * 80)

# Import libraries
try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers
    print("[OK] TensorFlow imported successfully")
except ImportError as e:
    print(f"[ERROR] TensorFlow not available: {e}")
    print("[ERROR] Please install: pip install tensorflow")
    exit(1)

# Load data
print("\nLoading data...")
import multi_omics_integration as moi

metadata_dir = 'data/sc_rnaseq/GSE122897'
try:
    metadata = pd.read_csv(metadata_dir + '/sample_metadata_with_clinical.tsv', sep='\t')
except:
    metadata = pd.read_csv(metadata_dir + '/sample_metadata.tsv', sep='\t')

expr_df = moi.load_tissue_expression()
labels = (metadata['Group'] == 'Ruptured').astype(int).values

if expr_df is None:
    print("[ERROR] Expression data not available")
    exit(1)

# Prepare data
expr_df = expr_df.T
expr_df = expr_df.reindex(metadata['Sample'])
expr_df = expr_df.dropna()

aligned_labels = labels[metadata['Sample'].isin(expr_df.index)]

# Select top genes
gene_vars = expr_df.var(axis=0)
top_genes = gene_vars.nlargest(5000).index.tolist()
expr_subset = expr_df[top_genes].values
expr_subset = np.log2(expr_subset + 1)

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_scaled = scaler.fit_transform(expr_subset)

print(f"[OK] Data prepared: {X_scaled.shape}")
print(f"  Samples: {X_scaled.shape[0]}")
print(f"  Genes: {X_scaled.shape[1]}")
print(f"  Ruptured: {aligned_labels.sum()}, Unruptured: {(aligned_labels == 0).sum()}")

# Test Transformer architecture
print("\nTesting Transformer architecture...")

n_genes = X_scaled.shape[1]
d_model = 128

try:
    # Build model
    inputs = keras.Input(shape=(n_genes,), name='gene_expression')
    
    # Reshape and project
    x = layers.Reshape((n_genes, 1))(inputs)
    x = layers.Dense(d_model, activation=None)(x)
    x = layers.Dropout(0.1)(x)
    
    # Positional encoding
    pos_emb = layers.Embedding(n_genes, d_model)(tf.range(n_genes))
    pos_emb = tf.expand_dims(pos_emb, 0)
    x = x + pos_emb
    
    # Transformer encoder
    def transformer_encoder(inputs, num_heads=4, ff_dim=256, dropout=0.1):
        attention_output = layers.MultiHeadAttention(
            num_heads=num_heads, key_dim=d_model, dropout=dropout
        )(inputs, inputs)
        attention_output = layers.Dropout(dropout)(attention_output)
        attention_output = layers.LayerNormalization(epsilon=1e-6)(inputs + attention_output)
        
        ffn_output = layers.Dense(ff_dim, activation='relu')(attention_output)
        ffn_output = layers.Dense(d_model)(ffn_output)
        ffn_output = layers.Dropout(dropout)(ffn_output)
        output = layers.LayerNormalization(epsilon=1e-6)(attention_output + ffn_output)
        
        return output
    
    x = transformer_encoder(x, num_heads=4, ff_dim=256, dropout=0.2)
    x = layers.GlobalAveragePooling1D()(x)
    
    # Classification head
    x = layers.Dense(64, activation='relu')(x)
    x = layers.Dropout(0.3)(x)
    outputs = layers.Dense(1, activation='sigmoid', name='classification')(x)
    
    model = keras.Model(inputs, outputs, name='transformer_test')
    
    model.compile(
        optimizer='adam',
        loss='binary_crossentropy',
        metrics=['accuracy', keras.metrics.AUC(name='auc')]
    )
    
    print(f"[OK] Transformer model built successfully")
    print(f"  Input shape: {model.input_shape}")
    print(f"  Output shape: {model.output_shape}")
    print(f"  Total parameters: {model.count_params():,}")
    
    # Test forward pass
    print("\nTesting forward pass...")
    test_input = X_scaled[:5]  # First 5 samples
    test_output = model.predict(test_input, verbose=0)
    
    print(f"[OK] Forward pass successful")
    print(f"  Input shape: {test_input.shape}")
    print(f"  Output shape: {test_output.shape}")
    print(f"  Output range: [{test_output.min():.4f}, {test_output.max():.4f}]")
    
    # Test training step
    print("\nTesting training step...")
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(
        X_scaled, aligned_labels, test_size=0.2, random_state=42, stratify=aligned_labels
    )
    
    # Single training step
    history = model.fit(
        X_train[:10], y_train[:10],  # Use small subset for test
        validation_data=(X_test[:5], y_test[:5]),
        epochs=1,
        batch_size=4,
        verbose=0
    )
    
    print(f"[OK] Training step successful")
    print(f"  Loss: {history.history['loss'][0]:.4f}")
    print(f"  Accuracy: {history.history['accuracy'][0]:.4f}")
    
    print("\n" + "=" * 80)
    print("TRANSFORMER MODEL TEST: SUCCESSFUL!")
    print("=" * 80)
    print("\nThe transformer model is working correctly with real data.")
    print("Ready for full training with advanced_ml_gnn_transformer.py")
    
except Exception as e:
    print(f"\n[ERROR] Transformer test failed: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

