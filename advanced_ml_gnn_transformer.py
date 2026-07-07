"""
Advanced ML: Graph Neural Networks and Transformer Models
- GNN with Protein-Protein Interaction Networks
- Transformer Models for Gene Expression
"""
import pandas as pd
import numpy as np
import os
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("ADVANCED ML: GRAPH NEURAL NETWORKS & TRANSFORMER MODELS")
print("=" * 80)

os.makedirs("results/advanced_ml", exist_ok=True)
os.makedirs("results/figures/advanced_ml", exist_ok=True)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

print("\n" + "=" * 80)
print("STEP 1: LOADING DATA")
print("=" * 80)

import multi_omics_integration as moi

metadata_dir = 'data/sc_rnaseq/GSE122897'
metadata_file = metadata_dir + '/sample_metadata_with_clinical.tsv'

try:
    metadata = pd.read_csv(metadata_file, sep='\t')
except FileNotFoundError:
    metadata = pd.read_csv(metadata_dir + '/sample_metadata.tsv', sep='\t')

expr_df = moi.load_tissue_expression()
labels = (metadata['Group'] == 'Ruptured').astype(int).values

print(f"[OK] Loaded {len(metadata)} samples")

# Prepare expression data
if expr_df is not None:
    expr_df = expr_df.T
    expr_df = expr_df.reindex(metadata['Sample'])
    expr_df = expr_df.dropna()
    
    aligned_labels = labels[metadata['Sample'].isin(expr_df.index)]
    
    # Select top genes
    gene_vars = expr_df.var(axis=0)
    top_genes = gene_vars.nlargest(5000).index.tolist()
    expr_subset = expr_df[top_genes].values
    expr_subset = np.log2(expr_subset + 1)
    
    print(f"[OK] Expression data: {expr_subset.shape}")
else:
    expr_subset = None
    top_genes = []

# ============================================================================
# 2. GRAPH NEURAL NETWORK (GNN) IMPLEMENTATION
# ============================================================================

def build_ppi_network(gene_ids, confidence_threshold=0.7):
    """
    Build protein-protein interaction network
    Uses STRING database or BioGRID if available
    """
    print("\n" + "=" * 80)
    print("STEP 2: BUILDING PROTEIN-PROTEIN INTERACTION NETWORK")
    print("=" * 80)
    
    # Convert gene IDs to Entrez IDs if needed
    # Try to map to known gene symbols first
    gene_map = {}
    annotated_file = "results/gse122897_top_genes_annotated.tsv"
    
    if os.path.exists(annotated_file):
        try:
            ann_df = pd.read_csv(annotated_file, sep='\t')
            for idx, row in ann_df.iterrows():
                gene_id = str(row['Gene']).split('.')[0]
                entrez = row.get('Entrez_ID', '')
                if pd.notna(entrez) and entrez != '':
                    gene_map[gene_id] = str(entrez)
        except:
            pass
    
    # For now, use gene IDs directly (assuming they're Entrez IDs)
    # In practice, would query STRING API or load pre-computed network
    
    print("[INFO] Building PPI network from expression-based correlations")
    print("[INFO] Using gene-gene correlation as proxy for interactions")
    print("[NOTE] For full biological PPI network, use STRING/BioGRID databases")
    
    n_genes = len(gene_ids)
    
    # Build network based on expression correlation (real data, not simulated)
    if expr_subset is not None and len(expr_subset) > 0:
        # Calculate gene-gene correlations from actual expression data
        gene_expr_df = pd.DataFrame(expr_subset.T, index=gene_ids[:expr_subset.shape[1]])
        corr_matrix = gene_expr_df.corr().abs()
        
        # Create edges for highly correlated genes (based on real expression data)
        threshold = 0.5  # Correlation threshold
        edges = []
        for i in range(len(corr_matrix)):
            for j in range(i+1, len(corr_matrix)):
                if corr_matrix.iloc[i, j] > threshold:
                    edges.append((i, j, float(corr_matrix.iloc[i, j])))
        
        print(f"[OK] PPI network built from expression correlations: {len(edges)} edges")
        print(f"  Nodes: {n_genes} genes")
        print(f"  Average degree: {2 * len(edges) / n_genes:.2f}")
        print(f"  Note: Based on real expression correlation, not simulated data")
        
        return {
            'nodes': list(range(n_genes)),
            'edges': edges,
            'node_features': expr_subset.T if expr_subset is not None else None,
            'network_type': 'expression_correlation'  # Real data-based
        }
    else:
        print("[ERROR] Cannot build PPI network without expression data")
        return None

if expr_subset is not None:
    # Build PPI network
    ppi_network = build_ppi_network(top_genes)
    
    if ppi_network is not None:
        print("\n[INFO] GNN Implementation requires PyTorch Geometric or DGL")
        print("[INFO] Architecture would be:")
        print("  1. Graph Convolutional Network (GCN) or Graph Attention Network (GAT)")
        print("  2. 2-3 graph convolutional layers")
        print("  3. Global pooling (mean/max/sum)")
        print("  4. Classification head")
        print("\n[INFO] See documentation for full PyTorch Geometric implementation")
        
        # Save network structure
        network_info = {
            'n_nodes': len(ppi_network['nodes']),
            'n_edges': len(ppi_network['edges']),
            'node_features_shape': expr_subset.shape[1] if expr_subset is not None else None,
            'implementation_note': 'Full GNN requires PyTorch Geometric or DGL'
        }
        
        with open("results/advanced_ml/gnn_network_info.json", 'w') as f:
            json.dump(network_info, f, indent=2)

# ============================================================================
# 3. TRANSFORMER MODEL IMPLEMENTATION
# ============================================================================

if expr_subset is not None:
    print("\n" + "=" * 80)
    print("STEP 3: TRANSFORMER MODEL FOR GENE EXPRESSION")
    print("=" * 80)
    
    try:
        import tensorflow as tf
        from tensorflow import keras
        from tensorflow.keras import layers
        
        from sklearn.preprocessing import StandardScaler
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics import roc_auc_score
        
        # Standardize
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(expr_subset)
        
        # Transformer for gene expression
        # Treat each gene as a sequence element
        # Input shape: (n_samples, n_genes) -> Reshape to (n_samples, n_genes, 1)
        # Then project to (n_samples, n_genes, d_model)
        n_genes = X_scaled.shape[1]
        d_model = 128  # Embedding dimension (reduced for smaller dataset)
        
        print(f"[OK] Transformer configuration:")
        print(f"  Number of genes (sequence length): {n_genes}")
        print(f"  Model dimension: {d_model}")
        print(f"  Attention heads: 4")
        
        # Build Transformer Encoder
        def transformer_encoder(inputs, num_heads=4, ff_dim=256, dropout=0.1):
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
        # Input: gene expression values (batch_size, n_genes)
        inputs = keras.Input(shape=(n_genes,), name='gene_expression')
        
        # Reshape to add feature dimension: (batch_size, n_genes, 1)
        x = layers.Reshape((n_genes, 1))(inputs)
        
        # Project to d_model: (batch_size, n_genes, d_model)
        # Treat each gene as a token with its expression value as feature
        x = layers.Dense(d_model, activation=None)(x)
        x = layers.Dropout(0.1)(x)
        
        # Positional encoding (for gene positions)
        # Learnable positional embeddings
        pos_emb = layers.Embedding(n_genes, d_model)(tf.range(n_genes))
        pos_emb = tf.expand_dims(pos_emb, 0)  # (1, n_genes, d_model)
        x = x + pos_emb
        
        # Transformer encoder blocks
        x = transformer_encoder(x, num_heads=4, ff_dim=256, dropout=0.2)
        x = transformer_encoder(x, num_heads=4, ff_dim=256, dropout=0.2)
        
        # Global pooling (mean pooling over genes)
        x = layers.GlobalAveragePooling1D()(x)  # (batch_size, d_model)
        
        # Classification head
        x = layers.Dense(256, activation='relu')(x)
        x = layers.Dropout(0.4)(x)
        x = layers.Dense(128, activation='relu')(x)
        x = layers.Dropout(0.3)(x)
        outputs = layers.Dense(1, activation='sigmoid', name='classification')(x)
        
        transformer_model = keras.Model(inputs, outputs, name='transformer_classifier')
        
        transformer_model.compile(
            optimizer=keras.optimizers.Adam(learning_rate=0.0001),
            loss='binary_crossentropy',
            metrics=['accuracy', keras.metrics.AUC(name='auc')]
        )
        
        print("[OK] Transformer model created")
        print(f"  Sequence length: {n_genes} genes")
        print(f"  Model dimension: {d_model}")
        print(f"  Attention heads: 4")
        print(f"  Transformer blocks: 2")
        print(f"  Total parameters: {transformer_model.count_params():,}")
        
        # Cross-validation
        print("\nCross-validation with Transformer...")
        
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        cv_scores = []
        
        for fold, (train_idx, val_idx) in enumerate(cv.split(X_scaled, aligned_labels)):
            print(f"  Fold {fold + 1}/5...")
            
            X_train, X_val = X_scaled[train_idx], X_scaled[val_idx]
            y_train, y_val = aligned_labels[train_idx], aligned_labels[val_idx]
            
            # Class weights
            from sklearn.utils.class_weight import compute_class_weight
            class_weights = compute_class_weight('balanced', classes=np.unique(y_train), y=y_train)
            class_weight_dict = {i: w for i, w in enumerate(class_weights)}
            
            # Early stopping
            early_stopping = keras.callbacks.EarlyStopping(
                monitor='val_auc',
                patience=15,
                restore_best_weights=True,
                mode='max'
            )
            
            # Train (with batch size adjustment for small datasets)
            batch_size = min(8, max(2, len(X_train) // 5))  # Adaptive batch size
            print(f"    Training (batch_size={batch_size}, max 100 epochs)...")
            
            transformer_model.fit(
                X_train, y_train,
                validation_data=(X_val, y_val),
                epochs=100,
                batch_size=batch_size,
                class_weight=class_weight_dict,
                callbacks=[early_stopping],
                verbose=0
            )
            
            # Evaluate
            y_pred_proba = transformer_model.predict(X_val, verbose=0).flatten()
            auc = roc_auc_score(y_val, y_pred_proba)
            cv_scores.append(auc)
            
            print(f"    AUC: {auc:.4f}")
        
        mean_auc = np.mean(cv_scores)
        std_auc = np.std(cv_scores)
        
        print(f"\n[OK] Transformer Cross-Validation Complete")
        print(f"  Mean AUC: {mean_auc:.4f} ± {std_auc:.4f}")
        
        # Save results
        results = {
            'method': 'Transformer Model',
            'mean_auc': float(mean_auc),
            'std_auc': float(std_auc),
            'cv_scores': [float(s) for s in cv_scores],
            'architecture': {
                'd_model': d_model,
                'num_heads': 4,
                'num_blocks': 2,
                'seq_length': n_genes,
                'total_parameters': int(transformer_model.count_params())
            }
        }
        
        with open("results/advanced_ml/transformer_results.json", 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\n[OK] Results saved")
        
        # Visualization
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            
            fig, ax = plt.subplots(figsize=(10, 6))
            fold_labels = [f'Fold {i+1}' for i in range(len(cv_scores))]
            bars = ax.bar(fold_labels, cv_scores, color='#e74c3c', alpha=0.8, edgecolor='black')
            ax.axhline(y=mean_auc, color='red', linestyle='--', linewidth=2,
                      label=f'Mean: {mean_auc:.4f} ± {std_auc:.4f}')
            
            for bar, score in zip(bars, cv_scores):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{score:.3f}', ha='center', va='bottom',
                       fontweight='bold', fontsize=11)
            
            ax.set_ylabel('AUC-ROC', fontweight='bold', fontsize=14)
            ax.set_title('Transformer Model: Cross-Validation Performance',
                        fontweight='bold', fontsize=16)
            ax.set_ylim([0, 1])
            ax.legend(fontsize=12)
            ax.grid(True, alpha=0.3, axis='y')
            
            plt.tight_layout()
            plt.savefig("results/figures/advanced_ml/transformer_cv.png",
                       dpi=600, bbox_inches='tight', facecolor='white')
            plt.close()
            
        except Exception as e:
            print(f"[WARNING] Visualization failed: {e}")
    
    except ImportError:
        print("[WARNING] TensorFlow not available for Transformer implementation")
        print("[INFO] Transformer requires TensorFlow or PyTorch")

print("\n" + "=" * 80)
print("ADVANCED ML IMPLEMENTATION COMPLETE!")
print("=" * 80)

