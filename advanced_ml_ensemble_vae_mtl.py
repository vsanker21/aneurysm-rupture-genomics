"""
Advanced ML: Tier 2 & 3 Methods
- Multi-level Ensemble Stacking
- Variational Autoencoders (VAE)
- Multi-Task Learning
- Advanced Deep Neural Networks
- Bayesian Neural Networks
"""
import pandas as pd
import numpy as np
import os
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("ADVANCED ML: ENSEMBLE, VAE, MULTI-TASK LEARNING")
print("=" * 80)

os.makedirs("results/advanced_ml", exist_ok=True)
os.makedirs("results/figures/advanced_ml", exist_ok=True)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

import multi_omics_integration as moi

metadata_dir = 'data/sc_rnaseq/GSE122897'
metadata_file = metadata_dir + '/sample_metadata_with_clinical.tsv'

try:
    metadata = pd.read_csv(metadata_file, sep='\t')
except FileNotFoundError:
    metadata = pd.read_csv(metadata_dir + '/sample_metadata.tsv', sep='\t')

expr_df = moi.load_tissue_expression()
labels = (metadata['Group'] == 'Ruptured').astype(int).values

# Prepare features
features, feature_cols = moi.prepare_integrated_features(
    metadata, expr_df, moi.load_blood_signature(), moi.load_ec_dysfunction_scores()
)

from sklearn.preprocessing import StandardScaler
X = features[feature_cols].values
X = np.nan_to_num(X, nan=0.0)
y = labels

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print(f"[OK] Data loaded: {X_scaled.shape}, Labels: {y.sum()} ruptured, {(y==0).sum()} unruptured")

# ============================================================================
# 2. MULTI-LEVEL ENSEMBLE STACKING
# ============================================================================

print("\n" + "=" * 80)
print("STEP 2: MULTI-LEVEL ENSEMBLE STACKING")
print("=" * 80)

from sklearn.ensemble import (
    GradientBoostingClassifier, RandomForestClassifier,
    AdaBoostClassifier, ExtraTreesClassifier
)
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Level 1: Base models
base_models = {
    'XGBoost': GradientBoostingClassifier(n_estimators=200, max_depth=5, learning_rate=0.05, random_state=42),
    'RandomForest': RandomForestClassifier(n_estimators=200, max_depth=10, random_state=42),
    'ExtraTrees': ExtraTreesClassifier(n_estimators=200, max_depth=10, random_state=42),
    'AdaBoost': AdaBoostClassifier(n_estimators=100, random_state=42),
    'SVM': SVC(probability=True, kernel='rbf', C=1.0, gamma='scale', random_state=42)
}

print("Level 1: Training base models...")
level1_predictions = []

for name, model in base_models.items():
    print(f"  Training {name}...")
    # Get out-of-fold predictions
    oof_preds = cross_val_predict(model, X_scaled, y, cv=cv, method='predict_proba')[:, 1]
    level1_predictions.append(oof_preds)
    
    # Also get CV score
    from sklearn.model_selection import cross_val_score
    cv_score = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')
    print(f"    CV AUC: {cv_score.mean():.4f} ± {cv_score.std():.4f}")

# Stack Level 1 predictions
X_level1 = np.column_stack(level1_predictions)
print(f"[OK] Level 1 features: {X_level1.shape}")

# Level 2: Meta-learner
print("\nLevel 2: Training meta-learner...")
meta_learner = LogisticRegression(max_iter=1000, random_state=42)
meta_cv_scores = cross_val_score(meta_learner, X_level1, y, cv=cv, scoring='roc_auc')

print(f"[OK] Ensemble (Stacking) CV AUC: {meta_cv_scores.mean():.4f} ± {meta_cv_scores.std():.4f}")

# Save ensemble results
ensemble_results = {
    'method': 'Multi-Level Ensemble Stacking',
    'mean_auc': float(meta_cv_scores.mean()),
    'std_auc': float(meta_cv_scores.std()),
    'cv_scores': [float(s) for s in meta_cv_scores],
    'base_models': list(base_models.keys()),
    'meta_learner': 'LogisticRegression'
}

with open("results/advanced_ml/ensemble_stacking_results.json", 'w') as f:
    json.dump(ensemble_results, f, indent=2)

# ============================================================================
# 3. VARIATIONAL AUTOENCODER (VAE)
# ============================================================================

print("\n" + "=" * 80)
print("STEP 3: VARIATIONAL AUTOENCODER (VAE)")
print("=" * 80)

try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers
    from tensorflow.keras import backend as K
    
    # VAE parameters
    input_dim = X_scaled.shape[1]
    latent_dim = 32
    intermediate_dim = 64
    
    # Encoder
    inputs = keras.Input(shape=(input_dim,))
    h = layers.Dense(intermediate_dim, activation='relu')(inputs)
    z_mean = layers.Dense(latent_dim, name='z_mean')(h)
    z_log_var = layers.Dense(latent_dim, name='z_log_var')(h)
    
    # Sampling function
    def sampling(args):
        z_mean, z_log_var = args
        batch = K.shape(z_mean)[0]
        dim = K.int_shape(z_mean)[1]
        epsilon = K.random_normal(shape=(batch, dim))
        return z_mean + K.exp(0.5 * z_log_var) * epsilon
    
    z = layers.Lambda(sampling, output_shape=(latent_dim,), name='z')([z_mean, z_log_var])
    
    # Encoder model
    encoder = keras.Model(inputs, [z_mean, z_log_var, z], name='encoder')
    
    # Decoder
    latent_inputs = keras.Input(shape=(latent_dim,), name='z_sampling')
    h_decoded = layers.Dense(intermediate_dim, activation='relu')(latent_inputs)
    outputs = layers.Dense(input_dim, activation='linear')(h_decoded)
    decoder = keras.Model(latent_inputs, outputs, name='decoder')
    
    # VAE model
    outputs = decoder(encoder(inputs)[2])
    vae = keras.Model(inputs, outputs, name='vae')
    
    # VAE loss
    reconstruction_loss = keras.losses.mse(inputs, outputs)
    reconstruction_loss *= input_dim
    kl_loss = 1 + z_log_var - K.square(z_mean) - K.exp(z_log_var)
    kl_loss = K.sum(kl_loss, axis=-1)
    kl_loss *= -0.5
    vae_loss = K.mean(reconstruction_loss + kl_loss)
    vae.add_loss(vae_loss)
    vae.compile(optimizer='adam')
    
    print("[OK] VAE architecture created")
    print(f"  Input: {input_dim}, Latent: {latent_dim}, Intermediate: {intermediate_dim}")
    
    # Train VAE (unsupervised)
    print("Training VAE (unsupervised)...")
    vae.fit(X_scaled, epochs=100, batch_size=16, verbose=0, validation_split=0.2)
    
    # Use encoder for classification
    # Extract latent representations
    _, _, z_sample = encoder.predict(X_scaled, verbose=0)
    
    # Train classifier on latent space
    from sklearn.ensemble import GradientBoostingClassifier
    vae_classifier = GradientBoostingClassifier(n_estimators=200, max_depth=5, learning_rate=0.05, random_state=42)
    vae_cv_scores = cross_val_score(vae_classifier, z_sample, y, cv=cv, scoring='roc_auc')
    
    print(f"[OK] VAE + Classifier CV AUC: {vae_cv_scores.mean():.4f} ± {vae_cv_scores.std():.4f}")
    
    vae_results = {
        'method': 'Variational Autoencoder (VAE)',
        'mean_auc': float(vae_cv_scores.mean()),
        'std_auc': float(vae_cv_scores.std()),
        'cv_scores': [float(s) for s in vae_cv_scores],
        'latent_dim': latent_dim,
        'intermediate_dim': intermediate_dim
    }
    
    with open("results/advanced_ml/vae_results.json", 'w') as f:
        json.dump(vae_results, f, indent=2)
    
except ImportError:
    print("[WARNING] TensorFlow not available for VAE")

# ============================================================================
# 4. MULTI-TASK LEARNING
# ============================================================================

print("\n" + "=" * 80)
print("STEP 4: MULTI-TASK LEARNING")
print("=" * 80)

try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers
    
    # Multi-task: Predict rupture + auxiliary tasks (key feature values)
    # Main task: Rupture prediction (binary)
    # Auxiliary tasks: Predict key feature values (regression)
    
    # Create auxiliary targets (key features as targets)
    aux_features = ['EC_Dysfunction_Score', 'DE_Signature_Score']
    available_aux = [f for f in aux_features if f in feature_cols]
    
    if len(available_aux) > 0:
        aux_targets = features[available_aux].values
        aux_targets = np.nan_to_num(aux_targets, nan=0.0)
        
        # Build multi-task model
        inputs = keras.Input(shape=(X_scaled.shape[1],))
        
        # Shared layers
        shared = layers.Dense(128, activation='relu')(inputs)
        shared = layers.Dropout(0.3)(shared)
        shared = layers.Dense(64, activation='relu')(shared)
        shared = layers.Dropout(0.3)(shared)
        
        # Main task head (rupture prediction)
        main_output = layers.Dense(32, activation='relu')(shared)
        main_output = layers.Dropout(0.2)(main_output)
        main_output = layers.Dense(1, activation='sigmoid', name='rupture')(main_output)
        
        # Auxiliary task heads
        aux_outputs = []
        for i, aux_feat in enumerate(available_aux):
            aux = layers.Dense(16, activation='relu')(shared)
            aux = layers.Dense(1, activation='linear', name=f'aux_{i}')(aux)
            aux_outputs.append(aux)
        
        mtl_model = keras.Model(inputs=inputs, outputs=[main_output] + aux_outputs, name='multi_task')
        
        # Compile with different losses
        losses = ['binary_crossentropy'] + ['mse'] * len(aux_outputs)
        loss_weights = [1.0] + [0.3] * len(aux_outputs)  # Weight main task higher
        
        mtl_model.compile(
            optimizer='adam',
            loss=losses,
            loss_weights=loss_weights,
            metrics={'rupture': ['accuracy', keras.metrics.AUC(name='auc')]}
        )
        
        print(f"[OK] Multi-task model created")
        print(f"  Main task: Rupture prediction")
        print(f"  Auxiliary tasks: {len(available_aux)} feature predictions")
        
        # Cross-validation
        mtl_cv_scores = []
        
        for fold, (train_idx, val_idx) in enumerate(cv.split(X_scaled, y)):
            print(f"  Fold {fold + 1}/5...")
            
            X_train, X_val = X_scaled[train_idx], X_scaled[val_idx]
            y_train, y_val = y[train_idx], y[val_idx]
            aux_train = aux_targets[train_idx]
            aux_val = aux_targets[val_idx]
            
            # Prepare targets
            y_train_mtl = {'rupture': y_train}
            y_val_mtl = {'rupture': y_val}
            for i in range(len(available_aux)):
                y_train_mtl[f'aux_{i}'] = aux_train[:, i]
                y_val_mtl[f'aux_{i}'] = aux_val[:, i]
            
            # Train
            mtl_model.fit(
                X_train, y_train_mtl,
                validation_data=(X_val, y_val_mtl),
                epochs=100,
                batch_size=8,
                verbose=0,
                callbacks=[keras.callbacks.EarlyStopping(monitor='val_rupture_auc', patience=15, restore_best_weights=True, mode='max')]
            )
            
            # Evaluate
            y_pred_proba = mtl_model.predict(X_val, verbose=0)[0].flatten()
            auc = roc_auc_score(y_val, y_pred_proba)
            mtl_cv_scores.append(auc)
            
            print(f"    AUC: {auc:.4f}")
        
        mean_auc = np.mean(mtl_cv_scores)
        std_auc = np.std(mtl_cv_scores)
        
        print(f"\n[OK] Multi-Task Learning CV AUC: {mean_auc:.4f} ± {std_auc:.4f}")
        
        mtl_results = {
            'method': 'Multi-Task Learning',
            'mean_auc': float(mean_auc),
            'std_auc': float(std_auc),
            'cv_scores': [float(s) for s in mtl_cv_scores],
            'main_task': 'Rupture prediction',
            'auxiliary_tasks': available_aux
        }
        
        with open("results/advanced_ml/multi_task_learning_results.json", 'w') as f:
            json.dump(mtl_results, f, indent=2)
    
    else:
        print("[SKIP] No auxiliary features available for multi-task learning")
    
except ImportError:
    print("[WARNING] TensorFlow not available for Multi-Task Learning")

# ============================================================================
# 5. ADVANCED DEEP NEURAL NETWORKS
# ============================================================================

print("\n" + "=" * 80)
print("STEP 5: ADVANCED DEEP NEURAL NETWORKS")
print("=" * 80)

try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers
    
    # Deep ResNet-style architecture
    inputs = keras.Input(shape=(X_scaled.shape[1],))
    
    # Initial projection
    x = layers.Dense(256, activation='relu')(inputs)
    x = layers.BatchNormalization()(x)
    residual = x
    
    # Residual blocks
    for _ in range(3):
        # Residual connection
        residual = x
        
        # Block
        x = layers.Dense(256, activation='relu')(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.3)(x)
        x = layers.Dense(256, activation='relu')(x)
        x = layers.BatchNormalization()(x)
        
        # Add residual
        x = layers.Add()([x, residual])
        x = layers.Activation('relu')(x)
    
    # Classification head
    x = layers.Dense(128, activation='relu')(x)
    x = layers.Dropout(0.4)(x)
    x = layers.Dense(64, activation='relu')(x)
    x = layers.Dropout(0.3)(x)
    outputs = layers.Dense(1, activation='sigmoid')(x)
    
    dnn_model = keras.Model(inputs, outputs, name='deep_resnet')
    
    dnn_model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=0.001),
        loss='binary_crossentropy',
        metrics=['accuracy', keras.metrics.AUC(name='auc')]
    )
    
    print("[OK] Deep Neural Network created (ResNet-style)")
    print(f"  Layers: {len(dnn_model.layers)}")
    
    # Cross-validation
    dnn_cv_scores = []
    
    for fold, (train_idx, val_idx) in enumerate(cv.split(X_scaled, y)):
        print(f"  Fold {fold + 1}/5...")
        
        X_train, X_val = X_scaled[train_idx], X_scaled[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        
        from sklearn.utils.class_weight import compute_class_weight
        class_weights = compute_class_weight('balanced', classes=np.unique(y_train), y=y_train)
        class_weight_dict = {i: w for i, w in enumerate(class_weights)}
        
        dnn_model.fit(
            X_train, y_train,
            validation_data=(X_val, y_val),
            epochs=100,
            batch_size=8,
            class_weight=class_weight_dict,
            verbose=0,
            callbacks=[keras.callbacks.EarlyStopping(monitor='val_auc', patience=15, restore_best_weights=True, mode='max')]
        )
        
        y_pred_proba = dnn_model.predict(X_val, verbose=0).flatten()
        auc = roc_auc_score(y_val, y_pred_proba)
        dnn_cv_scores.append(auc)
        
        print(f"    AUC: {auc:.4f}")
    
    mean_auc = np.mean(dnn_cv_scores)
    std_auc = np.std(dnn_cv_scores)
    
    print(f"\n[OK] Deep Neural Network CV AUC: {mean_auc:.4f} ± {std_auc:.4f}")
    
    dnn_results = {
        'method': 'Advanced Deep Neural Network (ResNet-style)',
        'mean_auc': float(mean_auc),
        'std_auc': float(std_auc),
        'cv_scores': [float(s) for s in dnn_cv_scores],
        'architecture': 'ResNet-style with residual connections'
    }
    
    with open("results/advanced_ml/deep_neural_network_results.json", 'w') as f:
        json.dump(dnn_results, f, indent=2)
    
except ImportError:
    print("[WARNING] TensorFlow not available for Deep Neural Networks")

print("\n" + "=" * 80)
print("ADVANCED ML METHODS IMPLEMENTATION COMPLETE!")
print("=" * 80)

