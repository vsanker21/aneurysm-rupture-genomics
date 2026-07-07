"""Run each analysis individually with error handling"""
import sys
import traceback
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve
import os

warnings.filterwarnings('ignore')

# Ultra HD settings
ULTRA_HD_DPI = 600
ULTRA_HD_SIZE = (16, 12)
plt.rcParams['figure.dpi'] = ULTRA_HD_DPI
plt.rcParams['savefig.dpi'] = ULTRA_HD_DPI
sns.set_style("whitegrid")

# Setup
output_dir = Path("results/advanced_analyses_ultra_hd")
output_dir.mkdir(parents=True, exist_ok=True)

print("Loading data...")
features_df = pd.read_csv('results/integrated_features.tsv', sep='\t')
features_df = features_df[features_df['Group'].isin(['Ruptured', 'Unruptured'])].copy()
labels = (features_df['Group'] == 'Ruptured').astype(int)

key_features = ['DE_Signature_Score', 'EC_Dysfunction_Score', 'ECM_Remodeling_Score']
X = features_df[key_features].fillna(0).values
y = labels.values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print(f"Loaded {len(features_df)} samples")
print(f"Ruptured: {y.sum()}, Unruptured: {(y==0).sum()}\n")

# Import the main script's analysis functions
# We'll run analyses 3-8 individually
analyses_to_run = [3, 4, 5, 6, 7, 8]

for analysis_num in analyses_to_run:
    print(f"\n{'='*80}")
    print(f"RUNNING ANALYSIS {analysis_num}")
    print(f"{'='*80}")
    try:
        # Read the main script and extract the analysis code
        with open('advanced_ml_ultra_hd_analyses.py', 'r', encoding='utf-8') as f:
            script_content = f.read()
        
        # Find the analysis section
        analysis_markers = {
            3: "ANALYSIS 3: TIME-TO-EVENT",
            4: "ANALYSIS 4: GRAPH NEURAL NETWORKS",
            5: "ANALYSIS 5: SHAP",
            6: "ANALYSIS 6: MULTI-TASK LEARNING",
            7: "ANALYSIS 7: LATENT FACTOR",
            8: "ANALYSIS 8: ROBUSTNESS"
        }
        
        marker = analysis_markers[analysis_num]
        start_idx = script_content.find(marker)
        
        if start_idx == -1:
            print(f"  [ERROR] Could not find {marker}")
            continue
        
        # Find the next analysis or end of script
        next_markers = [v for k, v in analysis_markers.items() if k > analysis_num]
        if next_markers:
            end_idx = min([script_content.find(m, start_idx + 1) for m in next_markers if script_content.find(m, start_idx + 1) != -1] + [len(script_content)])
        else:
            end_idx = len(script_content)
        
        # Extract and execute the analysis code
        analysis_code = script_content[start_idx:end_idx]
        
        # Create a namespace with all needed variables
        namespace = {
            'pd': pd, 'np': np, 'plt': plt, 'sns': sns, 'Path': Path,
            'warnings': warnings, 'stats': stats, 'StandardScaler': StandardScaler,
            'GradientBoostingClassifier': GradientBoostingClassifier,
            'RandomForestClassifier': RandomForestClassifier,
            'roc_auc_score': roc_auc_score, 'roc_curve': roc_curve,
            'os': os, 'output_dir': output_dir, 'X_scaled': X_scaled,
            'y': y, 'key_features': key_features, 'features_df': features_df,
            'scaler': scaler, 'ULTRA_HD_DPI': ULTRA_HD_DPI,
            'ULTRA_HD_SIZE': ULTRA_HD_SIZE
        }
        
        exec(analysis_code, namespace)
        print(f"  [OK] Analysis {analysis_num} completed successfully")
        
    except Exception as e:
        print(f"  [ERROR] Analysis {analysis_num} failed: {e}")
        traceback.print_exc()
        print()

print(f"\n{'='*80}")
print("ALL ANALYSES COMPLETE")
print(f"{'='*80}")

