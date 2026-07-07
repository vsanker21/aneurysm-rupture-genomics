"""
Complete remaining outputs - ensures all files are created
"""
import sys
import os

if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

print("Completing remaining outputs...")
sys.stdout.flush()

os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/gene_lists", exist_ok=True)
os.makedirs("results/risk_scores", exist_ok=True)

# 1. Risk Scores
print("\n[1] Creating risk scores...")
try:
    tissue_pred = pd.read_csv("results/ensemble_tissue_predictions.tsv", sep='\t')
    blood_pred = pd.read_csv("results/ensemble_blood_predictions.tsv", sep='\t')
    
    tissue_pred['Risk_Category'] = pd.cut(tissue_pred['Normalized_Score'], 
                                          bins=[0, 0.33, 0.67, 1.0], 
                                          labels=['Low', 'Medium', 'High'],
                                          include_lowest=True)
    blood_pred['Risk_Category'] = pd.cut(blood_pred['Normalized_Score'], 
                                         bins=[0, 0.33, 0.67, 1.0], 
                                         labels=['Low', 'Medium', 'High'],
                                         include_lowest=True)
    
    tissue_pred.to_csv('results/risk_scores/tissue_risk_scores_categorized.tsv', sep='\t', index=False)
    blood_pred.to_csv('results/risk_scores/blood_risk_scores_categorized.tsv', sep='\t', index=False)
    print("  ✓ Risk scores created")
    
    # Risk stratification plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    tissue_risk_counts = tissue_pred['Risk_Category'].value_counts().sort_index()
    axes[0].bar(tissue_risk_counts.index.astype(str), tissue_risk_counts.values, 
               color=['#2ca02c', '#ff7f0e', '#d62728'], alpha=0.8)
    axes[0].set_xlabel('Risk Category')
    axes[0].set_ylabel('Number of Samples')
    axes[0].set_title('Tissue Risk Stratification')
    
    blood_risk_counts = blood_pred['Risk_Category'].value_counts().sort_index()
    axes[1].bar(blood_risk_counts.index.astype(str), blood_risk_counts.values, 
               color=['#2ca02c', '#ff7f0e', '#d62728'], alpha=0.8)
    axes[1].set_xlabel('Risk Category')
    axes[1].set_ylabel('Number of Samples')
    axes[1].set_title('Blood Risk Stratification')
    
    plt.suptitle('Risk Score Stratification')
    plt.tight_layout()
    plt.savefig('results/figures/risk_stratification.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("  ✓ Risk stratification plot created")
except Exception as e:
    print(f"  ⚠ Error: {e}")

# 2. Gene Lists
print("\n[2] Creating gene lists...")
try:
    gene_lists = []
    
    if os.path.exists('results/multi_omics_tissue_signature_genes.tsv'):
        tissue_sig = pd.read_csv('results/multi_omics_tissue_signature_genes.tsv', sep='\t')
        for _, row in tissue_sig.iterrows():
            gene_lists.append({
                'Gene_ID': str(row.get('Gene', '')),
                'Dataset': 'Tissue (GSE122897)',
                'Log2FC': row.get('Log2FC', 0),
                'P_adj': row.get('P_adj_BH', row.get('P_adj', 1)),
                'Direction': 'Up' if row.get('Log2FC', 0) > 0 else 'Down'
            })
    
    if os.path.exists('results/multi_omics_blood_signature_genes.tsv'):
        blood_sig = pd.read_csv('results/multi_omics_blood_signature_genes.tsv', sep='\t')
        for _, row in blood_sig.iterrows():
            gene_lists.append({
                'Gene_ID': str(row.get('Gene', row.get('Gene_ID', ''))),
                'Dataset': 'Blood (GSE159610)',
                'Log2FC': row.get('Log2FC', 0),
                'P_adj': row.get('P_adj_BH', row.get('P_adj', 1)),
                'Direction': 'Up' if row.get('Log2FC', 0) > 0 else 'Down'
            })
    
    if len(gene_lists) > 0:
        all_genes_df = pd.DataFrame(gene_lists)
        all_genes_df.to_csv('results/gene_lists/comprehensive_signature_genes_annotated.tsv', sep='\t', index=False)
        
        tissue_genes_df = all_genes_df[all_genes_df['Dataset'] == 'Tissue (GSE122897)'].head(50)
        blood_genes_df = all_genes_df[all_genes_df['Dataset'] == 'Blood (GSE159610)'].head(50)
        
        if len(tissue_genes_df) > 0:
            tissue_genes_df.to_csv('results/gene_lists/tissue_top50_signature_genes.tsv', sep='\t', index=False)
        if len(blood_genes_df) > 0:
            blood_genes_df.to_csv('results/gene_lists/blood_top50_signature_genes.tsv', sep='\t', index=False)
        
        print("  ✓ Gene lists created")
    else:
        print("  ⚠ No signature genes found")
except Exception as e:
    print(f"  ⚠ Error: {e}")

# 3. Feature Importance (if available)
print("\n[3] Creating feature importance plot...")
try:
    if os.path.exists('results/feature_importance.tsv'):
        feature_imp = pd.read_csv('results/feature_importance.tsv', sep='\t')
        fig, axes = plt.subplots(2, 2, figsize=(16, 14))
        axes = axes.flatten()
        
        plot_idx = 0
        for model_type in ['Tissue', 'Blood']:
            for model_name in ['Random Forest', 'Gradient Boosting']:
                model_data = feature_imp[(feature_imp['Model'].str.contains(model_type, case=False)) & 
                                         (feature_imp['Model'].str.contains(model_name, case=False))]
                if len(model_data) > 0 and plot_idx < len(axes):
                    top_features = model_data.head(20)
                    if 'Importance' in top_features.columns:
                        feat_col = 'Importance'
                    else:
                        feat_col = top_features.columns[1] if len(top_features.columns) > 1 else top_features.columns[0]
                    
                    if 'Feature' in top_features.columns:
                        feat_names = top_features['Feature'].values
                    else:
                        feat_names = top_features.iloc[:, 0].values
                    
                    axes[plot_idx].barh(range(len(top_features)), top_features[feat_col].values)
                    axes[plot_idx].set_yticks(range(len(top_features)))
                    axes[plot_idx].set_yticklabels(feat_names, fontsize=9)
                    axes[plot_idx].set_title(f'{model_type} {model_name}')
                    axes[plot_idx].invert_yaxis()
                    plot_idx += 1
        
        plt.suptitle('Feature Importance Rankings')
        plt.tight_layout()
        plt.savefig('results/figures/feature_importance_all_models.png', bbox_inches='tight', dpi=300)
        plt.close()
        
        feature_imp.to_csv('results/feature_importance_all_models.tsv', sep='\t', index=False)
        print("  ✓ Feature importance plot created")
    else:
        print("  ⚠ Feature importance file not found")
except Exception as e:
    print(f"  ⚠ Error: {e}")

# 4. Pathway Enrichment
print("\n[4] Creating pathway enrichment plots...")
try:
    if os.path.exists('results/pathway_enrichment_prioritized_fdr_0.05.tsv'):
        tissue_pathways = pd.read_csv('results/pathway_enrichment_prioritized_fdr_0.05.tsv', sep='\t')
        if len(tissue_pathways) > 0:
            top_paths = tissue_pathways.head(10)
            p_col = 'P_adj' if 'P_adj' in top_paths.columns else ('P_value' if 'P_value' in top_paths.columns else None)
            if p_col:
                path_names = top_paths['Term'].values if 'Term' in top_paths.columns else top_paths.index[:10].astype(str)
                p_values = -np.log10(top_paths[p_col].values + 1e-10)
                
                fig, ax = plt.subplots(figsize=(12, 8))
                ax.barh(range(len(top_paths)), p_values, color='#1f77b4', alpha=0.8)
                ax.set_yticks(range(len(top_paths)))
                ax.set_yticklabels(path_names, fontsize=10)
                ax.set_xlabel('-log10(P-value)')
                ax.set_title('Tissue Pathway Enrichment - Top 10')
                ax.invert_yaxis()
                plt.tight_layout()
                plt.savefig('results/figures/pathway_enrichment_tissue.png', bbox_inches='tight', dpi=300)
                plt.close()
                print("  ✓ Tissue pathway plot created")
except Exception as e:
    print(f"  ⚠ Error: {e}")

try:
    if os.path.exists('results/pathway_enrichment_fdr_0.05.tsv'):
        blood_pathways = pd.read_csv('results/pathway_enrichment_fdr_0.05.tsv', sep='\t')
        if len(blood_pathways) > 0:
            top_paths = blood_pathways.head(10)
            p_col = 'P_adj' if 'P_adj' in top_paths.columns else ('P_value' if 'P_value' in top_paths.columns else None)
            if p_col:
                path_names = top_paths['Term'].values if 'Term' in top_paths.columns else top_paths.index[:10].astype(str)
                p_values = -np.log10(top_paths[p_col].values + 1e-10)
                
                fig, ax = plt.subplots(figsize=(12, 6))
                ax.barh(range(len(top_paths)), p_values, color='#ff7f0e', alpha=0.8)
                ax.set_yticks(range(len(top_paths)))
                ax.set_yticklabels(path_names, fontsize=10)
                ax.set_xlabel('-log10(P-value)')
                ax.set_title('Blood Pathway Enrichment')
                ax.invert_yaxis()
                plt.tight_layout()
                plt.savefig('results/figures/pathway_enrichment_blood.png', bbox_inches='tight', dpi=300)
                plt.close()
                print("  ✓ Blood pathway plot created")
except Exception as e:
    print(f"  ⚠ Error: {e}")

print("\n✅ All remaining outputs completed!")
sys.stdout.flush()

