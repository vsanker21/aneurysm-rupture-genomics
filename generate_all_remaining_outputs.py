"""
Generate ALL remaining outputs with explicit error reporting
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

print("=" * 80)
print("GENERATING ALL REMAINING OUTPUTS")
print("=" * 80)
sys.stdout.flush()

os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/gene_lists", exist_ok=True)
os.makedirs("results/risk_scores", exist_ok=True)

# 1. Risk Stratification Plot
print("\n[1] Creating risk stratification plot...")
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
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    tissue_risk_counts = tissue_pred['Risk_Category'].value_counts().sort_index()
    axes[0].bar(tissue_risk_counts.index.astype(str), tissue_risk_counts.values, 
               color=['#2ca02c', '#ff7f0e', '#d62728'], alpha=0.8)
    axes[0].set_xlabel('Risk Category', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
    axes[0].set_title('Tissue Risk Stratification', fontsize=13, fontweight='bold')
    axes[0].grid(axis='y', alpha=0.3)
    
    blood_risk_counts = blood_pred['Risk_Category'].value_counts().sort_index()
    axes[1].bar(blood_risk_counts.index.astype(str), blood_risk_counts.values, 
               color=['#2ca02c', '#ff7f0e', '#d62728'], alpha=0.8)
    axes[1].set_xlabel('Risk Category', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
    axes[1].set_title('Blood Risk Stratification', fontsize=13, fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)
    
    plt.suptitle('Risk Score Stratification', fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig('results/figures/risk_stratification.png', bbox_inches='tight', dpi=300)
    plt.close()
    print("  ✅ Saved: results/figures/risk_stratification.png")
    sys.stdout.flush()
except Exception as e:
    print(f"  ❌ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.stdout.flush()

# 2. Gene Lists
print("\n[2] Creating gene lists...")
try:
    gene_lists = []
    
    if os.path.exists('results/multi_omics_tissue_signature_genes.tsv'):
        tissue_sig = pd.read_csv('results/multi_omics_tissue_signature_genes.tsv', sep='\t')
        print(f"  Found {len(tissue_sig)} tissue signature genes")
        for _, row in tissue_sig.iterrows():
            gene_lists.append({
                'Gene_ID': str(row.get('Gene', '')),
                'Gene_Symbol': str(row.get('Gene_Symbol', row.get('Gene', ''))),
                'Dataset': 'Tissue (GSE122897)',
                'Log2FC': row.get('Log2FC', 0),
                'P_adj': row.get('P_adj_BH', row.get('P_adj', row.get('P_value', 1))),
                'Direction': 'Up' if row.get('Log2FC', 0) > 0 else 'Down'
            })
    
    if os.path.exists('results/multi_omics_blood_signature_genes.tsv'):
        blood_sig = pd.read_csv('results/multi_omics_blood_signature_genes.tsv', sep='\t')
        print(f"  Found {len(blood_sig)} blood signature genes")
        for _, row in blood_sig.iterrows():
            gene_lists.append({
                'Gene_ID': str(row.get('Gene', row.get('Gene_ID', ''))),
                'Gene_Symbol': str(row.get('Gene_Symbol', row.get('Gene', row.get('Gene_ID', '')))),
                'Dataset': 'Blood (GSE159610)',
                'Log2FC': row.get('Log2FC', 0),
                'P_adj': row.get('P_adj_BH', row.get('P_adj', row.get('P_value', 1))),
                'Direction': 'Up' if row.get('Log2FC', 0) > 0 else 'Down'
            })
    
    if len(gene_lists) > 0:
        all_genes_df = pd.DataFrame(gene_lists)
        all_genes_df.to_csv('results/gene_lists/comprehensive_signature_genes_annotated.tsv', sep='\t', index=False)
        print("  ✅ Saved: results/gene_lists/comprehensive_signature_genes_annotated.tsv")
        
        try:
            all_genes_df.to_excel('results/gene_lists/comprehensive_signature_genes_annotated.xlsx', index=False)
            print("  ✅ Saved: results/gene_lists/comprehensive_signature_genes_annotated.xlsx")
        except Exception as e:
            print(f"  ⚠ Excel export skipped: {e}")
        
        tissue_genes_df = all_genes_df[all_genes_df['Dataset'] == 'Tissue (GSE122897)'].head(50)
        blood_genes_df = all_genes_df[all_genes_df['Dataset'] == 'Blood (GSE159610)'].head(50)
        
        if len(tissue_genes_df) > 0:
            tissue_genes_df.to_csv('results/gene_lists/tissue_top50_signature_genes.tsv', sep='\t', index=False)
            print(f"  ✅ Saved: results/gene_lists/tissue_top50_signature_genes.tsv ({len(tissue_genes_df)} genes)")
        
        if len(blood_genes_df) > 0:
            blood_genes_df.to_csv('results/gene_lists/blood_top50_signature_genes.tsv', sep='\t', index=False)
            print(f"  ✅ Saved: results/gene_lists/blood_top50_signature_genes.tsv ({len(blood_genes_df)} genes)")
    else:
        print("  ⚠ No signature genes found")
    sys.stdout.flush()
except Exception as e:
    print(f"  ❌ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.stdout.flush()

# 3. Feature Importance
print("\n[3] Creating feature importance plot...")
try:
    if os.path.exists('results/feature_importance.tsv'):
        feature_imp = pd.read_csv('results/feature_importance.tsv', sep='\t')
        print(f"  Loaded feature importance file with {len(feature_imp)} entries")
        print(f"  Columns: {list(feature_imp.columns)}")
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 14))
        axes = axes.flatten()
        
        plot_idx = 0
        for model_type in ['Tissue', 'Blood']:
            for model_name in ['Random Forest', 'Gradient Boosting']:
                model_data = feature_imp[(feature_imp['Model'].str.contains(model_type, case=False, na=False)) & 
                                         (feature_imp['Model'].str.contains(model_name, case=False, na=False))]
                if len(model_data) > 0 and plot_idx < len(axes):
                    top_features = model_data.head(20)
                    
                    # Find importance column
                    if 'Importance' in top_features.columns:
                        feat_col = 'Importance'
                    elif 'importance' in top_features.columns:
                        feat_col = 'importance'
                    elif 'Score' in top_features.columns:
                        feat_col = 'Score'
                    else:
                        feat_col = top_features.columns[1] if len(top_features.columns) > 1 else top_features.columns[0]
                    
                    # Find feature name column
                    if 'Feature' in top_features.columns:
                        feat_names = top_features['Feature'].values
                    elif 'feature' in top_features.columns:
                        feat_names = top_features['feature'].values
                    elif 'Gene' in top_features.columns:
                        feat_names = top_features['Gene'].values
                    else:
                        feat_names = top_features.iloc[:, 0].values
                    
                    axes[plot_idx].barh(range(len(top_features)), top_features[feat_col].values, 
                                       color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'][plot_idx % 4])
                    axes[plot_idx].set_yticks(range(len(top_features)))
                    axes[plot_idx].set_yticklabels(feat_names, fontsize=9)
                    axes[plot_idx].set_xlabel('Importance', fontsize=11, fontweight='bold')
                    axes[plot_idx].set_title(f'{model_type} {model_name} - Top Features', fontsize=12, fontweight='bold')
                    axes[plot_idx].invert_yaxis()
                    axes[plot_idx].grid(axis='x', alpha=0.3)
                    plot_idx += 1
                else:
                    print(f"  ⚠ No data found for {model_type} {model_name}")
        
        plt.suptitle('Feature Importance Rankings - Top 20 Features', fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()
        plt.savefig('results/figures/feature_importance_all_models.png', bbox_inches='tight', dpi=300)
        plt.close()
        
        feature_imp.to_csv('results/feature_importance_all_models.tsv', sep='\t', index=False)
        print("  ✅ Saved: results/figures/feature_importance_all_models.png")
        print("  ✅ Saved: results/feature_importance_all_models.tsv")
    else:
        print("  ⚠ Feature importance file not found")
    sys.stdout.flush()
except Exception as e:
    print(f"  ❌ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.stdout.flush()

# 4. Pathway Enrichment - Tissue
print("\n[4] Creating tissue pathway enrichment plot...")
try:
    if os.path.exists('results/pathway_enrichment_prioritized_fdr_0.05.tsv'):
        tissue_pathways = pd.read_csv('results/pathway_enrichment_prioritized_fdr_0.05.tsv', sep='\t')
        print(f"  Loaded tissue pathways file with {len(tissue_pathways)} pathways")
        print(f"  Columns: {list(tissue_pathways.columns)}")
        
        if len(tissue_pathways) > 0:
            top_paths = tissue_pathways.head(10)
            
            # Find p-value column
            p_col = None
            for col in ['P_adj', 'P_value', 'P_adj_BH', 'p_adj', 'p_value']:
                if col in top_paths.columns:
                    p_col = col
                    break
            
            # Find pathway name column
            path_name_col = None
            for col in ['Term', 'Pathway', 'pathway', 'term']:
                if col in top_paths.columns:
                    path_name_col = col
                    break
            
            if p_col and path_name_col:
                path_names = top_paths[path_name_col].values[:10]
                p_values = -np.log10(top_paths[p_col].values[:10] + 1e-10)
                
                fig, ax = plt.subplots(figsize=(12, 8))
                y_pos = range(len(path_names))
                ax.barh(y_pos, p_values, color='#1f77b4', alpha=0.8)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(path_names, fontsize=10)
                ax.set_xlabel('-log10(P-value)', fontsize=12, fontweight='bold')
                ax.set_title('Tissue Pathway Enrichment - Top 10', fontsize=14, fontweight='bold')
                ax.invert_yaxis()
                ax.grid(axis='x', alpha=0.3)
                plt.tight_layout()
                plt.savefig('results/figures/pathway_enrichment_tissue.png', bbox_inches='tight', dpi=300)
                plt.close()
                print("  ✅ Saved: results/figures/pathway_enrichment_tissue.png")
            else:
                print(f"  ⚠ Missing columns: p_col={p_col}, path_name_col={path_name_col}")
    else:
        print("  ⚠ Tissue pathways file not found")
    sys.stdout.flush()
except Exception as e:
    print(f"  ❌ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.stdout.flush()

# 5. Pathway Enrichment - Blood
print("\n[5] Creating blood pathway enrichment plot...")
try:
    if os.path.exists('results/pathway_enrichment_fdr_0.05.tsv'):
        blood_pathways = pd.read_csv('results/pathway_enrichment_fdr_0.05.tsv', sep='\t')
        print(f"  Loaded blood pathways file with {len(blood_pathways)} pathways")
        print(f"  Columns: {list(blood_pathways.columns)}")
        
        if len(blood_pathways) > 0:
            top_paths = blood_pathways.head(10)
            
            # Find p-value column
            p_col = None
            for col in ['P_adj', 'P_value', 'P_adj_BH', 'p_adj', 'p_value']:
                if col in top_paths.columns:
                    p_col = col
                    break
            
            # Find pathway name column
            path_name_col = None
            for col in ['Term', 'Pathway', 'pathway', 'term']:
                if col in top_paths.columns:
                    path_name_col = col
                    break
            
            if p_col and path_name_col:
                path_names = top_paths[path_name_col].values[:10]
                p_values = -np.log10(top_paths[p_col].values[:10] + 1e-10)
                
                fig, ax = plt.subplots(figsize=(12, 6))
                y_pos = range(len(path_names))
                ax.barh(y_pos, p_values, color='#ff7f0e', alpha=0.8)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(path_names, fontsize=10)
                ax.set_xlabel('-log10(P-value)', fontsize=12, fontweight='bold')
                ax.set_title('Blood Pathway Enrichment', fontsize=14, fontweight='bold')
                ax.invert_yaxis()
                ax.grid(axis='x', alpha=0.3)
                plt.tight_layout()
                plt.savefig('results/figures/pathway_enrichment_blood.png', bbox_inches='tight', dpi=300)
                plt.close()
                print("  ✅ Saved: results/figures/pathway_enrichment_blood.png")
            else:
                print(f"  ⚠ Missing columns: p_col={p_col}, path_name_col={path_name_col}")
    else:
        print("  ⚠ Blood pathways file not found")
    sys.stdout.flush()
except Exception as e:
    print(f"  ❌ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.stdout.flush()

# 6. Common Pathways
print("\n[6] Identifying common pathways...")
try:
    if os.path.exists('results/pathway_enrichment_prioritized_fdr_0.05.tsv') and os.path.exists('results/pathway_enrichment_fdr_0.05.tsv'):
        tissue_pathways = pd.read_csv('results/pathway_enrichment_prioritized_fdr_0.05.tsv', sep='\t')
        blood_pathways = pd.read_csv('results/pathway_enrichment_fdr_0.05.tsv', sep='\t')
        
        tissue_col = 'Term' if 'Term' in tissue_pathways.columns else ('Pathway' if 'Pathway' in tissue_pathways.columns else None)
        blood_col = 'Term' if 'Term' in blood_pathways.columns else ('Pathway' if 'Pathway' in blood_pathways.columns else None)
        
        if tissue_col and blood_col:
            tissue_path_names = set(tissue_pathways[tissue_col].str.lower().dropna())
            blood_path_names = set(blood_pathways[blood_col].str.lower().dropna())
            common_pathways = tissue_path_names.intersection(blood_path_names)
            
            if len(common_pathways) > 0:
                # Map back to original names
                common_paths_list = []
                for path_lower in common_pathways:
                    orig = tissue_pathways[tissue_pathways[tissue_col].str.lower() == path_lower][tissue_col].values
                    if len(orig) > 0:
                        common_paths_list.append(orig[0])
                
                if len(common_paths_list) > 0:
                    common_paths_df = pd.DataFrame({
                        'Pathway': common_paths_list,
                        'Source': 'Both_Datasets'
                    })
                    common_paths_df.to_csv('results/common_enriched_pathways.tsv', sep='\t', index=False)
                    print(f"  ✅ Found {len(common_pathways)} common pathways")
                    print("  ✅ Saved: results/common_enriched_pathways.tsv")
            else:
                print("  ℹ No common pathways found")
        else:
            print(f"  ⚠ Missing pathway name columns")
    else:
        print("  ⚠ Pathway files not found")
    sys.stdout.flush()
except Exception as e:
    print(f"  ❌ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.stdout.flush()

print("\n" + "=" * 80)
print("✅ ALL REMAINING OUTPUTS GENERATION COMPLETE")
print("=" * 80)
sys.stdout.flush()

