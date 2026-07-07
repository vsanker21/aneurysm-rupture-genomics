"""Force generation of all remaining outputs"""
import sys
import os
import traceback

if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

try:
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    print("=" * 80)
    print("FORCING GENERATION OF ALL OUTPUTS")
    print("=" * 80)
    sys.stdout.flush()
    
    os.makedirs("results/figures", exist_ok=True)
    os.makedirs("results/gene_lists", exist_ok=True)
    os.makedirs("results/risk_scores", exist_ok=True)
    
    # 1. Risk Stratification Plot
    print("\n[1] Risk Stratification Plot...")
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
        axes[0].bar([str(x) for x in tissue_risk_counts.index], tissue_risk_counts.values, 
                   color=['#2ca02c', '#ff7f0e', '#d62728'], alpha=0.8)
        axes[0].set_xlabel('Risk Category')
        axes[0].set_ylabel('Number of Samples')
        axes[0].set_title('Tissue Risk Stratification')
        
        blood_risk_counts = blood_pred['Risk_Category'].value_counts().sort_index()
        axes[1].bar([str(x) for x in blood_risk_counts.index], blood_risk_counts.values, 
                   color=['#2ca02c', '#ff7f0e', '#d62728'], alpha=0.8)
        axes[1].set_xlabel('Risk Category')
        axes[1].set_ylabel('Number of Samples')
        axes[1].set_title('Blood Risk Stratification')
        
        plt.suptitle('Risk Score Stratification')
        plt.tight_layout()
        plt.savefig('results/figures/risk_stratification.png', bbox_inches='tight', dpi=300)
        plt.close()
        print("  ✅ Created: risk_stratification.png")
        sys.stdout.flush()
    except Exception as e:
        print(f"  ❌ Error: {e}")
        traceback.print_exc()
        sys.stdout.flush()
    
    # 2. Gene Lists
    print("\n[2] Gene Lists...")
    try:
        gene_lists = []
        
        if os.path.exists('results/multi_omics_tissue_signature_genes.tsv'):
            tissue_sig = pd.read_csv('results/multi_omics_tissue_signature_genes.tsv', sep='\t')
            print(f"  Tissue genes: {len(tissue_sig)}")
            for _, row in tissue_sig.iterrows():
                gene_lists.append({
                    'Gene_ID': str(row.get('Gene', '')),
                    'Gene_Symbol': str(row.get('Gene_Symbol', '')),
                    'Dataset': 'Tissue (GSE122897)',
                    'Log2FC': float(row.get('Log2FC', 0)),
                    'P_adj': float(row.get('P_adj_BH', row.get('P_adj', row.get('P_value', 1)))),
                    'Direction': 'Up' if float(row.get('Log2FC', 0)) > 0 else 'Down'
                })
        
        if os.path.exists('results/multi_omics_blood_signature_genes.tsv'):
            blood_sig = pd.read_csv('results/multi_omics_blood_signature_genes.tsv', sep='\t')
            print(f"  Blood genes: {len(blood_sig)}")
            for _, row in blood_sig.iterrows():
                gene_lists.append({
                    'Gene_ID': str(row.get('Gene', row.get('Gene_ID', ''))),
                    'Gene_Symbol': str(row.get('Gene_Symbol', '')),
                    'Dataset': 'Blood (GSE159610)',
                    'Log2FC': float(row.get('Log2FC', 0)),
                    'P_adj': float(row.get('P_adj_BH', row.get('P_adj', row.get('P_value', 1)))),
                    'Direction': 'Up' if float(row.get('Log2FC', 0)) > 0 else 'Down'
                })
        
        if len(gene_lists) > 0:
            all_genes_df = pd.DataFrame(gene_lists)
            all_genes_df.to_csv('results/gene_lists/comprehensive_signature_genes_annotated.tsv', sep='\t', index=False)
            print(f"  ✅ Created: comprehensive_signature_genes_annotated.tsv ({len(all_genes_df)} genes)")
            
            tissue_genes_df = all_genes_df[all_genes_df['Dataset'] == 'Tissue (GSE122897)'].head(50)
            blood_genes_df = all_genes_df[all_genes_df['Dataset'] == 'Blood (GSE159610)'].head(50)
            
            if len(tissue_genes_df) > 0:
                tissue_genes_df.to_csv('results/gene_lists/tissue_top50_signature_genes.tsv', sep='\t', index=False)
                print(f"  ✅ Created: tissue_top50_signature_genes.tsv ({len(tissue_genes_df)} genes)")
            
            if len(blood_genes_df) > 0:
                blood_genes_df.to_csv('results/gene_lists/blood_top50_signature_genes.tsv', sep='\t', index=False)
                print(f"  ✅ Created: blood_top50_signature_genes.tsv ({len(blood_genes_df)} genes)")
        sys.stdout.flush()
    except Exception as e:
        print(f"  ❌ Error: {e}")
        traceback.print_exc()
        sys.stdout.flush()
    
    # 3. Feature Importance
    print("\n[3] Feature Importance Plot...")
    try:
        if os.path.exists('results/feature_importance.tsv'):
            feature_imp = pd.read_csv('results/feature_importance.tsv', sep='\t')
            print(f"  Loaded {len(feature_imp)} features")
            
            fig, axes = plt.subplots(2, 2, figsize=(16, 14))
            axes = axes.flatten()
            
            # Create plots for each model type
            models = [
                ('Tissue', 'Random Forest'),
                ('Tissue', 'Gradient Boosting'),
                ('Blood', 'Random Forest'),
                ('Blood', 'Gradient Boosting')
            ]
            
            for idx, (model_type, model_name) in enumerate(models):
                model_data = feature_imp[
                    (feature_imp['Model'].str.contains(model_type, case=False, na=False)) & 
                    (feature_imp['Model'].str.contains(model_name, case=False, na=False))
                ]
                
                if len(model_data) > 0:
                    top_features = model_data.head(20).copy()
                    importance_col = 'Importance' if 'Importance' in top_features.columns else top_features.columns[1]
                    feature_col = 'Feature' if 'Feature' in top_features.columns else top_features.columns[0]
                    
                    top_features = top_features.sort_values(importance_col, ascending=True).tail(20)
                    
                    axes[idx].barh(range(len(top_features)), top_features[importance_col].values)
                    axes[idx].set_yticks(range(len(top_features)))
                    axes[idx].set_yticklabels(top_features[feature_col].values, fontsize=9)
                    axes[idx].set_title(f'{model_type} {model_name}')
                    axes[idx].invert_yaxis()
                else:
                    axes[idx].text(0.5, 0.5, f'No data for\n{model_type} {model_name}', 
                                 ha='center', va='center', transform=axes[idx].transAxes)
            
            plt.suptitle('Feature Importance Rankings')
            plt.tight_layout()
            plt.savefig('results/figures/feature_importance_all_models.png', bbox_inches='tight', dpi=300)
            plt.close()
            print("  ✅ Created: feature_importance_all_models.png")
            
            feature_imp.to_csv('results/feature_importance_all_models.tsv', sep='\t', index=False)
            print("  ✅ Created: feature_importance_all_models.tsv")
        else:
            print("  ⚠ File not found: results/feature_importance.tsv")
        sys.stdout.flush()
    except Exception as e:
        print(f"  ❌ Error: {e}")
        traceback.print_exc()
        sys.stdout.flush()
    
    # 4. Tissue Pathway Enrichment
    print("\n[4] Tissue Pathway Enrichment...")
    try:
        if os.path.exists('results/pathway_enrichment_prioritized_fdr_0.05.tsv'):
            tissue_pathways = pd.read_csv('results/pathway_enrichment_prioritized_fdr_0.05.tsv', sep='\t')
            print(f"  Loaded {len(tissue_pathways)} pathways")
            
            top_paths = tissue_pathways.head(10).copy()
            
            p_col = None
            for col in ['P_adj', 'P_value', 'P_adj_BH']:
                if col in top_paths.columns:
                    p_col = col
                    break
            
            path_col = None
            for col in ['Term', 'Pathway']:
                if col in top_paths.columns:
                    path_col = col
                    break
            
            if p_col and path_col:
                fig, ax = plt.subplots(figsize=(12, 8))
                path_names = top_paths[path_col].values
                p_values = -np.log10(top_paths[p_col].values + 1e-10)
                
                y_pos = range(len(path_names))
                ax.barh(y_pos, p_values, color='#1f77b4', alpha=0.8)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(path_names, fontsize=10)
                ax.set_xlabel('-log10(P-value)')
                ax.set_title('Tissue Pathway Enrichment - Top 10')
                ax.invert_yaxis()
                plt.tight_layout()
                plt.savefig('results/figures/pathway_enrichment_tissue.png', bbox_inches='tight', dpi=300)
                plt.close()
                print("  ✅ Created: pathway_enrichment_tissue.png")
            else:
                print(f"  ⚠ Missing columns: p_col={p_col}, path_col={path_col}")
        else:
            print("  ⚠ File not found")
        sys.stdout.flush()
    except Exception as e:
        print(f"  ❌ Error: {e}")
        traceback.print_exc()
        sys.stdout.flush()
    
    # 5. Blood Pathway Enrichment
    print("\n[5] Blood Pathway Enrichment...")
    try:
        if os.path.exists('results/pathway_enrichment_fdr_0.05.tsv'):
            blood_pathways = pd.read_csv('results/pathway_enrichment_fdr_0.05.tsv', sep='\t')
            print(f"  Loaded {len(blood_pathways)} pathways")
            
            top_paths = blood_pathways.head(10).copy()
            
            p_col = None
            for col in ['P_adj', 'P_value', 'P_adj_BH']:
                if col in top_paths.columns:
                    p_col = col
                    break
            
            path_col = None
            for col in ['Term', 'Pathway']:
                if col in top_paths.columns:
                    path_col = col
                    break
            
            if p_col and path_col:
                fig, ax = plt.subplots(figsize=(12, 6))
                path_names = top_paths[path_col].values
                p_values = -np.log10(top_paths[p_col].values + 1e-10)
                
                y_pos = range(len(path_names))
                ax.barh(y_pos, p_values, color='#ff7f0e', alpha=0.8)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(path_names, fontsize=10)
                ax.set_xlabel('-log10(P-value)')
                ax.set_title('Blood Pathway Enrichment')
                ax.invert_yaxis()
                plt.tight_layout()
                plt.savefig('results/figures/pathway_enrichment_blood.png', bbox_inches='tight', dpi=300)
                plt.close()
                print("  ✅ Created: pathway_enrichment_blood.png")
            else:
                print(f"  ⚠ Missing columns")
        else:
            print("  ⚠ File not found")
        sys.stdout.flush()
    except Exception as e:
        print(f"  ❌ Error: {e}")
        traceback.print_exc()
        sys.stdout.flush()
    
    print("\n" + "=" * 80)
    print("✅ ALL OUTPUTS GENERATION COMPLETE")
    print("=" * 80)
    sys.stdout.flush()
    
except Exception as e:
    print(f"FATAL ERROR: {e}")
    traceback.print_exc()
    sys.exit(1)

