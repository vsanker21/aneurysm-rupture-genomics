"""Generate Pathway Enrichment Plots - Individual Script"""
import sys
import os
import traceback

if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

print("=" * 80)
print("GENERATING PATHWAY ENRICHMENT PLOTS")
print("=" * 80)
sys.stdout.flush()

try:
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    print("✓ Imports successful")
    sys.stdout.flush()
    
    os.makedirs("results/figures", exist_ok=True)
    
    # Tissue pathways
    print("\n[1] Creating tissue pathway plot...")
    sys.stdout.flush()
    if os.path.exists('results/pathway_enrichment_prioritized_fdr_0.05.tsv'):
        tissue_pathways = pd.read_csv('results/pathway_enrichment_prioritized_fdr_0.05.tsv', sep='\t')
        print(f"✓ Loaded {len(tissue_pathways)} pathways")
        print(f"  Columns: {list(tissue_pathways.columns)}")
        sys.stdout.flush()
        
        if len(tissue_pathways) > 0:
            top_paths = tissue_pathways.head(10).copy()
            
            # Find columns
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
                path_names = top_paths[path_col].values[:10]
                p_values = -np.log10(top_paths[p_col].values[:10] + 1e-10)
                
                y_pos = range(len(path_names))
                ax.barh(y_pos, p_values, color='#1f77b4', alpha=0.8)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(path_names, fontsize=10)
                ax.set_xlabel('-log10(P-value)', fontsize=12, fontweight='bold')
                ax.set_title('Tissue Pathway Enrichment - Top 10', fontsize=14, fontweight='bold')
                ax.invert_yaxis()
                ax.grid(axis='x', alpha=0.3)
                plt.tight_layout()
                
                output_file = 'results/figures/pathway_enrichment_tissue.png'
                plt.savefig(output_file, bbox_inches='tight', dpi=300)
                plt.close()
                
                if os.path.exists(output_file):
                    file_size = os.path.getsize(output_file)
                    print(f"✓ SUCCESS: Created {output_file} ({file_size/1024:.2f} KB)")
                else:
                    print(f"✗ ERROR: File not created: {output_file}")
            else:
                print(f"✗ ERROR: Missing columns - p_col={p_col}, path_col={path_col}")
    else:
        print("✗ ERROR: File not found: results/pathway_enrichment_prioritized_fdr_0.05.tsv")
    
    # Blood pathways
    print("\n[2] Creating blood pathway plot...")
    sys.stdout.flush()
    if os.path.exists('results/pathway_enrichment_fdr_0.05.tsv'):
        blood_pathways = pd.read_csv('results/pathway_enrichment_fdr_0.05.tsv', sep='\t')
        print(f"✓ Loaded {len(blood_pathways)} pathways")
        print(f"  Columns: {list(blood_pathways.columns)}")
        sys.stdout.flush()
        
        if len(blood_pathways) > 0:
            top_paths = blood_pathways.head(10).copy()
            
            # Find columns
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
                path_names = top_paths[path_col].values[:10]
                p_values = -np.log10(top_paths[p_col].values[:10] + 1e-10)
                
                y_pos = range(len(path_names))
                ax.barh(y_pos, p_values, color='#ff7f0e', alpha=0.8)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(path_names, fontsize=10)
                ax.set_xlabel('-log10(P-value)', fontsize=12, fontweight='bold')
                ax.set_title('Blood Pathway Enrichment', fontsize=14, fontweight='bold')
                ax.invert_yaxis()
                ax.grid(axis='x', alpha=0.3)
                plt.tight_layout()
                
                output_file = 'results/figures/pathway_enrichment_blood.png'
                plt.savefig(output_file, bbox_inches='tight', dpi=300)
                plt.close()
                
                if os.path.exists(output_file):
                    file_size = os.path.getsize(output_file)
                    print(f"✓ SUCCESS: Created {output_file} ({file_size/1024:.2f} KB)")
                else:
                    print(f"✗ ERROR: File not created: {output_file}")
            else:
                print(f"✗ ERROR: Missing columns - p_col={p_col}, path_col={path_col}")
    else:
        print("✗ ERROR: File not found: results/pathway_enrichment_fdr_0.05.tsv")
    
    # Common pathways
    print("\n[3] Identifying common pathways...")
    sys.stdout.flush()
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
                        print(f"✓ SUCCESS: Found {len(common_pathways)} common pathways")
                        print("✓ Saved: results/common_enriched_pathways.tsv")
    except Exception as e:
        print(f"⚠ Warning: {e}")
    
    sys.stdout.flush()
    
except Exception as e:
    print(f"\n✗ ERROR: {e}")
    print("\nFull traceback:")
    traceback.print_exc()
    sys.stdout.flush()
    sys.exit(1)

print("\n" + "=" * 80)
print("COMPLETE")
print("=" * 80)

