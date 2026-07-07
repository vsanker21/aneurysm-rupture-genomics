"""
Generate Complete Ultra HD Visualizations with Gene Names
All figures use gene symbols/names instead of IDs, and only include available data
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import json
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality style
try:
    plt.style.use('seaborn-v0_8-paper')
except:
    try:
        plt.style.use('seaborn-paper')
    except:
        plt.style.use('default')
sns.set_palette("husl")
# Large, visible fonts for all figures (main and supplementary)
matplotlib.rcParams.update({
    'font.size': 18,
    'axes.titlesize': 22,
    'axes.labelsize': 20,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'figure.titlesize': 24,
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Liberation Sans', 'Arial'],
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.5,
    'grid.linewidth': 0.8,
    'lines.linewidth': 2.5,
    'patch.linewidth': 1.5,
    'axes.spines.left': True,
    'axes.spines.bottom': True,
    'axes.spines.top': False,
    'axes.spines.right': False
})
warnings.filterwarnings('ignore', category=UserWarning, message='.*Glyph.*')

# Create output directory
output_dir = Path("results/figures/ultra_hd")
output_dir.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("GENERATING COMPLETE ULTRA HD VISUALIZATIONS WITH GENE NAMES")
print("=" * 80)

# ============================================================================
# GENE MAPPING FUNCTIONS
# ============================================================================

def load_gene_mapping():
    """Load gene ID to symbol/name mapping"""
    mapping = {}
    
    # Try to load annotated DE results first
    annotated_file = "results/gse122897_top_genes_annotated.tsv"
    if os.path.exists(annotated_file):
        try:
            annotated_df = pd.read_csv(annotated_file, sep='\t')
            for idx, row in annotated_df.iterrows():
                gene_id = str(row['Gene']).split('.')[0]  # Remove .0
                symbol = row.get('Gene_Symbol', '')
                name = row.get('Gene_Name', '')
                if pd.notna(symbol) and str(symbol).strip():
                    mapping[gene_id] = {'symbol': str(symbol), 'name': str(name) if pd.notna(name) else ''}
        except:
            pass
    
    # Also try integrated file
    integrated_file = "results/gse122897_integrated_de.tsv"
    if os.path.exists(integrated_file):
        try:
            integrated_df = pd.read_csv(integrated_file, sep='\t')
            if 'Gene' in integrated_df.columns and 'Gene_Symbol' in integrated_df.columns:
                for idx, row in integrated_df.iterrows():
                    gene_id = str(row['Gene']).split('.')[0]
                    symbol = row.get('Gene_Symbol', '')
                    name = row.get('Gene_Name', '')
                    if pd.notna(symbol) and str(symbol).strip() and gene_id not in mapping:
                        mapping[gene_id] = {'symbol': str(symbol), 'name': str(name) if pd.notna(name) else ''}
        except:
            pass
    
    return mapping

# Load mapping once
GENE_MAPPING = load_gene_mapping()
print(f"\nLoaded gene mapping for {len(GENE_MAPPING)} genes")

def get_gene_name(gene_id, prefer_symbol=True):
    """Get gene symbol or name from ID"""
    gene_id_str = str(gene_id).split('.')[0]
    if gene_id_str in GENE_MAPPING:
        info = GENE_MAPPING[gene_id_str]
        if prefer_symbol and info.get('symbol'):
            return info['symbol']
        elif info.get('name'):
            return info['name']
    # If no mapping found, try to return as-is if it looks like a symbol
    if isinstance(gene_id, str) and not gene_id.isdigit():
        return gene_id
    return f"Gene_{gene_id_str}"

# ============================================================================
# 1. VOLCANO PLOT WITH GENE NAMES
# ============================================================================

def plot_volcano_ultra_hd():
    """Create ultra HD volcano plot with gene names"""
    print("\n[1/15] Creating Ultra HD Volcano Plot with Gene Names...")
    
    de_file = "results/gse122897_differential_expression.tsv"
    if not os.path.exists(de_file):
        print(f"  [SKIP] {de_file} not found")
        return
    
    de_df = pd.read_csv(de_file, sep='\t')
    
    # Try to get gene symbols if available
    annotated_file = "results/gse122897_top_genes_annotated.tsv"
    if os.path.exists(annotated_file):
        try:
            annotated_df = pd.read_csv(annotated_file, sep='\t')
            de_df = de_df.merge(
                annotated_df[['Gene', 'Gene_Symbol', 'Gene_Name']].drop_duplicates('Gene'),
                on='Gene', how='left'
            )
            de_df['Gene_Display'] = de_df.apply(
                lambda row: row['Gene_Symbol'] if pd.notna(row.get('Gene_Symbol')) 
                else get_gene_name(row['Gene']), axis=1
            )
        except:
            de_df['Gene_Display'] = de_df['Gene'].apply(get_gene_name)
    else:
        de_df['Gene_Display'] = de_df['Gene'].apply(get_gene_name)
    
    # Prepare data
    de_df['-log10_padj'] = -np.log10(de_df['P_adj_BH'].clip(lower=1e-300))
    de_df['abs_log2fc'] = np.abs(de_df['Log2FC'])
    
    # Categorize genes
    de_df['Significance'] = 'Not Significant'
    de_df.loc[(de_df['P_adj_BH'] < 0.05) & (de_df['Log2FC'] > 1), 'Significance'] = 'Upregulated (P<0.05, FC>2)'
    de_df.loc[(de_df['P_adj_BH'] < 0.05) & (de_df['Log2FC'] < -1), 'Significance'] = 'Downregulated (P<0.05, FC>2)'
    de_df.loc[(de_df['P_adj_BH'] < 0.05) & (de_df['abs_log2fc'] <= 1), 'Significance'] = 'Significant (P<0.05)'
    
    fig, ax = plt.subplots(figsize=(20, 14))
    
    # Plot points: darker, more visible colors (significant = red/blue, non-significant = dark grey)
    colors = {'Not Significant': '#4d4d4d',           # dark grey, clearly visible
              'Significant (P<0.05)': '#2d2d2d',       # darker grey for significant but |FC|<=1
              'Upregulated (P<0.05, FC>2)': '#b71c1c', # dark red, highly visible
              'Downregulated (P<0.05, FC>2)': '#1565c0'} # dark blue, highly visible
    
    for sig_type in de_df['Significance'].unique():
        subset = de_df[de_df['Significance'] == sig_type]
        ax.scatter(subset['Log2FC'], subset['-log10_padj'], 
                  c=colors.get(sig_type, '#4d4d4d'), s=120, alpha=0.85, 
                  label=sig_type, edgecolors='none')
    
    # Highlight top genes with names
    top_genes = de_df.nsmallest(30, 'P_adj_BH')
    for idx, row in top_genes.iterrows():
        if row['P_adj_BH'] < 0.1:
            gene_label = row['Gene_Display']
            # Truncate long names
            if len(gene_label) > 15:
                gene_label = gene_label[:12] + '...'
            ax.annotate(gene_label, 
                       xy=(row['Log2FC'], row['-log10_padj']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=14, alpha=0.9, 
                       bbox=dict(boxstyle='round,pad=0.3', 
                               facecolor='white', alpha=0.8, edgecolor='gray', linewidth=0.5))
    
    # Add thresholds
    ax.axhline(y=-np.log10(0.05), color='gray', linestyle='--', 
              linewidth=1.5, alpha=0.7, label='FDR = 0.05')
    ax.axvline(x=1, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axvline(x=-1, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    
    ax.set_xlabel('Log2 Fold Change (Ruptured vs Unruptured)', 
                 fontweight='bold', fontsize=16)
    ax.set_ylabel('-Log10 (FDR-adjusted P-value)', 
                 fontweight='bold', fontsize=16)
    ax.set_title('Volcano Plot: Differential Expression Analysis\n' + 
                'Tissue RNA-seq (GSE122897)', fontweight='bold', fontsize=18)
    ax.legend(loc='best', frameon=True, fancybox=True, shadow=True, fontsize=12)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', labelsize=12)
    
    plt.tight_layout()
    plt.savefig(output_dir / "volcano_plot_ultra_hd.png", 
               dpi=600, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  [OK] Saved: {output_dir}/volcano_plot_ultra_hd.png")

# ============================================================================
# 2. TOP GENES HEATMAP WITH GENE NAMES
# ============================================================================

def plot_top_genes_heatmap_ultra_hd():
    """Create ultra HD heatmap with gene names"""
    print("\n[2/15] Creating Ultra HD Top Genes Heatmap with Gene Names...")
    
    de_file = "results/gse122897_differential_expression.tsv"
    expr_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    metadata_file = "data/sc_rnaseq/GSE122897/sample_metadata_with_clinical.tsv"
    
    if not all(os.path.exists(f) for f in [de_file, expr_file]):
        print(f"  [SKIP] Required files not found")
        return
    
    de_df = pd.read_csv(de_file, sep='\t')
    
    # Get gene symbols
    annotated_file = "results/gse122897_top_genes_annotated.tsv"
    if os.path.exists(annotated_file):
        try:
            annotated_df = pd.read_csv(annotated_file, sep='\t')
            de_df = de_df.merge(
                annotated_df[['Gene', 'Gene_Symbol', 'Gene_Name']].drop_duplicates('Gene'),
                on='Gene', how='left'
            )
            de_df['Gene_Display'] = de_df.apply(
                lambda row: row['Gene_Symbol'] if pd.notna(row.get('Gene_Symbol')) 
                else get_gene_name(row['Gene']), axis=1
            )
        except:
            de_df['Gene_Display'] = de_df['Gene'].apply(get_gene_name)
    else:
        de_df['Gene_Display'] = de_df['Gene'].apply(get_gene_name)
    
    # Get top 50 genes
    top_genes = de_df.nsmallest(50, 'P_adj_BH')
    gene_ids = top_genes['Gene'].astype(str).tolist()
    gene_labels = top_genes['Gene_Display'].tolist()
    
    # Load expression data
    expr_df = pd.read_csv(expr_file, sep='\t', compression='gzip', index_col=0)
    expr_df.index = expr_df.index.astype(str)
    
    # Map gene IDs
    top_expr = expr_df.loc[expr_df.index.intersection(gene_ids)]
    
    if len(top_expr) < 10:
        print(f"  [SKIP] Too few genes matched ({len(top_expr)})")
        return
    
    # Create mapping for row labels
    id_to_label = dict(zip(gene_ids, gene_labels))
    row_labels = [id_to_label.get(str(idx), str(idx)) for idx in top_expr.index]
    
    # Load metadata and sort samples
    if os.path.exists(metadata_file):
        metadata = pd.read_csv(metadata_file, sep='\t')
        if 'Sample' in metadata.columns and 'Group' in metadata.columns:
            samples = top_expr.columns.tolist()
            metadata_aligned = metadata[metadata['Sample'].isin(samples)].set_index('Sample')
            group_order = ['Unruptured', 'Ruptured']
            metadata_aligned['Group_Order'] = metadata_aligned['Group'].map(
                {g: i for i, g in enumerate(group_order)})
            sample_order = metadata_aligned.sort_values('Group_Order').index.tolist()
            available_samples = [s for s in sample_order if s in top_expr.columns]
            top_expr = top_expr[available_samples]
            row_labels = [id_to_label.get(str(idx), str(idx)) for idx in top_expr.index]
    
    # Log transform and z-score normalize
    top_expr_log = np.log2(top_expr + 1)
    top_expr_z = top_expr_log.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(18, 16))
    
    # Create heatmap with gene names
    sns.heatmap(top_expr_z, cmap='RdBu_r', center=0, 
               vmin=-2, vmax=2, cbar_kws={'label': 'Z-score (Log2 TPM)'},
               xticklabels=False, yticklabels=row_labels, 
               linewidths=0.1, linecolor='white', ax=ax)
    
    ax.set_xlabel('Samples (Unruptured | Ruptured)', fontweight='bold', fontsize=20)
    ax.set_ylabel('Top 50 Differentially Expressed Genes', fontweight='bold', fontsize=20)
    ax.set_title('Heatmap: Top 50 Differentially Expressed Genes\n' + 
                'Z-score Normalized Expression (Log2 TPM + 1)', 
                fontweight='bold', fontsize=24)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=14, rotation=0)
    ax.tick_params(axis='both', labelsize=16)
    
    plt.tight_layout()
    plt.savefig(output_dir / "top_genes_heatmap_ultra_hd.png", 
               dpi=600, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  [OK] Saved: {output_dir}/top_genes_heatmap_ultra_hd.png")

# ============================================================================
# 3. MA PLOT (simplified, uses available data only)
# ============================================================================

def plot_ma_plot_ultra_hd():
    """Create ultra HD MA plot with gene names"""
    print("\n[3/15] Creating Ultra HD MA Plot with Gene Names...")
    
    de_file = "results/gse122897_differential_expression.tsv"
    if not os.path.exists(de_file):
        print(f"  [SKIP] {de_file} not found")
        return
    
    de_df = pd.read_csv(de_file, sep='\t')
    
    # Get gene symbols
    annotated_file = "results/gse122897_top_genes_annotated.tsv"
    if os.path.exists(annotated_file):
        try:
            annotated_df = pd.read_csv(annotated_file, sep='\t')
            de_df = de_df.merge(
                annotated_df[['Gene', 'Gene_Symbol']].drop_duplicates('Gene'),
                on='Gene', how='left'
            )
            de_df['Gene_Display'] = de_df.apply(
                lambda row: row['Gene_Symbol'] if pd.notna(row.get('Gene_Symbol')) 
                else get_gene_name(row['Gene']), axis=1
            )
        except:
            de_df['Gene_Display'] = de_df['Gene'].apply(get_gene_name)
    else:
        de_df['Gene_Display'] = de_df['Gene'].apply(get_gene_name)
    
    # Calculate A and M
    de_df['A'] = np.log2((de_df['Mean_Ruptured'] + de_df['Mean_Unruptured'] + 1e-6) / 2)
    de_df['M'] = de_df['Log2FC']
    de_df['-log10_padj'] = -np.log10(de_df['P_adj_BH'].clip(lower=1e-300))
    
    # Categorize
    de_df['Significance'] = 'Not Significant'
    de_df.loc[de_df['P_adj_BH'] < 0.05, 'Significance'] = 'Significant (FDR<0.05)'
    
    fig, ax = plt.subplots(figsize=(20, 14))
    
    # Plot
    sig_df = de_df[de_df['P_adj_BH'] < 0.05]
    nsig_df = de_df[de_df['P_adj_BH'] >= 0.05]
    
    ax.scatter(nsig_df['A'], nsig_df['M'], c='#999999', 
              s=60, alpha=0.5, label='Not Significant', edgecolors='none')
    ax.scatter(sig_df['A'], sig_df['M'], c='#d62728', 
              s=80, alpha=0.7, label='Significant (FDR<0.05)', edgecolors='none')
    
    # Label top genes
    top_genes = de_df.nsmallest(20, 'P_adj_BH')
    for idx, row in top_genes.iterrows():
        if row['P_adj_BH'] < 0.1:
            gene_label = row['Gene_Display']
            if len(gene_label) > 12:
                gene_label = gene_label[:10] + '...'
            ax.annotate(gene_label, 
                       xy=(row['A'], row['M']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=14, alpha=0.9,
                       bbox=dict(boxstyle='round,pad=0.3', 
                               facecolor='white', alpha=0.8, edgecolor='gray', linewidth=0.5))
    
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1.5, alpha=0.5)
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axhline(y=-1, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    
    ax.set_xlabel('Average Expression (Log2 TPM)', fontweight='bold', fontsize=20)
    ax.set_ylabel('Log2 Fold Change', fontweight='bold', fontsize=20)
    ax.set_title('MA Plot: Differential Expression\n' + 
                'Mean Expression vs Fold Change', fontweight='bold', fontsize=24)
    ax.legend(loc='best', frameon=True, fancybox=True, shadow=True, fontsize=16)
    ax.tick_params(axis='both', labelsize=16)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(output_dir / "ma_plot_ultra_hd.png", 
               dpi=600, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  [OK] Saved: {output_dir}/ma_plot_ultra_hd.png")

# Continue with other plots... (using existing code but ensuring gene names are used)
# Due to length, I'll create a function that imports and runs the rest

def run_all_visualizations_with_gene_names():
    """Run all visualizations with gene names"""
    
    # Core DE plots
    plot_volcano_ultra_hd()
    plot_ma_plot_ultra_hd()
    plot_top_genes_heatmap_ultra_hd()
    
    # Import and run other visualization functions from the original script
    # (These don't need gene name changes as they use aggregated features)
    try:
        from generate_ultra_hd_visualizations import (
            plot_model_performance_ultra_hd, plot_feature_importance_ultra_hd,
            plot_risk_stratification_ultra_hd, plot_pathway_enrichment_ultra_hd,
            plot_blood_signature_ultra_hd, plot_ec_dysfunction_ultra_hd,
            plot_integration_summary_ultra_hd, plot_gene_prioritization_ultra_hd,
            plot_model_comparison_ultra_hd, plot_feature_distributions_ultra_hd,
            plot_roc_curve_ultra_hd, create_master_figure
        )
        
        plot_model_performance_ultra_hd()
        plot_feature_importance_ultra_hd()
        plot_risk_stratification_ultra_hd()
        plot_pathway_enrichment_ultra_hd()
        plot_blood_signature_ultra_hd()
        plot_ec_dysfunction_ultra_hd()
        plot_integration_summary_ultra_hd()
        plot_gene_prioritization_ultra_hd()
        plot_model_comparison_ultra_hd()
        plot_feature_distributions_ultra_hd()
        plot_roc_curve_ultra_hd()
        create_master_figure()
        
    except Exception as e:
        print(f"[WARNING] Some visualizations failed: {e}")
        print("  Continuing with available figures...")

if __name__ == "__main__":
    run_all_visualizations_with_gene_names()
    print("\n" + "=" * 80)
    print("VISUALIZATION GENERATION COMPLETE!")
    print("=" * 80)

