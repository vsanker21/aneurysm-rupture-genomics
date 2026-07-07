"""
Create visualizations for multi-omics integration results
Including volcano plots, Venn diagrams, and summary plots
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Circle
import os

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

def create_volcano_plot(df, fc_col='Log2FC', pval_col='P_value', padj_col='P_adj_BH', 
                       title='Volcano Plot', output_file=None, highlight_genes=None):
    """Create volcano plot"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Calculate -log10(p)
    df = df.copy()
    df['neg_log10_p'] = -np.log10(df[pval_col] + 1e-300)
    
    # Define significance thresholds
    fc_threshold = 1.0
    p_threshold = 0.05
    if padj_col in df.columns:
        sig_mask = (df[padj_col] < 0.05) & (np.abs(df[fc_col]) > fc_threshold)
        highly_sig_mask = (df[padj_col] < 0.01) & (np.abs(df[fc_col]) > fc_threshold * 1.5)
    else:
        sig_mask = (df[pval_col] < 0.05) & (np.abs(df[fc_col]) > fc_threshold)
        highly_sig_mask = (df[pval_col] < 0.01) & (np.abs(df[fc_col]) > fc_threshold * 1.5)
    
    # Plot points
    ax.scatter(df[~sig_mask][fc_col], df[~sig_mask]['neg_log10_p'], 
              alpha=0.5, s=30, c='gray', label='Not significant')
    ax.scatter(df[sig_mask & ~highly_sig_mask][fc_col], 
              df[sig_mask & ~highly_sig_mask]['neg_log10_p'],
              alpha=0.7, s=50, c='blue', label='Significant (FDR < 0.05)')
    ax.scatter(df[highly_sig_mask][fc_col], 
              df[highly_sig_mask]['neg_log10_p'],
              alpha=0.8, s=70, c='red', label='Highly significant (FDR < 0.01)')
    
    # Highlight specific genes if provided
    if highlight_genes is not None:
        for gene_symbol in highlight_genes:
            gene_row = df[df.get('Gene_Symbol', '') == gene_symbol]
            if len(gene_row) > 0:
                row = gene_row.iloc[0]
                ax.scatter(row[fc_col], row['neg_log10_p'], 
                          s=200, c='orange', marker='*', 
                          edgecolors='black', linewidths=1.5, zorder=5)
                ax.annotate(gene_symbol, (row[fc_col], row['neg_log10_p']),
                           xytext=(5, 5), textcoords='offset points', fontsize=9)
    
    # Add threshold lines
    ax.axhline(-np.log10(p_threshold), color='red', linestyle='--', alpha=0.5, label=f'P = {p_threshold}')
    ax.axvline(fc_threshold, color='blue', linestyle='--', alpha=0.5, label=f'Log2FC = {fc_threshold}')
    ax.axvline(-fc_threshold, color='blue', linestyle='--', alpha=0.5)
    
    ax.set_xlabel('Log2 Fold Change', fontsize=12, fontweight='bold')
    ax.set_ylabel('-Log10 P-value', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"[OK] Saved volcano plot to {output_file}")
    
    return fig

def create_venn_diagram(sets_dict, title='Venn Diagram', output_file=None):
    """Create Venn diagram using matplotlib"""
    from matplotlib_venn import venn2, venn3
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if len(sets_dict) == 2:
        labels = list(sets_dict.keys())
        venn2([sets_dict[labels[0]], sets_dict[labels[1]]], 
              set_labels=labels, ax=ax)
    elif len(sets_dict) == 3:
        labels = list(sets_dict.keys())
        venn3([sets_dict[labels[0]], sets_dict[labels[1]], sets_dict[labels[2]]],
              set_labels=labels, ax=ax)
    else:
        print(f"[WARNING] Venn diagram supports 2-3 sets, got {len(sets_dict)}")
        return None
    
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"[OK] Saved Venn diagram to {output_file}")
    
    return fig

def create_venn_simple(sets_dict, title='Venn Diagram', output_file=None):
    """Create simple Venn diagram using matplotlib circles"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if len(sets_dict) == 2:
        labels = list(sets_dict.keys())
        set1, set2 = sets_dict[labels[0]], sets_dict[labels[1]]
        
        # Calculate overlaps
        overlap = set1 & set2
        only1 = set1 - set2
        only2 = set2 - set1
        
        # Draw circles
        circle1 = Circle((0.5, 0.5), 0.3, alpha=0.5, color='blue', label=labels[0])
        circle2 = Circle((0.6, 0.5), 0.3, alpha=0.5, color='red', label=labels[1])
        ax.add_patch(circle1)
        ax.add_patch(circle2)
        
        # Add text
        ax.text(0.5, 0.5, f'{len(only1)}\n{labels[0]} only', 
               ha='center', va='center', fontsize=10, fontweight='bold')
        ax.text(0.6, 0.5, f'{len(only2)}\n{labels[1]} only',
               ha='center', va='center', fontsize=10, fontweight='bold')
        ax.text(0.55, 0.5, f'{len(overlap)}\nOverlap',
               ha='center', va='center', fontsize=11, fontweight='bold', color='white')
        
        ax.set_xlim(0.2, 0.9)
        ax.set_ylim(0.2, 0.8)
        ax.set_aspect('equal')
        ax.axis('off')
    
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"[OK] Saved Venn diagram to {output_file}")
    
    return fig

def create_priority_score_distribution(df, score_col='Priority_Score', 
                                      category_col=None, output_file=None):
    """Create distribution plot of priority scores"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram
    axes[0].hist(df[score_col], bins=30, edgecolor='black', alpha=0.7)
    axes[0].axvline(df[score_col].median(), color='red', linestyle='--', 
                   label=f'Median: {df[score_col].median():.1f}')
    axes[0].axvline(df[score_col].mean(), color='blue', linestyle='--',
                   label=f'Mean: {df[score_col].mean():.1f}')
    axes[0].set_xlabel('Priority Score', fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Frequency', fontsize=11, fontweight='bold')
    axes[0].set_title('Priority Score Distribution', fontsize=12, fontweight='bold')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Box plot if category column provided
    if category_col and category_col in df.columns:
        categories = df[category_col].unique()
        data_by_category = [df[df[category_col] == cat][score_col].values 
                           for cat in categories]
        axes[1].boxplot(data_by_category, labels=categories)
        axes[1].set_ylabel('Priority Score', fontsize=11, fontweight='bold')
        axes[1].set_xlabel(category_col, fontsize=11, fontweight='bold')
        axes[1].set_title('Priority Score by Category', fontsize=12, fontweight='bold')
        axes[1].grid(True, alpha=0.3, axis='y')
        plt.setp(axes[1].xaxis.get_majorticklabels(), rotation=45, ha='right')
    else:
        axes[1].boxplot([df[score_col].values])
        axes[1].set_ylabel('Priority Score', fontsize=11, fontweight='bold')
        axes[1].set_title('Priority Score Box Plot', fontsize=12, fontweight='bold')
        axes[1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"[OK] Saved distribution plot to {output_file}")
    
    return fig

def create_fold_change_comparison(df_blood, df_tissue, output_file=None):
    """Compare fold changes between blood and tissue datasets"""
    # Get cross-validated genes
    if 'Gene_Symbol' not in df_blood.columns or 'Gene_Symbol' not in df_tissue.columns:
        print("[WARNING] Gene_Symbol column not found, skipping comparison")
        return None
    
    # Merge on gene symbol
    merged = pd.merge(
        df_blood[['Gene_Symbol', 'Log2FC']].rename(columns={'Log2FC': 'Blood_Log2FC'}),
        df_tissue[['Gene_Symbol', 'Log2FC']].rename(columns={'Log2FC': 'Tissue_Log2FC'}),
        on='Gene_Symbol',
        how='inner'
    )
    
    if len(merged) == 0:
        print("[WARNING] No common genes found for comparison")
        return None
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Scatter plot
    axes[0].scatter(merged['Blood_Log2FC'], merged['Tissue_Log2FC'], 
                   alpha=0.6, s=50, edgecolors='black', linewidths=0.5)
    axes[0].axhline(0, color='gray', linestyle='--', alpha=0.5)
    axes[0].axvline(0, color='gray', linestyle='--', alpha=0.5)
    
    # Add correlation
    corr = merged['Blood_Log2FC'].corr(merged['Tissue_Log2FC'])
    axes[0].text(0.05, 0.95, f'r = {corr:.3f}', transform=axes[0].transAxes,
                fontsize=12, fontweight='bold', verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    axes[0].set_xlabel('Blood Log2FC (GSE159610)', fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Tissue Log2FC (GSE122897)', fontsize=11, fontweight='bold')
    axes[0].set_title('Fold Change Correlation', fontsize=12, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    
    # Distribution comparison
    axes[1].hist(merged['Blood_Log2FC'], bins=20, alpha=0.6, label='Blood', 
                edgecolor='black')
    axes[1].hist(merged['Tissue_Log2FC'], bins=20, alpha=0.6, label='Tissue',
                edgecolor='black')
    axes[1].set_xlabel('Log2 Fold Change', fontsize=11, fontweight='bold')
    axes[1].set_ylabel('Frequency', fontsize=11, fontweight='bold')
    axes[1].set_title('Fold Change Distribution', fontsize=12, fontweight='bold')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"[OK] Saved fold change comparison to {output_file}")
    
    return fig

def main():
    """Create all visualizations"""
    print("=" * 80)
    print("CREATING INTEGRATION VISUALIZATIONS")
    print("=" * 80)
    
    os.makedirs("results/figures", exist_ok=True)
    
    # Load datasets
    print("\nLoading datasets...")
    
    # Blood dataset
    blood_file = "results/gse159610_top_candidate_genes.tsv"
    blood_df = None
    if os.path.exists(blood_file):
        blood_df = pd.read_csv(blood_file, sep='\t')
        print(f"Loaded blood data: {len(blood_df)} genes")
    
    # Tissue dataset
    tissue_file = "results/gse122897_differential_expression.tsv"
    tissue_df = None
    if os.path.exists(tissue_file):
        tissue_df = pd.read_csv(tissue_file, sep='\t')
        # Filter significant
        tissue_df = tissue_df[tissue_df['P_adj_BH'] < 0.1].copy()
        print(f"Loaded tissue data: {len(tissue_df)} significant genes")
    
    # Multi-omics integrated
    integrated_file = "results/multi_omics_comprehensive_prioritized_genes.tsv"
    integrated_df = None
    if os.path.exists(integrated_file):
        integrated_df = pd.read_csv(integrated_file, sep='\t')
        print(f"Loaded integrated data: {len(integrated_df)} genes")
    
    # Cross-validated dataset
    cross_file = "results/integrated_cross_dataset_genes.tsv"
    cross_df = None
    if os.path.exists(cross_file):
        cross_df = pd.read_csv(cross_file, sep='\t')
        print(f"Loaded cross-validated data: {len(cross_df)} genes")
    
    # 1. Blood volcano plot
    if blood_df is not None:
        print("\n1. Creating blood volcano plot...")
        create_volcano_plot(
            blood_df, 
            fc_col='Log2FC', 
            pval_col='P_value', 
            padj_col='P_adj_BH',
            title='Blood Transcriptomics - GSE159610\n(Aneurysm vs Control)',
            output_file='results/figures/volcano_plot_blood.png'
        )
        plt.close()
    
    # 2. Tissue volcano plot
    if tissue_df is not None:
        print("\n2. Creating tissue volcano plot...")
        create_volcano_plot(
            tissue_df,
            fc_col='Log2FC',
            pval_col='P_value',
            padj_col='P_adj_BH',
            title='Tissue Transcriptomics - GSE122897\n(Ruptured vs Unruptured)',
            output_file='results/figures/volcano_plot_tissue.png'
        )
        plt.close()
    
    # 3. Venn diagram
    if blood_df is not None and tissue_df is not None:
        print("\n3. Creating Venn diagram...")
        # Try to get gene sets
        blood_genes = set()
        tissue_genes = set()
        
        if 'Gene_Symbol' in blood_df.columns:
            blood_genes = set(blood_df['Gene_Symbol'].dropna().astype(str))
        if 'Gene' in tissue_df.columns:
            tissue_genes = set(tissue_df['Gene'].dropna().astype(str))
        elif 'Gene_ID' in tissue_df.columns:
            tissue_genes = set(tissue_df['Gene_ID'].dropna().astype(str))
        
        if blood_genes and tissue_genes:
            create_venn_simple(
                {'Blood (GSE159610)': blood_genes, 
                 'Tissue (GSE122897)': tissue_genes},
                title='Overlap Between Blood and Tissue Datasets',
                output_file='results/figures/venn_blood_tissue.png'
            )
            plt.close()
    
    # 4. Priority score distribution
    if integrated_df is not None and 'Multi_Omics_Priority_Score' in integrated_df.columns:
        print("\n4. Creating priority score distribution...")
        create_priority_score_distribution(
            integrated_df,
            score_col='Multi_Omics_Priority_Score',
            category_col='In_Blood' if 'In_Blood' in integrated_df.columns else None,
            output_file='results/figures/priority_score_distribution.png'
        )
        plt.close()
    
    # 5. Fold change comparison (if cross-validated data available)
    if cross_df is not None and len(cross_df) > 0:
        print("\n5. Creating fold change comparison...")
        # Need both blood and tissue FC columns
        if 'Blood_Log2FC' in cross_df.columns and 'Tissue_Log2FC' in cross_df.columns:
            create_fold_change_comparison(
                cross_df.rename(columns={'Blood_Log2FC': 'Log2FC'}),
                cross_df.rename(columns={'Tissue_Log2FC': 'Log2FC'}),
                output_file='results/figures/fold_change_comparison.png'
            )
            plt.close()
    
    print("\n" + "=" * 80)
    print("VISUALIZATION COMPLETE!")
    print("=" * 80)
    print("\nAll figures saved to results/figures/")

if __name__ == "__main__":
    main()






