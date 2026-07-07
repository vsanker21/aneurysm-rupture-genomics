"""
Create Visualizations for Multi-Omics Integration Results
Generates publication-quality figures for prioritized genes
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'Arial'

def create_output_dir():
    """Create output directory for figures"""
    os.makedirs("results/figures", exist_ok=True)
    return "results/figures"

def plot_priority_score_distribution(df, output_dir):
    """Plot distribution of multi-omics priority scores"""
    print("\nCreating priority score distribution plot...")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create histogram
    ax.hist(df['Multi_Omics_Priority_Score'], bins=30, edgecolor='black', alpha=0.7, color='steelblue')
    ax.axvline(df['Multi_Omics_Priority_Score'].mean(), color='red', linestyle='--', 
               linewidth=2, label=f'Mean: {df["Multi_Omics_Priority_Score"].mean():.2f}')
    ax.axvline(df['Multi_Omics_Priority_Score'].median(), color='green', linestyle='--', 
               linewidth=2, label=f'Median: {df["Multi_Omics_Priority_Score"].median():.2f}')
    
    ax.set_xlabel('Multi-Omics Priority Score', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of Multi-Omics Priority Scores', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'priority_score_distribution.png')
    plt.savefig(output_file, bbox_inches='tight')
    print(f"  [OK] Saved to {output_file}")
    plt.close()

def plot_dataset_coverage(df, output_dir):
    """Plot dataset coverage (blood, tissue, GWAS)"""
    print("\nCreating dataset coverage plot...")
    
    # Calculate coverage
    in_blood = df['In_Blood'].sum()
    in_tissue = df['In_Tissue'].sum()
    in_both = ((df['In_Blood'] == True) & (df['In_Tissue'] == True)).sum()
    in_gwas = df['In_GWAS'].sum()
    total = len(df)
    
    # Create data
    categories = ['Blood\n(GSE159610)', 'Tissue\n(GSE122897)', 'Both\nDatasets', 'GWAS\nSupport']
    counts = [in_blood, in_tissue, in_both, in_gwas]
    percentages = [c/total*100 for c in counts]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))
    
    # Bar plot
    bars = ax1.bar(categories, counts, color=['#2E86AB', '#A23B72', '#F18F01', '#C73E1D'], alpha=0.8)
    ax1.set_ylabel('Number of Genes', fontsize=20, fontweight='bold')
    ax1.set_title('Dataset Coverage', fontsize=22, fontweight='bold')
    ax1.tick_params(axis='both', labelsize=16)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar, count, pct in zip(bars, counts, percentages):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{count}\n({pct:.1f}%)',
                ha='center', va='bottom', fontsize=16, fontweight='bold')
    
    # Pie chart
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']
    wedges, texts, autotexts = ax2.pie(counts, labels=categories, autopct='%1.1f%%',
                                       colors=colors, startangle=90, textprops={'fontsize': 16})
    for t in autotexts:
        t.set_fontsize(14)
    ax2.set_title('Dataset Coverage Distribution', fontsize=22, fontweight='bold')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'dataset_coverage.png')
    plt.savefig(output_file, bbox_inches='tight')
    print(f"  [OK] Saved to {output_file}")
    plt.close()

def plot_top_genes_bar(df, top_n=20, output_dir=None):
    """Plot top N prioritized genes as bar chart"""
    print(f"\nCreating top {top_n} genes bar chart...")
    
    top_genes = df.head(top_n).copy()
    
    # Get gene symbols (use Gene_Symbol or Gene_ID)
    top_genes['Display_Name'] = top_genes['Gene_Symbol'].fillna(top_genes['Gene_ID'])
    top_genes = top_genes.sort_values('Multi_Omics_Priority_Score', ascending=True)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create horizontal bar chart
    y_pos = np.arange(len(top_genes))
    colors = ['#C73E1D' if row['In_GWAS'] else '#2E86AB' for _, row in top_genes.iterrows()]
    
    bars = ax.barh(y_pos, top_genes['Multi_Omics_Priority_Score'], color=colors, alpha=0.8)
    
    # Add evidence labels
    evidence_labels = []
    for _, row in top_genes.iterrows():
        evidence = []
        if row['In_Blood']:
            evidence.append('B')
        if row['In_Tissue']:
            evidence.append('T')
        if row['In_GWAS']:
            evidence.append('G')
        evidence_labels.append('+'.join(evidence) if evidence else 'None')
    
    # Set labels
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{name} ({ev})" for name, ev in zip(top_genes['Display_Name'], evidence_labels)], fontsize=9)
    ax.set_xlabel('Multi-Omics Priority Score', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {top_n} Multi-Omics Prioritized Genes', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#C73E1D', alpha=0.8, label='With GWAS Support'),
        Patch(facecolor='#2E86AB', alpha=0.8, label='No GWAS Support')
    ]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, f'top_{top_n}_genes_bar.png')
    plt.savefig(output_file, bbox_inches='tight')
    print(f"  [OK] Saved to {output_file}")
    plt.close()

def plot_evidence_venn(df, output_dir):
    """Create Venn diagram showing overlap between datasets"""
    print("\nCreating evidence overlap plot...")
    
    try:
        from matplotlib_venn import venn2, venn3
    except ImportError:
        print("  [WARNING] matplotlib_venn not installed, skipping Venn diagram")
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Calculate overlaps
    in_blood = set(df[df['In_Blood'] == True].index)
    in_tissue = set(df[df['In_Tissue'] == True].index)
    in_gwas = set(df[df['In_GWAS'] == True].index)
    
    # Create 3-way Venn
    venn3([in_blood, in_tissue, in_gwas],
          set_labels=('Blood\n(GSE159610)', 'Tissue\n(GSE122897)', 'GWAS'),
          ax=ax, alpha=0.7,
          set_colors=('#2E86AB', '#A23B72', '#F18F01'))
    
    ax.set_title('Multi-Omics Evidence Overlap', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'evidence_overlap_venn.png')
    plt.savefig(output_file, bbox_inches='tight')
    print(f"  [OK] Saved to {output_file}")
    plt.close()

def plot_log2fc_comparison(df, output_dir):
    """Plot Log2FC comparison between blood and tissue"""
    print("\nCreating Log2FC comparison plot...")
    
    # Get genes present in both datasets
    both_datasets = df[(df['In_Blood'] == True) & (df['In_Tissue'] == True)]
    
    if len(both_datasets) == 0:
        print("  [WARNING] No genes present in both datasets")
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Scatter plot
    ax.scatter(both_datasets['Blood_Log2FC'], both_datasets['Tissue_Log2FC'],
              s=100, alpha=0.6, c=both_datasets['Multi_Omics_Priority_Score'],
              cmap='viridis', edgecolors='black', linewidth=0.5)
    
    # Add diagonal line
    max_val = max(abs(both_datasets['Blood_Log2FC'].max()),
                  abs(both_datasets['Tissue_Log2FC'].max()))
    ax.plot([-max_val, max_val], [-max_val, max_val], 'r--', alpha=0.5, label='y=x')
    
    ax.set_xlabel('Blood Log2FC (GSE159610)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Tissue Log2FC (GSE122897)', fontsize=12, fontweight='bold')
    ax.set_title('Log2FC Comparison: Blood vs Tissue', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Add colorbar
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label('Priority Score', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'log2fc_comparison.png')
    plt.savefig(output_file, bbox_inches='tight')
    print(f"  [OK] Saved to {output_file}")
    plt.close()

def plot_pathway_enrichment(pathway_file, output_dir, top_n=15):
    """Plot top enriched pathways"""
    print("\nCreating pathway enrichment plot...")
    
    if not os.path.exists(pathway_file):
        print(f"  [WARNING] Pathway file not found: {pathway_file}")
        return
    
    pathways_df = pd.read_csv(pathway_file, sep='\t')
    top_pathways = pathways_df.head(top_n).copy()
    
    if len(top_pathways) == 0:
        print("  [WARNING] No pathways to plot")
        return
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create horizontal bar chart
    y_pos = np.arange(len(top_pathways))
    top_pathways = top_pathways.sort_values('P_value', ascending=True)
    
    # Truncate long pathway names
    pathway_names = [name[:60] + '...' if len(name) > 60 else name 
                     for name in top_pathways['Term']]
    
    bars = ax.barh(y_pos, -np.log10(top_pathways['P_value'] + 1e-300),
                   color='steelblue', alpha=0.8)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(pathway_names, fontsize=9)
    ax.set_xlabel('-log10(P-value)', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {len(top_pathways)} Enriched Pathways', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    # Add intersection size as text
    for i, (idx, row) in enumerate(top_pathways.iterrows()):
        ax.text(-np.log10(row['P_value'] + 1e-300) + 0.1, i,
                f"({row['Intersection_size']} genes)", va='center', fontsize=8)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'pathway_enrichment.png')
    plt.savefig(output_file, bbox_inches='tight')
    print(f"  [OK] Saved to {output_file}")
    plt.close()

def main():
    """Main visualization pipeline"""
    
    print("=" * 80)
    print("MULTI-OMICS INTEGRATION VISUALIZATIONS")
    print("=" * 80)
    
    # Create output directory
    output_dir = create_output_dir()
    
    # Load prioritized genes
    prioritized_file = "results/multi_omics_comprehensive_prioritized_genes.tsv"
    if not os.path.exists(prioritized_file):
        print(f"[ERROR] Prioritized genes file not found: {prioritized_file}")
        return
    
    df = pd.read_csv(prioritized_file, sep='\t')
    print(f"\nLoaded {len(df):,} prioritized genes")
    
    # Create visualizations
    plot_priority_score_distribution(df, output_dir)
    plot_dataset_coverage(df, output_dir)
    plot_top_genes_bar(df, top_n=20, output_dir=output_dir)
    plot_evidence_venn(df, output_dir)
    plot_log2fc_comparison(df, output_dir)
    
    # Plot pathway enrichment if available
    pathway_file = "results/pathway_enrichment_prioritized_fdr_0.05.tsv"
    if os.path.exists(pathway_file):
        plot_pathway_enrichment(pathway_file, output_dir, top_n=10)
    
    print("\n" + "=" * 80)
    print("VISUALIZATION COMPLETE")
    print("=" * 80)
    print(f"\nAll figures saved to: {output_dir}")

if __name__ == "__main__":
    main()


