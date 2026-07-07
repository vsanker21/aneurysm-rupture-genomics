"""
Generate Visualizations for Analysis Results
Creates publication-quality figures and summaries
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

def create_output_dir():
    """Create output directory for figures"""
    os.makedirs("results/figures", exist_ok=True)
    return "results/figures"

def plot_de_volcano(de_df, output_dir):
    """Create volcano plot for differential expression"""
    print("\nCreating volcano plot...")
    
    if 'Log2FC' not in de_df.columns or 'P_adj_BH' not in de_df.columns:
        print("  [WARNING] Required columns not found in DE data")
        return None
    
    # Prepare data
    de_df = de_df.copy()
    de_df['-log10_P_adj'] = -np.log10(de_df['P_adj_BH'] + 1e-300)
    
    # Identify significant genes
    sig_mask = de_df['P_adj_BH'] < 0.05
    highly_sig_mask = de_df['P_adj_BH'] < 0.01
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot non-significant
    ax.scatter(de_df[~sig_mask]['Log2FC'], 
               de_df[~sig_mask]['-log10_P_adj'],
               c='gray', alpha=0.5, s=20, label='Not significant')
    
    # Plot significant
    ax.scatter(de_df[sig_mask & ~highly_sig_mask]['Log2FC'],
               de_df[sig_mask & ~highly_sig_mask]['-log10_P_adj'],
               c='blue', alpha=0.6, s=30, label='P_adj < 0.05')
    
    # Plot highly significant
    ax.scatter(de_df[highly_sig_mask]['Log2FC'],
               de_df[highly_sig_mask]['-log10_P_adj'],
               c='red', alpha=0.7, s=40, label='P_adj < 0.01')
    
    # Add threshold lines
    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
    
    ax.set_xlabel('Log2 Fold Change', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log10(Adjusted P-value)', fontsize=12, fontweight='bold')
    ax.set_title('Volcano Plot: Ruptured vs Unruptured Aneurysms', fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, "volcano_plot_de.png")
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved to {output_file}")
    return output_file

def plot_ml_performance(ml_file, output_dir):
    """Plot ML model performance comparison"""
    print("\nCreating ML performance plot...")
    
    if not os.path.exists(ml_file):
        print(f"  [WARNING] ML results file not found: {ml_file}")
        return None
    
    ml_df = pd.read_csv(ml_file, sep='\t')
    
    if 'Model' not in ml_df.columns or 'Mean_AUC' not in ml_df.columns:
        print("  [WARNING] Required columns not found")
        return None
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    models = ml_df['Model'].values
    aucs = ml_df['Mean_AUC'].values
    stds = ml_df['Std_AUC'].values if 'Std_AUC' in ml_df.columns else np.zeros(len(aucs))
    
    x_pos = np.arange(len(models))
    bars = ax.bar(x_pos, aucs, yerr=stds, capsize=5, alpha=0.7, 
                  color=['#1f77b4', '#ff7f0e', '#2ca02c'])
    
    ax.set_xlabel('Model', fontsize=12, fontweight='bold')
    ax.set_ylabel('Cross-Validation AUC', fontsize=12, fontweight='bold')
    ax.set_title('Machine Learning Model Performance', fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(models, rotation=45, ha='right')
    ax.set_ylim([0, 1.0])
    ax.axhline(y=0.5, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Random (AUC=0.5)')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for i, (bar, auc, std) in enumerate(zip(bars, aucs, stds)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + std + 0.02,
                f'{auc:.3f}±{std:.3f}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, "ml_performance_comparison.png")
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved to {output_file}")
    return output_file

def plot_ec_dysfunction_comparison(ec_file, metadata_file, output_dir):
    """Plot EC dysfunction scores by group"""
    print("\nCreating EC dysfunction comparison plot...")
    
    if not os.path.exists(ec_file):
        print(f"  [WARNING] EC dysfunction file not found: {ec_file}")
        return None
    
    ec_df = pd.read_csv(ec_file, sep='\t')
    
    if 'EC_Dysfunction_Score' not in ec_df.columns:
        print("  [WARNING] EC_Dysfunction_Score column not found")
        return None
    
    # Try to load metadata
    metadata = None
    if os.path.exists(metadata_file):
        metadata = pd.read_csv(metadata_file, sep='\t')
        if 'Sample' in metadata.columns and 'Group' in metadata.columns:
            ec_df = ec_df.merge(metadata[['Sample', 'Group']], on='Sample', how='left')
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if metadata is not None and 'Group' in ec_df.columns:
        # Box plot by group
        groups = ec_df['Group'].dropna().unique()
        if len(groups) >= 2:
            data_to_plot = [ec_df[ec_df['Group'] == group]['EC_Dysfunction_Score'].dropna() 
                          for group in groups]
            
            bp = ax.boxplot(data_to_plot, labels=groups, patch_artist=True)
            
            # Color the boxes
            colors = ['#ff7f0e', '#2ca02c', '#1f77b4']
            for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            
            # Add statistical test
            if len(groups) == 2:
                group1_data = ec_df[ec_df['Group'] == groups[0]]['EC_Dysfunction_Score'].dropna()
                group2_data = ec_df[ec_df['Group'] == groups[1]]['EC_Dysfunction_Score'].dropna()
                if len(group1_data) > 0 and len(group2_data) > 0:
                    t_stat, p_val = stats.ttest_ind(group1_data, group2_data)
                    ax.text(0.5, 0.95, f'T-test: p={p_val:.4f}', 
                           transform=ax.transAxes, ha='center', 
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    else:
        # Histogram if no groups
        ax.hist(ec_df['EC_Dysfunction_Score'].dropna(), bins=20, alpha=0.7, color='steelblue')
        ax.set_xlabel('EC Dysfunction Score', fontsize=12, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    
    ax.set_title('Endothelial Dysfunction Score Distribution', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, "ec_dysfunction_comparison.png")
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved to {output_file}")
    return output_file

def plot_pathway_enrichment(pathway_file, output_dir, top_n=15):
    """Plot top enriched pathways"""
    print("\nCreating pathway enrichment plot...")
    
    if not os.path.exists(pathway_file):
        print(f"  [WARNING] Pathway file not found: {pathway_file}")
        return None
    
    pathway_df = pd.read_csv(pathway_file, sep='\t')
    
    if len(pathway_df) == 0:
        print("  [WARNING] No pathways to plot")
        return None
    
    # Take top N pathways
    top_pathways = pathway_df.head(top_n).copy()
    top_pathways['-log10_P'] = -np.log10(top_pathways['P_value'] + 1e-300)
    
    # Truncate long pathway names
    top_pathways['Term_Short'] = top_pathways['Term'].apply(
        lambda x: x[:50] + '...' if len(str(x)) > 50 else str(x)
    )
    
    fig, ax = plt.subplots(figsize=(10, max(6, len(top_pathways) * 0.4)))
    
    y_pos = np.arange(len(top_pathways))
    bars = ax.barh(y_pos, top_pathways['-log10_P'], alpha=0.7, 
                   color='steelblue')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_pathways['Term_Short'])
    ax.set_xlabel('-log10(P-value)', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {len(top_pathways)} Enriched Pathways', fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')
    
    # Add value labels
    for i, (bar, pval) in enumerate(zip(bars, top_pathways['P_value'])):
        width = bar.get_width()
        ax.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                f'p={pval:.2e}', ha='left', va='center', fontsize=8)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, "pathway_enrichment_top.png")
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved to {output_file}")
    return output_file

def plot_risk_stratification(risk_file, output_dir):
    """Plot risk stratification results"""
    print("\nCreating risk stratification plot...")
    
    if not os.path.exists(risk_file):
        print(f"  [WARNING] Risk stratification file not found: {risk_file}")
        return None
    
    risk_df = pd.read_csv(risk_file, sep='\t')
    
    if 'Risk_Category' not in risk_df.columns:
        print("  [WARNING] Risk_Category column not found")
        return None
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Pie chart of risk categories
    risk_counts = risk_df['Risk_Category'].value_counts()
    colors = ['#2ca02c', '#ff7f0e', '#d62728']  # Green, Orange, Red
    ax1.pie(risk_counts.values, labels=risk_counts.index, autopct='%1.1f%%',
           colors=colors[:len(risk_counts)], startangle=90)
    ax1.set_title('Risk Category Distribution', fontsize=12, fontweight='bold')
    
    # Histogram of rupture probabilities
    if 'Rupture_Probability' in risk_df.columns:
        ax2.hist(risk_df['Rupture_Probability'], bins=20, alpha=0.7, color='steelblue', edgecolor='black')
        ax2.axvline(x=0.33, color='green', linestyle='--', linewidth=2, label='Low Risk Threshold')
        ax2.axvline(x=0.67, color='red', linestyle='--', linewidth=2, label='High Risk Threshold')
        ax2.set_xlabel('Rupture Probability', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Frequency', fontsize=12, fontweight='bold')
        ax2.set_title('Rupture Probability Distribution', fontsize=12, fontweight='bold')
        ax2.legend()
        ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, "risk_stratification.png")
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved to {output_file}")
    return output_file

def generate_summary_dashboard(output_dir):
    """Generate a summary dashboard with key metrics"""
    print("\nCreating summary dashboard...")
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Load and summarize data
    summary_text = []
    
    # DE Summary
    de_file = "results/gse122897_differential_expression.tsv"
    if os.path.exists(de_file):
        de_df = pd.read_csv(de_file, sep='\t')
        sig_genes = (de_df['P_adj_BH'] < 0.05).sum() if 'P_adj_BH' in de_df.columns else 0
        summary_text.append(f"DE Analysis: {sig_genes} significant genes")
    
    # ML Summary
    ml_file = "results/blood_ml_cv_results.tsv"
    if os.path.exists(ml_file):
        ml_df = pd.read_csv(ml_file, sep='\t')
        if 'Mean_AUC' in ml_df.columns:
            best_auc = ml_df['Mean_AUC'].max()
            best_model = ml_df.loc[ml_df['Mean_AUC'].idxmax(), 'Model'] if 'Model' in ml_df.columns else 'N/A'
            summary_text.append(f"ML Best Model: {best_model} (AUC={best_auc:.3f})")
    
    # Create text summary
    ax_text = fig.add_subplot(gs[0, :])
    ax_text.axis('off')
    summary_str = "\n".join(summary_text) if summary_text else "No summary data available"
    ax_text.text(0.1, 0.5, summary_str, fontsize=14, verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax_text.set_title('Analysis Summary', fontsize=16, fontweight='bold', pad=20)
    
    plt.suptitle('Cerebral Aneurysm Genomics Project - Analysis Dashboard', 
                fontsize=18, fontweight='bold', y=0.98)
    
    output_file = os.path.join(output_dir, "summary_dashboard.png")
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Saved to {output_file}")
    return output_file

def main():
    """Generate all visualizations"""
    print("=" * 80)
    print("GENERATING VISUALIZATIONS")
    print("=" * 80)
    
    output_dir = create_output_dir()
    generated_files = []
    
    # 1. Volcano plot
    de_file = "results/gse122897_differential_expression.tsv"
    if os.path.exists(de_file):
        de_df = pd.read_csv(de_file, sep='\t')
        if len(de_df) > 0:
            volcano_file = plot_de_volcano(de_df, output_dir)
            if volcano_file:
                generated_files.append(volcano_file)
    
    # 2. ML performance
    ml_file = "results/blood_ml_cv_results.tsv"
    if os.path.exists(ml_file):
        ml_file_path = plot_ml_performance(ml_file, output_dir)
        if ml_file_path:
            generated_files.append(ml_file_path)
    
    # 3. EC dysfunction
    ec_file = "results/ec_dysfunction_scores.tsv"
    metadata_file = "data/sc_rnaseq/GSE122897/sample_metadata.tsv"
    if os.path.exists(ec_file):
        ec_file_path = plot_ec_dysfunction_comparison(ec_file, metadata_file, output_dir)
        if ec_file_path:
            generated_files.append(ec_file_path)
    
    # 4. Pathway enrichment
    pathway_file = "results/pathway_enrichment_fdr_0.05.tsv"
    if os.path.exists(pathway_file):
        pathway_file_path = plot_pathway_enrichment(pathway_file, output_dir)
        if pathway_file_path:
            generated_files.append(pathway_file_path)
    
    # 5. Risk stratification
    risk_file = "results/risk_stratification.tsv"
    if os.path.exists(risk_file):
        risk_file_path = plot_risk_stratification(risk_file, output_dir)
        if risk_file_path:
            generated_files.append(risk_file_path)
    
    # 6. Summary dashboard
    dashboard_file = generate_summary_dashboard(output_dir)
    if dashboard_file:
        generated_files.append(dashboard_file)
    
    print("\n" + "=" * 80)
    print("VISUALIZATION GENERATION COMPLETE")
    print("=" * 80)
    print(f"\nGenerated {len(generated_files)} visualization(s):")
    for f in generated_files:
        print(f"  - {f}")

if __name__ == "__main__":
    main()






