"""
Phase 7: Generate Publication-Quality Visualizations
Creates comprehensive figures for validation and comparison results
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("PHASE 7: GENERATING PUBLICATION-QUALITY VISUALIZATIONS")
print("=" * 80)

# Set publication-quality style — large, visible fonts for all figures
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 18
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 22
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['figure.titlesize'] = 24

# Create output directory
output_dir = Path("results/figures/phase7")
output_dir.mkdir(parents=True, exist_ok=True)
print(f"\nOutput directory: {output_dir}")

# Load data
print("\nLoading Phase 7 results...")
with open('results/phase7_comprehensive_summary.json', 'r') as f:
    summary = json.load(f)

roc_data = pd.read_csv('results/phase7_roc_curve_optimized.tsv', sep='\t')
pr_data = pd.read_csv('results/phase7_pr_curve_optimized.tsv', sep='\t')
risk_strat = pd.read_csv('results/phase7_risk_stratification.tsv', sep='\t')
metrics_df = pd.read_csv('results/phase7_optimized_model_metrics.tsv', sep='\t')

# Load features data for baseline comparison
features_df = pd.read_csv('results/integrated_features.tsv', sep='\t')
features_df = features_df[features_df['Group'].isin(['Ruptured', 'Unruptured'])].copy()
labels = (features_df['Group'] == 'Ruptured').astype(int)

print(f"  Loaded {len(roc_data)} ROC curve points")
print(f"  Loaded {len(pr_data)} PR curve points")
print(f"  Loaded {len(risk_strat)} risk stratification samples")

# ============================================================================
# 1. ROC CURVE
# ============================================================================
print("\n1. Creating ROC Curve...")
fig, ax = plt.subplots(figsize=(12, 12))

# Optimized model ROC
auc_opt = summary['optimized_model']['auc']
ax.plot(roc_data['FPR'], roc_data['TPR'], 
        linewidth=2.5, color='#2E86AB', label=f'Optimized Model (AUC = {auc_opt:.3f})')

# Baseline ROC (if available)
baseline_auc = summary['baseline_comparison']['baseline_auc']
if baseline_auc > 0:
    # Create diagonal line for baseline
    ax.plot([0, 1], [0, 1], '--', linewidth=2, color='#A23B72', 
            alpha=0.7, label=f'Baseline (AUC = {baseline_auc:.3f})')

# Diagonal reference line
ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.3, label='Random (AUC = 0.500)')

ax.set_xlabel('False Positive Rate (1 - Specificity)', fontweight='bold', fontsize=20)
ax.set_ylabel('True Positive Rate (Sensitivity)', fontweight='bold', fontsize=20)
ax.set_title('ROC Curve: Optimized Model Performance', fontweight='bold', fontsize=24, pad=15)
ax.legend(loc='lower right', frameon=True, fancybox=True, shadow=True, fontsize=16)
ax.tick_params(axis='both', labelsize=16)
ax.grid(True, alpha=0.3)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])

plt.tight_layout()
roc_file = output_dir / "phase7_roc_curve.png"
plt.savefig(roc_file, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  [OK] Saved to {roc_file}")

# ============================================================================
# 2. PRECISION-RECALL CURVE
# ============================================================================
print("\n2. Creating Precision-Recall Curve...")
fig, ax = plt.subplots(figsize=(12, 12))

pr_auc = summary['optimized_model']['pr_auc']
ax.plot(pr_data['Recall'], pr_data['Precision'], 
        linewidth=2.5, color='#F18F01', label=f'Optimized Model (PR AUC = {pr_auc:.3f})')

# Baseline PR (no skill line)
baseline_precision = summary['dataset']['rupture_rate']
ax.axhline(y=baseline_precision, color='#A23B72', linestyle='--', 
           linewidth=2, alpha=0.7, label=f'Baseline (PR = {baseline_precision:.3f})')

ax.set_xlabel('Recall (Sensitivity)', fontweight='bold')
ax.set_ylabel('Precision (Positive Predictive Value)', fontweight='bold')
ax.set_title('Precision-Recall Curve: Optimized Model Performance', fontweight='bold', pad=15)
ax.legend(loc='lower left', frameon=True, fancybox=True, shadow=True)
ax.grid(True, alpha=0.3)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])

plt.tight_layout()
pr_file = output_dir / "phase7_pr_curve.png"
plt.savefig(pr_file, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  [OK] Saved to {pr_file}")

# ============================================================================
# 3. CONFUSION MATRIX
# ============================================================================
print("\n3. Creating Confusion Matrix...")
from sklearn.metrics import confusion_matrix

# Get binary predictions at optimal threshold
opt_threshold = summary['optimized_model']['optimal_threshold']
opt_pred_binary = (risk_strat['Optimized_Score'] >= opt_threshold).astype(int)
cm = confusion_matrix(risk_strat['True_Label'], opt_pred_binary)

fig, ax = plt.subplots(figsize=(7, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=True,
            xticklabels=['Unruptured', 'Ruptured'],
            yticklabels=['Unruptured', 'Ruptured'],
            ax=ax, annot_kws={'size': 14, 'weight': 'bold'})

ax.set_xlabel('Predicted Label', fontweight='bold', fontsize=20)
ax.set_ylabel('True Label', fontweight='bold', fontsize=20)
ax.set_title('Confusion Matrix (Optimal Threshold = {:.3f})'.format(opt_threshold), 
             fontweight='bold', fontsize=22, pad=15)
ax.tick_params(axis='both', labelsize=16)

# Add performance metrics as text
tn, fp, fn, tp = cm.ravel()
accuracy = (tn + tp) / (tn + fp + fn + tp)
precision = tp / (tp + fp) if (tp + fp) > 0 else 0
recall = tp / (tp + fn) if (tp + fn) > 0 else 0
f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

metrics_text = f'Accuracy: {accuracy:.3f}\nPrecision: {precision:.3f}\nRecall: {recall:.3f}\nF1-Score: {f1:.3f}'
ax.text(0.5, -0.15, metrics_text, transform=ax.transAxes,
        ha='center', fontsize=16, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
cm_file = output_dir / "phase7_confusion_matrix.png"
plt.savefig(cm_file, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  [OK] Saved to {cm_file}")

# ============================================================================
# 4. RISK STRATIFICATION PLOT
# ============================================================================
print("\n4. Creating Risk Stratification Plot...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

# Panel A: Risk score distribution by true label
ruptured_scores = risk_strat[risk_strat['True_Label'] == 1]['Optimized_Score']
unruptured_scores = risk_strat[risk_strat['True_Label'] == 0]['Optimized_Score']

ax1.hist(unruptured_scores, bins=15, alpha=0.6, color='#06A77D', 
         label='Unruptured', edgecolor='black', linewidth=1.2)
ax1.hist(ruptured_scores, bins=15, alpha=0.6, color='#D00000', 
         label='Ruptured', edgecolor='black', linewidth=1.2)
ax1.axvline(x=opt_threshold, color='black', linestyle='--', 
            linewidth=2, label=f'Optimal Threshold ({opt_threshold:.3f})')
ax1.set_xlabel('Predicted Rupture Risk Score', fontweight='bold')
ax1.set_ylabel('Number of Samples', fontweight='bold')
ax1.set_title('A. Risk Score Distribution', fontweight='bold', pad=10)
ax1.legend(frameon=True, fancybox=True, shadow=True)
ax1.grid(True, alpha=0.3, axis='y')

# Panel B: Rupture rate by risk category
risk_categories = ['Low', 'Medium', 'High']
rupture_rates = []
n_samples = []
for category in risk_categories:
    mask = risk_strat['Optimized_Risk'] == category
    if mask.sum() > 0:
        rupture_rate = risk_strat[mask]['True_Label'].mean()
        n = mask.sum()
        rupture_rates.append(rupture_rate)
        n_samples.append(n)
    else:
        rupture_rates.append(0)
        n_samples.append(0)

colors = ['#06A77D', '#F18F01', '#D00000']
bars = ax2.bar(risk_categories, rupture_rates, color=colors, alpha=0.7, 
               edgecolor='black', linewidth=1.5)
ax2.set_ylabel('Rupture Rate', fontweight='bold')
ax2.set_xlabel('Risk Category', fontweight='bold')
ax2.set_title('B. Rupture Rate by Risk Category', fontweight='bold', pad=10)
ax2.set_ylim([0, 1.1])
ax2.grid(True, alpha=0.3, axis='y')

# Add sample counts on bars
for i, (bar, rate, n) in enumerate(zip(bars, rupture_rates, n_samples)):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + 0.02,
             f'{rate:.1%}\n(n={n})', ha='center', va='bottom', 
             fontweight='bold', fontsize=16)

plt.tight_layout()
risk_file = output_dir / "phase7_risk_stratification.png"
plt.savefig(risk_file, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  [OK] Saved to {risk_file}")

# ============================================================================
# 5. PERFORMANCE METRICS COMPARISON
# ============================================================================
print("\n5. Creating Performance Metrics Comparison...")
fig, ax = plt.subplots(figsize=(14, 9))

# Extract metrics
metrics_to_plot = ['ROC_AUC', 'PR_AUC', 'Accuracy', 'Precision', 'Recall', 'F1_Score']
metric_values = []
for metric in metrics_to_plot:
    value = metrics_df[metrics_df['Metric'] == metric]['Optimized_Model'].values[0]
    metric_values.append(value)

# Create bar plot
colors_map = plt.cm.viridis(np.linspace(0.2, 0.8, len(metrics_to_plot)))
bars = ax.bar(metrics_to_plot, metric_values, color=colors_map, alpha=0.7,
              edgecolor='black', linewidth=1.2)

# Add value labels on bars
for bar, value in zip(bars, metric_values):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
            f'{value:.3f}', ha='center', va='bottom', fontweight='bold', fontsize=10)

ax.set_ylabel('Score', fontweight='bold')
ax.set_xlabel('Metric', fontweight='bold')
ax.set_title('Optimized Model Performance Metrics', fontweight='bold', pad=15)
ax.set_ylim([0, 1.1])
ax.grid(True, alpha=0.3, axis='y')
plt.xticks(rotation=45, ha='right')

plt.tight_layout()
metrics_file = output_dir / "phase7_performance_metrics.png"
plt.savefig(metrics_file, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  [OK] Saved to {metrics_file}")

# ============================================================================
# 6. BASELINE COMPARISON
# ============================================================================
print("\n6. Creating Baseline Comparison...")
fig, ax = plt.subplots(figsize=(12, 10))

models = ['Baseline\n(Clinical Score)', 'Optimized\nModel']
aucs = [summary['baseline_comparison']['baseline_auc'], 
        summary['baseline_comparison']['optimized_auc']]
improvement = summary['baseline_comparison']['improvement']
relative_improvement = summary['baseline_comparison']['relative_improvement_percent']

colors = ['#A23B72', '#2E86AB']
bars = ax.bar(models, aucs, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)

# Add value labels
for bar, auc_val in zip(bars, aucs):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
            f'{auc_val:.3f}', ha='center', va='bottom', fontweight='bold', fontsize=16)

# Add improvement annotation
ax.annotate(f'+{improvement:.3f} AUC\n(+{relative_improvement:.1f}%)',
            xy=(1, aucs[1]), xytext=(0.5, 0.85),
            arrowprops=dict(arrowstyle='->', color='black', lw=2),
            fontsize=16, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7),
            ha='center', transform=ax.transAxes)

ax.set_ylabel('ROC AUC', fontweight='bold')
ax.set_title('Model Performance: Baseline vs Optimized', fontweight='bold', pad=15)
ax.set_ylim([0, 1.0])
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
comparison_file = output_dir / "phase7_baseline_comparison.png"
plt.savefig(comparison_file, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  [OK] Saved to {comparison_file}")

# ============================================================================
# 7. COMBINED FIGURE (Publication-ready)
# ============================================================================
print("\n7. Creating Combined Publication Figure...")
fig = plt.figure(figsize=(22, 16))
gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

# Panel A: ROC Curve
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(roc_data['FPR'], roc_data['TPR'], linewidth=2.5, color='#2E86AB',
         label=f'Optimized Model (AUC = {auc_opt:.3f})')
ax1.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.3, label='Random')
ax1.set_xlabel('False Positive Rate', fontweight='bold')
ax1.set_ylabel('True Positive Rate', fontweight='bold')
ax1.set_title('A. ROC Curve', fontweight='bold', loc='left')
ax1.legend(loc='lower right', frameon=True)
ax1.grid(True, alpha=0.3)
ax1.set_xlim([0, 1])
ax1.set_ylim([0, 1])

# Panel B: PR Curve
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(pr_data['Recall'], pr_data['Precision'], linewidth=2.5, color='#F18F01',
         label=f'Optimized Model (PR AUC = {pr_auc:.3f})')
ax2.axhline(y=baseline_precision, color='#A23B72', linestyle='--', 
            linewidth=2, alpha=0.7, label=f'Baseline')
ax2.set_xlabel('Recall', fontweight='bold')
ax2.set_ylabel('Precision', fontweight='bold')
ax2.set_title('B. Precision-Recall Curve', fontweight='bold', loc='left')
ax2.legend(loc='lower left', frameon=True)
ax2.grid(True, alpha=0.3)
ax2.set_xlim([0, 1])
ax2.set_ylim([0, 1])

# Panel C: Confusion Matrix
ax3 = fig.add_subplot(gs[0, 2])
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=True,
            xticklabels=['Unruptured', 'Ruptured'],
            yticklabels=['Unruptured', 'Ruptured'],
            ax=ax3, annot_kws={'size': 12, 'weight': 'bold'})
ax3.set_xlabel('Predicted', fontweight='bold')
ax3.set_ylabel('True', fontweight='bold')
ax3.set_title('C. Confusion Matrix', fontweight='bold', loc='left')

# Panel D: Risk Score Distribution
ax4 = fig.add_subplot(gs[1, :2])
ax4.hist(unruptured_scores, bins=15, alpha=0.6, color='#06A77D', 
         label='Unruptured', edgecolor='black', linewidth=1.2)
ax4.hist(ruptured_scores, bins=15, alpha=0.6, color='#D00000', 
         label='Ruptured', edgecolor='black', linewidth=1.2)
ax4.axvline(x=opt_threshold, color='black', linestyle='--', 
            linewidth=2, label=f'Threshold ({opt_threshold:.3f})')
ax4.set_xlabel('Predicted Rupture Risk Score', fontweight='bold')
ax4.set_ylabel('Number of Samples', fontweight='bold')
ax4.set_title('D. Risk Score Distribution', fontweight='bold', loc='left')
ax4.legend(frameon=True)
ax4.grid(True, alpha=0.3, axis='y')

# Panel E: Rupture Rate by Category
ax5 = fig.add_subplot(gs[1, 2])
bars = ax5.bar(risk_categories, rupture_rates, color=colors, alpha=0.7, 
               edgecolor='black', linewidth=1.5)
ax5.set_ylabel('Rupture Rate', fontweight='bold')
ax5.set_xlabel('Risk Category', fontweight='bold')
ax5.set_title('E. Rupture Rate by Category', fontweight='bold', loc='left')
ax5.set_ylim([0, 1.1])
ax5.grid(True, alpha=0.3, axis='y')
for bar, rate, n in zip(bars, rupture_rates, n_samples):
    height = bar.get_height()
    ax5.text(bar.get_x() + bar.get_width()/2., height + 0.02,
             f'{rate:.1%}\n(n={n})', ha='center', va='bottom', 
             fontweight='bold', fontsize=14)

# Panel F: Performance Metrics
ax6 = fig.add_subplot(gs[2, :])
bars = ax6.bar(metrics_to_plot, metric_values, color=colors_map, alpha=0.7,
               edgecolor='black', linewidth=1.2)
for bar, value in zip(bars, metric_values):
    height = bar.get_height()
    ax6.text(bar.get_x() + bar.get_width()/2., height + 0.01,
            f'{value:.3f}', ha='center', va='bottom', fontweight='bold', fontsize=10)
ax6.set_ylabel('Score', fontweight='bold')
ax6.set_xlabel('Metric', fontweight='bold')
ax6.set_title('F. Performance Metrics Summary', fontweight='bold', loc='left')
ax6.set_ylim([0, 1.1])
ax6.grid(True, alpha=0.3, axis='y')
plt.setp(ax6.xaxis.get_majorticklabels(), rotation=45, ha='right')

# Add overall title
fig.suptitle('Phase 7: Validation and Comparison Results', 
             fontsize=24, fontweight='bold', y=0.995)

combined_file = output_dir / "phase7_combined_figure.png"
plt.savefig(combined_file, bbox_inches='tight', facecolor='white', dpi=300)
plt.close()
print(f"  [OK] Saved to {combined_file}")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("VISUALIZATION GENERATION COMPLETE")
print("=" * 80)
print(f"\nGenerated {7} publication-quality figures:")
print(f"  1. ROC Curve: {roc_file}")
print(f"  2. Precision-Recall Curve: {pr_file}")
print(f"  3. Confusion Matrix: {cm_file}")
print(f"  4. Risk Stratification: {risk_file}")
print(f"  5. Performance Metrics: {metrics_file}")
print(f"  6. Baseline Comparison: {comparison_file}")
print(f"  7. Combined Figure: {combined_file}")
print(f"\nAll figures saved to: {output_dir}")
print("\n" + "=" * 80)

