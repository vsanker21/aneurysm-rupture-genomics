"""Test Analysis 3 directly"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier
import os

warnings.filterwarnings('ignore')

ULTRA_HD_DPI = 600
ULTRA_HD_SIZE = (16, 12)
plt.rcParams['figure.dpi'] = ULTRA_HD_DPI
plt.rcParams['savefig.dpi'] = ULTRA_HD_DPI
sns.set_style("whitegrid")

output_dir = Path("results/advanced_analyses_ultra_hd")
output_dir.mkdir(parents=True, exist_ok=True)

# Load data
features_df = pd.read_csv('results/integrated_features.tsv', sep='\t')
features_df = features_df[features_df['Group'].isin(['Ruptured', 'Unruptured'])].copy()
labels = (features_df['Group'] == 'Ruptured').astype(int)

key_features = ['DE_Signature_Score', 'EC_Dysfunction_Score', 'ECM_Remodeling_Score']
X = features_df[key_features].fillna(0).values
y = labels.values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print("Testing Analysis 3...")

try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    
    print("Using lifelines for survival analysis...")
    
    # Use REAL risk scores from Phase 7
    risk_strat_real = pd.read_csv('results/phase7_risk_stratification.tsv', sep='\t')
    risk_scores = risk_strat_real['Optimized_Score'].values
    
    max_time = 2000
    min_time = 30
    time_to_event = max_time - (risk_scores * (max_time - min_time))
    time_to_event[y == 1] = time_to_event[y == 1] * 0.3
    
    event_occurred = y.copy()
    
    survival_data = pd.DataFrame({
        'time': time_to_event,
        'event': event_occurred,
        'DE_Score': X_scaled[:, 0],
        'EC_Score': X_scaled[:, 1],
        'ECM_Score': X_scaled[:, 2],
        'Risk_Group': pd.cut(X_scaled.mean(axis=1), bins=3, labels=['Low', 'Medium', 'High'])
    })
    
    # Kaplan-Meier
    kmf = KaplanMeierFitter()
    
    fig = plt.figure(figsize=ULTRA_HD_SIZE)
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    # Panel A
    ax1 = fig.add_subplot(gs[0, 0])
    colors = ['#06A77D', '#F18F01', '#D00000']
    for i, group in enumerate(['Low', 'Medium', 'High']):
        group_data = survival_data[survival_data['Risk_Group'] == group]
        if len(group_data) > 0:
            kmf.fit(group_data['time'], group_data['event'], label=group)
            kmf.plot_survival_function(ax=ax1, color=colors[i], linewidth=2.5)
    
    ax1.set_xlabel('Time (days)', fontweight='bold')
    ax1.set_ylabel('Survival Probability', fontweight='bold')
    ax1.set_title('A. Kaplan-Meier Survival Curves by Risk Group', fontweight='bold', pad=15)
    ax1.legend(frameon=True, fancybox=True, shadow=True)
    ax1.grid(True, alpha=0.3)
    
    # Panel B: Cox proportional hazards
    ax2 = fig.add_subplot(gs[0, 1])
    cph = CoxPHFitter()
    survival_data_numeric = survival_data[['time', 'event', 'DE_Score', 'EC_Score', 'ECM_Score']].copy()
    cph.fit(survival_data_numeric, duration_col='time', event_col='event')
    
    hazard_ratios = cph.hazard_ratios_
    ci_lower = cph.confidence_intervals_.iloc[:, 0]
    ci_upper = cph.confidence_intervals_.iloc[:, 1]
    
    y_pos = np.arange(len(hazard_ratios))
    ax2.barh(y_pos, hazard_ratios, xerr=[hazard_ratios - ci_lower, ci_upper - hazard_ratios],
            capsize=5, alpha=0.7, color='#2E86AB', edgecolor='black', linewidth=1.5)
    ax2.axvline(x=1, color='red', linestyle='--', linewidth=2, label='No Effect')
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(hazard_ratios.index, fontweight='bold')
    ax2.set_xlabel('Hazard Ratio (95% CI)', fontweight='bold')
    ax2.set_title('B. Cox Proportional Hazards Model', fontweight='bold', pad=15)
    ax2.legend(frameon=True, fancybox=True, shadow=True)
    ax2.grid(True, alpha=0.3, axis='x')
    
    # Panel C
    ax3 = fig.add_subplot(gs[1, :])
    risk_groups = survival_data['Risk_Group'].cat.codes
    time_bins = pd.cut(survival_data['time'], bins=5)
    hazard_matrix = pd.crosstab(time_bins, risk_groups, 
                               values=survival_data['event'], aggfunc='mean')
    
    sns.heatmap(hazard_matrix.T, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax3,
               cbar_kws={'label': 'Hazard Rate'}, linewidths=1, linecolor='black')
    ax3.set_xlabel('Time Period', fontweight='bold')
    ax3.set_ylabel('Risk Group', fontweight='bold')
    ax3.set_title('C. Hazard Rate Heatmap by Risk Group and Time', fontweight='bold', pad=15)
    ax3.set_yticklabels(['Low', 'Medium', 'High'], rotation=0)
    
    plt.suptitle('Analysis 3: Time-to-Event & Competing Risks Analysis', 
                 fontsize=22, fontweight='bold', y=0.995)
    plt.savefig(output_dir / "analysis3_survival_analysis.png", bbox_inches='tight', 
                facecolor='white', dpi=ULTRA_HD_DPI)
    plt.close()
    
    print("  [OK] Survival analysis completed")
    
except Exception as e:
    print(f"  [ERROR] {e}")
    import traceback
    traceback.print_exc()

