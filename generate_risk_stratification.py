"""Generate Risk Stratification Plot - Individual Script"""
import sys
import os
import traceback

if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

print("=" * 80)
print("GENERATING RISK STRATIFICATION PLOT")
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
    
    # Load predictions
    print("\nLoading predictions...")
    sys.stdout.flush()
    tissue_pred = pd.read_csv("results/ensemble_tissue_predictions.tsv", sep='\t')
    blood_pred = pd.read_csv("results/ensemble_blood_predictions.tsv", sep='\t')
    print(f"✓ Tissue: {len(tissue_pred)} samples")
    print(f"✓ Blood: {len(blood_pred)} samples")
    sys.stdout.flush()
    
    # Create risk categories
    print("\nCreating risk categories...")
    sys.stdout.flush()
    tissue_pred['Risk_Category'] = pd.cut(tissue_pred['Normalized_Score'], 
                                          bins=[0, 0.33, 0.67, 1.0], 
                                          labels=['Low', 'Medium', 'High'],
                                          include_lowest=True)
    blood_pred['Risk_Category'] = pd.cut(blood_pred['Normalized_Score'], 
                                         bins=[0, 0.33, 0.67, 1.0], 
                                         labels=['Low', 'Medium', 'High'],
                                         include_lowest=True)
    
    # Create plot
    print("\nCreating plot...")
    sys.stdout.flush()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    tissue_risk_counts = tissue_pred['Risk_Category'].value_counts().sort_index()
    axes[0].bar([str(x) for x in tissue_risk_counts.index], tissue_risk_counts.values, 
               color=['#2ca02c', '#ff7f0e', '#d62728'], alpha=0.8)
    axes[0].set_xlabel('Risk Category', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
    axes[0].set_title('Tissue Risk Stratification', fontsize=13, fontweight='bold')
    axes[0].grid(axis='y', alpha=0.3)
    
    blood_risk_counts = blood_pred['Risk_Category'].value_counts().sort_index()
    axes[1].bar([str(x) for x in blood_risk_counts.index], blood_risk_counts.values, 
               color=['#2ca02c', '#ff7f0e', '#d62728'], alpha=0.8)
    axes[1].set_xlabel('Risk Category', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
    axes[1].set_title('Blood Risk Stratification', fontsize=13, fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)
    
    plt.suptitle('Risk Score Stratification', fontsize=15, fontweight='bold')
    plt.tight_layout()
    
    output_file = 'results/figures/risk_stratification.png'
    plt.savefig(output_file, bbox_inches='tight', dpi=150)  # Reduced DPI to avoid memory error
    plt.close()
    
    if os.path.exists(output_file):
        file_size = os.path.getsize(output_file)
        print(f"✓ SUCCESS: Created {output_file} ({file_size/1024:.2f} KB)")
    else:
        print(f"✗ ERROR: File not created: {output_file}")
    
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

