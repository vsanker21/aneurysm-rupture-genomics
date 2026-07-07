"""Generate Feature Importance Plot - Individual Script"""
import sys
import os
import traceback

if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

print("=" * 80)
print("GENERATING FEATURE IMPORTANCE PLOT")
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
    
    # Load feature importance
    print("\nLoading feature importance...")
    sys.stdout.flush()
    if os.path.exists('results/feature_importance.tsv'):
        feature_imp = pd.read_csv('results/feature_importance.tsv', sep='\t')
        print(f"✓ Loaded {len(feature_imp)} features")
        print(f"  Columns: {list(feature_imp.columns)}")
        sys.stdout.flush()
        
        # Create plot
        print("\nCreating plot...")
        sys.stdout.flush()
        fig, axes = plt.subplots(2, 2, figsize=(16, 14))
        axes = axes.flatten()
        
        models = [
            ('Tissue', 'Random Forest'),
            ('Tissue', 'Gradient Boosting'),
            ('Blood', 'Random Forest'),
            ('Blood', 'Gradient Boosting')
        ]
        
        for idx, (model_type, model_name) in enumerate(models):
            # Check if Model column exists
            if 'Model' in feature_imp.columns:
                model_data = feature_imp[
                    (feature_imp['Model'].str.contains(model_type, case=False, na=False)) & 
                    (feature_imp['Model'].str.contains(model_name, case=False, na=False))
                ]
            else:
                # If no Model column, use all features for visualization
                model_data = feature_imp.copy()
            
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
                print(f"  ✓ Added {model_type} {model_name}: {len(top_features)} features")
            else:
                axes[idx].text(0.5, 0.5, f'No data for\n{model_type} {model_name}', 
                             ha='center', va='center', transform=axes[idx].transAxes)
                print(f"  ⚠ No data for {model_type} {model_name}")
        
        plt.suptitle('Feature Importance Rankings - Top 20 Features', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        output_file = 'results/figures/feature_importance_all_models.png'
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()
        
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file)
            print(f"✓ SUCCESS: Created {output_file} ({file_size/1024:.2f} KB)")
        else:
            print(f"✗ ERROR: File not created: {output_file}")
        
        # Save TSV
        feature_imp.to_csv('results/feature_importance_all_models.tsv', sep='\t', index=False)
        print("✓ Saved: results/feature_importance_all_models.tsv")
        
    else:
        print("✗ ERROR: File not found: results/feature_importance.tsv")
        sys.exit(1)
    
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

