"""
Run comprehensive analysis pipeline with validation and explicit output
Ensures all features are implemented successfully
"""
import sys
import os

# Force UTF-8 encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

print("=" * 80, flush=True)
print("RUNNING COMPREHENSIVE ANALYSIS PIPELINE", flush=True)
print("=" * 80, flush=True)

try:
    # Execute the main pipeline
    with open('comprehensive_analysis_pipeline.py', 'r', encoding='utf-8', errors='replace') as f:
        code = f.read()
        exec(code, globals())
    
    print("\n" + "=" * 80, flush=True)
    print("VALIDATING OUTPUTS...", flush=True)
    print("=" * 80, flush=True)
    
    # Check for expected output files
    expected_files = {
        'figures': [
            'results/figures/roc_curves_all_models.png',
            'results/figures/roc_curves_combined.png',
            'results/figures/feature_importance_all_models.png',
            'results/figures/pathway_enrichment_tissue.png',
            'results/figures/pathway_enrichment_blood.png',
            'results/figures/risk_stratification.png'
        ],
        'gene_lists': [
            'results/gene_lists/comprehensive_signature_genes_annotated.tsv',
            'results/gene_lists/comprehensive_signature_genes_annotated.xlsx',
            'results/gene_lists/tissue_top50_signature_genes.tsv',
            'results/gene_lists/blood_top50_signature_genes.tsv'
        ],
        'risk_scores': [
            'results/risk_scores/tissue_risk_scores_categorized.tsv',
            'results/risk_scores/blood_risk_scores_categorized.tsv'
        ],
        'other': [
            'results/feature_importance_all_models.tsv',
            'results/common_enriched_pathways.tsv',
            'results/ensemble_model_summary.tsv'
        ]
    }
    
    all_found = True
    for category, files in expected_files.items():
        print(f"\n{category.upper()}:", flush=True)
        for file_path in files:
            if os.path.exists(file_path):
                size = os.path.getsize(file_path)
                print(f"  ✓ {file_path} ({size:,} bytes)", flush=True)
            else:
                print(f"  ✗ {file_path} - MISSING", flush=True)
                all_found = False
    
    if all_found:
        print("\n" + "=" * 80, flush=True)
        print("✅ ALL OUTPUTS SUCCESSFULLY GENERATED!", flush=True)
        print("=" * 80, flush=True)
    else:
        print("\n" + "=" * 80, flush=True)
        print("⚠️  SOME OUTPUTS MISSING", flush=True)
        print("=" * 80, flush=True)
        
except Exception as e:
    import traceback
    print(f"\nERROR: {e}", file=sys.stderr, flush=True)
    traceback.print_exc()
    sys.exit(1)

