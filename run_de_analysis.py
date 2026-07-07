"""
Wrapper script to run GSE122897 differential expression analysis
and save output to a file
"""
import sys
import os

# Redirect output to file
log_file = open('results/gse122897_analysis_output.txt', 'w', encoding='utf-8')
sys.stdout = log_file
sys.stderr = log_file

try:
    print("Starting GSE122897 differential expression analysis...")
    print("=" * 80)
    
    # Import and run
    import analyze_gse122897
    analyze_gse122897.main()
    
    print("\nAnalysis script completed successfully!")
    
except Exception as e:
    print(f"\nERROR: {str(e)}")
    import traceback
    traceback.print_exc()
finally:
    log_file.close()
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    
    # Print summary
    print("\nAnalysis run complete. Check results/gse122897_analysis_output.txt for full output.")
    
    # Check for result files
    result_files = [
        'results/gse122897_differential_expression.tsv',
        'results/gse122897_top_genes.tsv'
    ]
    
    for f in result_files:
        if os.path.exists(f):
            size = os.path.getsize(f)
            print(f"  ✓ {f} ({size:,} bytes)")
        else:
            print(f"  ✗ {f} (not found)")






