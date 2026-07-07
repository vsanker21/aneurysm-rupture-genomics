"""
Execute comprehensive multi-omics integration pipeline
Real data only - Gene ID conversion, Pathway analysis, Ensemble model
"""
import sys
import os

# Force UTF-8 encoding  
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

# Import and run
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import comprehensive_multi_omics_integration as cmi
    
    print("=" * 80, flush=True)
    print("RUNNING COMPREHENSIVE MULTI-OMICS INTEGRATION", flush=True)
    print("=" * 80, flush=True)
    
    cmi.main()
    
except Exception as e:
    import traceback
    print(f"\nERROR: {e}", file=sys.stderr, flush=True)
    traceback.print_exc()
    sys.exit(1)

