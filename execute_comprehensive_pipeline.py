"""
Execute Comprehensive Analysis Pipeline - Step by step with explicit output
Ensures all features are implemented and executed successfully
"""
import sys
import os

# Force UTF-8
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

print("=" * 80)
print("EXECUTING COMPREHENSIVE ANALYSIS PIPELINE")
print("All Features: ROC Curves, Feature Importance, Pathway Visualization, Gene Lists, Risk Scores")
print("=" * 80, flush=True)

# Import and execute the main pipeline
try:
    print("\n[0] Importing pipeline script...", flush=True)
    
    # Read and execute the comprehensive pipeline
    with open('comprehensive_analysis_pipeline.py', 'r', encoding='utf-8', errors='replace') as f:
        pipeline_code = f.read()
    
    # Create a namespace for execution
    namespace = {
        '__name__': '__main__',
        '__file__': 'comprehensive_analysis_pipeline.py',
        '__builtins__': __builtins__
    }
    
    print("[1] Executing pipeline...", flush=True)
    exec(pipeline_code, namespace)
    
    print("\n" + "=" * 80)
    print("✅ PIPELINE EXECUTION COMPLETE")
    print("=" * 80)
    
except KeyboardInterrupt:
    print("\n⚠️  Pipeline interrupted by user")
    sys.exit(1)
except Exception as e:
    import traceback
    print(f"\n❌ ERROR: {e}")
    print("\nTraceback:")
    traceback.print_exc()
    sys.exit(1)

