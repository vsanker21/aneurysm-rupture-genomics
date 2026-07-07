"""Run integration with explicit logging to file"""
import sys
import os

# Force UTF-8 encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

# Redirect output to file
log_file = open('integration_run.log', 'w', encoding='utf-8')
sys.stdout = sys.stderr = log_file

try:
    from integrate_real_data_multi_omics import main
    print("=" * 80)
    print("STARTING MULTI-OMICS INTEGRATION")
    print("=" * 80)
    main()
    print("\n[SUCCESS] Integration complete!")
except Exception as e:
    print(f"\n[ERROR] Integration failed: {e}")
    import traceback
    traceback.print_exc()
finally:
    log_file.close()

# Read and display results
with open('integration_run.log', 'r', encoding='utf-8') as f:
    content = f.read()
    print(content)


