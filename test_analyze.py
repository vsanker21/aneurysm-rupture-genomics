"""Test script to debug analyze_gse122897.py"""
import sys
import traceback
import io

# Redirect stdout to capture output
old_stdout = sys.stdout
sys.stdout = mystdout = io.StringIO()

try:
    import analyze_gse122897
    analyze_gse122897.main()
    
    # Get the output
    output = mystdout.getvalue()
    print("="*80)
    print("CAPTURED OUTPUT:")
    print("="*80)
    print(output)
    print("="*80)
    
except Exception as e:
    print("="*80)
    print("ERROR OCCURRED:")
    print("="*80)
    print(str(e))
    traceback.print_exc()
finally:
    sys.stdout = old_stdout






