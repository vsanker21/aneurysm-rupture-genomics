"""
Run GSE122897 analysis with proper error handling
"""
import sys
import traceback

try:
    from analyze_gse122897 import main
    main()
except Exception as e:
    print(f"Error: {e}")
    traceback.print_exc()
    sys.exit(1)






