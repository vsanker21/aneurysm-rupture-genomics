"""
Wrapper script to run comprehensive multi-omics integration with proper output
"""
import sys
import os

# Force UTF-8 encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

# Run the comprehensive integration
if __name__ == "__main__":
    print("Starting comprehensive multi-omics integration...")
    print("=" * 80)
    
    try:
        from comprehensive_multi_omics_integration import main
        main()
    except Exception as e:
        import traceback
        print(f"ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)

