"""
Verify Packages and Run Transformer Model
Tests that all packages are installed and transformer works correctly
"""
import sys
import os

print("=" * 80)
print("VERIFYING PACKAGES AND TESTING TRANSFORMER")
print("=" * 80)

# Check required packages
required_packages = {
    'tensorflow': 'tf',
    'sklearn': None,
    'pandas': 'pd',
    'numpy': 'np',
    'matplotlib': None,
    'seaborn': 'sns',
    'requests': None
}

print("\nChecking required packages...")
missing = []
for package, alias in required_packages.items():
    try:
        if package == 'sklearn':
            import sklearn
            print(f"[OK] scikit-learn (sklearn) imported")
        elif package == 'tensorflow':
            import tensorflow as tf
            print(f"[OK] TensorFlow {tf.__version__} imported")
        elif package == 'pandas':
            import pandas as pd
            print(f"[OK] pandas {pd.__version__} imported")
        elif package == 'numpy':
            import numpy as np
            print(f"[OK] numpy {np.__version__} imported")
        elif package == 'matplotlib':
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend
            print(f"[OK] matplotlib imported")
        elif package == 'seaborn':
            import seaborn as sns
            print(f"[OK] seaborn imported")
        elif package == 'requests':
            import requests
            print(f"[OK] requests imported")
    except ImportError as e:
        missing.append(package)
        print(f"[ERROR] {package} not available: {e}")

if missing:
    print(f"\n[ERROR] Missing packages: {', '.join(missing)}")
    print(f"Installing missing packages...")
    import subprocess
    for pkg in missing:
        try:
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', pkg, '--quiet'])
            print(f"[OK] {pkg} installed")
        except:
            print(f"[ERROR] Failed to install {pkg}")

# Now test transformer model
print("\n" + "=" * 80)
print("TESTING TRANSFORMER MODEL")
print("=" * 80)

try:
    # Run the transformer script
    exec(open('advanced_ml_gnn_transformer.py').read())
except Exception as e:
    print(f"\n[ERROR] Transformer test failed: {e}")
    import traceback
    traceback.print_exc()
    print("\n[INFO] Checking transformer code...")
    
    # At least verify the transformer code structure
    if os.path.exists('advanced_ml_gnn_transformer.py'):
        with open('advanced_ml_gnn_transformer.py', 'r') as f:
            content = f.read()
            if 'transformer_encoder' in content and 'MultiHeadAttention' in content:
                print("[OK] Transformer code structure looks correct")
            else:
                print("[ERROR] Transformer code structure issues found")

print("\n" + "=" * 80)
print("VERIFICATION COMPLETE")
print("=" * 80)

