"""Test script to run analyses individually and capture errors"""
import sys
import traceback

# Load the main script's setup
exec(open('advanced_ml_ultra_hd_analyses.py').read()[:5000])  # Just get imports and setup

# Now try to run each analysis section
print("Testing Analysis 2...")
try:
    # Analysis 2 code here - let's import and test
    exec(open('advanced_ml_ultra_hd_analyses.py').read()[5000:10000])
except Exception as e:
    print(f"Error in Analysis 2: {e}")
    traceback.print_exc()

