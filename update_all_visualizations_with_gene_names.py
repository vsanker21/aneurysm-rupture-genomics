"""
Update All Visualizations to Use Gene Names Instead of IDs
This script patches the existing visualization script to ensure all figures use gene symbols
"""
import pandas as pd
import numpy as np
import os
import shutil
from pathlib import Path

print("=" * 80)
print("UPDATING ALL VISUALIZATIONS WITH GENE NAMES")
print("=" * 80)

# First, let's run the new visualization script that uses gene names
print("\nRunning updated visualization script with gene names...")

# Import and run the gene name visualization script
exec(open('generate_complete_visualizations_with_gene_names.py').read())

print("\n" + "=" * 80)
print("UPDATED VISUALIZATIONS COMPLETE!")
print("=" * 80)

