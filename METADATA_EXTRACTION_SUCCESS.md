# GSE122897 Metadata Extraction - SUCCESS! ✅

## Summary

Successfully extracted sample metadata from `GSE122897_series_matrix.txt.gz`!

## Results

### Sample Groups Identified:
- **22 Ruptured samples** - Ruptured intracranial aneurysm
- **21 Unruptured samples** - Unruptured intracranial aneurysm  
- **16 Control samples** - Intracranial cortical artery (normal tissue)
- **2 Unknown samples** - Unknown rupture status

**Total: 60 samples**

### Files Created:
1. `data/sc_rnaseq/GSE122897/sample_metadata.tsv` - Full metadata with all columns
2. `results/gse122897_sample_metadata.tsv` - Same metadata in results folder

### Metadata Columns:
- Sample (GSM IDs)
- Group (Ruptured/Unruptured/Control/Unknown)
- All other GEO metadata (characteristics, age, sex, etc.)

## Next Steps

The metadata is now ready to use! You can:

1. **Run the analysis script:**
   ```bash
   python analyze_gse122897.py
   ```
   This will:
   - Load the expression data
   - Load the metadata (automatically)
   - Perform quality control
   - Run differential expression analysis (Ruptured vs Unruptured)
   - Identify top rupture-associated genes

2. **The analysis will compare:**
   - 22 Ruptured samples vs 21 Unruptured samples
   - Control samples are excluded from rupture comparison (but available for other analyses)

## What Was Found

The series matrix file contained:
- Sample titles: CO1-CO16 (controls), IA1-IA44 (aneurysms)
- Tissue/status information: "Ruptured intracranial aneurysm", "Unruptured intracranial aneurysm", "Intracranial cortical artery"
- Additional metadata: age, sex, aneurysm location, etc.

The extraction script successfully:
- Parsed the GEO series matrix format
- Extracted all sample metadata
- Identified rupture status from "tissue/status" column
- Created proper Group labels (Ruptured/Unruptured/Control)

## Status

✅ **Metadata extraction: COMPLETE**
✅ **Files saved and ready for analysis**
✅ **Analysis script updated to use metadata automatically**

You're ready to run the differential expression analysis!






