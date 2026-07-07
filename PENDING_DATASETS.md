# Pending Datasets Summary

**Last Updated:** Current Session

## Overview

| Status | Count | Datasets |
|--------|-------|----------|
| ✅ Complete | 3 | Blood RNA-seq, GWAS, GTEx |
| ⚠️ Partial | 1 | Single-cell RNA-seq (GSE186927) |
| 📋 Optional | 1 | UK Biobank |

---

## ⚠️ PENDING DATASETS

### 1. Single-cell RNA-seq - Intracranial Aneurysm (NOT FOUND)

**Status:** ❌ **ERROR CORRECTED** - GSE186927 was incorrect

**⚠️ CRITICAL CORRECTION:**
- **GSE186927 is WRONG** - This is a mouse intestinal study, NOT intracranial aneurysm
- **No human intracranial aneurysm single-cell RNA-seq dataset found** in GEO after comprehensive search

**What Was Searched:**
- Comprehensive GEO database search
- Multiple candidate datasets verified
- Result: No matching human single-cell dataset found

**Alternative Sources:**

#### Option 1: European Genome-Phenome Archive (EGA) - RECOMMENDED
- **URL:** https://ega-archive.org/
- **Search:** "intracranial aneurysm" + "single cell"
- **Process:** Registration → Application → Approval → Download
- **Note:** May require data access approval (can take weeks)

#### Option 2: Sequence Read Archive (SRA)
- **URL:** https://www.ncbi.nlm.nih.gov/sra
- **Search:** "intracranial aneurysm" + "single cell RNA-seq"
- May contain raw sequencing data

#### Option 3: Publication Supplementary Materials
- Check recent publications (2020-2024) on intracranial aneurysm single-cell studies
- Download processed data files (h5ad, h5, mtx) from journal supplements
- Look for authors' GitHub repositories

#### Option 4: Contact Authors Directly
- Reach out to corresponding authors of recent aneurysm single-cell publications
- Request processed data files
- Check institutional repositories

#### Option 5: Alternative Datasets
- **GSE158558:** Human intracranial aneurysm bulk RNA-seq (alternative)
- **GSE193533:** Mouse cerebral aneurysm single-cell (proof-of-concept)
- Other vascular disease single-cell datasets

**Priority:** 🔴 **HIGH** - Critical for single-cell analysis, but dataset needs to be identified first

---

### 2. UK Biobank GWAS Summary Statistics (OPTIONAL)

**Status:** 📋 **DOCUMENTATION ONLY** - Not downloaded

**Location:** `data/uk_biobank/`

**Why It's Optional:**
- The Bakker et al. 2020 meta-analysis **already includes UK Biobank data**
- May be redundant if you have Bakker summary statistics
- Only needed if you want to analyze UK Biobank separately

**Download Sources:**

1. **Neale Lab (Easiest):**
   - URL: http://www.nealelab.is/uk-biobank/
   - Browse phenotypes
   - Search for "intracranial aneurysm" or vascular phenotypes
   - Download summary statistics

2. **Pan-UK Biobank:**
   - URL: https://pan.ukbb.broadinstitute.org/
   - Download GWAS summary statistics
   - Filter by phenotype

3. **UK Biobank Official:**
   - URL: https://www.ukbiobank.ac.uk/
   - Requires application and approval
   - More comprehensive but requires access approval

**Priority:** 🟡 **LOW** - Optional, may be redundant

---

## ✅ COMPLETED DATASETS

1. **Blood RNA-seq (GSE159610)** - ✅ Complete
   - 2 files downloaded
   - Ready for analysis

2. **GWAS (Bakker et al. 2020)** - ✅ Complete
   - 763.6 MB summary statistics
   - Ready for PRS development

3. **GTEx Vessel Tissues** - ✅ Complete
   - 6 files (3 tissues × 2 formats)
   - All required vessel tissues present

---

## Priority Actions

### 🔴 HIGH PRIORITY

1. **Download GSE186927 raw counts**
   - **Method:** R GEOquery (recommended)
   - **Time:** ~30-60 minutes (download + processing)
   - **Impact:** Critical for single-cell analysis
   - **Blocking:** Yes - needed for endothelial cell analysis

### 🟡 LOW PRIORITY

2. **UK Biobank data** (if needed)
   - **Method:** Neale Lab or Pan-UK Biobank
   - **Time:** ~1-2 hours (if downloading)
   - **Impact:** Optional, may be redundant
   - **Blocking:** No - Bakker data may be sufficient

---

## Overall Progress

**Completed:** 3/5 datasets (60%)
- Blood RNA-seq ✅
- GWAS ✅
- GTEx ✅

**Pending:** 2/5 datasets (40%)
- GSE186927 raw counts ⚠️ (Critical)
- UK Biobank 📋 (Optional)

**Critical Pending:** 1 dataset
- GSE186927 raw counts

---

## Next Steps

1. **Set up R environment** (if not already done)
   - Install R: https://www.r-project.org/
   - Install BiocManager and GEOquery

2. **Download GSE186927**
   - Run R script to download expression data
   - Save to `data/sc_rnaseq/GSE186927/`

3. **Verify data integrity**
   - Check file sizes
   - Verify expression matrix dimensions
   - Confirm cell metadata is present

4. **Optional: UK Biobank**
   - Only if separate UK Biobank analysis is needed
   - Otherwise, Bakker data is sufficient

---

## Helpful Scripts

- `download_remaining_datasets.py` - Attempts alternative download methods
- `check_gtex_downloads.py` - Verifies downloaded files
- R script can be created: `data/sc_rnaseq/GSE186927/download_gse186927.R`

---

## Notes

- GSE186927 is the **only critical pending dataset**
- UK Biobank is optional and may be redundant
- All other datasets are complete and ready for analysis
- R environment setup is the main requirement for GSE186927

