# Remaining Datasets to Download

**Last Updated:** November 13, 2025

## Executive Summary

**Status:** ✅ **Most critical datasets are complete!**

- **Completed:** 9 major datasets (GWAS, blood RNA-seq, tissue RNA-seq, microarray, GTEx)
- **Remaining:** 3 optional/additional datasets
- **Not Available:** 1 dataset (single-cell RNA-seq - not found in public repositories)

---

## ✅ **COMPLETED DATASETS** (Ready for Analysis)

1. ✅ **Blood RNA-seq** - GSE159610 (47 samples)
2. ✅ **Tissue RNA-seq** - GSE122897, GSE158558, GSE66238 (ruptured vs unruptured)
3. ✅ **GWAS Summary Statistics** - Bakker et al. 2020 (10,754 cases, 306,882 controls)
4. ✅ **GTEx Reference** - All 3 vessel tissues (aorta, coronary, tibial)
5. ✅ **Microarray Data** - GSE13353, GSE54083, GSE75436

**Total:** 9 complete datasets, ~1.5 GB of data

---

## ⚠️ **REMAINING/OPTIONAL DATASETS**

### 1. **Additional GWAS Summary Statistics** (OPTIONAL)

#### Status: ⚠️ **PARTIALLY EXTRACTED FROM PDFs**

We've extracted some data from PDFs, but full summary statistics may be available:

#### a) **Foroud et al. 2014 (Stroke Journal)**
- **Status:** ⚠️ **Partially extracted** (4 SNPs from main paper)
- **What we have:**
  - ✅ 4 SNPs extracted from main paper PDF
  - ✅ Study metadata and sample sizes extracted
  - ⚠️ Full summary statistics may be in supplementary files
- **Location:** `data/supplementary_pdfs/extracted_gwas/foroud_2014_main_paper_gwas.tsv`
- **Action needed:**
  - Check journal supplement for full summary statistics table
  - Or download from: https://stroke.ahajournals.org/content/45/11/3194
  - May require manual download from supplementary materials

#### b) **Bilguvar et al. 2008 (Nature Genetics)**
- **Status:** ⚠️ **Partially extracted** (10 SNPs from supplement)
- **What we have:**
  - ✅ 10 SNPs extracted from supplementary PDF
  - ✅ Study metadata extracted
  - ⚠️ Full summary statistics may be available elsewhere
- **Location:** `data/supplementary_pdfs/extracted_gwas/NIHMS4348-supplement-s1_gwas_parsed.tsv`
- **Action needed:**
  - Check Nature Genetics supplement for full summary statistics
  - Or download from: https://www.nature.com/articles/ng.240
  - May require manual download from supplementary materials

**Priority:** 🟡 **LOW** - We have key SNPs extracted, full summary stats are optional

---

### 2. **UK Biobank GWAS Data** (OPTIONAL)

#### Status: 📋 **DOCUMENTATION ONLY**

- **Location:** `data/uk_biobank/`
- **Status:** ⚠️ **OPTIONAL** - Not downloaded
- **Why optional:**
  - Bakker et al. 2020 meta-analysis **already includes UK Biobank data**
  - May be redundant if you have Bakker summary statistics
  - Only needed if you want to analyze UK Biobank separately

**Download Sources:**
1. **Neale Lab** (Easiest): http://www.nealelab.is/uk-biobank/
2. **Pan-UK Biobank**: https://pan.ukbb.broadinstitute.org/
3. **UK Biobank Official**: https://www.ukbiobank.ac.uk/ (requires approval)

**Priority:** 🟡 **LOW** - Optional, may be redundant

---

### 3. **Poppenberg et al. 2023 Growth Rate Data** (UNCLEAR)

#### Status: ❓ **NEEDS VERIFICATION**

- **Mentioned in project docs:** 66-68 blood RNA-seq samples from aneurysm patients
- **Categories:** "Growing" vs "Stable" aneurysms
- **Status:** ❓ **Not clear if this is:**
  - Part of GSE159610 (already downloaded)
  - A separate dataset that needs to be downloaded
  - Data from a different publication

**Action needed:**
- Verify if this is part of GSE159610 or a separate dataset
- Check if it's available in GEO or requires author contact
- Search for: "Poppenberg" + "aneurysm" + "growth rate" in GEO

**Priority:** 🟡 **MEDIUM** - Would be valuable if separate dataset

---

## ❌ **DATASETS NOT AVAILABLE**

### 1. **Single-Cell RNA-seq of Aneurysm Tissue**

#### Status: ❌ **NOT FOUND IN PUBLIC REPOSITORIES**

- **Original target:** GSE186927 (was incorrect - mouse intestinal study)
- **Search result:** No human intracranial aneurysm single-cell RNA-seq dataset found in GEO
- **Impact:** Endothelial cell-specific analysis may need alternative approach

**Alternatives:**
1. Use bulk RNA-seq with cell-type deconvolution (GSE122897, GSE158558)
2. Contact authors of recent scRNA-seq studies
3. Search EGA, SRA repositories
4. Use mouse model data (GSE193533) for proof-of-concept

**Priority:** 🔴 **HIGH** (but not available) - Use alternatives

---

## 📊 **SUMMARY TABLE**

| Dataset | Status | Priority | Action Needed |
|---------|--------|----------|---------------|
| **Foroud 2014 GWAS** | ⚠️ Partial (4 SNPs) | 🟡 Low | Check journal supplement |
| **Bilguvar 2008 GWAS** | ⚠️ Partial (10 SNPs) | 🟡 Low | Check journal supplement |
| **UK Biobank** | 📋 Optional | 🟡 Low | Download if needed |
| **Poppenberg Growth Data** | ❓ Unclear | 🟡 Medium | Verify if separate dataset |
| **Single-cell RNA-seq** | ❌ Not found | 🔴 High | Use alternatives |

---

## 🎯 **RECOMMENDED ACTIONS**

### Immediate (Optional):
1. **Verify Poppenberg dataset:**
   - Check if 66-68 sample growth rate data is separate from GSE159610
   - Search GEO for additional Poppenberg datasets
   - Contact authors if needed

2. **Download full GWAS summary statistics (if needed):**
   - Foroud 2014: Check Stroke journal supplement
   - Bilguvar 2008: Check Nature Genetics supplement
   - Note: We already have key SNPs extracted from PDFs

### Low Priority:
3. **UK Biobank (if needed):**
   - Only if separate UK Biobank analysis is required
   - Otherwise, Bakker 2020 includes UK Biobank data

### Not Available (Use Alternatives):
4. **Single-cell RNA-seq:**
   - Use bulk RNA-seq with deconvolution
   - Contact authors of recent studies
   - Search EGA/SRA repositories

---

## ✅ **CONCLUSION**

**Good News:** All **critical datasets are complete!**

- ✅ GWAS data (Bakker 2020) - Complete
- ✅ Blood transcriptomics (GSE159610) - Complete
- ✅ Tissue transcriptomics (multiple datasets) - Complete
- ✅ Rupture comparisons (5 datasets) - Complete
- ✅ Reference data (GTEx) - Complete

**Remaining items are optional enhancements:**
- Additional GWAS studies (we have key SNPs extracted)
- UK Biobank (redundant with Bakker 2020)
- Poppenberg growth data (needs verification)

**Recommendation:** ✅ **PROCEED WITH ANALYSIS**

You have sufficient data to:
- Build polygenic risk scores
- Identify rupture prediction signatures
- Integrate multi-omics data
- Develop clinical prediction models

The remaining datasets are optional and can be added later if needed.

---

**Status:** ✅ **Ready for Analysis** - All critical datasets acquired!






