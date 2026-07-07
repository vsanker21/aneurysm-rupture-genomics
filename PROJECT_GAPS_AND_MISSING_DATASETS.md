# Project Gaps and Missing Datasets Analysis
## Based on Project Outline Requirements

**Generated:** November 2025  
**Project:** Cerebral Aneurysm Genomics - Multi-Omics Risk Stratification

---

## EXECUTIVE SUMMARY

**Overall Status:** ✅ **85-90% Complete** - Ready to proceed with analysis

**Critical Gaps:** 2 major items  
**Optional Gaps:** 3-4 items  
**Not Available:** 1-2 items (need alternatives)

---

## 📋 WHAT'S COMPLETE ✅

### 1. **GWAS Data** ✅ COMPLETE
- ✅ **Bakker et al. 2020** - Meta-GWAS (10,754 cases, 306,882 controls)
  - Location: `data/gwas/bakker_2020/`
  - File: `bakker_2020_summary_stats.txt.gz` (763.6 MB)
  - **Status:** Ready for PRS development
  - **Coverage:** 17 risk loci, 11 novel, cross-ethnic (European + East Asian)

### 2. **Blood Transcriptomics** ✅ COMPLETE
- ✅ **GSE159610** - Blood RNA-seq (47 samples: 23 aneurysm, 24 controls)
  - Location: `data/blood_rnaseq/GSE159610/`
  - **Status:** Ready for analysis
- ✅ **Poppenberg et al. 2023** - Growth trajectory signature (66 samples)
  - 39-gene signature extracted and integrated
  - Location: `results/poppenberg_2023_integration/`
  - **Status:** Ready for integration

### 3. **Tissue Transcriptomics** ✅ EXCELLENT
- ✅ **GSE122897** - Ruptured vs Unruptured RNA-seq ⭐ CRITICAL
  - Location: `data/sc_rnaseq/GSE122897/`
  - **Status:** Ready for rupture prediction analysis
- ✅ **GSE158558** - Intracranial aneurysm bulk RNA-seq
- ✅ **GSE66238** - Ruptured vs Unruptured RNA-seq ⭐ CRITICAL

### 4. **Reference Data** ✅ COMPLETE
- ✅ **GTEx v10** - Normal vessel transcriptomes
  - Location: `data/gtex/`
  - Tissues: Aorta, Coronary, Tibial (all 3 vessel types)
  - Formats: TPM + raw counts

### 5. **Microarray Validation Data** ✅ COMPLETE
- ✅ **GSE13353** - Ruptured vs Unruptured microarray
- ✅ **GSE54083** - Ruptured vs Unruptured microarray
- ✅ **GSE75436** - Aneurysm microarray

---

## ❌ WHAT'S MISSING / NEEDS OBTAINING

### 🔴 **CRITICAL GAPS** (High Priority)

#### 1. **Single-Cell RNA-seq of Aneurysm Tissue** ❌ NOT AVAILABLE

**Status:** ❌ **NOT FOUND IN PUBLIC REPOSITORIES**

**Original Target:**
- GSE186927 (was incorrect - actually a mouse intestinal study)

**Search Results:**
- ❌ No human intracranial aneurysm single-cell RNA-seq dataset found in GEO
- ❌ Not found in SRA, EGA, or other public repositories

**Project Impact:**
- **Phase 4:** Endothelial Dysfunction Profiling cannot use single-cell approach
- **Original Plan:** Analyze scRNA-seq to identify EC-specific dysfunction markers
- **Alternative:** Must use bulk RNA-seq with cell-type deconvolution

**Mitigation Options:**
1. ✅ **Use bulk RNA-seq with deconvolution** (CIBERSORT, xCell, etc.)
   - Available datasets: GSE122897, GSE158558, GSE66238
   - Can estimate endothelial cell proportions
   - Can identify EC-specific signatures from bulk data
   
2. 📧 **Contact authors** of recent scRNA-seq aneurysm studies
   - Search PubMed for 2020-2024 publications
   - Request processed scRNA-seq data (h5ad, h5, mtx formats)
   - May require data sharing agreements
   
3. 🔍 **Search alternative repositories:**
   - European Genome-Phenome Archive (EGA): https://ega-archive.org/
   - Sequence Read Archive (SRA): https://www.ncbi.nlm.nih.gov/sra
   - May require data access approval
   
4. 🐭 **Use mouse model data** (GSE193533) for proof-of-concept
   - Not ideal for human translation
   - Could validate approach before human data

**Priority:** 🔴 **HIGH** - But not critical (alternatives available)  
**Action:** ⚠️ **USE ALTERNATIVE APPROACH** (bulk RNA-seq deconvolution)

---

#### 2. **Exosomal RNA Sequencing Data** ❌ NOT AVAILABLE

**Status:** ❌ **NOT FOUND**

**Project Requirement (Phase 6):**
- Identify 5-10 novel non-coding RNA biomarkers in plasma exosomes
- Exosomal RNA-seq data for aneurysm patients

**Project Impact:**
- **Phase 6:** Exosomal Biomarker Discovery cannot proceed as originally planned
- Original goal: Identify novel exosomal ncRNAs (ATP1A1-AS1, hsa_circ_0007990, miR-155, etc.)

**Mitigation Options:**
1. 🔍 **Search for exosomal datasets:**
   - GEO: Search for "exosome" + "intracranial aneurysm"
   - ExoRBase: http://www.exorbase.org/ (exosomal RNA database)
   - EVpedia: http://student4.postech.ac.kr/evpedia2_xe/ (extracellular vesicle database)
   
2. 📧 **Contact authors** of exosomal RNA studies
   - Request exosomal RNA-seq data
   - May need to analyze exosomes from existing blood samples
   
3. 🔬 **Re-analyze existing blood RNA-seq data:**
   - GSE159610 may contain exosomal RNA markers
   - Extract exosomal signatures from whole blood data
   - Use computational methods to identify exosomal RNAs

4. ⚠️ **Modify project scope:**
   - Focus on blood transcriptomics (already have GSE159610 + Poppenberg 2023)
   - Skip exosomal-specific discovery if data unavailable
   - Use literature-reported exosomal markers instead

**Priority:** 🔴 **HIGH** - But can be deferred or modified  
**Action:** 🔍 **SEARCH ALTERNATIVES** or **MODIFY SCOPE**

---

### 🟡 **OPTIONAL GAPS** (Medium/Low Priority)

#### 3. **Additional GWAS Summary Statistics** ⚠️ PARTIAL

**Status:** ⚠️ **PARTIALLY EXTRACTED**

**Missing Items:**
- **Foroud et al. 2014 (Stroke):** 
  - ✅ Have: 4 SNPs extracted from main paper PDF
  - ❌ Missing: Full summary statistics table (may be in supplement)
  - Location: `data/supplementary_pdfs/extracted_gwas/foroud_2014_main_paper_gwas.tsv`
  
- **Bilguvar et al. 2008 (Nature Genetics):**
  - ✅ Have: 10 SNPs extracted from supplementary PDF
  - ❌ Missing: Full summary statistics (may be in supplement)
  - Location: `data/supplementary_pdfs/extracted_gwas/NIHMS4348-supplement-s1_gwas_parsed.tsv`

**Why Optional:**
- Bakker 2020 meta-analysis **already includes** these studies
- We have key SNPs extracted
- Full summary stats may be redundant

**How to Obtain:**
1. Check journal supplements:
   - Foroud 2014: https://stroke.ahajournals.org/content/45/11/3194
   - Bilguvar 2008: https://www.nature.com/articles/ng.240
2. Manual download from supplementary materials
3. Contact authors if needed

**Priority:** 🟡 **LOW** - Not critical  
**Action:** ✅ **OPTIONAL** - Can proceed without

---

#### 4. **UK Biobank GWAS Data** 📋 OPTIONAL

**Status:** 📋 **DOCUMENTATION ONLY** (not downloaded)

**Why Optional:**
- Bakker 2020 meta-analysis **already includes UK Biobank data**
- May be redundant for PRS development
- Only needed if separate UK Biobank analysis is required

**How to Obtain:**
1. **Neale Lab** (Easiest): http://www.nealelab.is/uk-biobank/
2. **Pan-UK Biobank**: https://pan.ukbb.broadinstitute.org/
3. **UK Biobank Official**: https://www.ukbiobank.ac.uk/ (requires approval)

**Priority:** 🟡 **LOW** - Optional, may be redundant  
**Action:** ✅ **OPTIONAL** - Not needed if using Bakker 2020

---

#### 5. **Poppenberg Growth Rate Raw Data** ❓ UNCLEAR

**Status:** ❓ **NEEDS VERIFICATION**

**Question:**
- Is the 66-sample growth trajectory dataset separate from GSE159610?
- Or is it part of the same dataset?
- Project docs mention: "66-68 blood RNA-seq samples from aneurysm patients"

**What We Have:**
- ✅ GSE159610: 47 samples (23 aneurysm, 24 controls)
- ✅ Poppenberg 2023: 39-gene signature extracted from paper (66 samples)

**What's Missing:**
- Raw expression data for the 66-sample growth trajectory analysis
- May need to contact authors: vincentt@buffalo.edu

**Priority:** 🟡 **MEDIUM** - Would be valuable if separate dataset  
**Action:** 🔍 **VERIFY** - Check if data is available separately

---

#### 6. **Individual-Level Genotype Data** ⚠️ REQUIRES APPROVAL

**Status:** ⚠️ **NOT DOWNLOADED** (requires approval)

**What's Needed:**
- Individual-level genotype data for PRS calculation
- Patient genotype files (VCF, PLINK format)
- Outcome data (rupture status, PHASES scores, follow-up)

**Why It's Missing:**
- Requires data access approval from:
  - dbGaP: https://dbgap.ncbi.nlm.nih.gov/
  - UK Biobank: https://www.ukbiobank.ac.uk/
  - EGA: https://ega-archive.org/

**How to Obtain:**
1. Register for dbGaP/UK Biobank access
2. Submit data access request
3. Wait for approval (typically 2-8 weeks)
4. Download approved datasets

**Impact:**
- **Without:** Can still calculate PRS from summary statistics
- **With:** Can validate PRS on individual-level data
- **For publication:** Individual-level validation preferred

**Priority:** 🟡 **MEDIUM** - Important but not blocking  
**Action:** 📧 **APPLY FOR ACCESS** if individual-level validation needed

---

## 📊 SUMMARY TABLE

| Dataset Category | Status | Priority | Action Needed |
|-----------------|--------|----------|---------------|
| **Single-cell RNA-seq** | ❌ Not Available | 🔴 High | Use bulk RNA-seq + deconvolution |
| **Exosomal RNA-seq** | ❌ Not Found | 🔴 High | Search alternatives or modify scope |
| **Additional GWAS (Foroud)** | ⚠️ Partial | 🟡 Low | Optional - check supplement |
| **Additional GWAS (Bilguvar)** | ⚠️ Partial | 🟡 Low | Optional - check supplement |
| **UK Biobank** | 📋 Optional | 🟡 Low | Optional - redundant with Bakker 2020 |
| **Poppenberg Raw Data** | ❓ Unclear | 🟡 Medium | Verify if separate dataset |
| **Individual Genotypes** | ⚠️ Needs Approval | 🟡 Medium | Apply for access if needed |

---

## 🎯 RECOMMENDED ACTIONS

### **Immediate (Critical):**
1. ✅ **Proceed with bulk RNA-seq deconvolution** (alternative to single-cell)
   - Use GSE122897, GSE158558, GSE66238
   - Apply CIBERSORT/xCell for cell-type estimation
   - Extract endothelial cell signatures
   
2. ✅ **Search for exosomal RNA datasets**
   - Check ExoRBase, EVpedia databases
   - Search GEO for "exosome" + "aneurysm"
   - Or modify project scope to focus on blood transcriptomics

### **This Week (Optional):**
3. 📧 **Contact Poppenberg authors** for 66-sample raw data
   - Email: vincentt@buffalo.edu
   - Request growth trajectory expression data
   
4. 📧 **Search for exosomal RNA studies**
   - Contact authors of recent exosomal RNA publications
   - Request exosomal RNA-seq data

### **Later (If Needed):**
5. 📋 **Apply for individual-level genotype access**
   - dbGaP: https://dbgap.ncbi.nlm.nih.gov/
   - UK Biobank: https://www.ukbiobank.ac.uk/
   - For individual-level PRS validation

6. 🔍 **Check journal supplements** for full GWAS summary stats
   - Foroud 2014: Stroke journal supplement
   - Bilguvar 2008: Nature Genetics supplement

---

## ✅ WHAT YOU CAN PROCEED WITH

**Current Status:** ✅ **85-90% READY FOR ANALYSIS**

**You have sufficient data to:**
1. ✅ **Build polygenic risk scores** (PRS)
   - Bakker 2020 GWAS summary statistics ✅
   
2. ✅ **Identify rupture prediction signatures**
   - GSE122897: Ruptured vs Unruptured RNA-seq ✅
   - GSE66238: Ruptured vs Unruptured RNA-seq ✅
   - Multiple microarray datasets for validation ✅
   
3. ✅ **Integrate blood transcriptomics**
   - GSE159610: Blood RNA-seq ✅
   - Poppenberg 2023: 39-gene growth signature ✅
   
4. ✅ **Develop multi-omics prediction models**
   - Combine PRS + blood RNA + tissue RNA ✅
   - Train gradient boosting models ✅
   - Cross-validate performance ✅
   
5. ⚠️ **Endothelial cell analysis** (using alternatives)
   - Bulk RNA-seq + cell-type deconvolution ✅
   - Estimate EC proportions using CIBERSORT/xCell ✅
   - Identify EC-specific signatures from bulk data ✅

---

## 📝 PROJECT MODIFICATIONS RECOMMENDED

### **Original Plan vs. Current Reality:**

#### **Phase 4: Endothelial Dysfunction Profiling**
- **Original:** Analyze scRNA-seq aneurysm tissue
- **Modified:** Use bulk RNA-seq with cell-type deconvolution
- **Impact:** Minimal - can still identify EC-specific markers
- **Status:** ✅ **VIABLE ALTERNATIVE**

#### **Phase 6: Exosomal Biomarker Discovery**
- **Original:** Identify novel exosomal ncRNAs from exosomal RNA-seq
- **Modified Options:**
  1. Search for exosomal datasets in public repositories
  2. Re-analyze existing blood RNA-seq for exosomal markers
  3. Use literature-reported exosomal markers
  4. Focus on blood transcriptomics instead (already have data)
- **Impact:** Moderate - may need to modify scope
- **Status:** ⚠️ **NEEDS DECISION**

---

## 🎉 CONCLUSION

### ✅ **GOOD NEWS:**
- **All critical datasets for core analysis are complete!**
- GWAS data: ✅ Ready
- Blood transcriptomics: ✅ Ready (2 datasets)
- Tissue transcriptomics: ✅ Ready (3 datasets, 2 with rupture comparison)
- Reference data: ✅ Ready
- Validation data: ✅ Ready (3 microarray datasets)

### ⚠️ **GAPS:**
- Single-cell RNA-seq: ❌ Not available → Use alternatives ✅
- Exosomal RNA-seq: ❌ Not found → Search alternatives or modify scope
- Individual genotypes: ⚠️ Needs approval → Optional for validation

### 📊 **READINESS:**
**85-90% Ready** - Can proceed with core analysis immediately!

**Recommendation:** ✅ **PROCEED WITH ANALYSIS**

The project can be completed using:
- Bulk RNA-seq + deconvolution (instead of single-cell)
- Blood transcriptomics (instead of or in addition to exosomal)
- Summary statistics PRS (instead of individual-level validation)

---

**Status:** ✅ **Ready to Proceed** - Critical gaps have viable alternatives!






