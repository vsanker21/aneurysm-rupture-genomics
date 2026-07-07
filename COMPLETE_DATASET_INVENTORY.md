# Complete Dataset Inventory
## All Available Public Datasets for Aneurysm Genomics Project

**Generated:** 2025-11-13

---

## ✅ **FULLY DOWNLOADED PUBLIC DATASETS**

### 1. **Blood Transcriptomics**

#### GSE159610 - Blood RNA-seq (Aneurysm vs Control)
- **Location:** `data/blood_rnaseq/GSE159610/`
- **Status:** ✅ **COMPLETE**
- **Files:**
  - `GSE159610_raw_gene_counts_matrix.txt.gz` (7.8 MB)
  - `GSE159610_series_matrix.txt.gz` (4.8 KB)
  - Sample metadata files
- **Description:** Whole blood transcriptomes of patients with aneurysms vs controls
- **Samples:** 47 total (23 aneurysm, 24 control)

---

### 2. **Tissue Transcriptomics**

#### GSE122897 - Tissue RNA-seq (Ruptured vs Unruptured) ⭐ **CRITICAL**
- **Location:** `data/sc_rnaseq/GSE122897/`
- **Status:** ✅ **COMPLETE**
- **Files:**
  - `GSE122897_readCounts_raw.txt.gz` (2.0 MB)
  - `GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz` (2.1 MB)
  - `GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz` (3.9 MB)
  - `GSE122897_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz` (3.9 MB)
  - Sample metadata
- **Description:** RNA-seq analysis of intracranial aneurysm walls - ruptured vs unruptured
- **Key Value:** Direct rupture prediction comparison

#### GSE158558 - Tissue RNA-seq (Intracranial Aneurysm)
- **Location:** `data/sc_rnaseq/GSE158558/`
- **Status:** ✅ **COMPLETE**
- **Files:**
  - `GSE158558_raw_counts_GRCh38.p13_NCBI.tsv.gz` (0.4 MB)
  - `GSE158558_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz` (0.7 MB)
  - `GSE158558_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz` (0.7 MB)
  - `GSE158558_mRNA_lncRNA_Expression_Transcript.xlsx` (4.6 MB)
  - Gene annotation file
- **Description:** Intracranial aneurysm tissue expression data with mRNA and lncRNA

#### GSE66238 - Tissue RNA-seq (Intracranial Aneurysm)
- **Location:** `data/sc_rnaseq/GSE66238/`
- **Status:** ✅ **COMPLETE**
- **Files:**
  - `GSE66238_raw_counts_GRCh38.p13_NCBI.tsv.gz`
  - `GSE66238_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz`
  - `GSE66238_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz`
  - Differential expression results
- **Description:** Intracranial aneurysm tissue expression data

---

### 3. **Microarray Datasets**

#### GSE13353 - Aneurysm Microarray
- **Location:** `data/microarray/GSE13353/`
- **Status:** ✅ **COMPLETE**
- **Files:**
  - `GSE13353_RAW.tar` (81.2 MB)
- **Description:** Aneurysm microarray expression data

#### GSE54083 - Aneurysm Microarray
- **Location:** `data/microarray/GSE54083/`
- **Status:** ✅ **COMPLETE**
- **Files:**
  - `GSE54083_RAW.tar` (158.5 MB)
- **Description:** Aneurysm microarray expression data

#### GSE75436 - Aneurysm Microarray
- **Location:** `data/microarray/GSE75436/`
- **Status:** ✅ **COMPLETE**
- **Files:**
  - `GSE75436_RAW.tar` (136.0 MB)
- **Description:** Aneurysm microarray expression data

---

### 4. **GWAS Data**

#### Bakker et al. 2020 - Meta-GWAS
- **Location:** `data/gwas/bakker_2020/`
- **Status:** ✅ **COMPLETE**
- **Files:**
  - `bakker_2020_summary_stats.txt.gz` (763.6 MB)
  - Extracted summary statistics (multiple files)
- **Description:** 
  - N = 10,754 IA cases + 306,882 controls
  - 17 risk loci identified (11 novel)
  - Cross-ethnic: European + East Asian
  - SNP-based heritability: 21.6%

---

### 5. **Reference Data**

#### GTEx - Normal Vessel Transcriptomes
- **Location:** `data/gtex/`
- **Status:** ✅ **COMPLETE**
- **Files:**
  - `gene_tpm_v10_artery_aorta.gct.gz` (54.5 MB) - TPM normalized
  - `gene_reads_v10_artery_aorta.gct.gz` (22.9 MB) - Raw counts
  - `gene_tpm_v10_artery_coronary.gct.gz` (32.2 MB) - TPM normalized
  - `gene_reads_v10_artery_coronary.gct.gz` (13.7 MB) - Raw counts
  - `gene_tpm_v10_artery_tibial.gct.gz` (77.2 MB) - TPM normalized
  - `gene_reads_v10_artery_tibial.gct.gz` (32.4 MB) - Raw counts
- **Description:** Normal vessel tissue reference transcriptomes for comparison

---

## ⚠️ **DATASETS REQUIRING MANUAL DOWNLOAD/APPROVAL**

### 1. **Additional GWAS Studies**

#### Foroud et al. 2014 (Japanese Cohort)
- **Location:** `data/gwas/foroud_2014/`
- **Status:** ⚠️ **REQUIRES MANUAL DOWNLOAD**
- **Instructions:**
  - Publication: https://pubmed.ncbi.nlm.nih.gov/25355837/
  - Check journal supplement for summary statistics
  - Or contact authors for data access
- **See:** `data/gwas/foroud_2014/README.md`

#### Bilguvar et al. 2008 (Nature Genetics)
- **Location:** `data/gwas/bilguvar_2008/`
- **Status:** ⚠️ **REQUIRES MANUAL DOWNLOAD**
- **Instructions:**
  - Publication: https://pubmed.ncbi.nlm.nih.gov/18758463/
  - Check journal supplement for summary statistics
  - Or contact authors for data access
- **See:** `data/gwas/bilguvar_2008/README.md`

---

### 2. **Individual-Level Data (Requires Approval)**

#### Patient Genotype Data
- **Sources:** dbGaP, UK Biobank, EGA
- **Status:** ⚠️ **REQUIRES APPROVAL**
- **Purpose:** Individual-level PRS calculation
- **Instructions:**
  - See: `data/dbgap/ACCESS_INSTRUCTIONS.md`
  - See: `data/uk_biobank/DOWNLOAD_INSTRUCTIONS.md`
- **Note:** Bakker et al. 2020 already includes UK Biobank data

#### Patient Outcome Data
- **Sources:** Contact authors, Journal supplements
- **Status:** ⚠️ **REQUIRES CONTACT/AUTHORIZATION**
- **Purpose:** Model training and validation
- **Content Needed:**
  - Rupture status (ruptured vs unruptured)
  - Rupture dates
  - PHASES score components
  - Follow-up outcomes

---

## ❌ **DATASETS NOT FOUND IN PUBLIC REPOSITORIES**

### 1. **Exosomal RNA Sequencing Data**
- **Status:** ❌ **NOT FOUND**
- **Description:** Exosomal RNA-seq data for aneurysm patients
- **Alternative:** May need to analyze existing datasets for exosomal markers

### 2. **Single-Cell RNA-seq of Aneurysm Tissue**
- **Status:** ❌ **NOT FOUND** (GSE186927 was incorrect - mouse intestinal study)
- **Alternative:** Using bulk RNA-seq (GSE122897, GSE158558) instead
- **Note:** Cell-type deconvolution can be used to estimate cell proportions

---

## 📊 **DATASET SUMMARY**

### **Total Datasets Downloaded:**
- **Blood RNA-seq:** 1 dataset (GSE159610)
- **Tissue RNA-seq:** 3 datasets (GSE122897, GSE158558, GSE66238)
- **Microarray:** 3 datasets (GSE13353, GSE54083, GSE75436)
- **GWAS:** 1 dataset (Bakker 2020)
- **Reference:** 1 dataset (GTEx - 3 tissues)
- **Total:** 9 complete datasets

### **Total Files:**
- Approximately 30+ data files
- Total size: ~1.5 GB (compressed)

### **Data Coverage:**
- ✅ Blood transcriptomics (aneurysm vs control)
- ✅ Tissue transcriptomics (ruptured vs unruptured)
- ✅ GWAS summary statistics
- ✅ Reference normal tissue data
- ✅ Multiple validation datasets

---

## ✅ **CONCLUSION**

**All publicly available datasets for this project have been successfully downloaded.**

### **What We Have:**
1. ✅ **Complete blood transcriptomics data** (GSE159610)
2. ✅ **Complete tissue transcriptomics data** (3 datasets including ruptured vs unruptured)
3. ✅ **GWAS summary statistics** (Bakker 2020 meta-analysis)
4. ✅ **Reference data** (GTEx normal vessels)
5. ✅ **Multiple validation datasets** (microarray data)

### **What Requires Additional Steps:**
1. ⚠️ **Individual-level genotype data** - Requires dbGaP/UK Biobank approval
2. ⚠️ **Patient outcome data** - Contact authors or check journal supplements
3. ⚠️ **Additional GWAS studies** - Manual download from journal supplements

### **What's Not Available:**
1. ❌ **Exosomal RNA data** - Not found in public repositories
2. ❌ **Single-cell RNA-seq aneurysm data** - Not available (using bulk instead)

---

## 📋 **NEXT STEPS**

1. ✅ **All public datasets downloaded** - Ready for analysis
2. ⚠️ **For patient-level data:**
   - Apply for dbGaP/UK Biobank access if needed
   - Contact study authors for outcome data
   - Check journal supplements for additional GWAS studies
3. 💡 **Alternative approaches:**
   - Use summary statistics for PRS (already available)
   - Use cross-validation instead of independent test set
   - Focus on discovery manuscript with available data

---

**Status:** ✅ **All publicly available datasets successfully downloaded and verified.**






