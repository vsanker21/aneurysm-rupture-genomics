# Multi-Omics Dataset Summary for Integration
## Quick Reference Guide

**Generated:** 2025-01-XX

---

## ✅ Currently Available Datasets

### 1. Blood Transcriptomics
- **GSE159610**: 47 samples (23 aneurysm, 24 control)
- **Status**: ✅ Downloaded
- **Model**: AUC 0.968 (blood-only model)
- **Location**: `data/blood_rnaseq/GSE159610/`

### 2. Tissue Transcriptomics  
- **GSE122897**: 35 samples (17 ruptured, 18 unruptured) ⭐ **PRIMARY**
- **GSE66238**: Ruptured vs unruptured + microRNA
- **GSE158558**: mRNA and lncRNA expression
- **Status**: ✅ Downloaded
- **Location**: `data/sc_rnaseq/`

### 3. GWAS Data
- **Bakker 2020**: 10,754 cases + 306,882 controls
- **Status**: ✅ Downloaded (summary statistics only)
- **PRS Weights**: Available in `results/prs_weights.tsv`
- **Location**: `data/gwas/bakker_2020/`

---

## 📥 Datasets That Can Be Obtained

### Priority 1: Contact Authors (Highest Priority) 🔥

#### **Dr. Ynte Ruigrok** (GSE122897 - Tissue)
- **Email**: ij.m.ruigrok@umcutrecht.nl
- **Institution**: UMC Utrecht, Netherlands
- **Request**: Paired blood samples or genotype data from same patients
- **Benefit**: Same patients across omics layers = true multi-omics integration
- **Timeline**: 1-2 weeks for response

#### **Dr. Poppenberg** (GSE159610 - Blood)
- **Contact**: Check GEO metadata
- **Request**: Paired tissue samples or genotype data from same patients
- **Benefit**: Same patients across omics layers = true multi-omics integration
- **Timeline**: 1-2 weeks for response

### Priority 2: Apply for dbGaP Access 🔥

- **URL**: https://www.ncbi.nlm.nih.gov/gap/
- **Search**: "intracranial aneurysm" or "cerebral aneurysm"
- **What**: Individual-level genotype data
- **Benefit**: Patient-level PRS calculation
- **Timeline**: 2-4 weeks for approval
- **Requirements**: PI account, IRB approval, DUA

### Priority 3: Download Additional Datasets 📥

#### **E-GEOD-15629** (ArrayExpress)
- **Type**: Tissue microarray (19 samples)
- **Purpose**: External validation
- **URL**: https://www.omicsdi.org/dataset/biostudies-arrayexpress/E-GEOD-15629
- **Status**: ⚠️ Not yet downloaded
- **Action**: Download immediately for validation

---

## 🎯 Integration Strategies

### Current Capability (Available Now)
1. **Separate Models**:
   - Blood Model: Aneurysm vs Control (AUC 0.968)
   - Tissue Model: Ruptured vs Unruptured (being built)
   
2. **Feature-Level Integration**:
   - Common signature genes (blood + tissue)
   - Cross-validation across datasets

### With Genotype Data
- **Blood + PRS**: Expected AUC 0.95-0.98
- **Tissue + PRS**: Expected AUC 0.75-0.85

### With Paired Samples (True Multi-Omics)
- **Integrated Model**: Expected AUC 0.85-0.95
- **Features**: PRS + Blood RNA + Tissue RNA + Clinical
- **Superior**: To any single-omics model

---

## 📋 Action Items

### Immediate (This Week)
1. ✅ Build tissue model (GSE122897 ruptured vs unruptured)
2. 📥 Download E-GEOD-15629 for external validation

### Short-Term (1-4 Weeks)
1. 📧 Contact Ruigrok (GSE122897) - **HIGHEST PRIORITY**
2. 📧 Contact Poppenberg (GSE159610) - **HIGH PRIORITY**
3. 📋 Apply for dbGaP access - **HIGH PRIORITY**

### Long-Term (1-6 Months)
1. Prospective study design (if needed)
2. Multi-center collaboration

---

## 📊 Expected Outcomes

| Scenario | Blood Model | Tissue Model | Multi-Omics | Status |
|----------|------------|--------------|-------------|--------|
| **Current Data** | AUC 0.968 | AUC > 0.70 (expected) | Not directly possible | ✅ Available |
| **With Genotype** | AUC 0.95-0.98 | AUC 0.75-0.85 | Blood/Tissue + PRS | ⚠️ Requires genotype |
| **With Paired Samples** | - | - | AUC 0.85-0.95 | ⚠️ Requires contact |

---

## 📚 Key Files

- **Dataset Guide**: `DATASET_ACQUISITION_GUIDE_MULTI_OMICS.md` (comprehensive)
- **Tissue Model**: `build_enhanced_tissue_model_ruptured_vs_unruptured.py`
- **Integration Scripts**: `integrate_blood_tissue_comprehensive.py`

---

**Summary**: 
- ✅ **Good data coverage** - Can proceed with tissue model building
- 🔥 **Highest priority**: Contact authors for paired samples/genotype data
- 📥 **Immediate action**: Download E-GEOD-15629 for validation


