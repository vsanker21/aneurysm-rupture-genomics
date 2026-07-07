# Multi-Omics Data Sources for Integration

## Current Situation

### Available Datasets

**1. Blood Dataset (GSE159610)**
- **Type**: Blood RNA-seq
- **Task**: Aneurysm vs Control (binary classification)
- **Samples**: 47 samples (23 aneurysm, 24 control)
- **Data**: ✓ Available
- **Genotype**: ✗ Not available
- **Tissue**: ✗ Not available
- **Model Performance**: AUC 0.968 (blood-only model)

**2. Tissue Dataset (GSE122897)**
- **Type**: Tissue RNA-seq
- **Task**: Ruptured vs Unruptured (binary classification)
- **Samples**: ~50 samples (ruptured vs unruptured)
- **Data**: ✓ Available
- **Genotype**: ✗ Not available
- **Blood**: ✗ Not available
- **Model Performance**: Being built (tissue-only model)

**3. GWAS Data (Bakker 2020)**
- **Type**: GWAS summary statistics
- **Samples**: 10,754 cases + 306,882 controls
- **Data**: ✓ Available (summary statistics only)
- **Genotype**: ✗ Individual genotypes not available
- **PRS Weights**: ✓ Available (`results/prs_weights.tsv`)

### Problem

**Cannot directly integrate** because:
- **Different patients**: Blood and tissue datasets have different samples
- **Different phenotypes**: Blood = aneurysm vs control; Tissue = ruptured vs unruptured
- **Missing genotype data**: No individual-level genotype data available

---

## Potential Data Sources for Multi-Omics Integration

### 1. dbGaP (Database of Genotypes and Phenotypes)

**URL**: https://www.ncbi.nlm.nih.gov/gap/

**What It Provides**:
- Genotype data (SNP arrays, WGS)
- Paired samples (blood + tissue, if available)
- Clinical metadata
- Phenotype data

**Access Requirements**:
- Principal Investigator (PI) account
- Institutional affiliation
- IRB approval
- Data Use Agreement (DUA)

**How to Apply**:
1. Search dbGaP: https://www.ncbi.nlm.nih.gov/gap/
2. Search for "intracranial aneurysm" or "cerebral aneurysm"
3. Find relevant studies with genotype data
4. Click "Request Access"
5. Complete application process
6. Wait for approval (typically 2-4 weeks)

**Potential Benefits**:
- Genotype data for PRS calculation
- May have paired samples (blood + tissue)
- Complete clinical metadata
- May have rupture outcomes

**Status**: Requires application and approval

---

### 2. Contact Study Authors

**Poppenberg et al. (GSE159610 - Blood)**
- **Contact**: Check GEO metadata for contact email
- **Potential**: May have:
  - Tissue samples from same patients
  - Genotype data
  - Additional clinical metadata
- **Action**: Email to request paired samples or genotype data
- **Benefits**: Same patients across omics layers

**Ruigrok et al. (GSE122897 - Tissue)**
- **Contact**: Ynte Ruigrok (ij.m.ruigrok@umcutrecht.nl)
- **Institution**: UMC Utrecht, Netherlands
- **Potential**: May have:
  - Blood samples from same patients
  - Genotype data
  - Clinical follow-up data
- **Action**: Email to request paired blood samples or genotype data
- **Benefits**: Same patients across omics layers

**Bakker et al. (GWAS)**
- **Potential**: May have:
  - Individual-level genotype data (with approval)
  - Phenotype data
  - Clinical metadata
- **Action**: Contact for collaboration or data sharing
- **Benefits**: Genotype data for PRS calculation

**Status**: Requires collaboration

---

### 3. E-GEOD-15629 (ArrayExpress)

**URL**: https://www.omicsdi.org/dataset/biostudies-arrayexpress/E-GEOD-15629

**What It Provides**:
- Tissue gene expression (microarray)
- **Samples**: 8 ruptured + 6 unruptured + 5 controls
- **Task**: Ruptured vs Unruptured vs Control

**Limitations**:
- Microarray (not RNA-seq) - different technology
- No blood samples
- No genotype data
- Small sample size (19 samples)

**Use Cases**:
- External validation of tissue model
- Cross-platform validation
- Additional tissue samples for analysis

**Status**: Available (open access)

---

### 4. UK Biobank

**URL**: https://www.ukbiobank.ac.uk/

**What It Provides**:
- Genotype data (if approved)
- Clinical data
- May have blood samples
- Population-based cohort

**Limitations**:
- Requires approval
- May not have aneurysm-specific samples
- Tissue data unlikely
- Large cohort but may not be aneurysm-focused

**Use Cases**:
- Genotype data for PRS calculation
- Clinical risk factors
- Population-level validation

**Status**: Requires application and approval

---

### 5. Prospective Clinical Study

**What It Provides**:
- Complete control over data collection
- Paired samples (blood + tissue from same patients)
- Genotype data (SNP array or WGS)
- Clinical metadata
- Follow-up for rupture outcomes

**Requirements**:
- IRB approval
- Patient recruitment (~200-300 patients)
- Sample collection (blood + tissue)
- Genotyping (SNP array or WGS)
- Follow-up (2-5 years)

**Timeline**: 2-5 years

**Benefits**:
- Perfect data for multi-omics integration
- All omics layers from same patients
- Complete clinical metadata
- Longitudinal follow-up

**Status**: Future option (requires funding and IRB)

---

## Recommended Approach

### Immediate (Current Capability)

**1. Separate Models**:
- ✓ **Blood Model**: Aneurysm vs Control (AUC 0.968)
  - Use for: Screening/detection
- ⚠ **Tissue Model**: Ruptured vs Unruptured (being built)
  - Use for: Risk stratification after detection

**2. Feature-Level Integration**:
- Identify common signature genes (blood + tissue)
- Cross-validate features across datasets
- Use for biomarker panel development (qPCR)

**3. Sequential Clinical Approach**:
- Step 1: Blood test (aneurysm present?) → Blood model
- Step 2: If positive, tissue analysis (rupture risk?) → Tissue model
- Step 3: If high risk, intervention

---

### Short-Term (1-6 months)

**1. Apply for dbGaP Access**:
- Search for aneurysm studies with genotype data
- Apply for genotype + phenotype data
- May provide paired samples (blood + tissue)

**2. Contact Authors**:
- Email Poppenberg (blood): Request paired tissue samples
- Email Ruigrok (tissue): Request paired blood samples
- Email Bakker (GWAS): Request genotype data

**3. Download E-GEOD-15629**:
- Use for external validation
- Cross-platform validation (microarray vs RNA-seq)

**4. Integrate Available Data**:
- If genotype data obtained: Add PRS to blood model
- If paired samples obtained: True multi-omics integration

---

### Long-Term (1-2 years)

**1. Prospective Study**:
- Collect blood + tissue from same patients
- Genotype all patients (SNP array)
- Follow for rupture outcomes
- Complete multi-omics integration

**2. Multi-Center Collaboration**:
- Pool resources across institutions
- Share samples and data
- Larger sample size for validation

---

## Data Requirements for True Multi-Omics Integration

### Required:

1. **Same Patients**: 
   - ✓ Blood sample
   - ✓ Tissue sample
   - ✓ Genotype data
   - ✓ Clinical metadata

2. **Same Phenotype**:
   - Rupture prediction (ruptured vs unruptured)
   - OR: Rupture risk (high vs low risk)

3. **Sufficient Sample Size**:
   - Training: ~100-200 patients
   - Validation: ~50-100 patients
   - Total: ~200-300 patients

4. **Data Quality**:
   - RNA-seq (not just microarray)
   - Genotype call rate ≥95%
   - Clinical metadata complete

### Current Capability:

- ✓ Blood transcriptomics: Available
- ✓ Tissue transcriptomics: Available
- ✓ GWAS summary stats: Available
- ✗ Individual genotypes: Not available
- ✗ Paired samples: Not available

---

## Action Items

### Priority 1 (Immediate):

1. ✅ **Build tissue model** (ruptured vs unruptured)
   - Status: In progress
   - Expected: AUC > 0.70

2. ✅ **qPCR panel design** (top 10-20 genes)
   - Status: Complete
   - File: `results/qpcr_panel_top20.tsv`

3. ⚠ **External validation** (GSE66238)
   - Status: Pending (needs sample labels)

### Priority 2 (Short-term):

1. **Apply for dbGaP access**
   - Timeline: 2-4 weeks for approval
   - Goal: Genotype data for PRS

2. **Contact authors**
   - Timeline: Variable
   - Goal: Paired samples or genotype data

3. **Download E-GEOD-15629**
   - Timeline: Immediate
   - Goal: External validation

### Priority 3 (Long-term):

1. **Prospective study design**
   - Timeline: 1-2 years
   - Goal: Complete multi-omics integration

2. **Multi-center collaboration**
   - Timeline: Variable
   - Goal: Larger sample size

---

## Expected Outcomes

### With Current Data:

- **Blood Model**: AUC 0.968 (aneurysm vs control) ✓
- **Tissue Model**: AUC > 0.70 (ruptured vs unruptured) (expected)
- **Multi-Omics**: Not possible (different patients, different phenotypes)

### With Genotype Data (PRS Integration):

- **Blood + PRS**: AUC 0.95-0.98 (expected improvement)
- **Tissue + PRS**: AUC 0.75-0.85 (expected improvement)

### With Paired Samples (True Multi-Omics):

- **Integrated Model**: AUC 0.85-0.95 (expected)
- **Features**: PRS + Blood RNA + Tissue RNA + Clinical
- **Benefit**: Superior to any single-omics model

---

## Summary

**Current Status**:
- ✓ Blood model: Available (AUC 0.968)
- ⚠ Tissue model: Being built
- ✗ Multi-omics: Not possible with current data

**Next Steps**:
1. Complete tissue model (ruptured vs unruptured)
2. Apply for dbGaP access (genotype data)
3. Contact authors (paired samples)
4. Consider prospective study (long-term)

**Key Insight**:
True multi-omics integration requires **same patients across all omics layers**. Current datasets have different patients, so direct integration is not scientifically valid. However, separate models (blood + tissue) can be used sequentially in clinical practice.

---

**Last Updated**: [Date]
**Files Created**:
- `results/multi_omics_data_sources_checklist.tsv`
- `results/multi_omics_data_sources_summary.md`


