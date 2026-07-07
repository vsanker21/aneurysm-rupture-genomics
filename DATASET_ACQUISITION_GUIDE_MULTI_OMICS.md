# Dataset Acquisition Guide for Multi-Omics Integration
## Comprehensive Resource for Obtaining Data for Cerebral Aneurysm Genomics Project

**Last Updated:** 2025-01-XX

---

## Executive Summary

This guide outlines **all available datasets and data sources** that can be obtained to proceed with multi-omics integration for cerebral aneurysm rupture prediction. It includes immediate options, short-term opportunities, and long-term strategies.

---

## Part 1: Currently Available Datasets ✅

### 1.1 Blood Transcriptomics

#### **GSE159610** (Poppenberg et al. 2020)
- **Type**: Whole blood RNA-seq
- **Samples**: 47 (23 aneurysm, 24 control)
- **Status**: ✅ **FULLY DOWNLOADED**
- **Location**: `data/blood_rnaseq/GSE159610/`
- **Platform**: Illumina HiSeq 2500
- **Model Performance**: AUC 0.968 (blood-only model)
- **What it provides**: Blood-based aneurysm detection signatures
- **Limitations**: 
  - Different patients than tissue datasets
  - Aneurysm vs Control (not ruptured vs unruptured)
  - No genotype data
  - No paired tissue samples

**Action**: ✅ Already available - ready for analysis

---

### 1.2 Tissue Transcriptomics

#### **GSE122897** (Ruigrok et al. 2018) ⭐ **PRIMARY**
- **Type**: Tissue RNA-seq (bulk)
- **Samples**: 60 total (17 Ruptured, 18 Unruptured, 25 Control/Unknown)
- **Status**: ✅ **FULLY DOWNLOADED**
- **Location**: `data/sc_rnaseq/GSE122897/`
- **Platform**: Illumina HiSeq 2500
- **Task**: Ruptured vs Unruptured classification
- **What it provides**: Direct rupture prediction comparison
- **Key Finding**: Focus on lysosomes and immunoglobulins in rupture
- **Contact**: Ynte Ruigrok (ij.m.ruigrok@umcutrecht.nl) - UMC Utrecht
- **Potential for Integration**:
  - May have blood samples from same patients (need to contact)
  - May have genotype data (need to contact)

**Action**: 
- ✅ Data available - building tissue model
- 📧 **RECOMMENDED**: Contact authors for paired blood samples or genotype data

#### **GSE66238** ⭐ **SECONDARY**
- **Type**: Tissue RNA-seq (bulk)
- **Samples**: Ruptured vs Unruptured
- **Status**: ✅ **FULLY DOWNLOADED**
- **Location**: `data/sc_rnaseq/GSE66238/`
- **Additional Data**: Includes microRNA expression
- **What it provides**: 
  - Rupture comparison
  - microRNA signatures
  - Isoform-level differential expression
- **Potential for Integration**: Can be meta-analyzed with GSE122897

**Action**: ✅ Data available - can integrate with GSE122897

#### **GSE158558**
- **Type**: Tissue RNA-seq (bulk)
- **Status**: ✅ **FULLY DOWNLOADED**
- **Location**: `data/sc_rnaseq/GSE158558/`
- **Additional Data**: mRNA and lncRNA expression
- **What it provides**: Long non-coding RNA signatures in aneurysm tissue

**Action**: ✅ Data available - can be integrated

#### **E-GEOD-15629** (ArrayExpress)
- **Type**: Tissue microarray
- **Samples**: 19 total (8 ruptured + 6 unruptured + 5 controls)
- **Status**: ⚠️ **NOT YET DOWNLOADED**
- **Platform**: Microarray (different technology)
- **URL**: https://www.omicsdi.org/dataset/biostudies-arrayexpress/E-GEOD-15629
- **What it provides**: External validation dataset
- **Limitations**: 
  - Microarray (not RNA-seq) - cross-platform validation
  - Small sample size
  - No blood samples
  - No genotype data

**Action**: 📥 **RECOMMENDED**: Download for external validation

---

### 1.3 GWAS Data

#### **Bakker et al. 2020** (Meta-GWAS)
- **Type**: GWAS summary statistics
- **Samples**: 10,754 cases + 306,882 controls
- **Status**: ✅ **FULLY DOWNLOADED**
- **Location**: `data/gwas/bakker_2020/`
- **What it provides**: 
  - 17 risk loci (11 novel)
  - PRS weights (`results/prs_weights.tsv`)
  - SNP-based heritability: 21.6%
- **Limitations**: 
  - Summary statistics only (no individual genotypes)
  - Cannot calculate patient-level PRS without genotypes

**Action**: ✅ Available - can use for feature-level integration

**Contact for Additional Data**: 
- **Potential**: Individual-level genotype data (with approval)
- **Action**: 📧 Contact authors for collaboration/data sharing

---

### 1.4 Reference Data

#### **GTEx v10** (Normal Vessels)
- **Type**: Normal tissue transcriptomes
- **Tissues**: Artery (Aorta, Coronary, Tibial)
- **Status**: ✅ **FULLY DOWNLOADED**
- **Location**: `data/gtex/`
- **What it provides**: Reference normal vessel expression

**Action**: ✅ Available - used for normalization/reference

---

## Part 2: Datasets That Can Be Obtained 📥

### 2.1 dbGaP (Database of Genotypes and Phenotypes)

**URL**: https://www.ncbi.nlm.nih.gov/gap/

**What It Provides**:
- Individual-level genotype data (SNP arrays, WGS)
- Paired samples (blood + tissue, if available)
- Complete clinical metadata
- Phenotype data
- Rupture outcomes

**Access Requirements**:
- Principal Investigator (PI) account
- Institutional affiliation
- IRB approval
- Data Use Agreement (DUA)
- Application review (typically 2-4 weeks)

**How to Apply**:
1. Go to https://www.ncbi.nlm.nih.gov/gap/
2. Search for "intracranial aneurysm" or "cerebral aneurysm"
3. Filter for studies with genotype data
4. Click "Request Access" on relevant studies
5. Complete application form:
   - Research objectives
   - IRB approval documentation
   - Data use restrictions
   - Institutional signatory
6. Wait for approval (typically 2-4 weeks)

**Potential Studies to Search**:
- Search terms: "intracranial aneurysm", "cerebral aneurysm", "subarachnoid hemorrhage"
- Filter: Studies with genotype data
- Priority: Studies with paired samples (blood + tissue)

**Benefits**:
- ✅ Genotype data for PRS calculation
- ✅ May have paired samples (blood + tissue from same patients)
- ✅ Complete clinical metadata
- ✅ Rupture outcomes
- ✅ **True multi-omics integration possible**

**Status**: ⚠️ **REQUIRES APPLICATION AND APPROVAL**

**Priority**: 🔥 **HIGH** - Enables patient-level PRS and multi-omics integration

---

### 2.2 Contact Study Authors (Collaboration)

#### **Poppenberg et al. (GSE159610 - Blood)**

**Contact Information**:
- Check GEO metadata: `data/blood_rnaseq/GSE159610/GSE159610_series_matrix.txt.gz`
- Publication: Search PubMed for "Poppenberg" + "GSE159610"

**Potential Data Available**:
- Tissue samples from same patients
- Genotype data (if collected)
- Additional clinical metadata (follow-up, rupture outcomes)
- Paired samples for multi-omics integration

**Request Template**:
```
Subject: Request for Collaboration - Cerebral Aneurysm Genomics Study

Dear Dr. Poppenberg,

I am writing to request collaboration regarding your dataset GSE159610 
(whole blood transcriptomes of aneurysm patients). 

Our research group is developing a multi-omics model for cerebral 
aneurysm rupture prediction, integrating:
1. GWAS polygenic risk scores
2. Blood transcriptomics (using your GSE159610 data)
3. Tissue transcriptomics (from GSE122897)

We would be grateful if you could share:
- Tissue samples from the same patients (if available)
- Genotype data (if available)
- Additional clinical metadata (follow-up outcomes)

We would be happy to acknowledge your contribution and collaborate 
on publications.

Best regards,
[Your Name]
[Institution]
```

**Benefits**:
- Same patients across omics layers
- True multi-omics integration
- Better statistical power

**Status**: ⚠️ **REQUIRES CONTACT**

**Priority**: 🔥 **HIGH** - Most direct path to paired samples

---

#### **Ruigrok et al. (GSE122897 - Tissue)**

**Contact Information**:
- **Email**: ij.m.ruigrok@umcutrecht.nl
- **Institution**: UMC Utrecht, Netherlands
- **GEO Profile**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122897

**Potential Data Available**:
- Blood samples from same patients (if collected)
- Genotype data (if available)
- Clinical follow-up data
- Additional metadata (PHASES scores, imaging data)

**Request Template**:
```
Subject: Request for Collaboration - Cerebral Aneurysm Genomics Study

Dear Dr. Ruigrok,

I am writing to request collaboration regarding your dataset GSE122897 
(RNA-seq analysis of intracranial aneurysm walls - ruptured vs unruptured).

Our research group is developing a multi-omics model for cerebral 
aneurysm rupture prediction. We are using your tissue data as the 
foundation for rupture prediction.

We would be grateful if you could share:
- Blood samples from the same patients (if available)
- Genotype data (if available)
- Additional clinical metadata (PHASES scores, follow-up outcomes)

We would be happy to acknowledge your contribution and collaborate 
on publications.

Best regards,
[Your Name]
[Institution]
```

**Benefits**:
- Paired blood + tissue samples from same patients
- Complete clinical metadata
- True multi-omics integration

**Status**: ⚠️ **REQUIRES CONTACT**

**Priority**: 🔥 **HIGH** - Most likely to have paired samples

---

#### **Bakker et al. (GWAS 2020)**

**Contact Information**:
- Check publication: Bakker et al. 2020 (European Heart Journal)
- Corresponding author email from publication

**Potential Data Available**:
- Individual-level genotype data (with approval)
- Phenotype data
- Clinical metadata
- Subset of patients with expression data (if available)

**Benefits**:
- Patient-level PRS calculation
- Larger sample size
- Genetic risk stratification

**Status**: ⚠️ **REQUIRES CONTACT**

**Priority**: 🟡 **MEDIUM** - Useful but GWAS summary stats already available

---

### 2.3 Additional GEO/ArrayExpress Datasets

#### **E-GEOD-15629** (ArrayExpress)
- **Status**: ⚠️ **NOT YET DOWNLOADED**
- **Action**: 📥 Download immediately for external validation
- **Script**: Create download script similar to `download_gse122897.py`

#### **GSE13353, GSE54083, GSE75436** (Microarray)
- **Status**: ✅ **DOWNLOADED** (in `data/microarray/`)
- **What they provide**: 
  - Rupture comparison (GSE13353, GSE54083)
  - Additional validation datasets
- **Limitations**: Microarray (different platform), smaller sample sizes
- **Action**: ✅ Available - can use for meta-analysis

---

### 2.4 UK Biobank

**URL**: https://www.ukbiobank.ac.uk/

**What It Provides**:
- Genotype data (if approved)
- Clinical data
- Blood samples (if approved)
- Population-based cohort

**Access Requirements**:
- Research registration
- IRB approval
- Data Use Agreement
- Application review

**Limitations**:
- May not have aneurysm-specific samples
- Tissue data unlikely
- Large cohort but may not be aneurysm-focused

**Use Cases**:
- Genotype data for PRS calculation (population-level)
- Clinical risk factors
- Population-level validation

**Status**: ⚠️ **REQUIRES APPLICATION**

**Priority**: 🟡 **MEDIUM** - Useful but Bakker 2020 already includes UK Biobank data

**Alternative Sources**:
- **Neale Lab**: http://www.nealelab.is/uk-biobank/ (summary stats)
- **Pan-UK Biobank**: https://pan.ukbb.broadinstitute.org/ (summary stats)

---

### 2.5 European Genome-Phenome Archive (EGA)

**URL**: https://ega-archive.org/

**What It Provides**:
- Individual-level sequencing data
- Genotype data
- Controlled access datasets

**How to Search**:
1. Go to https://ega-archive.org/
2. Search for "intracranial aneurysm" or "cerebral aneurysm"
3. Filter for human studies
4. Apply for access if relevant datasets found

**Access Requirements**: Similar to dbGaP (IRB, DUA, approval)

**Status**: ⚠️ **REQUIRES SEARCH AND APPLICATION**

**Priority**: 🟡 **MEDIUM** - Worth checking but may overlap with GEO

---

## Part 3: Data Integration Strategies

### 3.1 Current Capability (Available Now)

**Separate Models Approach**:
1. ✅ **Blood Model**: Aneurysm vs Control (AUC 0.968)
   - Use for: Screening/detection
2. ✅ **Tissue Model**: Ruptured vs Unruptured (being built)
   - Use for: Risk stratification after detection

**Feature-Level Integration**:
- Identify common signature genes (blood + tissue)
- Cross-validate features across datasets
- Use for biomarker panel development (qPCR)

**Sequential Clinical Approach**:
- Step 1: Blood test (aneurysm present?) → Blood model
- Step 2: If positive, tissue analysis (rupture risk?) → Tissue model
- Step 3: If high risk, intervention

**Status**: ✅ **FEASIBLE NOW** - Can proceed without additional data

---

### 3.2 With Genotype Data (dbGaP/Author Contact)

**Blood + PRS Integration**:
- Add PRS scores to blood model
- Expected improvement: AUC 0.95-0.98

**Tissue + PRS Integration**:
- Add PRS scores to tissue model
- Expected improvement: AUC 0.75-0.85

**Status**: ⚠️ **REQUIRES GENOTYPE DATA** - Apply to dbGaP or contact authors

---

### 3.3 With Paired Samples (True Multi-Omics)

**Integrated Model**:
- Features: PRS + Blood RNA + Tissue RNA + Clinical
- Same patients across all omics layers
- Expected AUC: 0.85-0.95

**Benefits**:
- Superior to any single-omics model
- Complete molecular profile per patient
- Better statistical power

**Status**: ⚠️ **REQUIRES PAIRED SAMPLES** - Contact authors (Ruigrok, Poppenberg)

---

## Part 4: Action Plan by Priority

### Priority 1: Immediate Actions (This Week)

1. ✅ **Build tissue model** (GSE122897 ruptured vs unruptured)
   - Status: Script exists, needs to be run
   - Expected: AUC > 0.70

2. 📥 **Download E-GEOD-15629** (ArrayExpress)
   - Purpose: External validation
   - Timeline: Immediate download

3. ✅ **Integrate GSE66238** with GSE122897
   - Purpose: Meta-analysis across tissue datasets
   - Status: Data available, need integration script

---

### Priority 2: Short-Term Actions (1-4 Weeks)

1. 📧 **Contact Ruigrok** (GSE122897)
   - Request: Paired blood samples or genotype data
   - Timeline: 1-2 weeks for response
   - Priority: 🔥 **HIGHEST** - Most likely to have paired samples

2. 📧 **Contact Poppenberg** (GSE159610)
   - Request: Paired tissue samples or genotype data
   - Timeline: 1-2 weeks for response
   - Priority: 🔥 **HIGH** - Direct path to paired samples

3. 📋 **Apply for dbGaP Access**
   - Search for aneurysm studies with genotype data
   - Complete application
   - Timeline: 2-4 weeks for approval
   - Priority: 🔥 **HIGH** - Enables PRS integration

4. 🔍 **Search EGA**
   - Search for aneurysm datasets
   - Apply if relevant datasets found
   - Timeline: 2-4 weeks for approval
   - Priority: 🟡 **MEDIUM**

---

### Priority 3: Long-Term Actions (1-6 Months)

1. 📊 **Prospective Study Design**
   - Timeline: 1-2 years
   - Goal: Complete multi-omics integration
   - Requirements: IRB, funding, patient recruitment

2. 🤝 **Multi-Center Collaboration**
   - Pool resources across institutions
   - Share samples and data
   - Larger sample size for validation

---

## Part 5: Expected Outcomes

### With Current Data:
- ✅ **Blood Model**: AUC 0.968 (aneurysm vs control)
- ✅ **Tissue Model**: AUC > 0.70 (ruptured vs unruptured) (expected)
- ⚠️ **Multi-Omics**: Not directly possible (different patients)

### With Genotype Data:
- ✅ **Blood + PRS**: AUC 0.95-0.98 (expected improvement)
- ✅ **Tissue + PRS**: AUC 0.75-0.85 (expected improvement)

### With Paired Samples:
- ✅ **Integrated Model**: AUC 0.85-0.95 (expected)
- ✅ **Features**: PRS + Blood RNA + Tissue RNA + Clinical
- ✅ **Benefit**: Superior to any single-omics model

---

## Part 6: Summary and Recommendations

### ✅ **What We Have**:
1. Complete blood transcriptomics data (GSE159610)
2. Complete tissue transcriptomics data (GSE122897, GSE66238, GSE158558)
3. GWAS summary statistics (Bakker 2020)
4. Reference data (GTEx)
5. Multiple validation datasets

### 📥 **What Can Be Obtained**:
1. **Paired samples**: Contact Ruigrok (GSE122897) - 🔥 **HIGHEST PRIORITY**
2. **Genotype data**: Apply to dbGaP or contact authors - 🔥 **HIGH PRIORITY**
3. **External validation**: Download E-GEOD-15629 - 📥 **IMMEDIATE**

### 🎯 **Recommended Strategy**:
1. **Immediate**: Build tissue model with available data
2. **Short-term**: Contact authors for paired samples/genotype data
3. **Medium-term**: Apply for dbGaP access
4. **Long-term**: Consider prospective study if needed

### 📊 **Integration Feasibility**:
- **Current**: ✅ Feature-level integration (common genes)
- **With genotype**: ✅ Blood/Tissue + PRS integration
- **With paired samples**: ✅ True multi-omics integration

---

## Part 7: Contact Templates and Resources

### Author Contact Information

**Dr. Ynte Ruigrok** (GSE122897 - Tissue)
- Email: ij.m.ruigrok@umcutrecht.nl
- Institution: UMC Utrecht, Netherlands
- Dataset: GSE122897

**Dr. Poppenberg** (GSE159610 - Blood)
- Contact: Check GEO metadata
- Dataset: GSE159610

**Dr. Bakker** (GWAS 2020)
- Contact: Check publication for corresponding author
- Publication: Bakker et al. 2020, European Heart Journal

### Data Repository Links

- **dbGaP**: https://www.ncbi.nlm.nih.gov/gap/
- **GEO**: https://www.ncbi.nlm.nih.gov/geo/
- **ArrayExpress**: https://www.ebi.ac.uk/arrayexpress/
- **EGA**: https://ega-archive.org/
- **UK Biobank**: https://www.ukbiobank.ac.uk/

---

**Last Updated**: 2025-01-XX  
**Status**: ✅ **COMPREHENSIVE GUIDE COMPLETE**


