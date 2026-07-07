# Programmatic Access Summary for Multi-Omics Datasets

**Generated:** 2025-01-XX  
**Status Check**: ✅ Complete

---

## Executive Summary

I've checked which datasets can be obtained programmatically. Here's what I found:

| Dataset/Contact | Programmatic? | Status | Action |
|----------------|---------------|--------|--------|
| **E-GEOD-15629** (ArrayExpress) | ✅ **YES** | Downloadable | Script created |
| **GSE122897** (GEO) | ✅ **YES** | Downloadable | Already have scripts |
| **GSE159610** (GEO) | ✅ **YES** | Downloadable | Already have scripts |
| **dbGaP** | ❌ **NO** | Requires approval | Manual application |
| **Author Contacts** | ❌ **NO** | Requires email | Manual contact |

---

## ✅ Programmatically Accessible

### 1. **E-GEOD-15629** (ArrayExpress) ✅

**Status**: ✅ **Can download programmatically**

**Dataset Info**:
- **Platform**: ArrayExpress
- **Samples**: 19 total (8 ruptured + 6 unruptured + 5 controls)
- **Technology**: Microarray
- **Purpose**: External validation

**Download Methods**:
1. **ArrayExpress API** (v3): 
   - URL: `https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-GEOD-15629`
   - Status: ✅ API accessible
   - Note: Some files may require manual download from website

2. **Direct HTTP/FTP**: 
   - Can attempt direct file downloads
   - Pattern: `https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-15629/`

3. **BioStudies** (new platform):
   - ArrayExpress datasets migrated to BioStudies
   - URL: `https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEOD-15629`

**Script Created**: ✅ `download_e_geod_15629.py`
- Attempts automatic download via API
- Falls back to manual download instructions if API fails
- Creates README with download links

**Action**: ✅ **READY TO DOWNLOAD** - Run `python download_e_geod_15629.py`

---

### 2. **GSE122897** (GEO) ✅

**Status**: ✅ **Can download programmatically**

**Dataset Info**:
- **Platform**: NCBI GEO
- **Samples**: 43 (22 ruptured, 21 unruptured)
- **Technology**: RNA-seq
- **Already Downloaded**: ✅ Yes

**Download Methods**:
1. **GEO API/URL Pattern**:
   - Base URL: `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE122897&format=file`
   - Status: ✅ Accessible

2. **FTP Download**:
   - Pattern: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122897/`

**Existing Scripts**: ✅ Already have download scripts
- Check: `download_datasets.py`
- Status: ✅ Data already downloaded in `data/sc_rnaseq/GSE122897/`

**Action**: ✅ **ALREADY DOWNLOADED** - No action needed

---

### 3. **GSE159610** (GEO) ✅

**Status**: ✅ **Can download programmatically**

**Dataset Info**:
- **Platform**: NCBI GEO
- **Samples**: 47 (23 aneurysm, 24 control)
- **Technology**: Blood RNA-seq
- **Already Downloaded**: ✅ Yes

**Download Methods**:
1. **GEO API/URL Pattern**:
   - Base URL: `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE159610&format=file`
   - Status: ✅ Accessible

**Existing Scripts**: ✅ Already have download scripts
- Status: ✅ Data already downloaded in `data/blood_rnaseq/GSE159610/`

**Action**: ✅ **ALREADY DOWNLOADED** - No action needed

---

## ❌ NOT Programmatically Accessible

### 1. **dbGaP Access** ❌

**Status**: ❌ **NOT programmatically accessible without approval**

**Requirements**:
- ✅ Website accessible: https://www.ncbi.nlm.nih.gov/gap/
- ❌ Download requires:
  - Principal Investigator (PI) account
  - IRB approval
  - Data Use Agreement (DUA)
  - Application review (2-4 weeks)

**What CAN be done programmatically**:
- ✅ Search for studies: Can programmatically search dbGaP
- ✅ View study metadata: Can access study information
- ❌ Download data: Requires approval process

**What CANNOT be done programmatically**:
- ❌ Automatically download genotype data
- ❌ Bypass approval process
- ❌ Access individual-level data

**Action Required**:
1. 📋 **Manual Application**: 
   - Visit: https://www.ncbi.nlm.nih.gov/gap/
   - Search: "intracranial aneurysm" with genotype data
   - Apply for access
   - Wait for approval (2-4 weeks)

2. 🔍 **Search Script** (Optional):
   - Can create script to search for relevant studies
   - But download still requires manual approval

**Priority**: 🔥 **HIGH** - Worth applying, but requires manual process

---

### 2. **Author Contacts** ❌

**Status**: ❌ **NOT programmatically accessible**

**Contacts**:

#### **Dr. Ynte Ruigrok** (GSE122897)
- **Email**: ij.m.ruigrok@umcutrecht.nl
- **Dataset**: GSE122897 (Tissue)
- **Request**: Paired blood samples or genotype data
- **Institution**: UMC Utrecht, Netherlands

#### **Dr. Poppenberg** (GSE159610)
- **Email**: Check GEO metadata
- **Dataset**: GSE159610 (Blood)
- **Request**: Paired tissue samples or genotype data

**What CANNOT be done programmatically**:
- ❌ No API for requesting paired samples
- ❌ No automated email system
- ❌ No programmatic data request system

**Action Required**:
1. 📧 **Manual Email Contact**:
   - Compose email to authors
   - Request paired samples/genotype data
   - Wait for response (1-2 weeks typical)

2. 📝 **Email Template** (Can create template, but requires manual sending):
   - I can create email templates
   - But sending requires manual action

**Priority**: 🔥 **HIGHEST** - Most likely to have paired samples

---

## 📋 Scripts Created

### ✅ Available Scripts

1. **`check_programmatic_access.py`**
   - Checks programmatic access for all datasets
   - Tests APIs and download methods
   - Generates summary report

2. **`download_e_geod_15629.py`**
   - Downloads E-GEOD-15629 from ArrayExpress
   - Attempts automatic download
   - Creates README with manual download instructions if needed

3. **Existing GEO Download Scripts**:
   - `download_datasets.py` - General GEO downloader
   - Data already downloaded for GSE122897 and GSE159610

---

## 🎯 Recommended Actions

### Immediate (Programmatic) ✅

1. ✅ **Download E-GEOD-15629**:
   ```bash
   python download_e_geod_15629.py
   ```
   - Attempts automatic download
   - Falls back to manual instructions if needed

2. ✅ **Verify GEO Datasets**:
   - GSE122897: ✅ Already downloaded
   - GSE159610: ✅ Already downloaded

### Short-Term (Manual) 📧

1. 📧 **Contact Authors** (Highest Priority):
   - Dr. Ynte Ruigrok: ij.m.ruigrok@umcutrecht.nl
   - Dr. Poppenberg: Check GEO metadata
   - **Request**: Paired samples or genotype data
   - **Timeline**: 1-2 weeks for response

2. 📋 **Apply for dbGaP Access**:
   - Visit: https://www.ncbi.nlm.nih.gov/gap/
   - Search: "intracranial aneurysm"
   - Apply for genotype data access
   - **Timeline**: 2-4 weeks for approval

---

## 📊 Summary Table

| Item | Programmatic? | Status | Action |
|------|---------------|--------|--------|
| **E-GEOD-15629** | ✅ Yes | Ready | Run download script |
| **GSE122897** | ✅ Yes | Downloaded | ✅ Complete |
| **GSE159610** | ✅ Yes | Downloaded | ✅ Complete |
| **dbGaP Search** | ✅ Yes | Can search | Manual application for download |
| **dbGaP Download** | ❌ No | Requires approval | Manual application |
| **Author Contact** | ❌ No | Requires email | Manual email |

---

## ✅ Conclusion

**Programmatically Accessible**:
- ✅ E-GEOD-15629 (ArrayExpress) - Script ready
- ✅ GSE122897 (GEO) - Already downloaded
- ✅ GSE159610 (GEO) - Already downloaded

**Requires Manual Action**:
- ❌ dbGaP access - Requires application (2-4 weeks)
- ❌ Author contacts - Requires email (1-2 weeks response)

**Next Steps**:
1. ✅ Run `download_e_geod_15629.py` to download validation dataset
2. 📧 Contact authors for paired samples (highest priority)
3. 📋 Apply for dbGaP access for genotype data

---

**Files Created**:
- `check_programmatic_access.py` - Access checker
- `download_e_geod_15629.py` - ArrayExpress downloader
- `PROGRAMMATIC_ACCESS_SUMMARY.md` - This document


