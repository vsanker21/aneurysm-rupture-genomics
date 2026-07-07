# Supplementary PDFs Summary

## ✅ Downloaded Supplementary PDFs

Two supplementary PDFs have been downloaded and organized in the project:

### 1. NIHMS4348-supplement-s1.pdf
- **Size**: 0.96 MB
- **Location**: `data/supplementary_pdfs/NIHMS4348-supplement-s1.pdf`
- **Likely Source**: Nature Genetics paper (possibly Bilguvar et al. 2008)
- **Potential Study**: Bilguvar et al. 2008 - Susceptibility loci for intracranial aneurysm
- **PubMed ID**: 18758463
- **Expected Content**:
  - GWAS summary statistics
  - Risk loci data
  - Supplementary tables
- **Status**: ✅ Downloaded - Needs extraction

### 2. str_stroke-2014-006096_supp1.pdf
- **Size**: 1.01 MB
- **Location**: `data/supplementary_pdfs/str_stroke-2014-006096_supp1.pdf`
- **Likely Source**: Stroke journal, 2014
- **Manuscript ID**: 006096
- **Potential Study**: Aneurysm-related study from Stroke journal 2014
- **Expected Content**:
  - Clinical data
  - Statistical analyses
  - Supplementary results
- **Status**: ✅ Downloaded - Needs extraction

## 📋 Next Steps

### 1. Identify Content
- Review PDFs to confirm their source papers
- Identify what data they contain (GWAS stats, tables, etc.)
- Verify if they contain data already in downloaded datasets

### 2. Extract Data
- **For GWAS Summary Statistics** (if present):
  - Extract tables containing SNP data
  - Convert to tabular format (TSV/CSV)
  - Standardize column names
  - Save in `data/gwas/[study_name]/`
  
- **For Supplementary Tables**:
  - Extract relevant tables
  - Save in appropriate directories
  - Document table contents

### 3. Tools Needed
- PDF text extraction: PyPDF2, pdfplumber
- Table extraction: tabula-py, camelot-py
- See `data/supplementary_pdfs/EXTRACTION_GUIDE.md` for details

### 4. Integration
- If GWAS data is present, integrate with existing Bakker 2020 dataset
- Compare with existing datasets
- Update dataset inventory

## 📊 Current Status

- ✅ PDFs downloaded and organized
- ✅ Documentation created
- ✅ Inventory created
- ⚠️ Data extraction pending
- ⚠️ Content verification pending

## 🔗 Related Datasets

- See `data/gwas/bilguvar_2008/` for Bilguvar 2008 GWAS
- See `data/gwas/foroud_2014/` for Foroud 2014 GWAS
- See `data/gwas/bakker_2020/` for Bakker 2020 GWAS (already downloaded)

## 📝 Notes

- PDFs may contain controlled-access data - verify citation requirements
- Some data may already be available in existing datasets
- Manual review may be necessary for complex tables
- Extraction tools may need to be installed

---

**Last Updated**: 2025-11-13
**Status**: PDFs downloaded, extraction pending






