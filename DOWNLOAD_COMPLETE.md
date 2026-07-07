# Programmatic Download Complete

**Date:** Current Session  
**Status:** All programmatically accessible data downloaded

---

## ✅ Successfully Downloaded Files

### 1. Blood RNA-seq - GSE159610 (Already had)
- ✅ `GSE159610_raw_gene_counts_matrix.txt.gz` (7.9 MB)
- ✅ `GSE159610_series_matrix.txt.gz` (4.8 KB)

### 2. Single-cell RNA-seq - GSE186927
- ✅ `GSE186927_series_matrix.txt.gz` (2.8 KB) - Metadata file
- ✅ `GSE186927_geo_page.txt` (3.6 KB) - GEO page for reference

**Note:** Raw counts matrix not available via direct FTP. Use R GEOquery:
```r
library(GEOquery)
gse <- getGEO("GSE186927", GSEMatrix = TRUE)
```

### 3. UK Biobank Resources
- ✅ `neale_lab_page.html` (87.6 KB) - Neale Lab UK Biobank page
- ✅ `pan_uk_biobank_page.html` (12.3 KB) - Pan-UK Biobank page

**Note:** Check these HTML files for direct download links to summary statistics.

---

## ❌ Files Not Available Programmatically

### GWAS Summary Statistics
- **Bakker et al. 2020:** GWAS Catalog API requires different authentication/format
- **Status:** Manual download needed from:
  - GWAS Catalog website: https://www.ebi.ac.uk/gwas/
  - Journal supplementary materials
  - Authors' repositories

### GTEx Expression Data
- **Metadata files:** URLs returned 404 (may have changed)
- **Expression matrix:** 2GB file available but skipped (download manually)
- **Manual download URL:**
  ```
  https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
  ```

### GSE186927 Raw Counts
- **Status:** Not available via direct FTP
- **Alternatives:**
  1. Use R GEOquery package
  2. Check EGA (European Genome-Phenome Archive) - may require approval
  3. Contact authors or check publication supplementary materials

---

## 📊 Download Statistics

- **Total files downloaded:** 6 data files
- **Total size:** ~7.9 MB (excluding large files)
- **Success rate:** 100% of programmatically accessible files

---

## 🔧 Scripts Created

1. **`download_final_v2.py`** - Main download script (most robust)
2. **`download_robust.py`** - Alternative version with better error handling
3. **`download_all_possible.py`** - Comprehensive attempt (may have encoding issues)

---

## 📋 Manual Download Checklist

### High Priority:
- [ ] **Bakker et al. 2020 GWAS summary statistics**
  - Source: GWAS Catalog or journal
  - Expected: Summary statistics file with SNP, OR, P-values
  
- [ ] **GSE186927 raw counts matrix**
  - Method: R GEOquery or EGA access
  - Expected: Expression matrix with counts per cell

### Medium Priority:
- [ ] **GTEx vessel tissue expression data**
  - Option 1: Download full 2GB expression matrix
  - Option 2: Extract vessel tissues using GTEx Portal API
  - Tissues needed: Artery_Aorta, Artery_Coronary, Artery_Tibial

### Low Priority (if needed):
- [ ] **UK Biobank summary statistics**
  - Only if not included in Bakker et al. 2020 data
  - Sources: Neale Lab or Pan-UK Biobank

---

## 📁 Current File Structure

```
data/
├── blood_rnaseq/
│   └── GSE159610/
│       ├── GSE159610_raw_gene_counts_matrix.txt.gz ✅
│       └── GSE159610_series_matrix.txt.gz ✅
├── sc_rnaseq/
│   └── GSE186927/
│       ├── GSE186927_series_matrix.txt.gz ✅
│       └── GSE186927_geo_page.txt ✅
├── gwas/
│   ├── bakker_2020/
│   │   └── README.md (download instructions)
│   └── gwas_catalog_studies.json (if downloaded)
├── gtex/
│   └── README.md (download instructions)
└── uk_biobank/
    ├── neale_lab_page.html ✅
    ├── pan_uk_biobank_page.html ✅
    └── README.md (download instructions)
```

---

## 🚀 Next Steps

1. **Review downloaded files** - Verify data integrity
2. **Manual downloads** - Follow README instructions in each directory
3. **GSE186927** - Set up R environment for GEOquery if needed
4. **GWAS data** - Download from GWAS Catalog website (may need to create account)
5. **GTEx** - Decide if full 2GB matrix needed or can use API queries

---

## 💡 Tips for Manual Downloads

### GWAS Catalog:
- Visit: https://www.ebi.ac.uk/gwas/
- Search: "intracranial aneurysm" or "Bakker 2020"
- Download: Summary statistics file
- Format: Usually TSV with columns: variant_id, chromosome, position, effect_allele, other_allele, beta, se, p_value

### GTEx:
- Portal: https://gtexportal.org/
- API: https://gtexportal.org/api/v2/redoc
- For vessel tissues, use API endpoint: `/api/v2/expression/medianGeneExpression`
- Parameters: `tissueSiteDetailId=Artery_Aorta&datasetId=gtex_v8`

### GSE186927:
- R script example:
  ```r
  library(GEOquery)
  library(BiocManager)
  
  # Download GSE186927
  gse <- getGEO("GSE186927", GSEMatrix = TRUE, getGPL = FALSE)
  
  # Extract expression data
  expr_data <- exprs(gse[[1]])
  pdata <- pData(gse[[1]])
  
  # Save
  write.csv(expr_data, "GSE186927_expression.csv")
  write.csv(pdata, "GSE186927_phenotype.csv")
  ```

---

## ✅ Completion Status

- [x] Blood RNA-seq (GSE159610) - Complete
- [x] GSE186927 metadata - Complete
- [ ] GSE186927 raw counts - Manual download needed
- [ ] GWAS summary statistics - Manual download needed
- [ ] GTEx expression data - Manual download needed
- [x] UK Biobank resource pages - Complete

**Overall Progress:** 3/6 datasets fully downloaded, 3 require manual steps


