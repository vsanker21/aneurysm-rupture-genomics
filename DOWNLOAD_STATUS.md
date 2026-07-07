# Dataset Download Status

**Last Updated:** Current Session

## Summary

| Dataset | Status | Files Downloaded | Notes |
|---------|--------|------------------|-------|
| **Blood RNA-seq (GSE159610)** | ✅ Complete | 2 files | Already downloaded |
| **Intracranial Aneurysm RNA-seq (GSE158558 & GSE122897)** | ✅ Complete | 9 files | Two bulk RNA-seq datasets (ruptured vs unruptured) |
| **GWAS (Bakker et al. 2020)** | ✅ Complete | 1 file | Summary statistics downloaded (763.6 MB) |
| **GTEx** | ✅ Complete | 6 files | All 3 vessel tissues with TPM + raw counts |
| **UK Biobank** | 📋 Documentation | 0 files | README created with download instructions |

---

## Detailed Status

### ✅ 1. Blood RNA-seq - GSE159610 (COMPLETE)
**Location:** `data/blood_rnaseq/GSE159610/`

**Files:**
- ✅ `GSE159610_raw_gene_counts_matrix.txt.gz` (~7.8 MB)
- ✅ `GSE159610_series_matrix.txt.gz` (~4.8 KB)

**Status:** Fully downloaded and ready for analysis.

---

### ✅ 2. Intracranial Aneurysm RNA-seq - GSE158558 & GSE122897 (COMPLETE)

**Status:** ✅ **TWO DATASETS DOWNLOADED** - Both human intracranial aneurysm bulk RNA-seq

---

#### Dataset 1: GSE158558

**Location:** `data/sc_rnaseq/GSE158558/`

**Dataset Information:**
- **Accession:** GSE158558
- **Title:** "Intracranial aneurysm sequencing data"
- **Organism:** Homo sapiens (Human) ✅
- **Type:** Bulk RNA-seq (Expression profiling by high throughput sequencing)

**Files Downloaded (5 files, ~10.6 MB):**
- ✅ `GSE158558_raw_counts_GRCh38.p13_NCBI.tsv.gz` (0.43 MB) - Raw gene counts
- ✅ `GSE158558_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz` (0.67 MB) - TPM normalized expression
- ✅ `GSE158558_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz` (0.68 MB) - FPKM normalized expression
- ✅ `GSE158558_mRNA_lncRNA_Expression_Transcript.xlsx` (4.62 MB) - mRNA and lncRNA expression matrix
- ✅ `Human.GRCh38.p13.annot.tsv.gz` (4.21 MB) - Gene annotation file

**What It Provides:**
- Human intracranial aneurysm tissue expression data
- Both mRNA and lncRNA expression profiles
- Multiple normalization formats

---

#### Dataset 2: GSE122897 ⭐ **HIGHLY VALUABLE**

**Location:** `data/sc_rnaseq/GSE122897/`

**Dataset Information:**
- **Accession:** GSE122897
- **Title:** "RNA Sequencing Analysis of Intracranial Aneurysm Walls Reveals Involvement of Lysosomes and Immunoglobulins in Rupture"
- **Organism:** Homo sapiens (Human) ✅
- **Type:** Bulk RNA-seq (Expression profiling by high throughput sequencing)
- **⭐ Key Feature:** Compares **ruptured vs unruptured** aneurysms - directly relevant to rupture prediction!

**Files Downloaded (4 files, ~12.0 MB):**
- ✅ `GSE122897_readCounts_raw.txt.gz` (2.00 MB) - Raw read counts
- ✅ `GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz` (2.10 MB) - Raw gene counts
- ✅ `GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz` (3.89 MB) - TPM normalized expression
- ✅ `GSE122897_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz` (3.91 MB) - FPKM normalized expression

**What It Provides:**
- **Ruptured vs unruptured aneurysm comparison** - Critical for rupture prediction!
- Human intracranial aneurysm wall expression data
- Identifies rupture-associated genes and pathways
- Focus on lysosomes and immunoglobulins in rupture

---

#### Combined Value of Both Datasets

**Total:** 2 datasets, 9 files, ~22.6 MB

**Advantages of Having Both:**
1. **Validation across studies** - Cross-validate findings between datasets
2. **Increased sample size** - More statistical power
3. **Rupture-specific analysis** - GSE122897 provides ruptured vs unruptured comparison
4. **Meta-analysis** - Combine both datasets for robust findings
5. **Pathway analysis** - GSE122897 identifies specific rupture pathways (lysosomes, immunoglobulins)

**How They Fit the Project:**
- ✅ Both provide human aneurysm tissue expression data
- ✅ GSE122897 directly addresses rupture prediction (ruptured vs unruptured)
- ✅ Can identify endothelial dysfunction markers
- ✅ Complement blood RNA-seq (GSE159610) and GWAS data
- ✅ Provide tissue-specific expression signatures
- ✅ Enable validation of rupture risk signatures

**Note:** While these are bulk RNA-seq (not single-cell), they are highly valuable and provide human aneurysm tissue expression data. GSE122897 is particularly valuable as it directly compares ruptured vs unruptured aneurysms.

**Status:** ✅ Both datasets ready for analysis

---

### ✅ 3. GWAS Summary Statistics - Bakker et al. 2020 (COMPLETE)
**Location:** `data/gwas/bakker_2020/`

**Status:** Successfully downloaded from Knowledge Portal Network.

**Files Downloaded:**
- ✅ `bakker_2020_summary_stats.txt.gz` (763.6 MB) - Summary statistics file
- ✅ `kp4cd_page.html` - Reference page from Knowledge Portal Network

**Source:** Downloaded from https://www.kp4cd.org/node/451

**File Format:** 
- Compressed summary statistics file (likely contains SNP, CHR, BP, A1, A2, OR, SE, P-value, etc.)
- Note: File extension is `.txt.gz` but may be a ZIP archive - verify before extraction

**Status:** Ready for analysis. Extract and verify file format before use.

---

### ✅ 4. GTEx Normal Vessel Reference Transcriptomes (COMPLETE)
**Location:** `data/gtex/`

**Status:** ✅ Successfully downloaded from GTEx Portal. All required vessel tissues obtained.

**Files Downloaded (6 files):**

**Artery - Aorta:**
- ✅ `gene_tpm_v10_artery_aorta.gct.gz` (54.55 MB) - TPM normalized expression
- ✅ `gene_reads_v10_artery_aorta.gct.gz` (22.85 MB) - Raw read counts

**Artery - Coronary:**
- ✅ `gene_tpm_v10_artery_coronary.gct.gz` (32.15 MB) - TPM normalized expression
- ✅ `gene_reads_v10_artery_coronary.gct.gz` (13.70 MB) - Raw read counts

**Artery - Tibial:**
- ✅ `gene_tpm_v10_artery_tibial.gct.gz` (77.18 MB) - TPM normalized expression
- ✅ `gene_reads_v10_artery_tibial.gct.gz` (32.36 MB) - Raw read counts

**Total Size:** ~232.79 MB

**File Format:**
- **GCT (Gene Cluster Text)** format - compressed (.gz)
- **Version:** GTEx v10
- **Two formats per tissue:**
  - **TPM files:** Normalized expression (Transcripts Per Million) - recommended for most analyses
  - **Read count files:** Raw gene expression counts - useful for differential expression analysis

**Coverage:**
- ✅ All 3 required vessel tissues present
- ✅ Both normalized (TPM) and raw counts available
- ✅ Ready for analysis

**Status:** Complete and ready for use as reference transcriptomes.

---

### 📋 5. UK Biobank GWAS Summary Statistics
**Location:** `data/uk_biobank/`

**Status:** Documentation created. Manual download required.

**Download Instructions:** See `data/uk_biobank/README.md`

**Key Sources:**
1. **Neale Lab (Easiest):** http://www.nealelab.is/uk-biobank/
   - Browse available phenotypes
   - Search for "intracranial aneurysm" or vascular phenotypes
   - Download summary statistics

2. **Pan-UK Biobank:** https://pan.ukbb.broadinstitute.org/
   - Download GWAS summary statistics
   - Filter by phenotype

3. **UK Biobank Official:**
   - May require application and approval
   - Visit: https://www.ukbiobank.ac.uk/

**Note:** The Bakker et al. 2020 meta-analysis already includes UK Biobank data, so separate UK Biobank files may not be necessary if you have the Bakker summary statistics.

---

## Next Steps

### Immediate Actions:
1. ✅ **GSE159610:** Ready for analysis
2. ✅ **GWAS (Bakker 2020):** Downloaded - Extract and verify format
3. ✅ **Intracranial Aneurysm RNA-seq (GSE158558 & GSE122897):** Downloaded - Two datasets, includes ruptured vs unruptured comparison
4. ✅ **GTEx:** Complete - All 3 vessel tissues with TPM + raw counts
5. 📋 **UK Biobank:** Download from Neale Lab or Pan-UK Biobank (if needed)

### Recommended Download Order:
1. ✅ **GWAS (Bakker et al. 2020)** - COMPLETE - Extract and verify
2. ✅ **Intracranial Aneurysm RNA-seq (GSE158558 & GSE122897)** - COMPLETE - Two datasets, GSE122897 includes ruptured vs unruptured
3. ✅ **GTEx vessel data** - COMPLETE - All 3 tissues with TPM + raw counts
4. **UK Biobank** - Optional (may be redundant with Bakker data)

### Automated Download Scripts:
- `download_datasets_v2.py` - Main download script
- `download_gse186927.py` - GSE186927 specific downloader
- `download_gwas.py` - GWAS data downloader

---

## File Structure

```
data/
├── blood_rnaseq/
│   └── GSE159610/
│       ├── GSE159610_raw_gene_counts_matrix.txt.gz ✅
│       └── GSE159610_series_matrix.txt.gz ✅
├── sc_rnaseq/
│   ├── GSE158558/
│   │   ├── GSE158558_raw_counts_GRCh38.p13_NCBI.tsv.gz ✅
│   │   ├── GSE158558_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz ✅
│   │   ├── GSE158558_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz ✅
│   │   ├── GSE158558_mRNA_lncRNA_Expression_Transcript.xlsx ✅
│   │   └── Human.GRCh38.p13.annot.tsv.gz ✅
│   └── GSE122897/
│       ├── GSE122897_readCounts_raw.txt.gz ✅
│       ├── GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz ✅
│       ├── GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz ✅
│       └── GSE122897_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz ✅
├── gwas/
│   └── bakker_2020/
│       ├── bakker_2020_summary_stats.txt.gz ✅ (763.6 MB)
│       ├── kp4cd_page.html ✅
│       └── README.md 📋
├── gtex/
│   ├── gene_tpm_v10_artery_aorta.gct.gz ✅
│   ├── gene_reads_v10_artery_aorta.gct.gz ✅
│   ├── gene_tpm_v10_artery_coronary.gct.gz ✅
│   ├── gene_reads_v10_artery_coronary.gct.gz ✅
│   ├── gene_tpm_v10_artery_tibial.gct.gz ✅
│   ├── gene_reads_v10_artery_tibial.gct.gz ✅
│   └── README.md 📋
└── uk_biobank/
    └── README.md 📋
```

---

## Notes

- Some datasets require registration/approval (GTEx dbGaP, UK Biobank full access)
- GSE186927 raw counts may be in EGA (European Genome-Phenome Archive) requiring approval
- GWAS summary statistics are typically available through GWAS Catalog or journal supplements
- All README files contain detailed download instructions with direct links

---

## Contact & Support

For issues with:
- **GEO datasets:** Check NCBI GEO help or contact dataset authors
- **GWAS data:** Check GWAS Catalog or contact study authors
- **GTEx:** Check GTEx Portal documentation
- **UK Biobank:** Check UK Biobank or Neale Lab websites


