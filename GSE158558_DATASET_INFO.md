# GSE158558 Dataset Information

## ✅ Dataset Successfully Added to Project

**Accession:** GSE158558  
**Title:** Intracranial aneurysm sequencing data  
**Organism:** Homo sapiens (Human)  
**Type:** Bulk RNA-seq (Expression profiling by high throughput sequencing)

## 📁 Files Downloaded (5 files, ~10.6 MB)

### Expression Data Files:

1. **GSE158558_raw_counts_GRCh38.p13_NCBI.tsv.gz** (0.43 MB)
   - Raw gene counts
   - Use for: Differential expression analysis, statistical modeling

2. **GSE158558_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz** (0.67 MB)
   - TPM (Transcripts Per Million) normalized expression
   - Use for: Cross-sample comparisons, visualization

3. **GSE158558_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz** (0.68 MB)
   - FPKM (Fragments Per Kilobase Million) normalized expression
   - Use for: Alternative normalization for comparisons

4. **GSE158558_mRNA_lncRNA_Expression_Transcript.xlsx** (4.62 MB)
   - Combined mRNA and lncRNA expression matrix
   - Use for: Comprehensive expression analysis including non-coding RNAs

### Annotation File:

5. **Human.GRCh38.p13.annot.tsv.gz** (4.21 MB)
   - Gene annotation file (GRCh38.p13 reference)
   - Use for: Gene ID mapping, annotation, functional analysis

## 🎯 How This Dataset Fits the Project

### ✅ Advantages:
- **Human tissue data** - Correct organism for intracranial aneurysm study
- **Tissue-specific** - Direct aneurysm tissue expression (not blood)
- **Comprehensive** - Includes both mRNA and lncRNA
- **Multiple formats** - Raw counts and normalized data available
- **Ready to use** - Processed data ready for analysis

### ⚠️ Note:
- This is **bulk RNA-seq** (not single-cell)
- Cannot resolve individual cell types directly
- Still valuable for identifying tissue-level gene expression patterns
- Can identify endothelial dysfunction markers at tissue level

## 📊 Use Cases in Project

1. **Endothelial Dysfunction Analysis:**
   - Identify dysregulated endothelial genes in aneurysm tissue
   - Compare with GTEx normal vessel expression
   - Find aneurysm-specific expression signatures

2. **Integration with Other Datasets:**
   - Compare tissue expression with blood RNA-seq (GSE159610)
   - Validate GWAS-identified risk genes
   - Cross-reference with GTEx normal vessel data

3. **Biomarker Discovery:**
   - Identify aneurysm-associated genes
   - Find rupture risk signatures
   - Discover novel therapeutic targets

4. **lncRNA Analysis:**
   - Analyze long non-coding RNA expression
   - Identify regulatory mechanisms
   - Discover novel regulatory RNAs

## 🔄 Comparison with Original Plan

**Original Plan:** Single-cell RNA-seq (GSE186927 - which was incorrect)  
**Current Status:** Bulk RNA-seq (GSE158558 - correct human dataset)

**Impact:**
- ✅ Still provides aneurysm tissue expression data
- ✅ Can identify tissue-level gene signatures
- ⚠️ Cannot resolve individual cell types (endothelial, smooth muscle, etc.)
- ✅ Alternative approach: Use bulk data + cell type deconvolution methods
- ✅ Or focus on tissue-level endothelial dysfunction markers

## 📝 Next Steps

1. **Extract and explore data:**
   - Load expression matrices
   - Check sample metadata
   - Verify data quality

2. **Analysis:**
   - Differential expression analysis
   - Compare with GTEx normal vessels
   - Identify endothelial dysfunction markers
   - Integrate with blood RNA-seq data

3. **Documentation:**
   - Create analysis scripts
   - Document findings
   - Update project status

## 📍 File Location

All files are located in: `data/sc_rnaseq/GSE158558/`

## ✅ Status

**Dataset Status:** ✅ Complete and ready for analysis  
**Integration Status:** ✅ Added to project  
**Documentation Status:** ✅ Updated






