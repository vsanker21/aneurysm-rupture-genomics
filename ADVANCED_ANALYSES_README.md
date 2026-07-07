# Advanced Analyses - Setup and Usage Guide

## Overview

This guide covers four advanced analyses for the blood RNA-seq data:

1. **Gene ID Conversion** - Map Ensembl IDs to gene symbols
2. **DESeq2 Analysis** - Proper count-based differential expression
3. **Literature Review** - Annotate top candidate genes
4. **Machine Learning** - Build classification models

---

## 1. Gene ID Conversion ✅

### Purpose
Convert Ensembl gene IDs to gene symbols for better pathway enrichment and annotation.

### Script
```bash
python convert_ensembl_to_symbols.py
```

### What it does
- Converts Ensembl IDs to gene symbols using mygene API
- Adds Entrez IDs and gene names
- Creates files ready for pathway enrichment

### Output Files
- `results/gse159610_substantial_fc_genes_annotated.tsv` - All genes with symbols
- `results/gse159610_top_candidate_genes_annotated.tsv` - Top candidates with symbols
- `results/gse159610_genes_with_symbols.tsv` - Genes with successful conversion
- `results/gse159610_gene_symbols_list.txt` - List of gene symbols for pathway analysis

### Time Required
~5-10 minutes (depends on API rate limits)

---

## 2. DESeq2 Analysis ⚠️

### Purpose
Run proper count-based differential expression analysis (better than t-tests for RNA-seq).

### Prerequisites
- **R** must be installed: https://www.r-project.org/
- **DESeq2** R package (will be auto-installed if R is available)

### Script
```bash
python setup_deseq2_analysis.py
```

### What it does
1. Checks if R is installed
2. Checks if DESeq2 is installed (installs if needed)
3. Creates R script for DESeq2 analysis
4. Runs DESeq2 analysis

### Manual Installation (if auto-install fails)
```r
# In R or RStudio:
install.packages("BiocManager")
BiocManager::install("DESeq2")
```

### Output Files
- `results/gse159610_deseq2_results.tsv` - DESeq2 DE results
- `deseq2_blood_analysis.R` - R script (can be run manually)

### Time Required
- Installation: 5-10 minutes (first time only)
- Analysis: 5-10 minutes

### Alternative: Run R Script Manually
```bash
Rscript deseq2_blood_analysis.R
```

---

## 3. Literature Review Annotation ✅

### Purpose
Annotate top candidate genes with:
- Gene information (name, summary, aliases)
- PubMed search for aneurysm-related publications

### Prerequisites
- Gene ID conversion must be completed first

### Script
```bash
python annotate_top_genes_literature.py
```

### What it does
- Fetches gene information from mygene API
- Searches PubMed for aneurysm-related publications
- Creates annotated gene list

### Output Files
- `results/gse159610_top_genes_literature_annotated.tsv` - All annotated genes
- `results/gse159610_genes_with_pubmed_results.tsv` - Genes with PubMed matches

### Time Required
~10-15 minutes (depends on API rate limits)

---

## 4. Machine Learning Classification ✅

### Purpose
Build classification models to distinguish aneurysm vs control samples using top DE genes.

### Script
```bash
python ml_classification_blood_genes.py
```

### What it does
- Selects top 50 DE genes as features
- Builds 3 models: Random Forest, SVM, Logistic Regression
- Performs cross-validation
- Reports feature importance

### Output Files
- `results/blood_ml_feature_importance.tsv` - Feature importance (Random Forest)
- `results/blood_ml_cv_results.tsv` - Cross-validation results
- `results/blood_ml_summary.tsv` - Summary statistics

### Time Required
~2-5 minutes

---

## Running All Analyses

### Option 1: Run All at Once
```bash
python run_all_advanced_analyses.py
```

This will run all analyses in sequence (may take 30-60 minutes total).

### Option 2: Run Individually
Run each script separately as needed:
```bash
# 1. Convert gene IDs
python convert_ensembl_to_symbols.py

# 2. Pathway enrichment (after conversion)
python pathway_enrichment_blood_genes_symbols.py

# 3. Literature annotation (after conversion)
python annotate_top_genes_literature.py

# 4. Machine learning
python ml_classification_blood_genes.py

# 5. DESeq2 (if R is available)
python setup_deseq2_analysis.py
```

---

## Expected Results

### Gene ID Conversion
- **~80-90%** of genes should successfully convert to symbols
- Remaining genes may be pseudogenes, non-coding RNAs, or deprecated IDs

### DESeq2 Analysis
- Should find **0-50 genes** with padj < 0.1 (due to small sample size)
- May find more genes than t-test analysis
- Results will be in different format (log2FoldChange, padj)

### Literature Review
- **~10-30%** of top genes may have aneurysm-related PubMed publications
- Gene summaries provide functional context

### Machine Learning
- **AUC scores: 0.6-0.8** expected (moderate performance due to small sample size)
- Random Forest typically performs best
- Feature importance identifies most predictive genes

---

## Troubleshooting

### Gene ID Conversion Fails
- **Issue:** API rate limits or network issues
- **Solution:** Wait and retry, or reduce batch size in script

### DESeq2 Not Installing
- **Issue:** R not found or BiocManager installation fails
- **Solution:** Install R manually, then run installation commands in R

### Pathway Enrichment Returns No Results
- **Issue:** Gene symbols not recognized or no enrichment
- **Solution:** Check gene symbol list, try less stringent FDR threshold

### Machine Learning Low Performance
- **Issue:** Small sample size (n=47) limits model performance
- **Solution:** This is expected; focus on feature importance rather than absolute performance

---

## Next Steps After Analyses

1. **Review pathway enrichment results** - Identify biological pathways
2. **Check literature annotations** - Focus on genes with PubMed matches
3. **Compare DESeq2 vs t-test results** - See if DESeq2 finds more significant genes
4. **Use ML feature importance** - Identify most predictive genes for biomarkers
5. **Integrate all results** - Combine findings for comprehensive gene prioritization

---

## Files Created

### Scripts
- `convert_ensembl_to_symbols.py` - Gene ID conversion
- `pathway_enrichment_blood_genes_symbols.py` - Pathway enrichment
- `annotate_top_genes_literature.py` - Literature annotation
- `ml_classification_blood_genes.py` - Machine learning
- `setup_deseq2_analysis.py` - DESeq2 setup
- `run_all_advanced_analyses.py` - Master script
- `deseq2_blood_analysis.R` - R script for DESeq2

### Results (will be created when scripts run)
- Multiple annotated gene files
- Pathway enrichment results
- Literature annotations
- ML classification results
- DESeq2 results (if available)

---

## Support

If you encounter issues:
1. Check that required data files exist in `data/blood_rnaseq/GSE159610/`
2. Ensure internet connection for API calls
3. Check Python package requirements (pandas, numpy, scikit-learn, requests)
4. For DESeq2, ensure R is properly installed






