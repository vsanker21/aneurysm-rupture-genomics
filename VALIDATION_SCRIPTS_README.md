# Validation and Integration Scripts

This document describes the scripts created for pathway enrichment, gene validation, functional annotation, and integration with GWAS and blood transcriptomics data.

## Overview

Five main scripts have been created to complete the validation and integration pipeline:

1. **Pathway Enrichment Analysis** (`pathway_enrichment_analysis.py`)
2. **Gene Validation** (`validate_de_genes.py`)
3. **Functional Annotation** (`gene_function_annotation.py`) - Enhanced
4. **GWAS Integration** (`integrate_gwas.py`)
5. **Blood Transcriptomics Integration** (`integrate_blood_transcriptomics.py`)

Additionally, a master script (`run_validation_pipeline.py`) runs all analyses sequentially.

---

## Scripts Description

### 1. Pathway Enrichment Analysis

**File:** `pathway_enrichment_analysis.py`

**Purpose:** Performs pathway enrichment analysis on differentially expressed genes using the g:Profiler API.

**Features:**
- Analyzes pathways at multiple FDR thresholds (0.05, 0.1, 0.2)
- Uses multiple pathway databases (GO:BP, GO:MF, GO:CC, KEGG, Reactome, WikiPathways)
- Identifies enriched pathways for both upregulated and downregulated genes
- Generates pathway enrichment tables with statistics

**Input:**
- `results/gse122897_differential_expression.tsv`

**Output:**
- `results/pathway_enrichment_fdr_0.05.tsv`
- `results/pathway_enrichment_fdr_0.1.tsv`
- `results/pathway_enrichment_fdr_0.2.tsv`
- `results/pathway_enrichment_summary.tsv`

**Usage:**
```bash
python pathway_enrichment_analysis.py
```

---

### 2. Gene Validation

**File:** `validate_de_genes.py`

**Purpose:** Validates DE genes against public databases, checks pathway memberships, and creates comprehensive validation reports.

**Features:**
- Validates gene IDs using MyGene.info API
- Maps genes to enriched pathways
- Checks gene symbols, names, and database IDs
- Creates validation reports with gene annotations
- Identifies genes with pathway associations

**Input:**
- `results/gse122897_differential_expression.tsv`
- Pathway enrichment results (if available)

**Output:**
- `results/gene_validation_report_fdr_0.05.tsv`
- `results/gene_validation_report_fdr_0.1.tsv`
- `results/gene_validation_report_fdr_0.2.tsv`

**Output Columns:**
- Gene_ID, Gene_Symbol, Gene_Name, Entrez_ID, Ensembl_ID
- Gene_Type, Validated status
- Associated_Pathways, Pathway_Count
- Log2FC, P_value, P_adj_BH
- Summary

**Usage:**
```bash
python validate_de_genes.py
```

---

### 3. Functional Annotation (Enhanced)

**File:** `gene_function_annotation.py`

**Purpose:** Comprehensive functional annotation of DE genes with GO terms, pathways, disease associations, and protein information.

**Features:**
- Annotates genes using MyGene.info API with extended fields
- Retrieves GO terms (Biological Process, Molecular Function, Cellular Component)
- Identifies KEGG and Reactome pathways
- Includes OMIM disease associations
- Provides UniProt protein IDs
- Adds gene aliases and types

**Input:**
- `results/gse122897_differential_expression.tsv`

**Output:**
- `results/gse122897_annotated_fdr_0.05.tsv`
- `results/gse122897_annotated_fdr_0.1.tsv`
- `results/gse122897_annotated_fdr_0.2.tsv`

**Output Columns:**
- Gene, Gene_Symbol, Gene_Name, Entrez_ID, Ensembl_ID, UniProt_ID
- Gene_Type, Aliases
- GO_Terms, Pathways, OMIM_Disease
- Gene_Summary
- Log2FC, P_value, P_adj_BH
- Mean expression values

**Usage:**
```bash
python gene_function_annotation.py
```

---

### 4. GWAS Integration

**File:** `integrate_gwas.py`

**Purpose:** Maps differentially expressed genes to GWAS variants and identifies genes that are both differentially expressed and in GWAS risk loci.

**Features:**
- Loads GWAS summary statistics (Bakker et al. 2020)
- Maps DE genes to GWAS significant regions
- Identifies prioritized genes (DE + GWAS)
- Calculates priority scores based on effect size and significance
- Maps to known IA GWAS loci (9p21.3, 8q11.23, 4q31, etc.)

**Input:**
- `results/gse122897_differential_expression.tsv`
- `data/gwas/bakker_2020/*.txt.gz` (GWAS summary statistics)

**Output:**
- `results/gwas_integration_report.tsv`

**Output Columns:**
- Gene_ID
- Log2FC, P_adj_BH
- Priority_Score
- In_GWAS_Locus (boolean)
- GWAS_Association

**Usage:**
```bash
python integrate_gwas.py
```

**Note:** Full gene-GWAS coordinate mapping requires:
1. Gene annotation file (GTF/GFF) with chromosome coordinates
2. GWAS summary statistics with chromosome and position columns
3. LD reference panel for fine-mapping (optional)

---

### 5. Blood Transcriptomics Integration

**File:** `integrate_blood_transcriptomics.py`

**Purpose:** Compares tissue DE genes with blood transcriptomics signatures to identify biomarker candidates.

**Features:**
- Loads blood RNA-seq data from GSE159610
- Performs differential expression analysis in blood samples
- Integrates tissue and blood signatures
- Identifies genes dysregulated in both tissue and blood
- Calculates biomarker scores
- Flags concordant direction of change
- Creates biomarker candidate lists

**Input:**
- `results/gse122897_differential_expression.tsv`
- `data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz`
- `data/blood_rnaseq/GSE159610/GSE159610_series_matrix.txt.gz` (metadata)

**Output:**
- `results/blood_de_analysis_fdr_*.tsv`
- `results/blood_tissue_integration_fdr_*.tsv`
- `results/biomarker_candidates_fdr_*.tsv`

**Output Columns:**
- Gene_ID
- Log2FC_Tissue, P_adj_BH_Tissue
- Log2FC_Blood, P_adj_BH_Blood
- In_Blood_Dataset (boolean)
- Dysregulated_in_Both (boolean)
- Concordant_Direction (boolean)
- Biomarker_Score

**Usage:**
```bash
python integrate_blood_transcriptomics.py
```

---

## Master Pipeline Script

**File:** `run_validation_pipeline.py`

**Purpose:** Runs all validation and integration scripts sequentially.

**Usage:**
```bash
python run_validation_pipeline.py
```

**Features:**
- Runs all five analysis scripts in order
- Provides progress updates
- Generates summary report
- Handles errors gracefully
- 10-minute timeout per script

---

## Prerequisites

1. **Required Python packages:**
   - pandas
   - numpy
   - scipy
   - statsmodels
   - requests

2. **Required data files:**
   - `results/gse122897_differential_expression.tsv` (from `analyze_gse122897.py`)
   - GWAS summary statistics (if running GWAS integration)
   - Blood transcriptomics data (if running blood integration)

3. **Internet connection:**
   - Required for API calls to g:Profiler, MyGene.info, etc.
   - Scripts include rate limiting to avoid API throttling

---

## Running the Complete Pipeline

To run all analyses sequentially:

```bash
python run_validation_pipeline.py
```

To run individual scripts:

```bash
# Pathway enrichment
python pathway_enrichment_analysis.py

# Gene validation
python validate_de_genes.py

# Functional annotation
python gene_function_annotation.py

# GWAS integration
python integrate_gwas.py

# Blood transcriptomics integration
python integrate_blood_transcriptomics.py
```

---

## Output Files Summary

All output files are saved in the `results/` directory:

### Pathway Enrichment:
- `pathway_enrichment_fdr_0.05.tsv`
- `pathway_enrichment_fdr_0.1.tsv`
- `pathway_enrichment_fdr_0.2.tsv`
- `pathway_enrichment_summary.tsv`

### Gene Validation:
- `gene_validation_report_fdr_0.05.tsv`
- `gene_validation_report_fdr_0.1.tsv`
- `gene_validation_report_fdr_0.2.tsv`

### Functional Annotation:
- `gse122897_annotated_fdr_0.05.tsv`
- `gse122897_annotated_fdr_0.1.tsv`
- `gse122897_annotated_fdr_0.2.tsv`

### GWAS Integration:
- `gwas_integration_report.tsv`

### Blood Transcriptomics:
- `blood_de_analysis_fdr_*.tsv`
- `blood_tissue_integration_fdr_*.tsv`
- `biomarker_candidates_fdr_*.tsv`

---

## Notes and Limitations

1. **API Rate Limiting:** Scripts include delays to avoid API throttling, but may still be slow for large gene lists.

2. **Gene ID Mapping:** Some scripts assume Entrez IDs. If your data uses different IDs, you may need to adjust the mapping logic.

3. **GWAS Mapping:** Full coordinate-based mapping requires gene annotation files (GTF/GFF) which may need to be downloaded separately.

4. **Blood Data:** Blood transcriptomics integration requires proper sample metadata. If metadata is missing, the script will create a placeholder structure.

5. **Pathway Enrichment:** Results depend on g:Profiler API availability. If the API is down, enrichment will fail.

---

## Next Steps

After running all scripts, you can:

1. Review enriched pathways for biological insights
2. Prioritize genes based on validation scores
3. Validate biomarker candidates in blood transcriptomics
4. Cross-reference GWAS loci with DE genes
5. Integrate all results for manuscript preparation

---

## Troubleshooting

**Issue: API errors/timeouts**
- Solution: Check internet connection, wait a few minutes and retry (rate limiting)

**Issue: Missing input files**
- Solution: Run `analyze_gse122897.py` first to generate DE results

**Issue: GWAS file not found**
- Solution: Check that GWAS files are in `data/gwas/bakker_2020/` directory

**Issue: Blood data parsing errors**
- Solution: Check that GSE159610 files are properly downloaded and in correct format

---

## Contact

For questions or issues with these scripts, please refer to the main project documentation.






