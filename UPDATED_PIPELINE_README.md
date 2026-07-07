# Updated Analysis Pipeline - README

## Overview

This updated pipeline integrates multiple datasets for comprehensive aneurysm genomics analysis:
- **GSE122897**: RNA-seq data (ruptured vs unruptured aneurysms)
- **GSE26969**: Microarray data (IA samples vs STA controls)
- **Single-cell atlas** (Maas et al. 2025): Framework ready for when data is available

## Pipeline Structure

### 1. Preprocessing (`preprocess_gse26969.py`)

**Purpose**: Normalize and prepare GSE26969 microarray CEL files

**Requirements**:
- R with `affy` and `Biobase` packages installed
- Extracted CEL.gz files in `results/gse26969_extracted/`

**Steps**:
1. Reads CEL.gz files
2. Performs RMA normalization using R
3. Quality control filtering
4. Creates metadata template

**Output**:
- `results/gse26969_normalized_expression.tsv`: Normalized expression matrix
- `results/gse26969_expression_qc.tsv`: QC-filtered expression matrix
- `results/gse26969_metadata_template.tsv`: Metadata template (needs manual annotation)

**Usage**:
```bash
python preprocess_gse26969.py
```

**Note**: If R is not available, the script will create metadata template only. You can normalize CEL files manually using R or alternative tools.

### 2. Integrated Analysis (`analyze_integrated_datasets.py`)

**Purpose**: Perform differential expression analysis on both platforms and integrate results

**Features**:
- Loads both RNA-seq and microarray data
- Performs platform-specific QC
- Differential expression analysis for each dataset
- Cross-platform validation framework (requires gene mapping)
- Saves integrated results

**Output**:
- `results/gse122897_integrated_de.tsv`: RNA-seq DE results
- `results/gse26969_integrated_de.tsv`: Microarray DE results
- `results/*_top_genes_integrated.tsv`: Top significant genes

**Usage**:
```bash
python analyze_integrated_datasets.py
```

### 3. Single-Cell Integration Framework (`prepare_single_cell_integration.py`)

**Purpose**: Prepare framework for single-cell atlas integration (when data becomes available)

**Features**:
- Loads atlas information from extracted PDF
- Creates deconvolution framework code
- Integration checklist
- Preparation for bulk RNA-seq deconvolution

**Output**:
- `results/single_cell_deconvolution_framework.py`: Framework code
- `results/single_cell_integration_checklist.json`: Integration checklist

**Usage**:
```bash
python prepare_single_cell_integration.py
```

### 4. Main Pipeline Runner (`run_updated_pipeline.py`)

**Purpose**: Run the complete updated pipeline

**Steps**:
1. Preprocess GSE26969
2. Run integrated analysis
3. Set up single-cell framework

**Usage**:
```bash
python run_updated_pipeline.py
```

## Data Requirements

### GSE122897 (RNA-seq)
- **Location**: `data/sc_rnaseq/GSE122897/`
- **Required files**:
  - `GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz`: TPM-normalized expression
  - `sample_metadata.tsv`: Sample metadata with Group column (Ruptured/Unruptured)

### GSE26969 (Microarray)
- **Location**: `results/gse26969_extracted/`
- **Required files**:
  - `GSM663769.CEL.gz` through `GSM663774.CEL.gz`: Raw CEL files
- **Metadata**: `results/gse26969_metadata.tsv` (needs manual creation/update)

### Single-Cell Atlas
- **Status**: Awaiting publication (medRxiv 2025.11.06.25339724)
- **Expected**: GEO accession to be determined
- **Format**: h5ad or count matrix + metadata

## Workflow

### Initial Setup

1. **Extract GSE26969 data** (if not already done):
   ```bash
   # Already extracted to results/gse26969_extracted/
   ```

2. **Update GSE26969 metadata**:
   - Edit `results/gse26969_metadata_template.tsv`
   - Add Group column: 'IA' (intracranial aneurysm) or 'STA' (superficial temporal artery control)
   - Save as `results/gse26969_metadata.tsv`

3. **Install R packages** (for microarray normalization):
   ```r
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install(c("affy", "Biobase"))
   ```

### Running Analysis

**Option 1: Run complete pipeline**
```bash
python run_updated_pipeline.py
```

**Option 2: Run steps individually**
```bash
# Step 1: Preprocess microarray data
python preprocess_gse26969.py

# Step 2: Run integrated analysis
python analyze_integrated_datasets.py

# Step 3: Set up single-cell framework
python prepare_single_cell_integration.py
```

## Output Files

### Preprocessing Outputs
- `results/gse26969_normalized_expression.tsv`: Normalized microarray expression
- `results/gse26969_expression_qc.tsv`: QC-filtered expression
- `results/gse26969_metadata_template.tsv`: Metadata template

### Analysis Outputs
- `results/gse122897_integrated_de.tsv`: RNA-seq differential expression
- `results/gse122897_top_genes_integrated.tsv`: Top RNA-seq genes
- `results/gse26969_integrated_de.tsv`: Microarray differential expression
- `results/gse26969_top_genes_integrated.tsv`: Top microarray genes

### Framework Outputs
- `results/single_cell_deconvolution_framework.py`: Deconvolution code
- `results/single_cell_integration_checklist.json`: Integration checklist

## Key Features

### Cross-Platform Integration
- Handles both RNA-seq and microarray data
- Platform-specific QC and normalization
- Framework for cross-platform validation (requires gene mapping)

### Single-Cell Atlas Preparation
- Framework ready for when atlas data is available
- Expected findings to validate:
  - Ruptured IAs: Higher myeloid cell proportions
  - Sex-specific gene regulation
  - 7 distinct cell types identified

### Quality Control
- **RNA-seq**: Filters genes with TPM > 1 in >= 3 samples
- **Microarray**: Filters probes with expression > 5 in >= 2 samples
- Platform-appropriate normalization

## Troubleshooting

### R/Affy Not Available
- **Issue**: Microarray normalization requires R
- **Solution**: 
  - Install R and required packages
  - Or use alternative normalization tools
  - Pipeline will still work for RNA-seq analysis

### Missing Metadata
- **Issue**: Metadata files not found
- **Solution**: 
  - Create metadata files from templates
  - Ensure Group column has proper labels (Ruptured/Unruptured for GSE122897, IA/STA for GSE26969)

### Gene Mapping
- **Issue**: Cross-platform comparison requires gene mapping
- **Solution**: 
  - Annotate microarray probes with gene symbols
  - Map Ensembl IDs (RNA-seq) to gene symbols
  - Use common gene symbols for comparison

## Next Steps

1. **Complete GSE26969 metadata**: Annotate samples as IA or STA
2. **Run integrated analysis**: Compare findings across platforms
3. **Monitor for single-cell atlas**: Check GEO for data release
4. **Gene annotation**: Map probes/IDs to gene symbols for cross-platform validation
5. **Pathway analysis**: Integrate DE results with pathway enrichment
6. **Validation**: Compare findings with published literature

## Integration with Existing Pipeline

This updated pipeline:
- **Maintains compatibility** with existing `analyze_gse122897.py`
- **Extends functionality** to include microarray data
- **Prepares framework** for single-cell integration
- **Follows same structure** for consistency

The original `analyze_gse122897.py` can still be run independently for RNA-seq-only analysis.

## References

- **GSE122897**: Ruptured vs unruptured aneurysm RNA-seq
- **GSE26969**: IA vs STA microarray data
- **Maas et al. 2025**: Single-cell atlas (medRxiv 2025.11.06.25339724)
  - DOI: 10.1101/2025.11.06.25339724
  - Key finding: Ruptured IAs have higher myeloid cell proportions

---

*Last updated: After dataset extraction and pipeline update*






