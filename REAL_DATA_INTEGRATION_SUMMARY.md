# Multi-Omics Integration with Real Data - Summary

## Overview

Created comprehensive integration script to combine real patient data from multiple datasets for aneurysm rupture prediction.

## Available Real Datasets

### 1. Blood Transcriptomics (GSE159610)
- **Source**: Poppenberg et al. 2020
- **Samples**: 47 (23 aneurysm, 24 control)
- **Platform**: Illumina HiSeq 2500
- **Location**: `data/blood_rnaseq/GSE159610/`
- **Files**:
  - Raw counts: `GSE159610_raw_gene_counts_matrix.txt.gz`
  - Metadata: `sample_metadata.tsv`
  - DE Results: `results/gse159610_full_de_results.tsv`

### 2. Tissue Transcriptomics (GSE122897)
- **Source**: Ruptured vs Unruptured aneurysm tissue
- **Samples**: Multiple samples (ruptured vs unruptured)
- **Platform**: RNA-seq
- **Location**: `data/sc_rnaseq/GSE122897/`
- **Files**:
  - TPM counts: `GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz`
  - Metadata: `sample_metadata.tsv`
  - DE Results: `results/gse122897_differential_expression.tsv`

### 3. GWAS Data (Phase 2 PRS)
- **Source**: Bakker et al. 2020
- **PRS Weights**: `results/prs_weights.tsv`
- **Risk Loci**: `results/prs_risk_loci_top17.tsv`

## Integration Script: `integrate_real_data_multi_omics.py`

### Features:

1. **Loads Real Expression Data**
   - Blood RNA-seq (GSE159610)
   - Tissue RNA-seq (GSE122897)
   - Pre-computed differential expression results

2. **Identifies Signature Genes**
   - Blood signature: Aneurysm vs Control
   - Tissue signature: Ruptured vs Unruptured
   - Combines both signatures for multi-omics features

3. **Creates Feature Matrix**
   - Extracts expression of signature genes for each sample
   - Handles different gene ID formats (ENSG, numeric)
   - Creates separate features for blood and tissue datasets

4. **Builds Prediction Models**
   - Gradient Boosting Classifier
   - Random Forest Classifier
   - 5-fold cross-validation
   - Feature importance analysis

### Output Files:

1. **`results/multi_omics_signature_genes.tsv`**
   - Combined signature genes from blood and tissue
   - Includes source (Blood/Tissue) and comparison type

2. **`results/multi_omics_features_real_data.tsv`**
   - Feature matrix with expression values
   - Sample IDs, dataset source, labels

3. **`results/multi_omics_predictions_real_data.tsv`**
   - Model predictions for each sample
   - Predicted probabilities from both models

4. **`results/multi_omics_feature_importance_real_data.tsv`**
   - Feature importance scores from both models
   - Helps identify most predictive genes

## Usage

```bash
python integrate_real_data_multi_omics.py
```

The script will:
1. Load expression data from both datasets
2. Load pre-computed DE results
3. Identify signature genes
4. Create multi-omics feature matrix
5. Train and validate prediction models
6. Save all results to `results/` directory

## Next Steps

1. **Run the integration script** to generate predictions
2. **Evaluate model performance** (AUC, sensitivity, specificity)
3. **Compare with PRS-only predictions** (when genotype data available)
4. **Integrate clinical factors** (hypertension, smoking, age, sex)
5. **Validate on external cohorts** if available

## Notes

- The script handles different gene ID formats automatically
- Uses P-value threshold for blood data (more lenient due to small sample size)
- Uses FDR threshold for tissue data (more robust with larger sample size)
- Feature importance helps identify key predictive genes
- Models can be saved and reused for new predictions

---

**Status**: Script created and ready to run

**Target**: AUC 0.85-0.92 for integrated multi-omics model


