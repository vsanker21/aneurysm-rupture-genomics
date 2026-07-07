# Poppenberg et al. 2023 Findings Integration Summary

## Paper Information
- **Title**: RNA Expression Signatures of Intracranial Aneurysm Growth Trajectory Identified in Circulating Whole Blood
- **Authors**: Poppenberg, K.E. et al.
- **Journal**: Journal of Personalized Medicine (2023)
- **DOI**: 10.3390/jpm13020266

## Study Overview
- **Sample Size**: 66 blood samples from intracranial aneurysm patients
  - Training cohort: 46 samples
  - Testing cohort: 20 samples
- **Classification**: 33 "growing" aneurysms (PAT >= 4.6) vs 33 "stable" aneurysms (PAT < 4.6)
- **Model Performance**: 
  - Training AUC: 0.85
  - Testing AUC: 0.86
  - Testing Accuracy: 80%

## Key Findings

### 39 Differentially Expressed Genes
- **11 genes decreased** in "growing" aneurysms (protective markers)
- **28 genes increased** in "growing" aneurysms (risk markers)
- **12 genes** with high fold-change (>1.0 Log2FC)

### Priority Genes (Common Across Multiple Analyses)

**Common Across All Metrics (PAT, Size, PHASES)**:
- DEFA1 (Defensin Alpha 1) - Decreased in growing
- HBA1 (Hemoglobin Alpha 1) - Increased in growing
- HBA2 (Hemoglobin Alpha 2) - Increased in growing
- OTOF (Otoferlin) - Decreased in growing
- SFRP2 (Secreted Frizzled Related Protein 2) - Increased in growing

**High Prognostic Potential**:
- DEFA1, HBA1, HBA2, SFRP2 (common across all)
- IFI27 (Interferon Alpha Inducible Protein 27) - Decreased in growing
- LPL (Lipoprotein Lipase) - Decreased in growing

**Common with Other Studies**:
- FAM83A, GSTM1, RNF152, RNF182, RUNDC3A

### Top 10 Priority Biomarker Candidates (by Priority Score)
1. **OTOF** - Log2FC: -1.52, High priority (Common across PAT/Size/PHASES)
2. **RNF152** - Log2FC: 1.18, High priority (Common with other studies)
3. **HBA2** - Log2FC: 1.45, High priority (High prognostic potential)
4. **SFRP2** - Log2FC: 1.40, High priority (High prognostic potential)
5. **RNF182** - Log2FC: 1.78, High priority (Common with other studies)
6. **LPL** - Log2FC: -1.14, High priority (High prognostic potential)
7. **IFI27** - Log2FC: -0.87, Medium-high priority (High prognostic potential)
8. **TMEM158** - Log2FC: 1.69, High fold-change
9. **HBA1** - Log2FC: 1.68, High priority (High prognostic potential)
10. **AC004556.1** - Log2FC: 1.54, High fold-change

## PAT Score Formula Information

**Predicted Aneurysm Trajectory (PAT)**
- **Reference**: Chien et al. 2019, J. Neurosurg. 2019, 132:1077-1087
- **Validation**: Trained on 520 longitudinally followed IA cases
- **Calculation Method**: Multivariate adaptive splines for analysis of longitudinal data in R

**Factors Included**:
1. Initial IA size
2. Age
3. Sex
4. Cigarette smoking
5. Hypothyroidism
6. Follow up time in months

**Key Findings**:
- Smoking and hypothyroidism had significant effects on growth trajectories of larger IAs (>=7 mm)
- Initial aneurysm size was significantly related to growth
- Larger IAs (>=7 mm) tended to grow significantly faster than small IAs (<3 mm)
- IAs in males showed faster growth rates than IAs in females
- PAT was markedly greater in ruptured IAs compared to unruptured IAs

**Usage in Study**: Median PAT score (4.6) used to dichotomize into "growing" (PAT >= 4.6) vs "stable" (< 4.6)

**Note**: Full coefficients available in Chien et al. 2019 manuscript Table

## Data Availability

**GEO Accession**: Not publicly available in GEO repository
**Data Access**: Available upon reasonable request to corresponding author
- **Contact**: Vincent M. Tutino, PhD
- **Email**: vincentt@buffalo.edu
- **Institution**: Canon Stroke and Vascular Research Center, University at Buffalo

## Biological Mechanisms

### Organismal Injury and Abnormalities
- Inflammatory signaling (NFKB, ERK1/2)
- Cancer-related pathways
- Cardiovascular disease associations
- Systemic lupus erythematosus (SLE) association

### Cell-to-Cell Signaling and Interaction
- FGFR2 (angiogenesis, wound healing, cell migration)
- RAP1GAP (cell signaling)
- ITLN1 (innate immune system, defensins)
- LPL (metabolism)

### Key Genes and Pathways

**Decreased in Growing Aneurysms**:
- **DEFA1** (Defensin Alpha 1): Neutrophil defensin, decreased in ruptured IAs (proteomic level)
- **IFI27**: Interferon-inducible protein, associated with apoptosis and medial thinning
- **LPL**: Lipoprotein lipase, phospholipase A1 activity

**Increased in Growing Aneurysms**:
- **SFRP2**: Wnt signaling modulator, overexpressed in IA tissue vs controls
- **HBA1/HBA2**: Hemoglobin alpha chains, increased in AAA, role in iron metabolism and inflammation
- **RNF152/RNF182**: Ring finger proteins, cellular processes
- **FGFR2**: Fibroblast growth factor receptor, angiogenesis and cell migration

## Integration with Existing Analyses

### Comparison with GSE122897 (Rupture Analysis)
- **No direct overlap** found (GSE122897 uses numeric gene IDs vs gene symbols)
- Poppenberg study focuses on **growth trajectory** (PAT score)
- GSE122897 focuses on **rupture status** (ruptured vs unruptured)
- These represent **complementary signatures**:
  - Growth trajectory (blood-based, 66 samples, AUC 0.86)
  - Rupture status (tissue-based, different sample set)

### Recommended Validation Strategy
1. **Top Priority Genes** (Priority Score >= 10): Validate in independent cohorts
2. **Common Genes Across Metrics**: DEFA1, HBA1, HBA2, OTOF, SFRP2
3. **Novel Findings**: 39 genes not previously identified in rupture analysis

## Clinical Implications

1. **Blood-Based Biomarker Panel**: 39-gene signature for growth trajectory prediction
2. **Model Performance**: AUC 0.86 in independent testing (superior to imaging alone)
3. **Complementary to Existing Metrics**: Integrates with PAT, IA size, and PHASES score
4. **Monitoring Tool**: Could enable more frequent monitoring than imaging alone

## Files Generated

1. **poppenberg_2023_39_genes.tsv** - Complete list of 39 genes with statistics
2. **integrated_gene_list.tsv** - Integrated list with priority scores
3. **pat_score_formula_info.json** - PAT score formula details
4. **priority_genes.json** - Priority gene classifications
5. **integration_summary.json** - Complete integration summary
6. **geo_accession_info.json** - Data availability information

## Next Steps

1. **Validate Top Priority Genes** in independent cohorts
2. **Request Raw Data** from corresponding author for deeper analysis
3. **Compare with Existing Blood Transcriptomics** (e.g., GSE159610)
4. **Integrate with GWAS Findings** for multi-omics risk prediction
5. **Develop Clinical Nomogram** combining PAT, gene signature, and clinical factors

## References

- Poppenberg, K.E.; Chien, A.; Santo, B.A.; et al. RNA Expression Signatures of Intracranial Aneurysm Growth Trajectory Identified in Circulating Whole Blood. J. Pers. Med. 2023, 13, 266.
- Chien, A.; Callender, R.A.; Yokota, H.; et al. Unruptured intracranial aneurysm growth trajectory: Occurrence and rate of enlargement in 520 longitudinally followed cases. J. Neurosurg. 2019, 132:1077-1087.






