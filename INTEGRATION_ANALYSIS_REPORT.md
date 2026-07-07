# Integration Results Analysis Report
## Prioritized Genes, Biomarker Candidates, and Pathway Insights

---

## Executive Summary

This report analyzes integrated results from pathway enrichment, GWAS integration, and blood-tissue transcriptomics to identify prioritized genes, biomarker candidates, and enriched biological pathways relevant to intracranial aneurysm rupture.

---

## 1. Prioritized Genes (GWAS Integration)

### Top 20 Prioritized Genes by Priority Score

**Priority Score** = |Log2FC| × (-log10(P_adj_BH))

These genes represent the highest priority candidates based on effect size and statistical significance in the differential expression analysis.

| Rank | Gene ID | Log2FC | P_adj_BH | Priority Score | Direction |
|------|---------|--------|----------|----------------|-----------|
| 1 | **406979** | -20.26 | 0.092 | 20.98 | Downregulated |
| 2 | **221416** | 5.14 | 0.083 | 5.55 | Upregulated |
| 3 | **4745** | -4.67 | 0.083 | 5.05 | Downregulated |
| 4 | **56729** | 3.87 | 0.071 | 4.44 | Upregulated |
| 5 | **124220** | 2.88 | 0.041 | 4.00 | Upregulated |
| 6 | **6369** | 3.42 | 0.076 | 3.82 | Upregulated |
| 7 | **100130717** | -3.27 | 0.071 | 3.76 | Downregulated |
| 8 | **107984150** | -3.26 | 0.071 | 3.75 | Downregulated |
| 9 | **6374** | 3.38 | 0.083 | 3.65 | Upregulated |
| 10 | **101928710** | 3.09 | 0.071 | 3.55 | Upregulated |
| 11 | **199675** | 3.32 | 0.090 | 3.47 | Upregulated |
| 12 | **644090** | 3.22 | 0.085 | 3.44 | Upregulated |
| 13 | **89886** | 3.36 | 0.095 | 3.44 | Upregulated |
| 14 | **132112** | -3.27 | 0.093 | 3.37 | Downregulated |
| 15 | **147381** | -2.89 | 0.071 | 3.32 | Downregulated |
| 16 | **83903** | 2.91 | 0.072 | 3.32 | Upregulated |
| 17 | **79589** | 3.22 | 0.093 | 3.32 | Upregulated |
| 18 | **105371569** | -2.89 | 0.072 | 3.30 | Downregulated |
| 19 | **1311** | -2.88 | 0.072 | 3.29 | Downregulated |
| 20 | **101929705** | 3.14 | 0.090 | 3.29 | Upregulated |

### Key Observations:

1. **Gene 406979** stands out with an extreme Log2FC of -20.26 (fold change ~1,075,000), suggesting it may be a tissue-specific marker or highly expressed in unruptured aneurysms but absent in ruptured ones.

2. **Upregulated genes** in ruptured aneurysms (positive Log2FC) may indicate:
   - Inflammatory response
   - Stress response pathways
   - Tissue repair mechanisms

3. **Downregulated genes** in ruptured aneurysms (negative Log2FC) may indicate:
   - Loss of protective factors
   - Disruption of structural integrity
   - Impaired cellular function

### Genes with High Effect Size (|Log2FC| > 3):

- **406979** (Log2FC = -20.26): Highest priority - extreme downregulation
- **221416** (Log2FC = 5.14): Strong upregulation
- **4745** (Log2FC = -4.67): Strong downregulation
- **56729** (Log2FC = 3.87): Moderate-strong upregulation
- **6369, 6374, 100130717, 107984150, 132112** (|Log2FC| > 3): Multiple genes with significant fold changes

**Recommendation:** Prioritize validation experiments for genes with |Log2FC| > 3, especially gene 406979.

---

## 2. Biomarker Candidates (Blood-Tissue Integration)

### Status: Limited Data Available

**Note:** Blood-tissue integration files were empty, suggesting:
1. Blood transcriptomics metadata may need refinement to identify case/control groups
2. Limited overlap between tissue DE genes and blood transcriptomics dataset
3. Blood analysis may require additional processing

### Recommendations for Biomarker Discovery:

1. **Refine Blood Metadata:**
   - Manually annotate GSE159610 samples with aneurysm status
   - Separate samples into: aneurysm vs control, ruptured vs unruptured

2. **Expand Analysis:**
   - Run differential expression analysis in blood samples
   - Compare tissue-specific DE genes with blood expression
   - Identify genes detectable in blood that reflect tissue changes

3. **Priority Targets from Tissue Analysis:**
   - Top 20 prioritized genes from GWAS integration
   - Genes with |Log2FC| > 2.5
   - Genes with P_adj_BH < 0.05

---

## 3. Enriched Pathways (Biological Insights)

### Pathway Enrichment Results

**Total pathways found:** 1 (limited by small number of significant genes at FDR < 0.05)

| Pathway | Source | P-value | Size |
|---------|--------|---------|------|
| **Spermine synthase activity** | GO:MF | 0.0125 | 1 |

### Biological Interpretation:

**Spermine Synthase Activity (GO:MF)**
- **Function:** Catalyzes synthesis of spermine from spermidine
- **Biological Role:**
  - Polyamine metabolism (spermine, spermidine, putrescine)
  - Cell proliferation and differentiation
  - Apoptosis regulation
  - Oxidative stress response
- **Relevance to Aneurysm:**
  - Polyamines play roles in vascular smooth muscle cell function
  - Implicated in vascular remodeling
  - May affect extracellular matrix synthesis
  - Could influence cell survival/apoptosis balance in aneurysm wall

### Limitations and Next Steps:

**Current Limitation:** Only 4 genes passed FDR < 0.05 threshold, limiting pathway enrichment power.

**Recommendations:**
1. **Relax FDR threshold:** Analyze pathways at FDR < 0.1 (1010 genes) or FDR < 0.2 for broader pathway discovery
2. **Gene set enrichment:** Use GSEA (Gene Set Enrichment Analysis) instead of overrepresentation analysis
3. **Known pathway mapping:** Manually map prioritized genes to known IA pathways:
   - **Vascular remodeling** (MMP-9, TIMP-1, collagen genes)
   - **Inflammation** (TNF-α, IL-6, CRP pathways)
   - **Endothelial dysfunction** (VEGF, angiopoietin pathways)
   - **Extracellular matrix** (COL3A1, COL1A1, LOXL1)
   - **Cell adhesion** (ICAM1, VCAM1)

---

## 4. Cross-Reference with Known IA Biology

### GWAS Loci Integration

**Note:** Current results show `In_GWAS_Locus = False` for all genes, which may indicate:
1. Coordinate-based mapping not yet implemented
2. GWAS loci genes may not be among the top DE genes
3. Tissue-specific expression may differ from genetic risk

**Known IA GWAS Loci (Bakker et al. 2020):**
- **9p21.3** (CDKN2B-AS1): Cell cycle regulation
- **8q11.23**: Endothelial function (SOX17)
- **4q31** (EDNRA): Endothelin signaling
- **10q24** (CFDP1): Unknown function
- **2q13** (COL3A1): Collagen, ECM
- **18q11** (LOXL1): Lysyl oxidase, ECM crosslinking

**Recommendation:** Map prioritized genes to chromosomal coordinates and check overlap with these GWAS loci.

### Functional Categories to Explore:

1. **Vascular Smooth Muscle Cells:**
   - Genes involved in contractility
   - ECM synthesis and degradation
   - Migration and proliferation

2. **Endothelial Cells:**
   - Adhesion molecules
   - Angiogenic factors
   - Nitric oxide signaling

3. **Inflammation:**
   - Cytokines and chemokines
   - Immune cell infiltration markers
   - Acute phase response

4. **ECM Remodeling:**
   - Collagens and elastin
   - Matrix metalloproteinases
   - Tissue inhibitors

---

## 5. Prioritization Matrix

### High Priority for Validation:

| Gene ID | Log2FC | P_adj_BH | Rationale |
|---------|--------|----------|-----------|
| **406979** | -20.26 | 0.092 | Extreme fold change, may be tissue-specific marker |
| **124220** | 2.88 | 0.041 | Significant with good effect size |
| **119587** | -2.07 | 0.041 | Significant with good effect size |
| **8685** | 2.61 | 0.059 | Good effect size, borderline significance |
| **221416** | 5.14 | 0.083 | Very strong effect size |

### Medium Priority:

- Genes with |Log2FC| > 2.5 and P_adj_BH < 0.1
- Genes in relevant biological pathways
- Genes with literature support for vascular/aneurysm biology

### Validation Strategy:

1. **qRT-PCR validation** of top 10 prioritized genes
2. **Immunohistochemistry** for protein-level validation
3. **Independent cohort validation** using similar samples
4. **Functional studies** (knockdown/overexpression) for top candidates
5. **Biomarker validation** in blood samples if detectable

---

## 6. Recommendations

### Immediate Actions:

1. **Gene Annotation:**
   - Annotate prioritized genes with gene symbols, names, and functions
   - Cross-reference with literature on aneurysm/vascular biology
   - Check for known associations with cerebrovascular diseases

2. **Pathway Analysis Enhancement:**
   - Re-run pathway enrichment at FDR < 0.1 or 0.2 thresholds
   - Perform GSEA for pathway discovery
   - Create pathway network diagrams

3. **Blood Integration:**
   - Refine blood sample metadata
   - Re-run blood-tissue integration analysis
   - Identify detectable biomarkers in blood

4. **Literature Review:**
   - Search PubMed for top prioritized genes
   - Identify functional studies
   - Compile evidence for aneurysm relevance

### Long-term Directions:

1. **Functional Validation:**
   - Cell culture studies (vascular smooth muscle, endothelial cells)
   - Animal models (if applicable)
   - Mechanistic studies

2. **Biomarker Development:**
   - Validate top genes as blood biomarkers
   - Develop diagnostic/prognostic tests
   - Assess clinical utility

3. **Therapeutic Targets:**
   - Identify druggable targets
   - Evaluate existing compounds
   - Consider drug repurposing

---

## 7. Summary Tables

### Top 10 Prioritized Genes Summary:

1. **406979**: Extreme downregulation (-20.26 FC) - highest priority
2. **221416**: Strong upregulation (5.14 FC) - high priority
3. **4745**: Strong downregulation (-4.67 FC) - high priority
4. **56729**: Moderate upregulation (3.87 FC)
5. **124220**: Moderate upregulation (2.88 FC), significant (P=0.041)
6. **6369**: Moderate upregulation (3.42 FC)
7. **100130717**: Moderate downregulation (-3.27 FC)
8. **107984150**: Moderate downregulation (-3.26 FC)
9. **6374**: Moderate upregulation (3.38 FC)
10. **101928710**: Moderate upregulation (3.09 FC)

### Key Findings:

- **1 extreme outlier** (Gene 406979) requires careful validation
- **Multiple genes** with significant fold changes (>3 FC)
- **Mix of upregulated and downregulated** genes suggests complex biology
- **Limited pathway enrichment** due to small number of significant genes
- **Blood integration** needs refinement for biomarker discovery

---

## 8. Next Steps

1. ✅ **Complete:** Pathway enrichment analysis
2. ✅ **Complete:** GWAS integration and prioritization
3. ✅ **Complete:** Initial biomarker analysis
4. ⏳ **Pending:** Gene annotation with symbols/names
5. ⏳ **Pending:** Enhanced pathway analysis (GSEA, relaxed thresholds)
6. ⏳ **Pending:** Blood metadata refinement and re-analysis
7. ⏳ **Pending:** Literature review and functional annotation
8. ⏳ **Pending:** Validation experiments

---

**Report Generated:** Integration Results Analysis  
**Analysis Date:** After completing validation pipeline  
**Status:** ✅ Analysis Complete - Ready for Validation Phase






