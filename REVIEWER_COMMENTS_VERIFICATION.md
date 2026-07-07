# Reviewer Comments – Verification Checklist

## Figures (5 main, max 4 panels each)

| Figure | Content | Panels | Status |
|--------|---------|--------|--------|
| Fig 1 | Study design / dataset coverage | 1 | ✓ Compliant |
| Fig 2 | Volcano plot (Welch DE) | 1 | ✓ Compliant |
| Fig 3 | Feature distributions (DE, EC, ECM) | 3 | ✓ Compliant |
| Fig 4 | 4a Feature importance, 4b ROC | 2 | ✓ Compliant |
| Fig 5 | Risk stratification | 1 | ✓ Compliant |

**Figure regeneration**: None required. Current figures reflect Welch t-test DE, signature scores, and model results. DESeq2, nested CV, and pathway enrichment are reported in text and supplementary tables; no new main figures needed.

---

## Supplementary Tables – Main Text References

| Table | Content | Referenced in |
|-------|---------|---------------|
| S1 | Risk scores, stratification | Table 1 inference, Methods (risk stratification) |
| S2 | CV performance summary | Methods (Statistical model) |
| S3 | Top DESeq2 genes | Results (Differential expression), Methods (DE) |
| S4 | Welch vs DESeq2 comparison | Results (Differential expression) |
| S5 | Pathway enrichment | Results (Biological signature scores) |
| S6 | Sample metadata (Sample, Group, Sex) | Table 1 inference, Methods (Data acquisition) |

---

## Reviewer Comments Addressed

| Comment | Addressed | Location |
|---------|-----------|----------|
| **DE method** – Justify Welch t-test or use DESeq2/edgeR | ✓ | Methods: Welch justification; DESeq2 run on raw counts (S3); both reported |
| **Circularity** – DE genes from same data as model | ✓ | Methods: nested DE-within-CV; Results: nested AUC 0.64 ± 0.17; Limitations |
| **Confounding** – Pre-rupture vs peri-rupture; age/sex | ✓ | Limitations: confounding; Data acquisition + S6: age/sex not available |
| **Pathway enrichment** – GO/KEGG/Reactome, FDR | ✓ | S5, pathway_enrichment_gse122897_de.py, g:Profiler |
| **Bulk RNA-seq limitations** | ✓ | Discussion, Limitations |
| **Single-cell/spatial** – Future directions | ✓ | Discussion, Future directions |
| **Modeling tone** – Cautious; AUC n=43, 8 false negatives | ✓ | Results, Discussion, Figure 4b inference |
| **Introduction** – Biology vs prediction gap | ✓ | Introduction |
| **Metadata** – Age, sex, covariates | ✓ | Data acquisition, Limitations, Table 1, S6 |

---

## Metadata

- **Available**: Sample (GEO ID), Group (Ruptured/Unruptured), Sex (1/43 reported)
- **Not available**: Age (not in GEO records)
- **Limitation**: Cannot adjust for age or sex in DE or modelling
