# Literature Review Guide for Prioritized Genes
## Aneurysm and Vascular Associations

This guide provides search strategies and key publications for reviewing prioritized genes in the context of intracranial aneurysm and vascular biology.

---

## Search Strategy

### PubMed Search Terms:
For each prioritized gene, use these search patterns:

1. **Direct gene search:**
   ```
   "[Gene_Symbol] OR [Gene_Name] OR [Aliases]"
   ```

2. **Aneurysm-specific:**
   ```
   "([Gene_Symbol] OR [Gene_Name]) AND (aneurysm OR intracranial aneurysm OR cerebral aneurysm OR IA)"
   ```

3. **Vascular biology:**
   ```
   "([Gene_Symbol] OR [Gene_Name]) AND (vascular OR vessel OR artery OR endothelium OR smooth muscle)"
   ```

4. **Stroke/cerebrovascular:**
   ```
   "([Gene_Symbol] OR [Gene_Name]) AND (stroke OR cerebrovascular OR subarachnoid hemorrhage OR SAH)"
   ```

---

## Top 20 Prioritized Genes for Literature Review

### 1. Gene 406979 (Priority: 20.98, Log2FC: -20.26)
**Status:** EXTREME PRIORITY - Highest fold change

**Search Terms:**
- "406979 OR (Entrez:406979)"
- "406979 AND aneurysm"
- "406979 AND vascular"

**Key Questions:**
- What is this gene's function?
- Any known associations with vascular diseases?
- Expression in cerebrovascular tissue?

---

### 2. Gene 221416 (Priority: 5.55, Log2FC: 5.14)
**Status:** HIGH PRIORITY - Strong upregulation

**Search Terms:**
- "221416 OR (Entrez:221416)"
- "221416 AND aneurysm"

**Key Questions:**
- Function in vascular remodeling?
- Inflammatory response gene?
- Role in tissue repair?

---

### 3. Gene 4745 (Priority: 5.05, Log2FC: -4.67)
**Status:** HIGH PRIORITY - Strong downregulation

**Search Terms:**
- "4745 OR (Entrez:4745)"
- "4745 AND aneurysm"

**Key Questions:**
- Structural or functional gene?
- Protective factor?
- Why downregulated in rupture?

---

### 4. Gene 56729 (Priority: 4.44, Log2FC: 3.87)
**Status:** HIGH PRIORITY - Moderate-strong upregulation

**Search Terms:**
- "56729 OR (Entrez:56729)"
- "56729 AND vascular"

---

### 5. Gene 124220 (Priority: 4.00, Log2FC: 2.88, P_adj: 0.041)
**Status:** HIGH PRIORITY - Significant with good effect size

**Search Terms:**
- "124220 OR (Entrez:124220)"
- "124220 AND aneurysm"

**Key Questions:**
- Most significant gene among top priorities
- Validated biomarker candidate
- Functional studies available?

---

## Key Databases and Resources

### 1. PubMed (https://pubmed.ncbi.nlm.nih.gov/)
- Primary literature search
- Use MeSH terms: "Intracranial Aneurysm", "Cerebral Arteries", "Vascular Remodeling"
- Filter by: Review articles, Human studies, Recent publications (last 5 years)

### 2. GeneCards (https://www.genecards.org/)
- Gene function and expression
- Pathway information
- Disease associations

### 3. DisGeNET (https://www.disgenet.org/)
- Disease-gene associations
- Variant associations
- Evidence scores

### 4. STRING (https://string-db.org/)
- Protein-protein interactions
- Functional associations
- Pathway networks

### 5. GWAS Catalog (https://www.ebi.ac.uk/gwas/)
- Genetic associations
- Variant information
- Related GWAS studies

---

## Key Publications to Review

### Intracranial Aneurysm GWAS:
1. **Bakker et al. 2020** - "Genome-wide association study of intracranial aneurysms identifies 17 risk loci"
   - Nature Genetics 2020
   - Key loci: 9p21.3, 8q11.23, 4q31, 10q24, 2q13, 18q11

2. **Foroud et al. 2014** - "Genome-wide association study of intracranial aneurysm"
   - Neurology 2014

3. **Yasuno et al. 2010** - "Genome-wide association study of intracranial aneurysm"
   - Nature Genetics 2010

### Blood Biomarkers:
1. **Poppenberg et al. 2023** - Blood transcriptomics for aneurysm detection
   - GSE159610 dataset

### Pathway Reviews:
1. Review articles on vascular remodeling
2. Review articles on aneurysm pathophysiology
3. Review articles on endothelial dysfunction

---

## Functional Categories to Review

### 1. Vascular Remodeling
**Key genes to check:**
- MMPs (matrix metalloproteinases)
- TIMPs (tissue inhibitors)
- Collagens (COL1A1, COL3A1)
- LOX/LOXL1 (lysyl oxidases)

**Search:** "vascular remodeling AND aneurysm"
**Keywords:** extracellular matrix, elastin, collagen, degradation

### 2. Inflammation
**Key genes to check:**
- Cytokines (TNF, IL-6, IL-8)
- Chemokines
- Acute phase proteins

**Search:** "inflammation AND aneurysm"
**Keywords:** immune response, cytokines, chemokines

### 3. Endothelial Dysfunction
**Key genes to check:**
- Adhesion molecules (ICAM1, VCAM1)
- Angiogenic factors (VEGFA, ANGPT2)
- Nitric oxide pathway

**Search:** "endothelial dysfunction AND aneurysm"
**Keywords:** endothelium, angiogenesis, nitric oxide

### 4. Smooth Muscle Cell Function
**Key genes to check:**
- Contractile proteins
- ECM synthesis
- Migration factors

**Search:** "smooth muscle AND aneurysm"
**Keywords:** VSMC, contractility, migration

---

## Literature Review Checklist

For each prioritized gene, document:

- [ ] **Gene Symbol and Name**
- [ ] **Function** (from GeneCards, UniProt, literature)
- [ ] **Expression** (tissue-specific expression patterns)
- [ ] **Pathways** (KEGG, Reactome, GO terms)
- [ ] **Disease Associations** (DisGeNET, literature)
- [ ] **Aneurysm-Specific Studies** (PubMed search)
- [ ] **Vascular Biology Studies** (PubMed search)
- [ ] **Genetic Variants** (GWAS Catalog, ClinVar)
- [ ] **Functional Studies** (knockout, overexpression, siRNA)
- [ ] **Biomarker Potential** (detectable in blood, serum, CSF)

---

## Priority Ranking for Literature Review

1. **Immediate Priority:**
   - Gene 406979 (extreme fold change)
   - Gene 124220 (most significant, P_adj < 0.05)
   - Gene 119587 (significant, P_adj < 0.05)

2. **High Priority:**
   - Top 10 prioritized genes (Priority Score > 3.5)
   - Genes with |Log2FC| > 3.0

3. **Medium Priority:**
   - Top 20-50 prioritized genes
   - Genes with |Log2FC| > 2.5

---

## Next Steps

1. **Create gene annotation table** with symbols, names, functions
2. **Perform PubMed searches** for top 20 genes
3. **Compile literature findings** in summary table
4. **Identify functional studies** for validation
5. **Assess biomarker potential** for blood detection

---

## Example Search Query Template

```python
# Template for automated literature search
gene_symbols = ["406979", "221416", "4745", "56729", "124220"]

for gene in gene_symbols:
    queries = [
        f'"{gene}" AND (aneurysm OR intracranial aneurysm)',
        f'"{gene}" AND (vascular OR vessel OR artery)',
        f'"{gene}" AND (stroke OR cerebrovascular)',
        f'"{gene}" AND (biomarker OR serum OR blood)'
    ]
    
    for query in queries:
        # Search PubMed API or web interface
        # Document results: PMID, title, relevance
        pass
```

---

**Last Updated:** After completing prioritized gene identification  
**Status:** Ready for literature review implementation






