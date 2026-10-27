Stouffer meta-analysis (v1.1.1)
================================

GSE122897 Entrez IDs are mapped to HGNC approved symbols using
metadata/hgnc_entrez_to_symbol.tsv (HGNC complete set, approved entries
with an Entrez ID). The three-cohort inner join is then on gene symbol,
so BH-FDR is computed across the full overlap with GSE54083 and GSE13353
rather than the 16-gene subset produced by the sparse local map in v1.1.0.

On the current HGNC complete-set snapshot this overlap is 12,521 genes
(Jin’s 12,562 figure is the same join with a slightly different map).
The six-gene wall core is derived from that overlap: directionally
concordant in all three cohorts and Meta_FDR_BH < 1e-6
(TGFBI, ANPEP, MARCO, MPP1, CLEC5A, CTSL).

Run:
  python rerun_stouffer_meta_21v21.py
