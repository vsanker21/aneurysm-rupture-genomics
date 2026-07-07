# Template for formal colocalization with R package `coloc` (install: install.packages("coloc"))
# Prerequisite: harmonize alleles between gwas_region_*.tsv and eqtl_region_*.tsv (same effect allele).
#
# Example workflow:
#   gwas <- read.table(".../gwas_region_MARCO_rs2422391.tsv", header=TRUE, sep="\t")
#   eqtl <- read.table(".../eqtl_region_MARCO_rs2422391.tsv", header=TRUE, sep="\t")
#   # Build dataset1 (IA case-control) and dataset2 (blood eQTL quant) per coloc::check_dataset
#   # res <- coloc.abf(dataset1, dataset2)
#
# See: https://chr1swallis.github.io/coloc/

message("This file is a pointer template. Load regional TSVs from results/phases_g/stronger_genetics/")
