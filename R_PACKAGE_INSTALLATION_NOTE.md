# R Package Installation Note

## Status
R is not currently installed on this system. The pipeline will skip R-dependent preprocessing steps but will proceed with RNA-seq analysis.

## To Install R Packages (when R is available)

1. Install R from: https://cran.r-project.org/

2. Run the installation script:
   ```R
   source("install_r_packages.R")
   ```

   Or manually in R:
   ```R
   if (!requireNamespace("BiocManager", quietly = TRUE)) {
       install.packages("BiocManager")
   }
   BiocManager::install(c("affy", "Biobase"), update = FALSE, ask = FALSE)
   ```

## Required Packages
- **affy**: For Affymetrix microarray normalization
- **Biobase**: Base classes for Bioconductor

## Impact
- GSE26969 microarray preprocessing requires R/affy
- GSE122897 RNA-seq analysis can proceed without R
- The pipeline will automatically skip R-dependent steps if R is not available






