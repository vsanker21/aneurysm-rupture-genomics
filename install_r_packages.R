# R script to install required packages for microarray normalization
# Run this in R or RStudio before preprocessing GSE26969

# Check if BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install required Bioconductor packages
BiocManager::install(c("affy", "Biobase"), update = FALSE, ask = FALSE)

# Verify installation
if (requireNamespace("affy", quietly = TRUE) && requireNamespace("Biobase", quietly = TRUE)) {
    cat("SUCCESS: Required packages installed\n")
    cat("  - affy: Available\n")
    cat("  - Biobase: Available\n")
} else {
    cat("ERROR: Package installation failed\n")
}






