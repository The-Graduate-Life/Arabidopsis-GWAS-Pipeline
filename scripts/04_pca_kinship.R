#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT:      04_pca_kinship.R
# DESCRIPTION: Compute PCA and IBS kinship matrix from QC-filtered PLINK files
#              using the SNPRelate package. Called by 04_pca_kinship.sh.
#
# ARGUMENTS:
#   args[1]  PLINK prefix (e.g. data/plink/qc)
#   args[2]  Output directory (e.g. results/pca)
#
# OUTPUTS:
#   <outdir>/PCA.csv          Principal component scores (Taxa, PC1, PC2)
#   <outdir>/Kinship.csv      IBS kinship matrix (samples x samples)
#   <outdir>/arabidopsis.gds  Intermediate GDS genotype file
#
# DEPENDENCIES:
#   SNPRelate (Bioconductor), gdsfmt
# =============================================================================

suppressPackageStartupMessages({
  library(SNPRelate)
})

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
args   <- commandArgs(trailingOnly = TRUE)
bfile  <- args[1]
outdir <- args[2]

bed <- paste0(bfile, ".bed")
bim <- paste0(bfile, ".bim")
fam <- paste0(bfile, ".fam")
gds <- file.path(outdir, "arabidopsis.gds")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Convert PLINK → GDS
# -----------------------------------------------------------------------------
cat("Converting PLINK to GDS...\n")
snpgdsBED2GDS(bed, fam, bim, gds, verbose = TRUE)

geno <- snpgdsOpen(gds)

# -----------------------------------------------------------------------------
# PCA
# Version-safe call: snpgdsPCA dropped the nthread argument in some builds.
# We call without it to ensure compatibility across SNPRelate versions.
# -----------------------------------------------------------------------------
cat("Computing PCA...\n")
pca <- snpgdsPCA(geno, verbose = FALSE)

pcs <- data.frame(
  Taxa = pca$sample.id,
  PC1  = pca$eigenvect[, 1],
  PC2  = pca$eigenvect[, 2],
  PC3  = pca$eigenvect[, 3]
)
write.csv(pcs, file.path(outdir, "PCA.csv"), row.names = FALSE)
cat("  →", nrow(pcs), "samples written to PCA.csv\n")

# -----------------------------------------------------------------------------
# Kinship (IBS matrix)
# -----------------------------------------------------------------------------
cat("Computing kinship (IBS)...\n")
kin <- snpgdsIBS(geno, verbose = FALSE)
write.csv(kin$ibs, file.path(outdir, "Kinship.csv"), row.names = TRUE)
cat("  → Kinship matrix written to Kinship.csv\n")

# -----------------------------------------------------------------------------
# Cleanup
# -----------------------------------------------------------------------------
snpgdsClose(geno)
cat("✓ PCA & kinship complete\n")
