#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT:      05_gwas.R
# DESCRIPTION: Run GWAS using the FarmCPU model in GAPIT3. PLINK binary files
#              are read using SNPRelate (which we already have from step 04),
#              converted to the numeric GD matrix and GM map that GAPIT3
#              expects. PCA covariates computed in step 04 are passed via CV
#              to control for population structure.
#
#              GAPIT3 interface used:
#                Y   — phenotype data.frame (Taxa + trait columns)
#                GD  — numeric genotype matrix (samples x SNPs), coded 0/1/2
#                GM  — SNP map data.frame (SNP, Chr, Pos)
#                CV  — covariate data.frame (Taxa + PC columns)
#
# ARGUMENTS:
#   args[1]  Path to phenotype CSV (must contain a 'Taxa' column)
#   args[2]  PLINK prefix (e.g. data/plink/qc — without .bed/.bim/.fam)
#   args[3]  Output directory (e.g. results/gwas)
#
# OUTPUTS:
#   <outdir>/GAPIT.<model>.<trait>.GWAS.Results.csv
#   <outdir>/GAPIT.<model>.<trait>.Manhattan.Plot.pdf
#   <outdir>/GAPIT.<model>.<trait>.QQ.Plot.pdf
#
# DEPENDENCIES:
#   GAPIT3, SNPRelate, gdsfmt
# =============================================================================

suppressPackageStartupMessages({
  library(GAPIT)
  library(SNPRelate)
})

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
pheno_file  <- args[1]
geno_prefix <- args[2]
outdir      <- args[3]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Resolve all paths to absolute BEFORE setwd so relative paths stay valid.
# For geno_prefix we resolve via the .bed file (which exists) then strip
# the extension — normalizePath(mustWork=FALSE) does not resolve correctly
# for prefixes that don't exist as files.
pheno_file  <- normalizePath(pheno_file, mustWork = TRUE)
bed_file    <- normalizePath(paste0(geno_prefix, ".bed"), mustWork = TRUE)
geno_prefix <- sub("\\.bed$", "", bed_file)
outdir      <- normalizePath(outdir, mustWork = TRUE)

# GAPIT3 writes output files to the working directory
setwd(outdir)

# -----------------------------------------------------------------------------
# Load phenotype
# -----------------------------------------------------------------------------
cat("Loading phenotype:", pheno_file, "\n")
pheno <- read.csv(pheno_file, header = TRUE)
pheno$Taxa <- as.character(pheno$Taxa)

# Keep Taxa + numeric trait columns only
trait_cols <- names(pheno)[sapply(pheno, is.numeric)]
pheno      <- pheno[, c("Taxa", trait_cols), drop = FALSE]
cat("  →", nrow(pheno), "accessions,", length(trait_cols),
    "trait(s):", paste(trait_cols, collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# Load PCA covariates (from step 04)
# -----------------------------------------------------------------------------
pca_file <- file.path(dirname(outdir), "pca", "PCA.csv")
cat("Loading PCA:", pca_file, "\n")
CV        <- read.csv(pca_file, header = TRUE)
CV$Taxa   <- as.character(CV$Taxa)

# -----------------------------------------------------------------------------
# Read PLINK → GDS → GD + GM for GAPIT3
# SNPRelate is used to parse the binary .bed format correctly.
# GD: numeric matrix, rows = samples, cols = SNPs, values in {0, 1, 2}
# GM: data.frame with columns SNP, Chr, Pos
# -----------------------------------------------------------------------------
cat("Reading PLINK files:", geno_prefix, "\n")

bed <- paste0(geno_prefix, ".bed")
bim <- paste0(geno_prefix, ".bim")
fam <- paste0(geno_prefix, ".fam")
gds <- file.path(outdir, "gwas_tmp.gds")

# Convert PLINK → GDS (temporary file in outdir)
snpgdsBED2GDS(bed, fam, bim, gds, verbose = FALSE)
geno_gds <- snpgdsOpen(gds)

# Extract genotype matrix as {0,1,2} dosage
geno_mat <- snpgdsGetGeno(geno_gds, snpfirstdim = FALSE, verbose = FALSE)
sample_ids <- read.gdsn(index.gdsn(geno_gds, "sample.id"))
snp_ids    <- read.gdsn(index.gdsn(geno_gds, "snp.id"))
snp_chr    <- read.gdsn(index.gdsn(geno_gds, "snp.chromosome"))
snp_pos    <- read.gdsn(index.gdsn(geno_gds, "snp.position"))

snpgdsClose(geno_gds)

# Build GD (add Taxa column as first column, required by GAPIT3)
GD        <- as.data.frame(geno_mat)
GD        <- cbind(data.frame(Taxa = sample_ids), GD)
names(GD) <- c("Taxa", snp_ids)

# Build GM (SNP map)
GM <- data.frame(
  SNP  = snp_ids,
  Chr  = snp_chr,
  Pos  = snp_pos,
  stringsAsFactors = FALSE
)

cat("  →", nrow(GD), "samples,", nrow(GM), "SNPs loaded\n")

# -----------------------------------------------------------------------------
# Align phenotype, covariates, and genotype to the same samples
# -----------------------------------------------------------------------------
common_taxa <- Reduce(intersect, list(pheno$Taxa, CV$Taxa, GD$Taxa))
cat("  →", length(common_taxa), "samples in common across pheno/CV/geno\n")

pheno <- pheno[pheno$Taxa %in% common_taxa, ]
CV    <- CV[CV$Taxa %in% common_taxa, ]
GD    <- GD[GD$Taxa %in% common_taxa, ]

# Ensure chromosome is numeric — Arabidopsis uses integers 1-5 which
# SNPRelate may read as character, causing GAPIT.Genotype.View to fail
GM$Chr <- as.numeric(GM$Chr)

# -----------------------------------------------------------------------------
# Run GAPIT3
# Genotype.View is disabled (plot.type excludes "GD") because the
# chromosome axis plot fails with non-standard chromosome naming.
# GAPIT3 will still produce Manhattan and QQ plots via file.output = TRUE.
# -----------------------------------------------------------------------------
cat("Running GAPIT3 (FarmCPU)...\n")
GAPIT(
  Y            = pheno,
  GD           = GD,
  GM           = GM,
  CV           = CV,
  PCA.total    = 0,      # PCA already supplied via CV
  model        = "FarmCPU",
  file.output  = TRUE,
  Geno.View.output = FALSE  # Disable genotype view plot
)

# Clean up temporary GDS file
if (file.exists(gds)) file.remove(gds)

cat("✓ GWAS complete → results written to", outdir, "\n")