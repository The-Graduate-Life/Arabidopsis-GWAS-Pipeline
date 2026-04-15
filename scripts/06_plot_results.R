#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT:      06_plot_results.R
# DESCRIPTION: Generate a Manhattan plot and a QQ plot from GAPIT3 GWAS
#              association results using the qqman R package. Saves both
#              plots as PNG files in results/figures/.
#
# USAGE:
#   Rscript scripts/06_plot_results.R <gwas_file>
#
# ARGUMENTS:
#   gwas_file   Path to GAPIT3 results CSV containing columns:
#               Chr, Pos, SNP, P.value
#
# OUTPUTS:
#   results/figures/Manhattan.png   Manhattan plot (2000 x 1000 px)
#   results/figures/QQ.png          QQ plot (1000 x 1000 px)
#
# AUTHOR:      Fritzner Pierre
# DATE:        Spring 2026
# =============================================================================

suppressPackageStartupMessages(library(qqman))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 06_plot_results.R <gwas_file>")
}

gwas_file <- args[1]

# -----------------------------------------------------------------------------
# Load results
# -----------------------------------------------------------------------------
message("Loading GWAS results: ", gwas_file)
res <- read.csv(gwas_file, stringsAsFactors = FALSE)

# Drop rows with missing p-values
res <- res[!is.na(res$P.value), ]

# qqman expects columns: SNP, CHR (integer), BP, P
res$CHR <- as.integer(res$Chr)
res$BP  <- as.integer(res$Pos)
res$P   <- res$P.value

# Keep only standard chromosomes (1-5 for Arabidopsis)
res <- res[res$CHR %in% 1:5, ]

outdir <- "results/figures"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Derive trait name from filename for per-trait output files
# e.g. GAPIT.Association.GWAS_Results.FarmCPU.FT10(NYC).csv → FT10_NYC
trait <- basename(gwas_file)
trait <- sub("GAPIT\\.Association\\.GWAS_Results\\.FarmCPU\\.", "", trait)
trait <- sub("\\.csv$", "", trait)
trait <- gsub("[()]", "_", trait)
trait <- gsub("\\s+", "", trait)

# -----------------------------------------------------------------------------
# Manhattan plot
# -----------------------------------------------------------------------------
manhattan_out <- file.path(outdir, paste0("Manhattan_", trait, ".png"))
message("Saving Manhattan plot -> ", manhattan_out)

png(manhattan_out, width = 2000, height = 1000, res = 150)
manhattan(
    res,
    chr     = "CHR",
    bp      = "BP",
    snp     = "SNP",
    p       = "P",
    main    = "Manhattan Plot — Arabidopsis Flowering Time (FarmCPU)",
    col     = c("#4393C3", "#2166AC"),
    suggestiveline = -log10(1e-5),
    genomewideline = -log10(5e-8)
)
dev.off()

# -----------------------------------------------------------------------------
# QQ plot
# -----------------------------------------------------------------------------
qq_out <- file.path(outdir, paste0("QQ_", trait, ".png"))
message("Saving QQ plot -> ", qq_out)

png(qq_out, width = 1000, height = 1000, res = 150)
qq(res$P, main = "QQ Plot — Arabidopsis Flowering Time (FarmCPU)")
dev.off()

message("Done.")
