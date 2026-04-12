#!/usr/bin/env bash
# =============================================================================
# SCRIPT:      04_pca_kinship.sh
# DESCRIPTION: Compute principal components (PCA) and a kinship matrix from
#              QC-filtered PLINK files using the SNPRelate R package. PCA is
#              used to visualise population structure and as covariates in the
#              GWAS model. The kinship matrix captures pairwise genetic
#              relatedness and is used by GAPIT3 to control for population
#              stratification.
#
# USAGE:
#   bash scripts/04_pca_kinship.sh <bfile> <outdir>
#
# ARGUMENTS:
#   bfile    PLINK prefix for QC-filtered files (e.g. data/plink/qc)
#   outdir   Directory for PCA and kinship outputs (e.g. results/pca)
#
# OUTPUTS:
#   <outdir>/PCA.csv       Principal component scores (PC1, PC2 per sample)
#   <outdir>/Kinship.csv   IBS kinship matrix (samples x samples)
#   <outdir>/arabidopsis.gds  GDS format genotype file (intermediate)
#
# DEPENDENCIES:
#   R >= 4.0, SNPRelate (Bioconductor)
#
# EXAMPLE:
#   bash scripts/04_pca_kinship.sh \
#       data/plink/qc \
#       results/pca
#
# AUTHOR:      Fritzner Pierre
# DATE:        Spring 2026
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Usage
# -----------------------------------------------------------------------------
usage() {
    sed -n '/^# USAGE:/,/^# AUTHOR:/p' "$0" | sed 's/^# \?//'
    exit 1
}

[[ $# -lt 2 ]] && usage

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
readonly BFILE=$1
readonly OUTDIR=$2

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
[[ -f "${BFILE}.bed" ]] || { echo "ERROR: PLINK .bed not found: ${BFILE}.bed"; exit 1; }

mkdir -p "$OUTDIR"

# -----------------------------------------------------------------------------
# Environment
# -----------------------------------------------------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gwas_env

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
echo "Running PCA and kinship estimation..."
echo "  Input  : $BFILE"
echo "  Output : $OUTDIR"

Rscript "$(dirname "$0")/04_pca_kinship.R" "$BFILE" "$OUTDIR"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "✓ PCA & kinship complete → $OUTDIR"
echo ""
echo "Next step:"
echo "  bash scripts/05_gwas.sh data/subset/phenotype_subset.csv $BFILE results/gwas"
