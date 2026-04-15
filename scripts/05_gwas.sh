#!/usr/bin/env bash
# =============================================================================
# SCRIPT:      05_gwas.sh
# DESCRIPTION: Run a genome-wide association study (GWAS) using the FarmCPU
#              model implemented in the GAPIT3 R package. FarmCPU (Fixed and
#              random model Circulating Probability Unification) controls for
#              population structure and kinship while avoiding model
#              overfit, making it well suited to the Arabidopsis diversity
#              panel.
#
#              The phenotype file must have a 'Taxa' column with accession IDs
#              matching the PLINK .fam file. This is produced automatically by
#              00_subset_data.sh as phenotype_subset.csv.
#
# USAGE:
#   bash scripts/05_gwas.sh <phenotype> <geno> <outdir>
#
# ARGUMENTS:
#   phenotype   Path to the phenotype CSV (must contain a 'Taxa' column)
#   geno        PLINK prefix for QC-filtered genotype files
#   outdir      Directory for GWAS result files
#
# OUTPUTS:
#   <outdir>/GAPIT.FarmCPU.csv        Full association results
#   <outdir>/GAPIT.FarmCPU.*.pdf      QQ and Manhattan plots (from GAPIT3)
#
# DEPENDENCIES:
#   R >= 4.0, GAPIT3
#
# EXAMPLE:
#   bash scripts/05_gwas.sh \
#       data/subset/phenotype_subset.csv \
#       data/plink/qc \
#       results/gwas
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

[[ $# -lt 3 ]] && usage

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
readonly PHENO=$1
readonly GENO=$2
readonly OUTDIR=$3

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
[[ -f "$PHENO" ]]         || { echo "ERROR: Phenotype file not found: $PHENO";   exit 1; }
[[ -f "${GENO}.bed" ]]    || { echo "ERROR: PLINK .bed not found: ${GENO}.bed";  exit 1; }

mkdir -p "$OUTDIR"

# -----------------------------------------------------------------------------
# Environment
# -----------------------------------------------------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gwas_env

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
echo "Running GWAS (FarmCPU via GAPIT3)..."
echo "  Phenotype : $PHENO"
echo "  Genotype  : $GENO"
echo "  Output    : $OUTDIR"

Rscript "$(dirname "$0")/05_gwas.R" "$PHENO" "$GENO" "$OUTDIR"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "✓ GWAS complete → $OUTDIR"
echo ""
echo "Next step:"
echo "  bash scripts/06_plot_results.sh $OUTDIR/GAPIT.FarmCPU.csv"
