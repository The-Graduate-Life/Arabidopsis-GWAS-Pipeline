#!/usr/bin/env bash
# =============================================================================
# SCRIPT:      06_plot_results.sh
# DESCRIPTION: Generate a Manhattan plot and a QQ (quantile-quantile) plot
#              from GWAS association results using the qqman R package.
#              The Manhattan plot displays -log10(p-value) per SNP across
#              chromosomes, highlighting genome-wide significant loci. The QQ
#              plot compares the observed p-value distribution against the null
#              expectation to assess inflation or deflation of test statistics.
#
# USAGE:
#   bash scripts/06_plot_results.sh <gwas_file>
#
# ARGUMENTS:
#   gwas_file   Path to GWAS results CSV from GAPIT3, containing columns:
#               Chr, Pos, SNP, P
#
# OUTPUTS:
#   results/figures/Manhattan.png   Manhattan plot (2000 x 1000 px)
#   results/figures/QQ.png          QQ plot (1000 x 1000 px)
#
# DEPENDENCIES:
#   R >= 4.0, qqman
#
# EXAMPLE:
#   bash scripts/06_plot_results.sh \
#       results/gwas/GAPIT.FarmCPU.csv
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

[[ $# -lt 1 ]] && usage

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
readonly GWAS_FILE=$1

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
[[ -f "$GWAS_FILE" ]] || { echo "ERROR: GWAS results file not found: $GWAS_FILE"; exit 1; }

mkdir -p results/figures

# -----------------------------------------------------------------------------
# Environment
# -----------------------------------------------------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gwas_env

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
echo "Generating plots..."
echo "  Input  : $GWAS_FILE"
echo "  Output : results/figures/Manhattan.png"
echo "           results/figures/QQ.png"

Rscript "$(dirname "$0")/06_plot_results.R" "$GWAS_FILE"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
TRAIT=$(basename "$GWAS_FILE" .csv | sed 's/GAPIT.Association.GWAS_Results.FarmCPU.//; s/[()]/_/g')
echo ""
echo "✓ Plots saved → results/figures/Manhattan_${TRAIT}.png & QQ_${TRAIT}.png"
