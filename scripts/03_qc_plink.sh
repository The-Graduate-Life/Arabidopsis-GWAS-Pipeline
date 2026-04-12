#!/usr/bin/env bash
# =============================================================================
# SCRIPT:      03_qc_plink.sh
# DESCRIPTION: Apply standard GWAS quality control filters to PLINK binary
#              files. Filters on minor allele frequency (MAF >= 5%), per-SNP
#              missingness (< 10%), and per-sample missingness (< 10%).
#
#              NOTE: HWE filtering is intentionally disabled (--hwe 0).
#              Arabidopsis thaliana is a highly selfing species, meaning that
#              virtually all loci deviate from Hardy-Weinberg equilibrium by
#              design. Applying an HWE filter would incorrectly remove the
#              majority of real, informative SNPs.
#
# USAGE:
#   bash scripts/03_qc_plink.sh <bfile> <out>
#
# ARGUMENTS:
#   bfile   Input PLINK prefix (e.g. data/plink/raw)
#   out     Output PLINK prefix for QC-filtered files (e.g. data/plink/qc)
#
# OUTPUTS:
#   <out>.bed   QC-filtered binary genotype matrix
#   <out>.bim   QC-filtered variant information file
#   <out>.fam   QC-filtered sample information file
#
# DEPENDENCIES:
#   plink >= 1.9
#
# EXAMPLE:
#   bash scripts/03_qc_plink.sh \
#       data/plink/raw \
#       data/plink/qc
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
readonly OUT=$2

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
[[ -f "${BFILE}.bed" ]] || { echo "ERROR: PLINK .bed not found: ${BFILE}.bed"; exit 1; }

mkdir -p "$(dirname "$OUT")"
