#!/usr/bin/env bash
# =============================================================================
# SCRIPT:      02_vcf_to_plink.sh
# DESCRIPTION: Convert a filtered VCF file into PLINK binary format (.bed,
#              .bim, .fam). Uses --allow-extra-chr to handle Arabidopsis
#              chromosome naming (1–5 rather than chr1–chr5) and
#              --set-missing-var-ids to assign IDs to any unlabelled variants.
#
# USAGE:
#   bash scripts/02_vcf_to_plink.sh <vcf> <out>
#
# ARGUMENTS:
#   vcf   Path to filtered VCF (.vcf.gz)
#   out   Output prefix for PLINK files (e.g. data/plink/raw)
#
# OUTPUTS:
#   <out>.bed   Binary genotype matrix
#   <out>.bim   Variant information file
#   <out>.fam   Sample information file
#
# DEPENDENCIES:
#   plink >= 1.9
#
# EXAMPLE:
#   bash scripts/02_vcf_to_plink.sh \
#       data/subset/filtered.vcf.gz \
#       data/plink/raw
#
# AUTHOR:      Fritzner Pierre
# DATE:        Spring 2025
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
readonly VCF=$1
readonly OUT=$2
