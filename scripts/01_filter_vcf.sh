#!/usr/bin/env bash
# =============================================================================
# SCRIPT:      01_filter_vcf.sh
# DESCRIPTION: Filter a VCF file to retain only biallelic SNPs with less than
#              10% missing genotypes across samples. This is the first quality
#              control step after subsetting and prepares the VCF for
#              conversion to PLINK format.
#
# USAGE:
#   bash scripts/01_filter_vcf.sh <vcf> <out>
#
# ARGUMENTS:
#   vcf   Path to the input VCF (.vcf.gz), typically data/subset/subset.vcf.gz
#   out   Path for the filtered output VCF (.vcf.gz)
#
# OUTPUTS:
#   <out>       Filtered, bgzipped VCF
#   <out>.tbi   Tabix index
#
# DEPENDENCIES:
#   bcftools >= 1.19
#
# EXAMPLE:
#   bash scripts/01_filter_vcf.sh \
#       data/subset/subset.vcf.gz \
#       data/subset/filtered.vcf.gz
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
readonly VCF=$1
readonly OUT=$2

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
[[ -f "$VCF" ]] || { echo "ERROR: VCF not found: $VCF"; exit 1; }

mkdir -p "$(dirname "$OUT")"

# -----------------------------------------------------------------------------
# Filter
# -----------------------------------------------------------------------------
echo "Filtering VCF..."
echo "  Input  : $VCF"
echo "  Output : $OUT"
echo "  Filters: biallelic SNPs only, missingness < 10%"

bcftools view -m2 -M2 -v snps "$VCF" \
    | bcftools filter -i 'F_MISSING<0.1' \
    -Oz -o "$OUT"

bcftools index "$OUT"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
N_VARS=$(bcftools stats "$OUT" | grep "^SN" | grep "number of SNPs" | cut -f4)
echo ""
echo "✓ VCF filtered → $N_VARS SNPs retained: $OUT"
echo ""
echo "Next step:"
echo "  bash scripts/02_vcf_to_plink.sh $OUT data/plink/raw"
