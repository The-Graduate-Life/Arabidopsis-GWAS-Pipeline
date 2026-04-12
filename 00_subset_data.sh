#!/usr/bin/env bash
# =============================================================================
# SCRIPT:      00_subset_data.sh
# DESCRIPTION: Subset the 1001 Genomes VCF and AraPheno phenotype file to a
#              manageable number of accessions for classroom use. Performs a
#              geographically stratified random sample to preserve population
#              structure across Eurasia.
#
#              This script assumes the raw data files have already been
#              downloaded and placed in data/raw/ as instructed in the README:
#
#              Genotype (VCF):
#                wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/ \
#                     intersection_snp_short_indel_vcf/ \
#                     intersection.snp.short_indel.prioritized_pubSNP.192accessions.vcf.gz \
#                     -O data/raw/genotype.vcf.gz
#
#              Phenotype (AraPheno study 12 — FT16 + FT10):
#                wget https://arapheno.1001genomes.org/rest/study/12/values.csv \
#                     -O data/raw/phenotype_raw.csv
#
# USAGE:
#   bash scripts/00_subset_data.sh [n_accessions]
#
# ARGUMENTS:
#   n_accessions  Number of accessions to sample (default: 150)
#
# INPUTS:
#   data/raw/genotype.vcf.gz     Full 1001 Genomes VCF
#   data/raw/phenotype_raw.csv   AraPheno phenotype CSV (accession_id column)
#
# OUTPUTS:
#   data/subset/subset.vcf.gz         Subsetted and indexed VCF
#   data/subset/phenotype_subset.csv  Matching phenotype file (Taxa column)
#   data/subset/sample_ids.txt        List of selected accession IDs
#
# DEPENDENCIES:
#   bcftools >= 1.15, tabix, conda (pandas, numpy)
#
# EXAMPLE:
#   cd ~/Arabidopsis-GWAS-Pipeline
#   bash scripts/00_subset_data.sh 100
#
# AUTHOR:      Fritzner Pierre
# DATE:        Spring 2026
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Always run from the project root regardless of where the script is called from
# -----------------------------------------------------------------------------
cd "$(dirname "$0")/.."

# -----------------------------------------------------------------------------
# Usage
# -----------------------------------------------------------------------------
usage() {
    sed -n '/^# USAGE:/,/^# AUTHOR:/p' "$0" | sed 's/^# \?//'
    exit 1
}

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
readonly N_ACCESSIONS=${1:-150}

readonly VCF="data/raw/genotype.vcf.gz"
readonly PHENO="data/raw/phenotype_raw.csv"
readonly SUBSET_DIR="data/subset"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
[[ -f "$VCF" ]]   || { echo "ERROR: VCF not found: $VCF"; echo "See README for download instructions."; exit 1; }
[[ -f "$PHENO" ]] || { echo "ERROR: Phenotype not found: $PHENO"; echo "See README for download instructions."; exit 1; }

mkdir -p "$SUBSET_DIR"

# -----------------------------------------------------------------------------
# Environment setup
# Creates the conda environment with all system-level dependencies.
# R packages that are not available on conda (SNPRelate, GAPIT3, qqman)
# are installed via Bioconductor, GitHub, and CRAN respectively.
# The environment is only created once — conda skips creation if it
# already exists.
# -----------------------------------------------------------------------------
conda create -n gwas_env -c conda-forge -c bioconda \
    pandas numpy plink \
    r-base r-devtools r-biocmanager \
    r-lme4 r-nloptr r-fs \
    -y

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gwas_env

echo "Installing R packages..."
Rscript -e "BiocManager::install(c('gdsfmt', 'SNPRelate'), ask=FALSE, update=TRUE, force=TRUE)"
Rscript -e "devtools::install_github('jiabowang/GAPIT3', force=TRUE)"
Rscript -e "install.packages('qqman', repos='https://cloud.r-project.org')"
echo "✓ R packages installed"

# -----------------------------------------------------------------------------
# Step 1: Select accessions
# Extracts VCF sample IDs and uses Python to select a geographically
# stratified random sample of N_ACCESSIONS from the phenotype file,
# keeping only accessions present in both the VCF and phenotype file.
# -----------------------------------------------------------------------------
echo "============================================="
echo " Step 1: Select $N_ACCESSIONS accessions"
echo "============================================="

bcftools query -l "$VCF" > /tmp/vcf_ids.txt
echo "VCF contains $(wc -l < /tmp/vcf_ids.txt) samples"

python3 - <<PYEOF
import pandas as pd
import numpy as np

pheno   = pd.read_csv("$PHENO")
vcf_ids = set(open("/tmp/vcf_ids.txt").read().split())

# Filter to accessions present in the VCF
pheno["accession_id"] = pheno["accession_id"].astype(str)
pheno = pheno[pheno["accession_id"].isin(vcf_ids)]
print(f"Accessions overlapping with VCF: {len(pheno)}")

# Rename for GAPIT3 compatibility and drop unused columns
pheno.rename(columns={"accession_id": "Taxa"}, inplace=True)
if "replicate_id" in pheno.columns:
    pheno.drop(columns=["replicate_id"], inplace=True)

# Deduplicate: one row per accession
pheno = pheno.groupby("Taxa", as_index=False).first()
print(f"Unique accessions after dedup: {len(pheno)}")

# Geographically stratified sample if coordinates are available,
# otherwise fall back to simple random sample
if "latitude" in pheno.columns and "longitude" in pheno.columns:
    pheno = pheno.dropna(subset=["latitude", "longitude"])
    pheno["lat_bin"]  = pd.cut(pheno["latitude"],  bins=4, labels=False)
    pheno["lon_bin"]  = pd.cut(pheno["longitude"], bins=4, labels=False)
    pheno["geo_cell"] = pheno["lat_bin"].astype(str) + "_" + pheno["lon_bin"].astype(str)

    n_cells  = pheno["geo_cell"].nunique()
    per_cell = max(1, $N_ACCESSIONS // n_cells)
    sampled  = (pheno.groupby("geo_cell", group_keys=False)
                     .apply(lambda g: g.sample(min(len(g), per_cell), random_state=42)))

    # Top up to exactly N_ACCESSIONS if stratification fell short
    if len(sampled) < $N_ACCESSIONS:
        remaining = pheno[~pheno["Taxa"].isin(sampled["Taxa"])]
        extra     = remaining.sample(
                        min($N_ACCESSIONS - len(sampled), len(remaining)),
                        random_state=42)
        sampled = pd.concat([sampled, extra])

    sampled = sampled.head($N_ACCESSIONS)
    print(f"Stratified sample across {n_cells} geographic cells")
else:
    print("No lat/lon columns found — falling back to random sample")
    sampled = pheno.sample($N_ACCESSIONS, random_state=42)

# Write outputs
sampled["Taxa"].astype(str).to_csv("$SUBSET_DIR/sample_ids.txt", index=False, header=False)

drop_cols = [c for c in ["lat_bin", "lon_bin", "geo_cell"] if c in sampled.columns]
sampled.drop(columns=drop_cols, inplace=True)
sampled.to_csv("$SUBSET_DIR/phenotype_subset.csv", index=False)

print(f"Selected {len(sampled)} accessions  →  $SUBSET_DIR/sample_ids.txt")
print(f"Phenotype subset          →  $SUBSET_DIR/phenotype_subset.csv")
PYEOF

# -----------------------------------------------------------------------------
# Step 2: Index VCF with tabix
# -----------------------------------------------------------------------------
echo ""
echo "============================================="
echo " Step 2: Index VCF with tabix"
echo "============================================="

if [[ ! -f "${VCF}.tbi" ]]; then
    echo "Indexing $VCF ..."
    tabix -p vcf "$VCF"
else
    echo "Index already exists, skipping."
fi

# -----------------------------------------------------------------------------
# Step 3: Subset VCF with bcftools
# -----------------------------------------------------------------------------
echo ""
echo "============================================="
echo " Step 3: Subset VCF with bcftools"
echo "============================================="

bcftools view \
    --samples-file "$SUBSET_DIR/sample_ids.txt" \
    --force-samples \
    -m2 -M2 -v snps \
    --min-ac 1:minor \
    -Oz -o "$SUBSET_DIR/subset.vcf.gz" \
    "$VCF"

bcftools index "$SUBSET_DIR/subset.vcf.gz"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
N_SAMPLES=$(bcftools query -l "$SUBSET_DIR/subset.vcf.gz" | wc -l)
N_VARS=$(bcftools stats "$SUBSET_DIR/subset.vcf.gz" | grep "^SN" | grep "number of SNPs" | cut -f4)

echo ""
echo "============================================="
echo " Summary"
echo "============================================="
echo "  Samples in subsetted VCF : $N_SAMPLES"
echo "  SNPs retained            : $N_VARS"
echo "  VCF output               : $SUBSET_DIR/subset.vcf.gz"
echo "  Phenotype output         : $SUBSET_DIR/phenotype_subset.csv"
echo ""
echo "✓ Done. Run the pipeline with:"
echo "  bash scripts/run_pipeline.sh"
