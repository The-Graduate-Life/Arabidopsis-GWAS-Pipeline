#!/usr/bin/env bash
# =============================================================================
# SCRIPT:      00b_subset_vcf.sh
# DESCRIPTION: Subset the full 1001 Genomes VCF to the accessions selected
#              by 00_subset_data.sh. This step is separated from the main
#              setup script because it is memory intensive (~32 GB) and must
#              be run as a batch job on memory-constrained systems such as
#              HPC login nodes.
#
#              On a local machine with sufficient RAM (>= 16 GB), this script
#              can be run directly with bash. On an HPC, submit it as a job:
#
#                qsub 00b_subset_vcf.sh         # PBS/Torque (ASC)
#                sbatch 00b_subset_vcf.sh       # SLURM
#                bash 00b_subset_vcf.sh         # local machine
#
#              Prerequisites: 00_subset_data.sh must have been run first to
#              produce data/subset/sample_ids.txt.
#
# PBS DIRECTIVES (ignored when run with bash locally):
#PBS -N subset_vcf
#PBS -q class
#PBS -W group_list=classq
#PBS -j oe
#PBS -l walltime=02:00:00
#PBS -l select=ncpus=4:mpiprocs=4:mem=32gb
#
# USAGE:
#   bash scripts/00b_subset_vcf.sh               # local
#   qsub scripts/00b_subset_vcf.sh               # HPC (PBS)
#
# ARGUMENTS:
#   None. All paths are defined as variables below.
#
# INPUTS:
#   data/raw/1001genomes_snp-short-indel_only_ACGTN.vcf.gz  Full VCF
#   data/subset/sample_ids.txt                               Sample list
#
# OUTPUTS:
#   data/subset/subset.vcf.gz      Subsetted and indexed VCF
#   data/subset/subset.vcf.gz.csi  Tabix index
#
# DEPENDENCIES:
#   bcftools >= 1.15, conda (gwas_env)
#
# EXAMPLE:
#   # On ASC HPC:
#   qsub scripts/00b_subset_vcf.sh
#
#   # Locally:
#   bash scripts/00b_subset_vcf.sh
#
# AUTHOR:      Fritzner Pierre
# DATE:        Spring 2026
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Working directory
# On HPC: PBS sets $PBS_O_WORKDIR to where qsub was called from
# Locally: resolve from the script's location
# -----------------------------------------------------------------------------
cd "${PBS_O_WORKDIR:-$(cd "$(dirname "$0")/.." && pwd)}"

# -----------------------------------------------------------------------------
# Environment
# module load is silently skipped if not on an HPC
# -----------------------------------------------------------------------------
module load anaconda 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gwas_env

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
readonly VCF="data/raw/1001genomes_snp-short-indel_only_ACGTN.vcf.gz"
readonly SAMPLE_IDS="data/subset/sample_ids.txt"
readonly SUBSET_VCF="data/subset/subset.vcf.gz"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
[[ -f "$VCF" ]]        || { echo "ERROR: VCF not found: $VCF";                    exit 1; }
[[ -f "$SAMPLE_IDS" ]] || { echo "ERROR: Sample IDs not found: $SAMPLE_IDS";      
                             echo "Run scripts/00_subset_data.sh first.";          exit 1; }

N_SAMPLES=$(wc -l < "$SAMPLE_IDS")
echo "Subsetting VCF to $N_SAMPLES accessions..."
echo "  Input  : $VCF"
echo "  Samples: $SAMPLE_IDS"
echo "  Output : $SUBSET_VCF"

# -----------------------------------------------------------------------------
# Subset VCF with bcftools
# Keeps only biallelic SNPs present in at least one of the selected samples
# -----------------------------------------------------------------------------
bcftools view \
    --samples-file "$SAMPLE_IDS" \
    --force-samples \
    -m2 -M2 -v snps \
    --min-ac 1:minor \
    -Oz -o "$SUBSET_VCF" \
    "$VCF"

bcftools index "$SUBSET_VCF"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
N_VARS=$(bcftools stats "$SUBSET_VCF" | grep "^SN" | grep "number of SNPs" | cut -f4)

echo ""
echo "============================================="
echo " Summary"
echo "============================================="
echo "  Samples : $N_SAMPLES"
echo "  SNPs    : $N_VARS"
echo "  Output  : $SUBSET_VCF"
echo ""
echo "✓ Done. Run the pipeline with:"
echo "  bash scripts/run_pipeline.sh               # local"
echo "  qsub scripts/run_pipeline.sh               # HPC"
