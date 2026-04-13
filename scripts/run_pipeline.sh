#!/usr/bin/env bash
# =============================================================================
# SCRIPT:      run_pipeline.sh
# DESCRIPTION: Master script to run the full Arabidopsis thaliana GWAS
#              pipeline end-to-end. Executes all steps in order, from VCF
#              filtering through to Manhattan and QQ plot generation.
#
#              Before running this script, the following setup scripts must
#              have been run in order:
#                1. bash scripts/00_subset_data.sh 100
#                   Sets up the conda environment, downloads phenotype data,
#                   and selects accessions.
#                2. bash scripts/00b_subset_vcf.sh   (local)
#                   qsub scripts/00b_subset_vcf.sh   (HPC)
#                   Subsets the VCF to selected accessions. Run as a PBS job
#                   on HPC due to memory requirements (~32 GB).
#
#              This script can be run locally or submitted as a PBS job on HPC:
#                bash scripts/run_pipeline.sh         # local
#                qsub scripts/run_pipeline.sh         # HPC (PBS)
#
# PBS DIRECTIVES (ignored when run with bash locally):
#PBS -N gwas_pipeline
#PBS -q class
#PBS -W group_list=classq
#PBS -j oe
#PBS -m abe
#PBS -l walltime=12:00:00
#PBS -l select=ncpus=8:mpiprocs=8:mem=32gb
#
# USAGE:
#   bash run_pipeline.sh
#
# ARGUMENTS:
#   None. All paths are defined as variables below — edit them to match
#   your local setup before running.
#
# OUTPUTS:
#   data/subset/filtered.vcf.gz        Filtered VCF
#   data/plink/raw.*                   PLINK binary files
#   data/plink/qc.*                    QC-filtered PLINK files
#   results/pca/PCA.csv                PCA scores
#   results/pca/Kinship.csv            Kinship matrix
#   results/gwas/GAPIT.Association.GWAS_Results.FarmCPU.*.csv  GWAS results (one per trait)
#   results/figures/Manhattan_<trait>.png  Manhattan plot per trait
#   results/figures/QQ_<trait>.png         QQ plot per trait
#
# DEPENDENCIES:
#   bcftools >= 1.19, tabix, plink >= 1.9, R >= 4.0
#   R packages: SNPRelate, GAPIT3, qqman
#
# EXAMPLE:
#   cd ~/Arabidopsis-GWAS-Pipeline
#   bash run_pipeline.sh
#
# AUTHOR:      Fritzner Pierre
# DATE:        Spring 2026
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Always run from the project root regardless of where the script is called from
# -----------------------------------------------------------------------------
# Working directory
# On HPC: PBS sets $PBS_O_WORKDIR to where qsub was called from
# Locally: resolve from the script's location
# -----------------------------------------------------------------------------
cd "${PBS_O_WORKDIR:-$(cd "$(dirname "$0")/.." && pwd)}"

# -----------------------------------------------------------------------------
# Paths — edit these if your directory structure differs
# -----------------------------------------------------------------------------
readonly SCRIPTS_DIR="scripts"
readonly VCF_IN="data/subset/subset.vcf.gz"
readonly VCF_FILTERED="data/subset/filtered.vcf.gz"
readonly PLINK_RAW="data/plink/raw"
readonly PLINK_QC="data/plink/qc"
readonly PHENO="data/subset/phenotype_subset.csv"
readonly OUTDIR_PCA="results/pca"
readonly OUTDIR_GWAS="results/gwas"

# -----------------------------------------------------------------------------
# Validation — fail early if key inputs are missing
# -----------------------------------------------------------------------------
[[ -f "$VCF_IN" ]] || {
    echo "ERROR: Subset VCF not found: $VCF_IN"
    echo "Run scripts/00_subset_data.sh first, then scripts/00b_subset_vcf.sh"
    exit 1
}
[[ -f "$PHENO" ]] || {
    echo "ERROR: Phenotype not found: $PHENO"
    echo "Run scripts/00_subset_data.sh first."
    exit 1
}

# -----------------------------------------------------------------------------
# Environment setup
# -----------------------------------------------------------------------------
module load anaconda 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gwas_env
export PATH="$(conda info --base)/envs/gwas_env/bin:$PATH"

# -----------------------------------------------------------------------------
# Pipeline
# -----------------------------------------------------------------------------
echo "====================================="
echo "   Arabidopsis GWAS Pipeline"
echo "====================================="
echo ""

echo "[1/6] Filtering VCF..."
bash "$SCRIPTS_DIR/01_filter_vcf.sh" "$VCF_IN" "$VCF_FILTERED"

echo "[2/6] Converting VCF to PLINK..."
bash "$SCRIPTS_DIR/02_vcf_to_plink.sh" "$VCF_FILTERED" "$PLINK_RAW"

echo "[3/6] Running QC..."
bash "$SCRIPTS_DIR/03_qc_plink.sh" "$PLINK_RAW" "$PLINK_QC"

echo "[4/6] Computing PCA and kinship..."
bash "$SCRIPTS_DIR/04_pca_kinship.sh" "$PLINK_QC" "$OUTDIR_PCA"

echo "[5/6] Running GWAS..."
bash "$SCRIPTS_DIR/05_gwas.sh" "$PHENO" "$PLINK_QC" "$OUTDIR_GWAS"

echo "[6/6] Plotting results..."
# GAPIT3 produces one results CSV per trait — loop over all of them
shopt -s nullglob
gwas_files=("$OUTDIR_GWAS"/GAPIT.Association.GWAS_Results.FarmCPU.*.csv)
if [[ ${#gwas_files[@]} -eq 0 ]]; then
    echo "ERROR: No GAPIT3 results CSV found in $OUTDIR_GWAS"
    exit 1
fi
for gwas_file in "${gwas_files[@]}"; do
    echo "  Plotting: $(basename "$gwas_file")"
    bash "$SCRIPTS_DIR/06_plot_results.sh" "$gwas_file"
done

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "====================================="
echo "   ✓ Pipeline complete"
echo "====================================="
echo ""
echo "  Filtered VCF    : $VCF_FILTERED"
echo "  PLINK QC files  : ${PLINK_QC}.bed / .bim / .fam"
echo "  PCA & kinship   : $OUTDIR_PCA"
echo "  GWAS results    : $OUTDIR_GWAS/GAPIT.Association.GWAS_Results.FarmCPU.*.csv"
echo "  Figures         : results/figures/Manhattan_<trait>.png & QQ_<trait>.png"
