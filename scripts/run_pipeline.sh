#!/usr/bin/env bash
# =============================================================================
# SCRIPT:      run_pipeline.sh
# DESCRIPTION: Master script to run the full Arabidopsis thaliana GWAS
#              pipeline end-to-end. Executes all steps in order, from VCF
#              filtering through to Manhattan and QQ plot generation.
#
#              Assumes data has already been subsetted using 00_subset_data.sh.
#              All scripts are expected to live in the scripts/ directory and
#              all data in the data/ directory relative to the project root.
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
#   bcftools >= 1.15, tabix, plink >= 1.9, R >= 4.0
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
cd "$(dirname "$0")/.."

# -----------------------------------------------------------------------------
# Paths — edit these if your directory structure differs
# -----------------------------------------------------------------------------
readonly SCRIPTS_DIR="scripts"
readonly VCF_IN="data/raw/1001genomes_snp-short-indel_only_ACGTN.vcf.gz"
readonly VCF_FILTERED="data/subset/filtered.vcf.gz"
readonly PLINK_RAW="data/plink/raw"
readonly PLINK_QC="data/plink/qc"
readonly PHENO="data/subset/phenotype_subset.csv"
readonly OUTDIR_PCA="results/pca"
readonly OUTDIR_GWAS="results/gwas"

# -----------------------------------------------------------------------------
# Validation — fail early if key inputs are missing
# -----------------------------------------------------------------------------
[[ -f "$VCF_IN" ]] || { echo "ERROR: VCF not found: $VCF_IN"; echo "Place the VCF in data/raw/ as instructed in the README."; exit 1; }
[[ -f "$PHENO" ]]  || { echo "ERROR: Phenotype not found: $PHENO"; echo "Place the phenotype CSV in data/raw/ as instructed in the README."; exit 1; }

# -----------------------------------------------------------------------------
# Environment setup
# -----------------------------------------------------------------------------
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
