**Author:** Fritzner Pierre  
**Course:** Scripting for Biologists (BIOL-7180)  
**Semester:** Spring 2026    
**Institution/Affiliation**: Auburn University

---

## Project: **GWAS Analysis Pipeline for *Arabidopsis thaliana* Flowering Time Using Public Data**

### Overview:
This tutorial walks you through a complete genome-wide association study (GWAS) pipeline using publicly available data from *Arabidopsis thaliana*. By the end, you will have run every step from raw VCF filtering to Manhattan and QQ plots using a series of modular shell scripts.

---

## 1. Background

### What is GWAS?

A **genome-wide association study (GWAS)** tests hundreds of thousands to millions of single-nucleotide polymorphisms (SNPs) across the genome for statistical association with a trait of interest. For each SNP, we ask: do individuals carrying one allele tend to have a systematically different phenotype than individuals carrying the other allele?

### Why *Arabidopsis thaliana* flowering time?

*A. thaliana* is a model plant with a well-characterized genome and extensive publicly available data. **Flowering time** (the number of days from germination to first flower) is a classic quantitative trait with well-known genetic architecture, making it ideal for teaching GWAS concepts. We use two phenotypes:

|Column|Meaning|
|-|-|
|`FT16`|Flowering time measured at 16 °C (long day)|
|`FT10`|Flowering time measured at 10 °C (short day, vernalization)|

### Why FarmCPU?

*A. thaliana* has strong geographic population structure. If ignored, variants that differ in frequency between populations — but are not causally linked to the trait — will produce false-positive associations. **FarmCPU** (Fixed and Random Model Circulating Probability Unification) is a mixed-model method that simultaneously controls for population structure (via PCA covariates) and kinship (via a random effect), while avoiding model overfitting.

> *Note on Hardy-Wenberg Equilibrium (HWE) filtering: *A. thaliana* reproduces almost exclusively by self-fertilization, meaning the vast majority of loci deviate strongly from Hardy–Weinberg equilibrium by design. The QC script disables HWE filtering for this reason.*

---

## 2. Pipeline Overview

```
Raw VCF + Phenotype CSV
        │
        ▼
   [00] Subset phenotype data          ← stratified random sample of ~100 accessions
        │
        ▼
   [00b] Subset VCF data               ← stratified random sample of ~100 accessions
        │
        ▼
   [01] Filter VCF                     ← biallelic SNPs only, < 10% missingness
        │
        ▼
   [02] VCF → PLINK                    ← .bed / .bim / .fam binary format
        │
        ▼
   [03] QC (PLINK)                     ← MAF ≥ 5%, SNP/sample missingness < 10%
        │
        ▼
   [04] PCA & Kinship (R)              ← SNPRelate: population structure + relatedness
        │
        ▼
   [05] GWAS (R / GAPIT3)               ← FarmCPU model, p-values per SNP
        │
        ▼
   [06] Plot results (R)                ← Manhattan plot + QQ plot
```

> *Each numbered step is a standalone shell script. A master script (`run_pipeline.sh`) chains steps 01-06 together. Steps 00 must run first, following by 00b before runing the master script.*

---

## 3. Repository Structure

```
Arabidopsis-GWAS-Pipeline/
├── data/
│   ├── raw/
│   │   ├── 1001genomes_snp-short-indel_only_ACGTN.vcf.gz           # Full 1001G VCF
│   │   ├── 1001genomes_snp-short-indel_only_ACGTN.vcf.gz.tbi
│   │   └── FT_Field_phenotype.csv                                   # AraPheno phenotype data
│   ├── subset/                                                       # Created by Step 00b
│   │   ├── subset.vcf.gz
│   │   ├── subset.vcf.gz.csi
│   │   ├── sample_ids.txt
│   │   └── phenotype_subset.csv                                        # Created by Step 00b
│   └── plink/                                                          # Created by Steps 2–3
│       ├── raw.bed / raw.bim / raw.fam
│       └── qc.bed  / qc.bim  / qc.fam
├── results/                                                             # Created by Steps 4–6
│   ├── pca/
│   │   ├── PCA.csv
│   │   └── Kinship.csv
│   ├── gwas/
│   │   └── GAPIT.FarmCPU.csv
│   └── figures/
│       ├── Manhattan.png
│       └── QQ.png
├── scripts/
│   ├── 00_subset_data.sh
|   ├── 00b_subset_vcf.sh
│   ├── 01_filter_vcf.sh
│   ├── 02_vcf_to_plink.sh
│   ├── 03_qc_plink.sh
│   ├── 04_pca_kinship.sh
│   ├── 04_pca_kinship.R
│   ├── 05_gwas.sh
│   ├── 05_gwas.R
│   ├── 06_plot_results.sh
│   ├── 06_plot_results.R
│   └── run_pipeline.sh
└── README.md
└── TUTORIAL.md
```
---
## 4. Prerequisites

### System tools

|Tool|Version|Install|
|-|-|-|
|`bcftools`|≥ 1.15|conda install -c bioconda bcftools|
|`tabix`|≥ 1.15|included with bcftools|
|`plink`|≥ 1.9|conda install -c bioconda plink|
|`R`|≥ 4.0|conda install -c conda-forge r-base|
|`conda`|any|[Miniconda installer](https://docs.conda.io/en/latest/miniconda.html)|

### R packages

Install these once inside an R session:

```r
# SNPRelate (Bioconductor)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SNPRelate")

# GAPIT3
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3", force = TRUE)

# qqman
install.packages("qqman")
```

### Python packages (for Step 0 only)

```bash
conda create -n gwas_env -c conda-forge pandas numpy -y
```

> The subsetting script activates `gwas_env` automatically. All other steps use only shell tools and R.

### Verify installation

```bash
bcftools --version | head -1   # bcftools 1.x
tabix --version | head -1      # tabix (htslib) 1.x
plink --version                # PLINK v1.9x
Rscript --version              # R scripting front-end version 4.x.x
```
---

## 5. Data

Before running the pipeline, the raw input files must be downloaded into `~/Arabidopsis-GWAS-Pipeline/data/raw/`. Run the following commands from the **project root**:

```bash
# Download the genotype VCF (~19 GB — this will take a while)
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz \
    -P data/raw/

# Download the phenotype data
wget "https://arapheno.1001genomes.org/phenotype/6/values.csv" \
    -O data/raw/FT_Field_phenotype.csv
```

After downloading, `~/Arabidopsis-GWAS-Pipeline/data/raw/` should contain:

```
Arabidopsis-GWAS-Pipeline/data/raw/
├── 1001genomes_snp-short-indel_only_ACGTN.vcf.gz
└── FT_Field_phenotype.csv
```

### Genotype data — 1001 Genomes Project

The full VCF (`1001genomes_snp-short-indel_only_ACGTN.vcf.gz`) contains SNPs and short indels for **1,135 natural accessions** of *A. thaliana* across five chromosomes.

* Source: [https://1001genomes.org/data/GMI-MPI/releases/v3.1/](https://1001genomes.org/data/GMI-MPI/releases/v3.1/)
* Format: bgzipped VCF with tabix index (`.tbi`)
* Chromosome naming: `1`, `2`, `3`, `4`, `5` (no "chr" prefix)

### Phenotype data — AraPheno

The phenotype file (`FT_Field_phenotype.csv`) contains field-measured flowering time for hundreds of accessions.

```
accession_id,replicate_id,FT16,FT10
10000,6520,54.5,61.0
100000,6521,48.25,71.67
...
```

* Source: [https://arapheno.1001genomes.org/study/12/](https://arapheno.1001genomes.org/study/12/)
* `FT16`: days to flower at 16 °C
* `FT10`: days to flower at 10 °C (vernalization treatment)

### The pre-built subset

For classroom use, a pre-built subset of **100 accessions** is already in `data/subset/`. You can skip Step 0 and start directly at Step 1 if you want to use this subset.

---
