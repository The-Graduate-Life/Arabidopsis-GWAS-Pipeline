**Author:** Fritzner Pierre  
**Course:** Scripting for Biologists (BIOL-7180)  
**Semester:** Spring 2026    
**Institution/Affiliation**: Auburn University

---

## Table of Contents

- [1. Background](#1-background)
- [2. Pipeline Overview](#2-pipeline-overview)
- [3. Repository Structure](#3-repository-structure)
- [4. Prerequisites](#4-prerequisites)
- [5. Data](#5-data)
- [6. How to Cite](#6-how-to-cite)
- [7. AI Disclosure](#7-ai-disclosure)
- [8. References](#8-references)

---

## Project: **GWAS Analysis Pipeline for *Arabidopsis thaliana* Flowering Time Using Public Data**

## 1\. Background

### What is GWAS?

A **genome-wide association study (GWAS)** is a genetic technique to test hundreds of thousands to millions of single-nucleotide polymorphisms (SNPs) across the genome for statistical association with a trait of interest. For each SNP, the GWAS inform if the individuals carrying one allele tend to have a systematically different phenotype than individuals carrying the other allele.

### Why *Arabidopsis thaliana* flowering time?

*Arabidopsis thaliana* is a model plant most worldwidly studied and with a well-characterized genome and extensive publicly available data. The **Flowering time** (the number of days from germination to first flower) in *A. thaliana* is a classic quantitative trait with well-known genetic architecture, making it ideal for teaching GWAS concepts. We use two phenotypes:

|Traits|Meaning|
|-|-|
|`FT16`|Flowering time measured at 16 °C (long day)|
|`FT10`|Flowering time measured at 10 °C (short day)|

### Why FarmCPU?

**FarmCPU** (Fixed and Random Model Circulating Probability Unification) is a statistical mixed-model method that simultaneously controls for population structure (via PCA covariates) and kinship (via a random effect), while avoiding model overfitting. Since *A. thaliana* has strong geographic population structure, if it is ignored, variants that differ in frequency between populations (but are not causally linked to the trait) will produce false-positive associations. Based on this, **FarmCPU** is one of thest mixed-models to handle these false-positive associations.

> *Note on Hardy-Wenberg Equilibrium (HWE) filtering: *A. thaliana* reproduces almost exclusively by self-fertilization, meaning the vast majority of loci deviate strongly from HWE by design. To handle that, the quality control (QC) helps to disable HWE filtering for this reason in the script.*


## 2\. Pipeline Overview

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

> *Each numbered step is a standalone shell script. A master script (`run_pipeline.sh`) chains steps 01-06 together. Script `00_subset_data.sh` must run first, following by `00b_subset_vcf.sh` before runing the master script.*


## 3\. Repository Structure

Here is a physical structure well-organized for a clean analysis in this project:
```
Arabidopsis-GWAS-Pipeline/
├── data/
│   ├── raw/
│   │   ├── 1001genomes_snp-short-indel_only_ACGTN.vcf.gz   # Full 1001G VCF
│   │   ├── 1001genomes_snp-short-indel_only_ACGTN.vcf.gz.tbi
│   │   └── FT_Field_phenotype.csv                           # AraPheno phenotype data
│   ├── subset/                                               # Created by steps 00 & 00b
│   │   ├── subset.vcf.gz
│   │   ├── subset.vcf.gz.csi
│   │   ├── sample_ids.txt
│   │   └── phenotype_subset.csv
│   └── plink/                                               # Created by steps 2–3
│       ├── raw.bed / raw.bim / raw.fam
│       └── qc.bed  / qc.bim  / qc.fam
├── results/                                                  # Created by steps 4–6
│   ├── pca/
│   │   ├── PCA.csv
│   │   └── Kinship.csv
│   ├── gwas/
│   │   └── GAPIT.Association.GWAS_Results.FarmCPU.*.csv
│   └── figures/
│       ├── Manhattan_<trait>.png
│       └── QQ_<trait>.png
├── scripts/
│   ├── 00_subset_data.sh
│   ├── 00b_subset_vcf.sh
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
├── README.md
├── FOLLOWME.md
└── PROJECT_DETAILS.md
```


## 4\. Prerequisites

For a smooth computational journey during this project, the following tools and packages must be installed prior running `01_filter_vcf.sh`. Luky you 😊, the script  `00_subset_data.sh` will smoothly install all the required tools, packages, and dependencies.

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
BiocManager::install(c("gdsfmt", "SNPRelate"))

# GAPIT3
install.packages("pak")
pak::pak("jiabowang/GAPIT3")

# qqman
install.packages("qqman")
```

### Python packages (for Step 0 only)

```bash
conda create -n gwas_env -c conda-forge -c bioconda \
    pandas numpy plink \
    r-base r-devtools r-biocmanager \
    r-lme4 r-nloptr r-fs r-matrix r-rcpp \
    zlib libcurl openssl \
    compilers make -y
```

> The subsetting script activates `gwas_env` automatically. All other steps use only shell tools and R.

### Verify installation

To make sure these tools are installed and check the active version, run the following lines of code in the command line:

```bash
bcftools --version | head -1   # bcftools 1.x
tabix --version | head -1      # tabix (htslib) 1.x
plink --version                # PLINK v1.9x
Rscript --version              # R scripting front-end version 4.x.x
```


## 5\. Data

Before running the pipeline, the raw input files must be downloaded into `~/Arabidopsis-GWAS-Pipeline/data/raw/`. Run the following commands from the **project root directory**:

```bash
# Download the genotype VCF (~19 GB — this will take a while)
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz \
    -P data/raw/

# Download the phenotype data
wget "https://arapheno.1001genomes.org/phenotype/6/values.csv" \
    -O data/raw/FT_Field_phenotype.csv
```

After downloading, the project root directory (`~/Arabidopsis-GWAS-Pipeline/data/raw/`) should contain at least these raw files:

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

The phenotype file (`FT_Field_phenotype.csv`) contains field-measured flowering time for hundreds of accessions and has the following structure:

```
accession_id,replicate_id,FT16,FT10
10000,6520,54.5,61.0
100000,6521,48.25,71.67
...
```

* Source: [https://arapheno.1001genomes.org/study/12/](https://arapheno.1001genomes.org/study/12/)
* `FT16`: days to flower at 16 °C
* `FT10`: days to flower at 10 °C (vernalization treatment)



## 6\. How to Cite

If you use this pipeline in your work, please cite:

> Pierre, F. (2026). *GWAS Analysis Pipeline for Arabidopsis thaliana Flowering Time Using Public Data*. GitHub. https://github.com/The-Graduate-Life/Arabidopsis-GWAS-Pipeline

Please also cite the underlying tools (see [References](#8-references)).



## 7\. AI Disclosure

This pipeline was completed with the assistance of **[Claude](https://claude.ai)** (Anthropic).

AI assistance was used for:
- Debugging shell scripts and R code
- Refining docstring-style script headers
- Troubleshooting conda environment and R package installation
- Adapting the pipeline for HPC submission (PBS)

All biological decisions, pipeline design, data selection, and interpretation of results were made by the author. AI-generated code and documentation were reviewed, tested, and modified before inclusion in the final pipeline.

---

## 8\. References

* The 1001 Genomes Consortium. 1001 Genomes Project. [https://1001genomes.org/](https://1001genomes.org/)
* AraPheno: A public database of Arabidopsis thaliana phenotypes. [https://arapheno.1001genomes.org/](https://arapheno.1001genomes.org/)
* Wang, J., & Zhang, Z. (2021). GAPIT version 3: boosting power and accuracy for genomic association and prediction. *Genomics, proteomics & bioinformatics, 19(4)*, 629-640. [https://zzlab.net/GAPIT/](https://zzlab.net/GAPIT/)
* Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A., Bender, D., ... & Sham, P. C. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. *The American journal of human genetics, 81(3)*, 559-575. [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/)
* Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., ... & Li, H. (2021). Twelve years of SAMtools and BCFtools. *Gigascience, 10(2)*, giab008. [http://www.htslib.org/](http://www.htslib.org/)
* Zheng, X., Levine, D., Shen, J., Gogarten, S. M., Laurie, C., & Weir, B. S. (2012). A high-performance computing toolset for relatedness and principal component analysis of SNP data. *Bioinformatics, 28(24)*, 3326-3328. SNPRelate: [https://bioconductor.org/packages/SNPRelate/](https://bioconductor.org/packages/SNPRelate/)
* Turner, S. D. (2014). qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots. *bioRxiv*. [https://cran.r-project.org/web/packages/qqman/](https://cran.r-project.org/web/packages/qqman/)
