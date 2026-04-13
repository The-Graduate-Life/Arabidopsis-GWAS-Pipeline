### Overview:
This tutorial walks you through a complete genome-wide association study (GWAS) pipeline using publicly available data from *Arabidopsis thaliana*. By the end, you will have run every step from raw VCF filtering to Manhattan and QQ plots using a series of modular shell scripts.

---
## Step-by-Step Tutorial
### Getting Started

**1. Fork [this repository](https://github.com/The-Graduate-Life/Arabidopsis-GWAS-Pipeline)**
Click the **Fork** button in the upper right corner of this page. After a few seconds you will have your own copy of the repository in your GitHub account.

**2. Clone your fork**
On your fork's page, click the green **Code** button and copy the URL. Then in your terminal:

```bash
git clone <the-url-you-just-copied>
cd Arabidopsis-GWAS-Pipeline
```

**3. Verify**
You should now be in your local copy of the repository:

```bash
ls scripts/
```

You are ready to proceed with the setup. 

**4. Download the raw data
```bash
# Download the genotype VCF (~19 GB — this will take a while)
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz \
    -P data/raw/

# Download the phenotype data
wget "https://arapheno.1001genomes.org/phenotype/6/values.csv" \
    -O data/raw/FT_Field_phenotype.csv
```

>All commands below assume you are in the **project root**.

---

### Step 0 — Subset the Data

> *Skip this step if you are using the pre-built subset in `data/subset/`. Proceed to Step 1.

This step selects a manageable number of accessions from the full 1001 Genomes VCF using a **geographically stratified random sample**, which preserves population structure across Eurasia. It then extracts those accessions from the VCF and writes a matching phenotype file.

**Script:** `scripts/00_subset_data.sh`

**Usage:**

```bash
bash scripts/00_subset_data.sh <vcf> <phenotype> <outdir> [n_accessions]
```

**Arguments:**

|Argument|Description|Example|
|-|-|-|
|`vcf`|Path to the full 1001G VCF|`data/raw/1001genomes_snp-short-indel_only_ACGTN.vcf.gz`|
|`phenotype`|Path to the AraPheno CSV|`data/raw/FT_Field_phenotype.csv`|
|`outdir`|Output directory|`data/subset`|
|`n_accessions`|Number of accessions to sample (default: 150)|`100`|

**Example run:**

```bash
bash scripts/00_subset_data.sh \
    data/raw/1001genomes_snp-short-indel_only_ACGTN.vcf.gz \
    data/raw/FT_Field_phenotype.csv \
    data/subset 100
```

**What it does:**

1. Reads all sample IDs from the VCF header
2. Joins with the phenotype file to find overlapping accessions
3. Performs a geographically stratified random sample (bins by latitude/longitude)
4. Writes `sample_ids.txt` (one ID per line) and `phenotype_subset.csv`
5. Calls `bcftools view --samples-file` to extract the chosen accessions into a new VCF

**Outputs:**

|File|Description|
|-|-|
|`data/subset/subset.vcf.gz`|Subsetted, bgzipped VCF|
|`data/subset/subset.vcf.gz.csi`|CSI index|
|`data/subset/sample_ids.txt`|Accession IDs, one per line|
|`data/subset/phenotype_subset.csv`|Phenotype file for selected accessions|

**Expected terminal output:**

```
=============================================
 Step 1: Select 100 accessions
=============================================
VCF contains 1135 samples
Accessions overlapping with VCF: 312
Unique accessions after dedup: 312
Stratified sample across 16 geographic cells
Selected 100 accessions  →  data/subset/sample_ids.txt
...
=============================================
 Summary
=============================================
  Samples in subsetted VCF : 100
  SNPs retained            : 4679817
  VCF output               : data/subset/subset.vcf.gz
  Phenotype output         : data/subset/phenotype_subset.csv
```

---
