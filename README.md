# Arabidopsis GWAS Pipeline — Class Project Proposal

**Student:** Fritzner Pierre  
**Course:** Scripting for Biologists (BIOL-7180)  
**Date:** February 2026  

---

## 1. Project Title

**Shell-Scripting–Driven GWAS Pipeline for Arabidopsis Flowering Time Using Public Data**

---

## 2. Project Overview

Genome-wide association studies (GWAS) are a fundamental tool in plant genomics, allowing researchers to link genetic variation to phenotypic traits. In this project, I will develop a **fully shell-script–based pipeline** for conducting GWAS on **Arabidopsis thaliana flowering time**, leveraging publicly available data from the **1001 Genomes Project**.  

The pipeline will integrate **VCF processing, quality control, PCA, kinship estimation, GWAS modeling, and result visualization**, while emphasizing **shell scripting** as the main method for automation. R scripts will be called via shell scripts for statistical computation and plotting, allowing students to learn both shell scripting and workflow management in genomics.

---

## 3. Objectives

The goals of this project are to:

1. Develop a **reproducible shell-script–driven workflow** for GWAS.  
2. Process publicly available Arabidopsis SNP and phenotype data.  
3. Perform standard GWAS steps:
   - SNP filtering and quality control  
   - Population structure analysis (PCA)  
   - Kinship estimation  
   - GWAS using mixed linear models (MLM / FarmCPU via GAPIT)  
   - Visualization (Manhattan and QQ plots)  
4. Produce a **self-contained tutorial** that allows other students to run the pipeline without instructor intervention.  
5. Emphasize **automation, modularity, and reproducibility**, demonstrating how shell scripting can manage complex bioinformatics workflows.

---

## 4. Data Source

- **SNP Data:** 1001 Genomes Project VCF files ([link](https://1001genomes.org/data/GMI-MPI/releases/v3.1/))  
- **Phenotype Data:** Flowering time for Arabidopsis accessions ([link](https://arapheno.1001genomes.org/study/12/?utm_source=chatgpt.com))  

> The dataset will be reduced to a manageable subset (~100–200 accessions) for class use.

---

## 5. Methodology

The pipeline is divided into **six main stages**, each implemented as a **standalone shell script**:

1. **VCF Filtering**  
   - Keep biallelic SNPs  
   - Remove SNPs with high missingness (>10%)  

2. **VCF to PLINK Conversion**  
   - Convert VCF to PLINK binary files (.bed, .bim, .fam)  

3. **Quality Control (QC)**  
   - Minor allele frequency filtering (MAF > 0.05)  
   - Individual and SNP missingness filtering  
   - Hardy–Weinberg equilibrium filter (p > 1e-6)  

4. **PCA & Kinship Estimation**  
   - Compute principal components to correct for population structure  
   - Compute kinship matrix to account for relatedness  

5. **GWAS Analysis**  
   - Run GWAS using FarmCPU / MLM models via **GAPIT** in R  
   - Correct for population structure (PCA) and kinship  

6. **Visualization**  
   - Generate Manhattan and QQ plots for interpretation of results  

All steps are automated in `run_pipeline.sh` and designed for **reproducibility and modularity**.

---

## 6. Software and Tools

- **Shell scripting:** Primary workflow automation  
- **PLINK:** SNP QC and conversion  
- **bcftools:** VCF filtering  
- **R & GAPIT3:** PCA, kinship, GWAS modeling, and plotting  
- **qqman (R):** Manhattan and QQ plots  

---

## 7. Expected Deliverables

1. **GitHub Repository** containing:
   - All shell and R scripts  
   - Sample data placeholders  
   - Results directories for QC, PCA/kinship, GWAS, and plots  
   - `run_pipeline.sh` master script  

2. **Tutorial documentation** (`tutorial.md`)  
   - Step-by-step instructions  
   - Explanation of each script and workflow stage  

3. **Final presentation / workshop**  
   - Demonstration of pipeline usage  
   - Interpretation of results  

---

## 8. Learning Outcomes

By completing this project, students (including myself) will learn:

- Practical shell scripting for bioinformatics workflows  
- How to integrate R scripts into automated pipelines  
- Quality control and preprocessing of real genomic data  
- Population structure correction and kinship in GWAS  
- Interpretation of GWAS results (Manhattan & QQ plots)  
- Reproducibility and workflow management in bioinformatics  

---

## 9. Timeline

| Week | Task |
|------|------|
| 1    | Gather data and subset for class use |
| 2    | Develop shell scripts for VCF filtering & PLINK conversion |
| 3    | Implement QC & PCA/Kinship scripts |
| 4    | Implement GWAS & visualization scripts |
| 5    | Test full pipeline, debug, and optimize |
| 6    | Write tutorial and prepare presentation |
| 7    | Submit GitHub repository & present in class |

---

## 10. References

- The 1001 Genomes Consortium. 1001 Genomes Project Data. [https://1001genomes.org/](https://1001genomes.org/)  
- GAPIT3: Genome Association and Prediction Integrated Tool. [https://zzlab.net/GAPIT/](https://zzlab.net/GAPIT/)  
- PLINK: Whole genome association analysis toolset. [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/)  
- qqman: R package for Manhattan and QQ plots. [https://cran.r-project.org/web/packages/qqman/](https://cran.r-project.org/web/packages/qqman/)

---

**This project is designed to emphasize scripting skills while providing a complete end-to-end workflow for a real genomic analysis.**
