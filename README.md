**Author:** Fritzner Pierre  
**Course:** Scripting for Biologists (BIOL-7180)  
**Semester:** Spring 2026    
**Institution/Affiliation**: Auburn University

---

## Project: **GWAS Analysis Pipeline for *Arabidopsis thaliana* Flowering Time Using Public Data**

---

## 1. Project Overview

Genome-wide association studies (GWAS) are a fundamental tool in plant genomics that allows researchers to link genetic variation to phenotypic traits. In this project, I will develop a **fully shell-script–based pipeline** for conducting GWAS on ***Arabidopsis thaliana* flowering time**, leveraging publicly available data from the **[1001 Genomes Project](https://1001genomes.org/)**.  

The pipeline will combine **VCF processing, quality control, PCA, GWAS modeling, and result visualization**, with the emphasis on **shell scripting** as the main method for automation. Some R scripts might be called via shell scripts for statistical computation and plotting if needed. This project focuses on both shell scripting and workflow management in genomics.



## 2. Objectives

Throughout this project, I aim to:

1. Develop a **reproducible shell-script–driven workflow** for GWAS analysis.  
2. Process publicly available Arabidopsis SNP and phenotype data.  
3. Perform standard GWAS steps:
   - SNP filtering and quality control  
   - Population structure analysis (PCA)  
   - GWAS using mixed linear models (MLM / FarmCPU via GAPIT)  
   - Visualization (Manhattan and QQ plots)  
4. Produce a **complete tutorial** that allows other researchers to run the pipeline without instructor intervention.  
5. Emphasize **automation, modularity, and reproducibility** to demonstrate how shell scripting can manage complex bioinformatics workflows.



## 3. Data Source

- **SNP Data:** 1001 Genomes Project VCF files ([link](https://1001genomes.org/data/GMI-MPI/releases/v3.1/))  
- **Phenotype Data:** Flowering time for Arabidopsis accessions ([link](https://arapheno.1001genomes.org/study/12/?utm_source=chatgpt.com))  

> The dataset will be reduced to a manageable subset (~100–200 accessions) for class use.



## 4. Methodology

The pipeline is divided into **six main stages**, each will be implemented as a **standalone shell script**:

1. **VCF Filtering**  
   - To keep biallelic SNPs  
   - To remove SNPs with high missingness (>10%)  

2. **VCF to PLINK Conversion**  
   - To convert VCF to PLINK binary files (.bed, .bim, .fam)  

3. **Quality Control (QC)**  
   - Minor allele frequency filtering (MAF > 0.05)  
   - Individual and SNP missingness filtering  
   - Hardy–Weinberg equilibrium filter (p > $1 \times 10^{-6}$)  

4. **PCA & Kinship Estimation**  
   - To compute principal components to correct for population structure  
   - To compute kinship matrix to account for relatedness  

5. **GWAS Analysis**  
   - Run GWAS using FarmCPU / MLM models via **GAPIT** in R  
   - Correct for population structure (PCA) and kinship  

6. **Visualization**  
   - Generate Manhattan and QQ plots for interpretation of results  

All steps will be automated using shell script `run_pipeline.sh` as master script and designed for **reproducibility and modularity**.



## 5. Software and Tools that might be used in this project

- **Shell scripting:** Primary workflow automation  
- **PLINK:** SNP QC and conversion  
- **bcftools:** VCF filtering  
- **R & GAPIT3:** PCA, kinship, GWAS modeling, and plotting  
- **qqman package (R):** Manhattan and QQ plots  



## 6. Expected Deliverables

1. **GitHub Repository** containing:
   - All shell and R scripts  
   - Sample data placeholders  
   - Results directories for QC, PCA/kinship, GWAS, and plots  
   - The master script `run_pipeline.sh`  

2. **Tutorial documentation** (`tutorial.md`)  
   - Step-by-step instructions  
   - Explanation of each script and workflow stage  



## 7. Learning Outcomes

By completing this project, I will learn:

- Practical shell scripting for bioinformatics workflows  
- How to integrate R scripts into automated pipelines  
- Quality control and preprocessing of real genomic data   
- Interpretation of GWAS results (Manhattan & QQ plots)  
- Reproducibility and workflow management in scripting  



## 8. Step-by-Step Flow

| Steps | Tasks |
|------|------|
| 1    | Gather data and subset for class use |
| 2    | Develop shell scripts for VCF filtering & PLINK conversion |
| 3    | Implement QC & PCA scripts |
| 4    | Implement GWAS & visualization scripts |
| 5    | Test full pipeline, debug, and optimize |
| 6    | Write tutorial and prepare presentation |
| 7    | Submit GitHub repository |



## 9. References

- The 1001 Genomes Consortium. 1001 Genomes Project Data. [https://1001genomes.org/](https://1001genomes.org/)  
- GAPIT3: Genome Association and Prediction Integrated Tool. [https://zzlab.net/GAPIT/](https://zzlab.net/GAPIT/)  
- PLINK: Whole genome association analysis toolset. [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/)  
- qqman: R package for Manhattan and QQ plots. [https://cran.r-project.org/web/packages/qqman/](https://cran.r-project.org/web/packages/qqman/)

---

## Feedback

This project idea sounds good.
I like the goal of making your code modular and extensible; try to make each
component (e.g., function or script) do one thing and do it well.
You are welcome to use Bash and R for your project.
[This Bash style guide](https://google.github.io/styleguide/shellguide.html)
will be helpful for translating the best practices we learn in class with
Python to Bash.
For example,
[the section about comments](https://google.github.io/styleguide/shellguide.html#comments)
shows how to follow best practices for documenting Bash code (i.e., the
equivalent of docstrings in Python).
Similarly,
[roxygen2](https://roxygen2.r-lib.org/)
provides a way to follow documentation best practices for R code.
All the best practices we learn in class using Python will be applicable to
your Bash and R coding.
