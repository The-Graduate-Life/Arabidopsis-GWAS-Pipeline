**Author:** Fritzner Pierre  
**Course:** Scripting for Biologists (BIOL-7180)  
**Semester:** Spring 2026    
**Institution/Affiliation**: Auburn University

---

## Project: **GWAS Analysis Pipeline for *Arabidopsis thaliana* Flowering Time Using Public Data**

### Overview:
This tutorial walks you through a complete genome-wide association study (GWAS) pipeline using publicly available data from *Arabidopsis thaliana*. By the end, you will have run every step from raw VCF filtering to Manhattan and QQ plots using a series of modular shell scripts.

---

## 1\. Background

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
