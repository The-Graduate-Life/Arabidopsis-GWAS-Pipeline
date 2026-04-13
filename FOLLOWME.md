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
>All commands below assume you are in the **project root**.

---

