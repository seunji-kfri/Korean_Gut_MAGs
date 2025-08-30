# Korean Gut MAGs

This repository contains scripts and workflows used in the study:  
**"Bridging the Amplicon-Metagenome Gap: ASV-Level Pan-genomes via Hybrid MAGs for Studying the Korean Gut Microbiome in Health and Disease"**

---

## 📂 Repository Structure
- `01_QC/` : Short- and long-read quality control scripts
- `02_Assembly/` : Hybrid and single-technology assembly scripts
- `03_Binning/` : Binning workflows (MetaBAT, MaxBin, etc.)
- `04_Annotation/` : Functional annotation (Prokka, EggNOG-mapper, etc.)
- `05_Taxonomy_analysis/` : KO, GO, EC, CAZy functional pipelines

---

## ⚙️ Environment
We used conda-based environments.  
Recreate the environment as:
```bash
conda env create -f environment.yml
conda activate KoreanGutMAGs
