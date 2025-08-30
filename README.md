# Korean Gut MAGs

This repository contains scripts and workflows used in the study:  
**"Bridging the Amplicon-Metagenome Gap: ASV-Level Pan-genomes via Hybrid MAGs for Studying the Korean Gut Microbiome in Health and Disease"**

---

## 📂 Repository Structure
### 01_QC/
- `preprocess_shortreads.py`  
  Illumina short-read preprocessing pipeline:  
  - Trimming (Trimmomatic)  
  - Human read removal (Bowtie2)  
  - Sequence statistics (seqkit)

- `preprocess_longreads.py`  
  ONT long-read preprocessing pipeline:  
  - Adapter trimming (Porechop)  
  - Quality filtering (Filtlong)  
  - Human read removal (Minimap2 + seqkit)  
  - Sequence statistics (seqkit)

- `02_Assembly/` : Hybrid and single-technology assembly scripts
- `03_Binning/` : Binning workflows (MetaBAT, MaxBin, etc.)
- `04_Annotation/` : Functional annotation (Prokka, EggNOG-mapper, etc.)
- `05_Taxonomy_Analysis/` : Taxonomic classification (GTDB-Tk) and summary scripts

---

## ⚙️ Environment
We used conda-based environments.  
Recreate the environment as:
```bash
conda env create -f environment.yml
conda activate KoreanGutMAGs
