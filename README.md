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
  - Human read removal  
    - **Default**: keep unmapped reads only (samtools -f 4)  
    - **Option**: PAF-based filtering (identity ≥80%, coverage ≥30%)  
  - Sequence statistics (seqkit)

#### Example usage
Short-reads:
```bash
python 01_QC/preprocess_shortreads.py \
  --sample SAMPLE01 \
  --r1 SAMPLE01_R1.fastq.gz --r2 SAMPLE01_R2.fastq.gz \
  --outdir preprocess_out \
  --trimmomatic_adapters adapters/TruSeq3-PE.fa \
  --human_bowtie2_index /path/to/human_index
```
Long-reads (default: unmapped only):
```bash
python 01_QC/preprocess_longreads.py \
  --sample SAMPLE01 \
  --in_fastq SAMPLE01.fastq.gz \
  --outdir preprocess_out \
  --human_mmi /path/to/human.mmi
```
Long-reads with PAF filter (identity ≥80%, coverage ≥30%):
```bash
python 01_QC/preprocess_longreads.py \
  --sample SAMPLE01 \
  --in_fastq SAMPLE01.fastq.gz \
  --outdir preprocess_out \
  --human_mmi /path/to/human.mmi \
  --paf_filter 80 30 \
  --paf_parser 01_QC/paf_parser.py
```
---

### 02_Assembly/
- Hybrid and single-technology assembly scripts
  
---

- `03_Binning/` : Binning workflows (MetaBAT, MaxBin, etc.)
- `04_Annotation/` : Functional annotation (Prokka, EggNOG-mapper, etc.)
- `05_Taxonomy_Analysis/` : Taxonomic classification (GTDB-Tk) and summary scripts

## 🔧 Required resources (not included in this repo)
1) Trimmomatic adapters (TruSeq3-PE.fa)  
   - provided with Trimmomatic installation  

2) Human reference indexes (GRCh38 or T2T-CHM13 FASTA)  
   - build Bowtie2 index (`bowtie2-build`) and Minimap2 index (`minimap2 -d`)

---

## ⚙️ Environment
We used conda-based environments.  
Recreate the environment as:
```bash
conda env create -f environment.yml
conda activate KoreanGutMAGs
