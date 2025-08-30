# Korean Gut MAGs

This repository contains scripts and workflows used in the study:  
**"Bridging the Amplicon-Metagenome Gap: ASV-Level Pan-genomes via Hybrid MAGs for Studying the Korean Gut Microbiome in Health and Disease"**

---

## 📂 Repository Structure

### 01_QC/
Illumina short-read and ONT long-read preprocessing pipelines.

- **`preprocess_shortreads.py`**  
  Short-read preprocessing pipeline for Illumina reads:  
  - Trimming (Trimmomatic)  
  - Human read removal (Bowtie2)  
  - Sequence statistics (seqkit)

- **`preprocess_longreads.py`**  
  Long-read preprocessing pipeline for ONT reads:  
  - Adapter trimming (Porechop)  
  - Quality filtering (Filtlong)  
  - Human read removal  
    - **Default**: keep unmapped reads only (samtools -f 4)  
    - **Option**: PAF-based filtering (identity ≥80%, coverage ≥30%)  
  - Sequence statistics (seqkit)

#### Example usage

**Short-reads:**
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
Assembly workflows for short-read, long-read, and hybrid assemblies.

This section contains scripts for performing assembly of short and long reads using different tools.

- **`run_flye.sh`**: Long-read assembly using Flye
- **`run_operams.sh`**: Short-read and Hybrid assembly using OPERA-MS (short-read / short-read + long-read)

#### Example usage

**Long-read assembly (Flye)**  
Run Flye to assemble ONT long-reads.
```bash
bash 02_Assembly/run_flye.sh SAMPLE01
```
- Input: SAMPLE01.nohuman.fastq (Long-read QC output)
- Output: SAMPLE01/contigs.fasta (Flye assembly output, long-read only)

**Short-read and Hybrid assembly (OPERA-MS with SPAdes output)** 
Run OPERA-MS to perform hybrid assembly using short and long reads. 
The short-read assembly output comes from SPAdes (stored in spades_assembly/contigs.fasta).
```bash
bash 02_Assembly/run_operams.sh \
  /path/to/SAMPLE01.nohuman.fastq \
  /path/to/SAMPLE01_1.nohuman.fq.gz \
  /path/to/SAMPLE01_2.nohuman.fq.gz \
  02_Assembly/out/SAMPLE01_hybrid
```
- Input:
  - SAMPLE01.nohuman.fastq (long-read)
  - SAMPLE01_1.nohuman.fq.gz, SAMPLE01_2.nohuman.fq.gz (short-read)

- Output: 
  - SAMPLE01_hybrid/intermediate_files/spades_assembly/contigs.fasta (SPAdes output, short-read only)
  - SAMPLE01_hybrid/contigs.polished.fasta (OPERA-MS hybrid assembly output, combining short- and long-reads)

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
