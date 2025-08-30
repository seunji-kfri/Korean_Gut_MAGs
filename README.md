# Korean Gut MAGs

This repository contains scripts and workflows used in the study:  
**"Bridging the Amplicon-Metagenome Gap: ASV-Level Pan-genomes via Hybrid MAGs for Studying the Korean Gut Microbiome in Health and Disease"**

---

## 📂 Repository Structure

### 01_QC/
Preprocessing pipelines for Illumina short-read and ONT long-read data.

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
  --outdir 01_QC/out \
  --trimmomatic_adapters adapters/TruSeq3-PE.fa \
  --human_bowtie2_index /path/to/human_index
```
**Long-reads** (default: unmapped only):
```bash
python 01_QC/preprocess_longreads.py \
  --sample SAMPLE01 \
  --in_fastq SAMPLE01.fastq.gz \
  --outdir 01_QC/out \
  --human_mmi /path/to/human.mmi
```
**Long-reads with PAF filter** (identity ≥80%, coverage ≥30%):
```bash
python 01_QC/preprocess_longreads.py \
  --sample SAMPLE01 \
  --in_fastq SAMPLE01.fastq.gz \
  --outdir 01_QC/out \
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

**Long-read assembly:**  
```bash
bash 02_Assembly/run_flye.sh SAMPLE01 ./02_Assembly/out/SAMPLE01_long
```
- Input:
  -SAMPLE01.nohuman.fastq (Long-read QC output)
- Output:
  -02_Assembly/out/SAMPLE01_long/contigs_polished.fasta (Flye assembly output, long-read only)

**Short-read and Hybrid assembly:** 
```bash
bash 02_Assembly/run_operams.sh \
  01_QC/out/SAMPLE01.nohuman.fastq \
  01_QC/out/SAMPLE01_1.nohuman.fq.gz \
  01_QC/out/SAMPLE01_2.nohuman.fq.gz \
  02_Assembly/out/SAMPLE01_hybrid
```
- Input:
  - SAMPLE01.nohuman.fastq (long-read)
  - SAMPLE01_1.nohuman.fq.gz, SAMPLE01_2.nohuman.fq.gz (short-read)

- Output: 
  - 02_Assembly/out/SAMPLE01_hybrid/spades_assembly/contigs.fasta (SPAdes output, short-read only)
  - 02_Assembly/out/SAMPLE01_hybrid/contigs_polished.fasta (OPERA-MS hybrid assembly output, combining short- and long-reads)

---

### 03_Binning/
Binning workflows using MetaWRAP.

- **`run_metawrap.sh`**: Perform binning and refinement using MetaWRAP.

#### Example usage
```bash
bash 03_Binning/run_metawrap.sh SAMPLE01 short_srr_1 ./03_Binning/out
```
- Input:
  - SAMPLE01_1.nohuman.fq.gz, SAMPLE01_2.nohuman.fq.gz (short-read files)
  - 02_Assembly/out/SAMPLE01_hybrid/contigs_polished.fasta (hybrid assembly)

- Output: 
  - Binning and refinement results saved to 03_Binning/out/SAMPLE01.bins

Note:
The second argument, short_srr_1, refers to the prefix of your short-read files.
For example, if your files are named SAMPLE01_1.nohuman.fq.gz and SAMPLE01_2.nohuman.fq.gz, use SAMPLE01 as the short_srr_1 argument.

---

### 04_Annotation/
Functional annotation of genome bins using Prokka and EggNOG.

- **`run_annotation.sh`**: A script that performs both **Prokka** and **EggNOG** annotation for genome bins.

#### Example usage:
```bash
bash run_annotation.sh SAMPLE01 SAMPLE01 ./04_Annotation/out
```
- Input:
  - SAMPLE01_1.nohuman.fq.gz, SAMPLE01_2.nohuman.fq.gz (short-read files)
  - 03_Binning/out/SAMPLE01/SAMPLE01.bins/refined/metawrap_50_10_bins/bin.*.fa (Binning result files from MetaWRAP)

- Output: 
  - Prokka and EggNOG annotation results are saved to ./04_Annotation/out/SAMPLE01/

---

### 05_Taxonomy_Analysis/
Taxonomic classification using GTDB-Tk.

This section uses **GTDB-Tk** to perform taxonomic classification of genome bins.

#### Example usage
```bash
gtdbtk classify_wf \
  --genome_dir ./03_Binning/out/SAMPLE01.bins/refined/metawrap_50_10_bins \
  --out_dir ./04_Taxonomy/out \
  -x fna \
  --cpus 40 \
  --pplacer_cpus 10 \
  --keep_intermediates \
  --write_single_copy_genes \
  --skip_ani_screen
```
- Input:
  - ./03_Binning/out/SAMPLE01.bins/refined/metawrap_50_10_bins/bin.*.fa (Binned genome files)

- Output: 
  - Taxonomic classification results saved to ./04_Taxonomy/out
---

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
