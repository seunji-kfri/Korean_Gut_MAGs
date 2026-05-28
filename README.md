# Korean Gut MAGs

This repository contains the scripts and workflows used in the study:

**"Population-specific MAG-guided machine learning enables functional inference
of the Korean gut microbiome from 16S rRNA sequencing"**

The pipeline (i) builds a Korean-specific gut metagenome-assembled genome (MAG)
catalog through hybrid (short- + long-read) assembly, and (ii) uses an
XGBoost-based machine learning framework to predict MAG cluster abundances from
16S rRNA amplicon sequence variants (ASVs) and convert them into KEGG Orthology
(KO) functional profiles. This lets functional profiles be inferred from
amplicon-only cohorts without paired shotgun sequencing.

> **Citation:** Manuscript in preparation / under review. Citation details and
> DOI will be added here upon publication.

---

## 📂 Repository Structure

```
Korean_Gut_MAGs/
├── 01_QC/                 Read QC (Illumina short reads, ONT long reads)
├── 02_Assembly/           Short-read, long-read, and hybrid assembly
├── 03_Binning/            Binning and refinement (MetaWRAP)
├── 04_Annotation/         Functional annotation (Prokka, eggNOG-mapper)
├── 05_Taxonomy/           Taxonomic classification (GTDB-Tk)
├── 06_Linking/            ASV–MAG linking (reference-based baseline)
├── 07_ML_Prediction/      Two-step ML framework (main method)
│   ├── config.R           Shared paths & parameters (edit this one file)
│   ├── predict.R          Standalone prediction tool (apply trained model)
│   ├── 04_external_apply.R  Train full-cohort model + apply to a new cohort
│   └── models/            Pre-trained model + carrier table (reusable assets)
├── environment.yml        Conda environment (assembly/annotation tools)
└── README.md
```

The pipeline has two parts. **Sections 01–05** reconstruct and annotate the MAG
catalog. **Sections 06–07** infer functional profiles from amplicon data:
`06_Linking` is the reference-matching baseline, and `07_ML_Prediction` is the
machine-learning framework that is the main contribution of the study.

---

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

- **`run_flye.sh`**: Long-read assembly using Flye
- **`run_operams.sh`**: Short-read and hybrid assembly using OPERA-MS

#### Example usage

**Long-read assembly:**
```bash
bash 02_Assembly/run_flye.sh SAMPLE01 ./02_Assembly/out/SAMPLE01_long
```
- Input: `SAMPLE01.nohuman.fastq` (long-read QC output)
- Output: `02_Assembly/out/SAMPLE01_long/contigs_polished.fasta`

**Short-read and hybrid assembly:**
```bash
bash 02_Assembly/run_operams.sh \
  01_QC/out/SAMPLE01.nohuman.fastq \
  01_QC/out/SAMPLE01_1.nohuman.fq.gz \
  01_QC/out/SAMPLE01_2.nohuman.fq.gz \
  02_Assembly/out/SAMPLE01_hybrid
```
- Input: `SAMPLE01.nohuman.fastq` (long-read);
  `SAMPLE01_1.nohuman.fq.gz`, `SAMPLE01_2.nohuman.fq.gz` (short-read)
- Output:
  - `02_Assembly/out/SAMPLE01_hybrid/spades_assembly/contigs.fasta` (short-read only)
  - `02_Assembly/out/SAMPLE01_hybrid/contigs_polished.fasta` (hybrid)

---

### 03_Binning/
Binning workflows using MetaWRAP.

- **`run_metawrap.sh`**: Perform binning and refinement using MetaWRAP.

#### Example usage
```bash
bash 03_Binning/run_metawrap.sh SAMPLE01 short_srr_1 ./03_Binning/out
```
- Input: `SAMPLE01_1.nohuman.fq.gz`, `SAMPLE01_2.nohuman.fq.gz` (short-read);
  `02_Assembly/out/SAMPLE01_hybrid/contigs_polished.fasta` (hybrid assembly)
- Output: binning and refinement results in `03_Binning/out/SAMPLE01.bins`

> Note: the second argument `short_srr_1` is the prefix of your short-read
> files. For `SAMPLE01_1.nohuman.fq.gz` / `SAMPLE01_2.nohuman.fq.gz`, use
> `SAMPLE01`.

---

### 04_Annotation/
Functional annotation of genome bins using Prokka and eggNOG-mapper.

- **`run_annotation.sh`**: Runs both **Prokka** and **eggNOG** annotation for
  genome bins.

#### Example usage
```bash
bash 04_Annotation/run_annotation.sh SAMPLE01 SAMPLE01 ./04_Annotation/out
```
- Input: short-read files; `bin.*.fa` from MetaWRAP refinement
- Output: Prokka and eggNOG annotation in `./04_Annotation/out/SAMPLE01/`

---

### 05_Taxonomy/
Taxonomic classification of genome bins using GTDB-Tk.

#### Example usage
```bash
gtdbtk classify_wf \
  --genome_dir ./03_Binning/out/SAMPLE01.bins/refined/metawrap_50_10_bins \
  --out_dir ./05_Taxonomy/out \
  -x fna \
  --cpus 40 \
  --pplacer_cpus 10 \
  --keep_intermediates \
  --write_single_copy_genes \
  --skip_ani_screen
```
- Input: binned genome files (`bin.*.fa`)
- Output: taxonomic classification in `./05_Taxonomy/out`

---

### 06_Linking/
Reference-based (link) functional inference: ASVs are linked to MAGs via their
16S rRNA genes, and KO content is aggregated across linked MAGs. This serves as
the baseline against which the machine-learning framework is compared
(manuscript Fig. S2).

- **`build_asv_features.py`**: builds per-ASV KO presence features from
  per-bin annotations (completeness-weighted or binary), used as input to the
  linking analysis.

#### Example usage
```bash
python 06_Linking/build_asv_features.py \
  --basedir /path/to/per_ASV_folders \
  --annotation_base /path/to/annotation_root \
  --mode binary
```
- Input: per-ASV folders containing `binsummary.txt` (and optional Roary output)
- Output: `asv_ko_*.csv`, `asv_keggmodule_*.csv` (ASV × KO/module tables)

---

### 07_ML_Prediction/
The two-step machine-learning framework (main method). An XGBoost model predicts
dRep95 cluster abundances from CLR-transformed ASV abundances (Step 1), and a
binary carrier matrix converts predicted cluster abundances into KO profiles
(Step 2). The model is trained once on the paired internal cohort and can then
be applied to any amplicon-only cohort.

```
07_ML_Prediction/
├── config.R               Shared paths & parameters — edit ONLY this file
├── predict.R              Standalone tool: apply the trained model to new ASVs
├── 04_external_apply.R    Train full-cohort model and apply to a new cohort
├── models/
│   ├── full_xgb_models.rds    Pre-trained XGBoost models (trained on 224 samples)
│   └── ko_cluster_table.tsv   Cluster × KO carrier table (carrier_frac column)
├── input/                 User-supplied data (not tracked; see "Input data")
└── output/                Generated results (not tracked)
```

#### Configuration
All paths and parameters live in **`config.R`**. To run the pipeline elsewhere,
either run R from the repository root, or set the root explicitly:
```bash
export KGM_ROOT=/path/to/Korean_Gut_MAGs
```
No other file needs editing unless your input file names differ (edit the
`FILES` list in `config.R`).

#### Quick start — apply the trained model to your own amplicon data
`predict.R` is the external-facing entry point. Given a 16S ASV
relative-abundance table, it outputs an inferred sample × KO profile using the
pre-trained model and carrier table in `models/` — no retraining required.
```bash
Rscript 07_ML_Prediction/predict.R \
  --asv     your_asv_relabund.csv \
  --models  07_ML_Prediction/models/full_xgb_models.rds \
  --carrier 07_ML_Prediction/models/ko_cluster_table.tsv \
  --out     predicted_sample_x_ko.tsv
```
Optional arguments: `--prev prev10|prev20` (model set, default `prev10`),
`--thr 0.25` (carrier-fraction threshold, default 0.25).

**Input ASV table format:** CSV/TSV with the first column = ASV id and the
remaining columns = samples; values are per-sample relative abundance. ASVs
absent from the user table are zero-padded; ASVs not seen during training are
ignored.

#### Reproducing the cross-cohort analysis
`04_external_apply.R` trains the full-cohort model on the internal paired data,
saves it to `models/`, applies it to an independent amplicon-only cohort, and
computes the cross-cohort consistency statistics (manuscript Fig. 4).
```bash
Rscript 07_ML_Prediction/04_external_apply.R
```
This requires the input data described below to be present in
`07_ML_Prediction/input/`.

#### Input data (not included in this repository)
The trained model and carrier table in `models/` are sufficient for `predict.R`.
Reproducing the full analysis additionally requires the following per-sample
tables in `07_ML_Prediction/input/` (cohort/participant-level data, available
from the authors on reasonable request):

| File | Description |
|------|-------------|
| `ASV_RA_224.csv` | Internal cohort ASV × sample relative abundance |
| `drep_95_bigliver_hybrid.csv` | Genome → dRep95 cluster map |
| `taxonomy_cpm_pivot_gtdb.csv` | Per-sample MAG CPM abundance |
| `sample_meta.csv` | Internal sample metadata (`sample_id`, `group`) |
| `all_bins_KEGG_ko_FPKM_pivot.csv` | Per-bin KO FPKM (for carrier profile) |
| `fpkm_gold_long.csv` | Shotgun-derived gold KO abundance |
| `ASV_RA_external.csv` | External cohort ASV table |
| `meta_external.csv` | External sample metadata |
| `pic_rel.csv`, `picrust2_ext_ko.tsv` | PICRUSt2 profiles (for benchmarking) |

---

## 🔧 Required resources (not included in this repo)
1) Trimmomatic adapters (`TruSeq3-PE.fa`) — provided with Trimmomatic.
2) Human reference indexes (GRCh38 or T2T-CHM13 FASTA) — build a Bowtie2 index
   (`bowtie2-build`) and a Minimap2 index (`minimap2 -d`).

---

## ⚙️ Environment

**Assembly / annotation tools (conda):**
```bash
conda env create -f environment.yml
conda activate KoreanGutMAGs
```

**R (sections 06–07):** R ≥ 4.3 with `data.table`, `tidyverse`, `xgboost`,
`glmnet`, and `Maaslin2`. For example:
```r
install.packages(c("data.table", "tidyverse", "xgboost", "glmnet"))
# Maaslin2 from Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Maaslin2")
```

---

## 📄 License
See the `LICENSE` file in this repository.
