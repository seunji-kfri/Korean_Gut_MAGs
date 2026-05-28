# Korean Gut MAGs

This repository contains the scripts and workflows used in the study:

**"Population-specific MAG-guided machine learning enables functional inference
of the Korean gut microbiome from 16S rRNA sequencing"**

The pipeline (i) builds a Korean-specific gut metagenome-assembled genome (MAG)
catalog through hybrid (short- + long-read) assembly, and (ii) uses an
XGBoost-based machine learning framework to predict MAG cluster abundances from
16S rRNA amplicon sequence variants (ASVs) and convert them into KEGG Orthology
(KO) functional profiles. This allows functional profiles to be inferred from
amplicon-only cohorts without paired shotgun sequencing.

> **Citation:** Manuscript in preparation / under review. Citation details and
> DOI will be added here upon publication.

---

## Quick start: apply the trained model to your own amplicon data

If you only want to infer KO functional profiles from your own 16S ASV table,
you do **not** need any of the cohort data below. You need three things, two of
which are already in this repository:

1. your ASV relative-abundance table (first column = ASV id, other columns = samples),
2. `07_ML_Prediction/models/full_xgb_models.rds` (the trained model, included), and
3. `07_ML_Prediction/models/ko_cluster_table.tsv` (the carrier table, included).

```bash
Rscript 07_ML_Prediction/predict.R \
  --asv     your_asv_relabund.csv \
  --models  07_ML_Prediction/models/full_xgb_models.rds \
  --carrier 07_ML_Prediction/models/ko_cluster_table.tsv \
  --out     predicted_sample_x_ko.tsv
```

The output is a sample × KO relative-abundance table. ASVs that are not in the
training feature set are ignored; training ASVs absent from your table are
zero-padded. See the `07_ML_Prediction/` section for details and options.

Reproducing the full analysis (training, benchmarking, biology) additionally
requires cohort-level input data — see "Input data for full reproduction" below.

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
│   ├── 01_train_step1_cluster.R   Step 1: ASV → cluster (XGBoost, CV)
│   ├── 02_predict_step2_ko.R      Step 2: cluster → KO (carrier matrix)
│   ├── 03_benchmark.R             Benchmark vs link/PICRUSt2 (Fig. 2/4/S5)
│   ├── 04_external_apply.R        Train full model + apply to new cohort (Fig. 4)
│   ├── 05_biology.R               MASLD pathway interpretation (Fig. 5)
│   ├── predict.R          Standalone prediction tool (apply trained model)
│   ├── design_choices/    Design-choice analyses (Fig. 3)
│   └── models/            Pre-trained model + carrier table (reusable)
├── environment.yml        Conda environment (assembly/annotation tools)
└── README.md
```

Sections 01–05 reconstruct and annotate the MAG catalog. Sections 06–07 infer
functional profiles from amplicon data: `06_Linking` is the reference-matching
baseline, and `07_ML_Prediction` is the machine-learning framework that is the
main contribution of the study.

---

### 01_QC/
Preprocessing pipelines for Illumina short-read and ONT long-read data.

- **`preprocess_shortreads.py`** — Illumina reads: trimming (Trimmomatic),
  human read removal (Bowtie2), sequence statistics (seqkit).
- **`preprocess_longreads.py`** — ONT reads: adapter trimming (Porechop),
  quality filtering (Filtlong), human read removal (default: keep unmapped
  reads; option: PAF-based filter, identity ≥80%, coverage ≥30%), seqkit stats.

```bash
python 01_QC/preprocess_shortreads.py \
  --sample SAMPLE01 \
  --r1 SAMPLE01_R1.fastq.gz --r2 SAMPLE01_R2.fastq.gz \
  --outdir 01_QC/out \
  --trimmomatic_adapters adapters/TruSeq3-PE.fa \
  --human_bowtie2_index /path/to/human_index

python 01_QC/preprocess_longreads.py \
  --sample SAMPLE01 --in_fastq SAMPLE01.fastq.gz \
  --outdir 01_QC/out --human_mmi /path/to/human.mmi
```

---

### 02_Assembly/
- **`run_flye.sh`** — long-read assembly (Flye)
- **`run_operams.sh`** — short-read and hybrid assembly (OPERA-MS)

```bash
bash 02_Assembly/run_flye.sh SAMPLE01 ./02_Assembly/out/SAMPLE01_long
bash 02_Assembly/run_operams.sh \
  01_QC/out/SAMPLE01.nohuman.fastq \
  01_QC/out/SAMPLE01_1.nohuman.fq.gz \
  01_QC/out/SAMPLE01_2.nohuman.fq.gz \
  02_Assembly/out/SAMPLE01_hybrid
```

---

### 03_Binning/
- **`run_metawrap.sh`** — binning and refinement (MetaWRAP).

```bash
bash 03_Binning/run_metawrap.sh SAMPLE01 short_srr_1 ./03_Binning/out
```
The second argument is the prefix of your short-read files
(`SAMPLE01_1.nohuman.fq.gz` / `SAMPLE01_2.nohuman.fq.gz` → `SAMPLE01`).

---

### 04_Annotation/
- **`run_annotation.sh`** — Prokka + eggNOG annotation for genome bins.

```bash
bash 04_Annotation/run_annotation.sh SAMPLE01 SAMPLE01 ./04_Annotation/out
```

---

### 05_Taxonomy/
Taxonomic classification of genome bins using GTDB-Tk (GTDB r226 database).

- **`run_gtdbtk.sh`** — wraps `gtdbtk classify_wf`.

```bash
export GTDBTK_DATA_PATH=/path/to/gtdbtk/db
bash 05_Taxonomy/run_gtdbtk.sh <genome_dir> 05_Taxonomy/out allbins
```
- Input: binned genome files (`bin.*.fna`)
- Output: GTDB-Tk classification (`allbins.bac120.summary.tsv`, etc.)

---

### 06_Linking/
Reference-based (link) functional inference: ASVs are linked to MAGs via their
16S rRNA genes, and KO content is aggregated across linked MAGs. This is the
baseline against which the ML framework is compared (Fig. S2; benchmarked in
`07_ML_Prediction/03_benchmark.R`).

- **`build_asv_features.py`** — builds per-ASV KO presence features from per-bin
  annotations (completeness-weighted or binary).
- **`01_build_linked_tables.R`** — sample × KO tables for six linking configs.
- **`02_run_maaslin2.R`** — MaAsLin2 differential abundance per config + gold.
- **`03_run_wilcoxon.R`** — per-KO Wilcoxon test per config + gold.
- **`04_compute_metrics.R`** — F1 vs gold and per-sample Spearman similarity.

Run in order (`01` → `02` → `03` → `04`) from the repository root, after placing
the required inputs in `06_Linking/input/` (see below).

---

### 07_ML_Prediction/
The two-step machine-learning framework (main method). An XGBoost model predicts
dRep95 cluster abundances from CLR-transformed ASV abundances (Step 1), and a
binary carrier matrix converts predicted cluster abundances into KO profiles
(Step 2). The model is trained once on the paired internal cohort and can then
be applied to any amplicon-only cohort.

#### Configuration
All paths and parameters live in **`config.R`**. Run R from the repository root,
or set the root explicitly:
```bash
export KGM_ROOT=/path/to/Korean_Gut_MAGs
```
Only edit the `FILES` list in `config.R` if your input file names differ.

#### Scripts and the figures they produce
| Script | Purpose | Manuscript |
|--------|---------|-----------|
| `01_train_step1_cluster.R` | ASV → cluster, 5-fold CV (XGBoost) | Step 1 |
| `02_predict_step2_ko.R` | cluster → KO via carrier matrix; F1 vs gold | Step 2 |
| `03_benchmark.R` | XGB vs link vs PICRUSt2 (MaAsLin2, Age/Sex/BMI adj.) | Fig. 2/4/S5 |
| `04_external_apply.R` | train full model, apply to external cohort | Fig. 4 |
| `05_biology.R` | MASLD pathway enrichment & interpretation | Fig. 5 |
| `design_choices/claim1_ko_specificity.R` | ML benefit by KO specificity | Fig. 3G |
| `design_choices/claim2_nonlinearity.R` | XGB vs Ridge (non-linearity) | Fig. 3B-E |
| `design_choices/claim3_cluster_resolution.R` | dRep95 vs dRep99 | Fig. 3A |
| `design_choices/claim4_aggregation.R` | binary carrier vs FPKM | Fig. 3F |
| `predict.R` | standalone tool to apply the trained model | — |

Typical reproduction order: `01` → `02` → `03` → `05`, then the
`design_choices/` scripts; `04` trains the full-cohort model and applies it to
an external cohort.

#### Input data for full reproduction (not included)
The trained model and carrier table in `models/` are sufficient for `predict.R`.
Full reproduction additionally requires the following per-sample tables in
`07_ML_Prediction/input/` (and the analogous files in `06_Linking/input/`).
These contain cohort/participant-level data and are available from the authors
on reasonable request.

| File | Format (key columns) |
|------|----------------------|
| `ASV_RA_224.csv` | first column `ASV`; remaining columns = sample IDs; values = relative abundance |
| `ASV_RA_external.csv` | same layout, external cohort |
| `drep_95_bigliver_hybrid.csv` | `genome` (e.g. `<bin>.fna`), `secondary_cluster` |
| `taxonomy_cpm_pivot_gtdb.csv` | `sample_id`, `bin_id`, `CPM_RA` (long) |
| `all_bins_KEGG_ko_FPKM_pivot.csv` | first column `bin_id`; remaining columns = KO ids (FPKM) |
| `fpkm_gold_long.csv` | `sample_id`, `KO`, `fpkm_sum` (shotgun gold, long) |
| `sample_meta.csv` | `sample_id`, `group` (Healthy/MASLD), `Age`, `Sex`, `BMI` |
| `meta_external.csv` | `sample_id`, `group`, `Age`/`AgeGroup`, `Sex`, `BMI` |
| `pic_rel.csv`, `picrust2_ext_ko.tsv` | first column = sample id; KO columns (PICRUSt2) |
| `ko_pathway_hierarchy_merged.csv` | `feature` (KO), `Pathway`, `Level1`, `Level2`, `Level3` |

Generated outputs and intermediate tables are written under each section's
`output/` directory and are git-ignored (see `.gitignore`).

---

## 🔧 Required resources (not included in this repo)
1) Trimmomatic adapters (`TruSeq3-PE.fa`) — provided with Trimmomatic.
2) Human reference indexes (GRCh38 or T2T-CHM13 FASTA) — build a Bowtie2 index
   (`bowtie2-build`) and a Minimap2 index (`minimap2 -d`).
3) GTDB-Tk release database (set `GTDBTK_DATA_PATH`).

---

## ⚙️ Environment

**Assembly / annotation tools (conda):**
```bash
conda env create -f environment.yml
conda activate KoreanGutMAGs
```

**R (sections 06–07):** R ≥ 4.3 with `data.table`, `tidyverse`, `xgboost`,
`glmnet`, and `Maaslin2`:
```r
install.packages(c("data.table", "tidyverse", "xgboost", "glmnet"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Maaslin2")
```

---

## 📄 License
See the `LICENSE` file in this repository.
