#!/bin/bash
# Usage: run_metawrap.sh <sample_id> <short_srr_1>
# <sample_id> : The ID of the sample (e.g., SAMPLE01)
# <short_srr_1> : The prefix of your short-read files (e.g., SAMPLE01)

# Activate conda environment (ensure the environment created with 'environment.yml' is available)
conda activate KoreanGutMAGs  # Ensure this is the environment created with 'environment.yml'

# Define input/output paths based on the directory structure in the project
ASSEMBLY_DIR=./02_Assembly/out   # Path to assembly results from 02_Assembly
BINNING_DIR=./03_Binning/out  # Final output directory for binning results
FASTQ_DIR=./01_QC/out  # Preprocessed data directory from 01_QC

# Define input/output files
SAMPLE=$1
SHORT_SRR=$2  # short_srr_1 refers to the prefix of your short-read files

# Correct path to the assembly output
ASSEMBLY=${ASSEMBLY_DIR}/${SAMPLE}_hybrid/contigs_polished.fasta  # Assembly from 02_Assembly (hybrid output)

# Input files (should be user-provided)
FASTQ1_ORIG=${FASTQ_DIR}/${SHORT_SRR}_1.nohuman.fq.gz
FASTQ2_ORIG=${FASTQ_DIR}/${SHORT_SRR}_2.nohuman.fq.gz

# Create necessary directories for final output
mkdir -p $BINNING_DIR/${SAMPLE}

# Unzip the fastq files
gunzip -c $FASTQ1_ORIG > ./reads_1.fastq
gunzip -c $FASTQ2_ORIG > ./reads_2.fastq

# Step 1: Binning
metawrap binning \
  -o $BINNING_DIR/${SAMPLE}.bins \
  -a $ASSEMBLY \
  --metabat2 --maxbin2 --concoct \
  ./reads_1.fastq ./reads_2.fastq

# Step 2: Refinement
metawrap bin_refinement \
  -o $BINNING_DIR/${SAMPLE}.bins/refined \
  -c 50 -x 10 \
  -A $BINNING_DIR/${SAMPLE}.bins/metabat2_bins \
  -B $BINNING_DIR/${SAMPLE}.bins/maxbin2_bins \
  -C $BINNING_DIR/${SAMPLE}.bins/concoct_bins

# Clean up intermediate files
rm -f ./reads_1.fastq ./reads_2.fastq $ASSEMBLY

echo "[INFO] Binning and refinement for ${SAMPLE} completed and saved to ${BINNING_DIR}/${SAMPLE}.bins"
