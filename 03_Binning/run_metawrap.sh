#!/bin/bash
# Usage: run_metawrap.sh <sample_id> <short_srr>

# Activate conda environment (ensure the environment created with 'environment.yml' is available)
conda activate KoreanGutMAGs  # Ensure this is the environment created with 'environment.yml'

# Define input/output paths based on the directory structure in the project
ASSEMBLY_DIR=./02_Assembly/out   # Path to assembly results from 02_Assembly
FINALDIR=./operams_results  # Final output directory for binning results
FASTQ_DIR=./preprocessed_data  # Preprocessed data directory from 01_QC

# Define input/output files
SAMPLE=$1
SHORT_SRR=$2

# Correct path to the assembly output
ASSEMBLY=${ASSEMBLY_DIR}/${SAMPLE}_hybrid/contigs_polished.fasta  # Assembly from 02_Assembly (hybrid output)

# Input files (should be user-provided)
FASTQ1_ORIG=${FASTQ_DIR}/${SHORT_SRR}_1.nohuman.fq.gz
FASTQ2_ORIG=${FASTQ_DIR}/${SHORT_SRR}_2.nohuman.fq.gz

LOGDIR=./03_Binning/logs  # Directory to store logs

# Create necessary directories
mkdir -p $FINALDIR
mkdir -p $LOGDIR

# Copy input files to the working directory
cp $FASTQ1_ORIG ./reads_1.fq.gz
cp $FASTQ2_ORIG ./reads_2.fq.gz
cp $ASSEMBLY ./assembly.fasta

# Unzip the fastq files
gunzip reads_1.fq.gz
gunzip reads_2.fq.gz

# Rename fastq files to .fastq
mv reads_1.fq reads_1.fastq
mv reads_2.fq reads_2.fastq

# Step 1: Binning
metawrap binning \
  -o binning \
  -t 40 -m 200 \
  -a assembly.fasta \
  --metabat2 --maxbin2 --concoct \
  reads_1.fastq reads_2.fastq \
  &> $LOGDIR/${SAMPLE}.binning.log

# Step 2: Refinement
metawrap bin_refinement \
  -o binning/refined \
  -t 40 -m 200 \
  -c 50 -x 10 \
  -A binning/metabat2_bins \
  -B binning/maxbin2_bins \
  -C binning/concoct_bins \
  &> $LOGDIR/${SAMPLE}.refine.log

# Clean up intermediate files
rm -f reads_1.fastq reads_2.fastq assembly.fasta

# Move results to the final output directory
mkdir -p $FINALDIR/${SAMPLE}
mv $WORKDIR/binning/refined/* $FINALDIR/${SAMPLE}.bins

# Clean up working directory
rm -rf binning

echo "[INFO] Binning and refinement for ${SAMPLE} completed and moved to ${FINALDIR}/${SAMPLE}.bins"
