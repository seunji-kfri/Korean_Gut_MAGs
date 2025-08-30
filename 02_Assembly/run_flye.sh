#!/bin/bash

# Number of threads (adjust as necessary)
CPU=20

# Activate conda environment
. ~/enable_conda.rc
conda activate kfri_genome_pipeline

# Define temporary directories
T=/tmp/flye_${1##*/}
C=/tmp/${1##*/}_contigs.fasta

# Run Flye assembly (ONT long-read)
flye --meta --nano-raw /home/caefs/microbiome/projects/metagenome_analysis/KMAGs/assembly_long/reads/${1}.nohuman.fastq --out-dir $T --threads $CPU

# Move output to designated folder
mkdir -p $1 && mv -f /tmp/flye_${1##*/} $1

# Deactivate conda environment
conda deactivate
