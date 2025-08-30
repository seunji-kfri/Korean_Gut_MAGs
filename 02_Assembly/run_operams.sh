#!/bin/bash
# Usage: run_operams.sh <long_read_file> <short_srr_1> <short_srr_2> <output_directory>
# <long_read_file> : The long-read file (e.g., SAMPLE01.nohuman.fastq)
# <short_srr_1> : The prefix for short-read file 1 (e.g., SAMPLE01_1.nohuman.fq.gz)
# <short_srr_2> : The prefix for short-read file 2 (e.g., SAMPLE01_2.nohuman.fq.gz)
# <output_directory> : The directory to store the output (e.g., 02_Assembly/out/SAMPLE01_hybrid)

# Number of threads (adjust as necessary)
CPU=40

# Activate conda environment (ensure the environment created with 'environment.yml' is available)
conda activate KoreanGutMAGs  # Ensure this is the environment created with 'environment.yml'

# Define output directory
OUTPUT_DIR=$4   # The final output directory (user-defined)

# Run OPERA-MS assembly (SPAdes + long-read)
operams --short-read-assembler spades --long-read $1 --short-read1 $2 --short-read2 $3 --out-dir $OUTPUT_DIR --num-processors $CPU &> ${OUTPUT_DIR##*/}.log

# Deactivate conda environment
conda deactivate
