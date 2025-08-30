#!/bin/bash

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
