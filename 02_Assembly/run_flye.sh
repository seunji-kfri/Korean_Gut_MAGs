#!/bin/bash

# Number of threads (adjust as necessary)
CPU=20

# Activate conda environment (ensure the environment created via 'environment.yml' is available)
conda activate KoreanGutMAGs  # Ensure this is the environment created with 'environment.yml'

# Define output directory
OUTPUT_DIR=$1   # Use the sample ID as the output directory

# Run Flye assembly (ONT long-read)
flye --meta --nano-raw /path/to/your/input/${1}.nohuman.fastq --out-dir $OUTPUT_DIR --threads $CPU

# Deactivate conda environment
conda deactivate
