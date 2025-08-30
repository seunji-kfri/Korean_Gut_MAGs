#!/bin/bash
# Usage: run_flye.sh <sample_id> <output_folder>
# <sample_id> : The ID of the sample (e.g., SAMPLE01)
# <output_folder> : The directory where Flye results will be saved (e.g., 02_Assembly/out/SAMPLE01)

# Number of threads (adjust as necessary)
CPU=20

# Activate conda environment (ensure the environment created via 'environment.yml' is available)
conda activate KoreanGutMAGs  # Ensure this is the environment created with 'environment.yml'

# Define output directory
OUTPUT_DIR=$2   # Use the second argument as the output directory

# Run Flye assembly (ONT long-read)
flye --meta --nano-raw /path/to/your/input/${1}.nohuman.fastq --out-dir $OUTPUT_DIR --threads $CPU

# Deactivate conda environment
conda deactivate
