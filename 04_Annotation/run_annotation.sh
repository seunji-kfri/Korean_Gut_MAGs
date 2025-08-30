#!/bin/bash
# Usage: run_annotation.sh <sample_id> <short_srr_1> <output_directory>
# <sample_id> : The ID of the sample (e.g., SAMPLE01)
# <short_srr_1> : The prefix for short-read files (e.g., SAMPLE01_1.nohuman.fq.gz)
# <output_directory> : The directory where the annotation results will be saved (e.g., annotation_results)

# Activate conda environment for annotation (using the correct environment name)
conda activate KoreanGutMAGs  # Ensure this is the environment created with 'environment.yml'

# Directories
BINNING_DIR=./03_Binning/out  # Path to binning results
ANNOTATION_DIR=$3  # Output directory for annotations

# Input files (e.g., binning results)
BINS=$BINNING_DIR/$1/$1.bins/refined/metawrap_50_10_bins/bin.*.fa  # Binning files for the sample

# Create output directory if it doesn't exist
mkdir -p $ANNOTATION_DIR

# Loop through each bin file for annotation
for bin_file in $BINS; do
    # Extract bin ID from the file name
    bin_id=$(basename $bin_file | cut -d'.' -f2)

    # Create output directories for Prokka and EggNOG
    bin_output_dir=$ANNOTATION_DIR/$1/bin_$bin_id
    mkdir -p $bin_output_dir

    # Run Prokka annotation
    prokka --cpus 20 --outdir $bin_output_dir/prokka --prefix $bin_id --locustag $bin_id $bin_file

    # Run EggNOG annotation
    emapper.py --cpu 20 --dbmem -i $bin_output_dir/prokka/$bin_id.faa --output $bin_output_dir/eggnog --temp_dir /tmp

done

# Deactivate conda environment
conda deactivate
