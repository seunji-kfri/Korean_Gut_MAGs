#!/bin/bash
# Usage: run_metawrap.sh <sample_id> <short_srr>

. ~/enable_conda.rc
conda activate metawrap

J=$1
short_srr=$2


WORKDIR=/tmp/${J}.bins
FINALDIR=/home/caefs/microbiome/projects/metagenome_analysis/China/operams
FASTQ_DIR=/home/caefs/microbiome/projects/metagenome_analysis/China/preprocess
BASE=/home/caefs/microbiome/projects/metagenome_analysis/China/operams/${J}
FASTQ1_ORIG=${FASTQ_DIR}/${short_srr}_1.nohuman.fq.gz
FASTQ2_ORIG=${FASTQ_DIR}/${short_srr}_2.nohuman.fq.gz
ASSEMBLY_ORIG=${BASE}/contigs.polished.fasta


LOGDIR=/home/caefs/microbiome/projects/metagenome_analysis/China


mkdir -p $WORKDIR
cd $WORKDIR


cp $FASTQ1_ORIG ${WORKDIR}/reads_1.fq.gz
cp $FASTQ2_ORIG ${WORKDIR}/reads_2.fq.gz
cp $ASSEMBLY_ORIG ${WORKDIR}/assembly.fasta

gunzip reads_1.fq.gz
gunzip reads_2.fq.gz

mv reads_1.fq reads_1.fastq
mv reads_2.fq reads_2.fastq

# Step 1: Binning
metawrap binning \
  -o binning \
  -t 40 -m 200 \
  -a assembly.fasta \
  --metabat2 --maxbin2 --concoct \
  reads_1.fastq reads_2.fastq \
  &> ${LOGDIR}/${J}.binning.log

# Step 2: Refinement
metawrap bin_refinement \
  -o binning/refined \
  -t 40 -m 200 \
  -c 50 -x 10 \
  -A binning/metabat2_bins \
  -B binning/maxbin2_bins \
  -C binning/concoct_bins \
  &> ${LOGDIR}/${J}.refine.log

rm -f reads_1.fastq reads_2.fastq assembly.fasta

mv ${WORKDIR}/* ${FINALDIR}/${J}.bins
rm -rf ${WORKDIR}


echo "[INFO] Binning and refinement for ${J} completed and moved to ${FINALDIR}/${J}.bins"
