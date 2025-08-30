#!/usr/bin/env python3
import os
import sys
import shutil
import subprocess
import os.path as path

# Arguments
if len(sys.argv) < 2:
    print("Usage: python preprocess_shortreads.py <sample_id>")
    sys.exit(1)

sample_id = sys.argv[1].strip()

# Paths
rawdata_dir = '/home/caefs/microbiome/projects/metagenome_analysis/China/rawread'
shortreads_dir = f'{rawdata_dir}'
shortreads_1 = f'{shortreads_dir}/{sample_id}_1.fastq'
shortreads_2 = f'{shortreads_dir}/{sample_id}_2.fastq'

workdir = '/home/caefs/microbiome/projects/metagenome_analysis/China/preprocess'
os.makedirs(workdir, exist_ok=True)

trimmomatic_path = '/home/caefs/microbiome/jupyconda/envs/metagenome_preprocess/share/trimmomatic-0.39-2'
human_bowtie_index = '/home/caefs/microbiome/projects/metagenome_analysis/China/resource/human.bowtie2/human'

# Output files
trimmed_1 = f'{workdir}/{sample_id}_1.fq.gz'
trimmed_2 = f'{workdir}/{sample_id}_2.fq.gz'
unpair_1 = f'{workdir}/{sample_id}_unpair1.fq.gz'
unpair_2 = f'{workdir}/{sample_id}_unpair2.fq.gz'
nohuman_prefix = f'{workdir}/{sample_id}_nohuman'
nohuman_1 = f'{workdir}/{sample_id}_1.nohuman.fq.gz'
nohuman_2 = f'{workdir}/{sample_id}_2.nohuman.fq.gz'

# Step 1: Trimming
if not (os.path.exists(trimmed_1) and os.path.exists(trimmed_2)):
    cmd = f"trimmomatic PE -threads 8 -phred33 {shortreads_1} {shortreads_2} " \
          f"{trimmed_1} {unpair_1} {trimmed_2} {unpair_2} " \
          f"ILLUMINACLIP:{trimmomatic_path}/adapters/TruSeq3-PE.fa:2:30:10 " \
          f"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

# Step 2: Remove human reads
if not (os.path.exists(nohuman_1) and os.path.exists(nohuman_2)):
    bam_out = f'{workdir}/{sample_id}.human.bam'
    cmd = f"bowtie2 -p 20 -x {human_bowtie_index} -1 {trimmed_1} -2 {trimmed_2} " \
          f"--un-conc-gz {nohuman_prefix} | samtools view -b -o {bam_out}"
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

    # Rename paired output
    if os.path.exists(f'{nohuman_prefix}.1'):
        shutil.move(f'{nohuman_prefix}.1', nohuman_1)
    if os.path.exists(f'{nohuman_prefix}.2'):
        shutil.move(f'{nohuman_prefix}.2', nohuman_2)

# Step 3: Stats for filtered reads
for fq, outstat in [(nohuman_1, f"{workdir}/{sample_id}_1.nohuman.seqstat.txt"),
                    (nohuman_2, f"{workdir}/{sample_id}_2.nohuman.seqstat.txt")]:
    if not os.path.exists(outstat):
        cmd = f"seqkit stats -a -T {fq} > {outstat}"
        print(f"Running: {cmd}")
        subprocess.run(cmd, shell=True, check=True)

print(f"[DONE] Preprocessing complete for {sample_id}")
