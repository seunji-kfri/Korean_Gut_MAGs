#!/usr/bin/env python3
import os
import sys
from os import path


if len(sys.argv) != 2:
    print("Usage: python preprocess_longreads.py <sample_id>")
    sys.exit(1)

sample_id = sys.argv[1]


rawdata_dir = '/home/caefs/microbiome/projects/metagenome_analysis/China/rawread'
longreads_dir = f'{rawdata_dir}'
longreads_fastq = f'{longreads_dir}/{sample_id}.fastq'


workdir = '/home/caefs/microbiome/projects/metagenome_analysis/China/preprocess'
os.makedirs(workdir, exist_ok=True)

# --- Step 1: Porechop ---
porechop_fastq = f'{workdir}/{sample_id}.porechop.fastq'
if not path.exists(porechop_fastq):
    cmd = f"porechop -i {longreads_fastq} -o {porechop_fastq} -t 20 &> {workdir}/{sample_id}.porechop.log"
    print(f"[RUN] {cmd}")
    os.system(cmd)
else:
    print(f"[SKIP] Porechop output exists: {porechop_fastq}")

# --- Step 2: Filtlong ---
filtlong_fastq = f'{workdir}/{sample_id}.filtlong.fastq'
if not path.exists(filtlong_fastq):
    cmd = f"filtlong --min_length 1000 --keep_percent 95 {porechop_fastq} > {filtlong_fastq} 2> {workdir}/{sample_id}.filtlong.log"
    print(f"[RUN] {cmd}")
    os.system(cmd)
else:
    print(f"[SKIP] Filtlong output exists: {filtlong_fastq}")

# --- Step 3: Human contamination check ---
human_mmi = '/home/caefs/microbiome/projects/metagenome_analysis/China/resource/human.ont.mmi'
human_paf = f'{workdir}/{sample_id}.human.paf'
if not path.exists(human_paf):
    cmd = f"minimap2 -t 26 -c -x map-ont {human_mmi} {filtlong_fastq} > {human_paf}"
    print(f"[RUN] {cmd}")
    os.system(cmd)
else:
    print(f"[SKIP] Human PAF exists: {human_paf}")

# --- Step 4: Parse PAF ---
progdir = path.abspath(os.getcwd())
human_paf_x = f'{human_paf}x'
if not path.exists(human_paf_x):
    cmd = f"{progdir}/paf_parser.py {human_paf} > {human_paf_x}"
    print(f"[RUN] {cmd}")
    os.system(cmd)
else:
    print(f"[SKIP] Parsed PAF exists: {human_paf_x}")

# --- Step 5: Filter by identity (80%) and coverage (30%) ---
human_filtered = f"{workdir}/{sample_id}.humanI80C30"
pafx_file = f"{human_filtered}.pafx"
if not path.exists(pafx_file):
    cmd = f"awk -F'\t' '$13 > 80 && $15 > 30' {human_paf_x} > {pafx_file}"
    print(f"[RUN] {cmd}")
    os.system(cmd)
else:
    print(f"[SKIP] Filtered PAF exists: {pafx_file}")

# --- Step 6: Extract human read IDs ---
ids_file = f"{human_filtered}.ids"
if not path.exists(ids_file):
    cmd = f"cut -f1 {pafx_file} | uniq > {ids_file}"
    print(f"[RUN] {cmd}")
    os.system(cmd)
else:
    print(f"[SKIP] ID list exists: {ids_file}")

# --- Step 7: Remove human reads ---
nohuman_fastq = f"{workdir}/{sample_id}.nohuman.fastq"
if not path.exists(nohuman_fastq):
    cmd = f"seqkit grep -v -f {ids_file} {filtlong_fastq} > {nohuman_fastq}"
    print(f"[RUN] {cmd}")
    os.system(cmd)
else:
    print(f"[SKIP] No-human FASTQ exists: {nohuman_fastq}")

# --- Step 8: Sequence statistics ---
seqstat_file = f"{workdir}/{sample_id}.nohuman.seqstat.txt"
if not path.exists(seqstat_file):
    cmd = f"seqkit stats -a -T {nohuman_fastq} > {seqstat_file}"
    print(f"[RUN] {cmd}")
    os.system(cmd)
else:
    print(f"[SKIP] Seqkit stats exists: {seqstat_file}")

print("Pipeline completed successfully.")
