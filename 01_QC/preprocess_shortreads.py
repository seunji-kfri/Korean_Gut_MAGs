#!/usr/bin/env python3
# Short-read preprocessing (Illumina PE)
# Trimming, Human read removal, Stats generation
# Usage: preprocess_shortreads.py --sample <sample_id> --r1 <R1_fastq> --r2 <R2_fastq> --outdir <output_dir> --trimmomatic_adapters <adapters_file> --human_bowtie2_index <index_path>

import argparse
import subprocess
from pathlib import Path

def run(cmd, log=None):
    res = subprocess.run(cmd, shell=True)
    if res.returncode != 0:
        raise RuntimeError(f"Command failed (exit {res.returncode}): {cmd}")
    if log:
        with open(log, "a") as f:
            f.write(cmd + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--trimmomatic_adapters", required=True)
    ap.add_argument("--human_bowtie2_index", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / f"{args.sample}.preprocess.log"

    # Trimming step
    trimmed_1 = outdir / f"{args.sample}_1.trim.fq.gz"
    trimmed_2 = outdir / f"{args.sample}_2.trim.fq.gz"
    if not (trimmed_1.exists() and trimmed_2.exists()):
        cmd_trim = f"trimmomatic PE -threads {args.threads} {args.r1} {args.r2} {trimmed_1} {trimmed_2} ILLUMINACLIP:{args.trimmomatic_adapters}:2:30:10"
        run(cmd_trim, log=str(log_file))

    # Human read removal (Bowtie2)
    nohuman_1 = outdir / f"{args.sample}_1.nohuman.fq.gz"
    nohuman_2 = outdir / f"{args.sample}_2.nohuman.fq.gz"
    if not (nohuman_1.exists() and nohuman_2.exists()):
        cmd_bt2 = f"bowtie2 -p {args.threads} -x {args.human_bowtie2_index} -1 {trimmed_1} -2 {trimmed_2} --un-conc-gz {outdir}/{args.sample}_nohuman"
        run(cmd_bt2, log=str(log_file))

    # Stats generation
    stat_out = outdir / f"{args.sample}.nohuman.seqstat.txt"
    if not stat_out.exists():
        run(f"seqkit stats -a -T {nohuman_1} > {stat_out}", log=str(log_file))

if __name__ == "__main__":
    main()
