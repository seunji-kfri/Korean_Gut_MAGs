#!/usr/bin/env python3
# Long-read preprocessing (ONT)
# Adapter trimming, Quality filtering, Human read removal, Stats generation

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
    ap.add_argument("--in_fastq", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--human_mmi", required=True)
    ap.add_argument("--threads", type=int, default=20)
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / f"{args.sample}.preprocess.log"

    # Adapter trimming (optional)
    porechop_out = outdir / f"{args.sample}.porechop.fastq.gz"
    if not porechop_out.exists():
        run(f"porechop -i {args.in_fastq} -o {porechop_out} -t {args.threads}", log=str(log_file))

    # Quality filtering (filtlong)
    filtlong_out = outdir / f"{args.sample}.filtlong.fastq.gz"
    run(f"filtlong --min_length 1000 --keep_percent 95 {porechop_out} | gzip > {filtlong_out}", log=str(log_file))

    # Human read removal (Minimap2)
    bam_out = outdir / f"{args.sample}.human.bam"
    run(f"minimap2 -t {args.threads} -x map-ont {args.human_mmi} {filtlong_out} | samtools view -b -o {bam_out}", log=str(log_file))

if __name__ == "__main__":
    main()
