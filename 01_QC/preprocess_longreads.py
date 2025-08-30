#!/usr/bin/env python3
"""
Long-read preprocessing (ONT)

Steps:
  1) Adapter trimming (porechop)  [--no-porechop to skip]
  2) Length/quality filtering (filtlong)
  3) Remove human reads (minimap2 vs human.mmi) -> keep UNMAPPED reads (samtools fastq -f 4)
  4) Basic stats (seqkit)

Usage:
  python preprocess_longreads.py \
    --sample SAMPLE01 \
    --in_fastq /path/to/SAMPLE01.ont.fastq.gz \
    --outdir ./preprocess_out \
    --human_mmi /path/to/human.mmi \
    --threads 20 \
    --min_length 1000 \
    --keep_percent 95

Requirements:
  - porechop (optional), filtlong, minimap2, samtools, seqkit
  - (optional) pigz

Notes:
  - We avoid PAF post-processing by directly extracting unmapped reads from BAM.
  - If you prefer identity/coverage thresholds, we can add an optional PAF-based mode.
"""

import argparse
import subprocess
from pathlib import Path

def run(cmd, log=None):
    print(f"[RUN] {cmd}")
    res = subprocess.run(cmd, shell=True)
    if res.returncode != 0:
        raise RuntimeError(f"Command failed (exit {res.returncode}): {cmd}")
    if log:
        with open(log, "a") as f:
            f.write(cmd + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, help="Sample ID")
    ap.add_argument("--in_fastq", required=True, help="Input ONT FASTQ(.gz)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--human_mmi", required=True, help="Minimap2 human index (.mmi)")
    ap.add_argument("--threads", type=int, default=20)
    ap.add_argument("--min_length", type=int, default=1000, help="filtlong --min_length")
    ap.add_argument("--keep_percent", type=int, default=95, help="filtlong --keep_percent")
    ap.add_argument("--no-porechop", action="store_true", help="Skip porechop")
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / f"{args.sample}.preprocess.log"

    # 1) porechop (optional)
    porechop_out = outdir / f"{args.sample}.porechop.fastq.gz"
    if args.no_porechop:
        # If skipping, just symlink/copy input to porechop_out-like name
        if not porechop_out.exists():
            run(f"ln -s {args.in_fastq} {porechop_out}", log=str(log_file))
    else:
        if not porechop_out.exists():
            run(
                f"porechop -i {args.in_fastq} -o {porechop_out} -t {args.threads} 2> {outdir}/{args.sample}.porechop.log",
                log=str(log_file),
            )
        else:
            print(f"[SKIP] porechop exists: {porechop_out.name}")

    # 2) filtlong
    filtlong_out = outdir / f"{args.sample}.filtlong.fastq.gz"
    if not filtlong_out.exists():
        run(
            f"filtlong --min_length {args.min_length} --keep_percent {args.keep_percent} "
            f"{porechop_out} | gzip > {filtlong_out} 2> {outdir}/{args.sample}.filtlong.log",
            log=str(log_file),
        )
    else:
        print(f"[SKIP] filtlong exists: {filtlong_out.name}")

    # 3) remove human (minimap2 ¡æ BAM ¡æ unmapped FASTQ)
    bam_out = outdir / f"{args.sample}.human.bam"
    if not bam_out.exists():
        run(
            f"minimap2 -t {args.threads} -x map-ont {args.human_mmi} {filtlong_out} "
            f"| samtools view -b -o {bam_out}",
            log=str(log_file),
        )
    else:
        print(f"[SKIP] BAM exists: {bam_out.name}")

    nohuman_fastq = outdir / f"{args.sample}.nohuman.fastq.gz"
    if not nohuman_fastq.exists():
        # -f 4: unmapped
        run(
            f"samtools fastq -f 4 {bam_out} | gzip > {nohuman_fastq}",
            log=str(log_file),
        )
    else:
        print(f"[SKIP] nohuman FASTQ exists: {nohuman_fastq.name}")

    # 4) stats
    stat_out = outdir / f"{args.sample}.nohuman.seqstat.txt"
    if not stat_out.exists():
        run(f"seqkit stats -a -T {nohuman_fastq} > {stat_out}", log=str(log_file))
    else:
        print(f"[SKIP] stats exists: {stat_out.name}")

    print("[DONE] Long-read preprocessing complete.")

if __name__ == "__main__":
    main()
