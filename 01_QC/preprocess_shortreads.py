#!/usr/bin/env python3
"""
Short-read preprocessing (Illumina PE)

Steps:
  1) Adapter/quality trimming (Trimmomatic)
  2) Remove human reads (Bowtie2 vs human reference) -> keep UNMAPPED pairs
  3) Basic stats (seqkit)

Usage:
  python preprocess_shortreads.py \
    --sample SAMPLE01 \
    --r1 /path/to/SAMPLE01_R1.fastq.gz \
    --r2 /path/to/SAMPLE01_R2.fastq.gz \
    --outdir ./preprocess_out \
    --trimmomatic_adapters /path/to/TruSeq3-PE.fa \
    --human_bowtie2_index /path/to/human_bowtie2_index_prefix \
    --threads 16

Requirements:
  - trimmomatic, bowtie2, samtools, seqkit
  - (optional) pigz

Notes:
  - Bowtie2 produces unmapped paired FASTQs via --un-conc-gz <prefix>
  - We then rename them to <sample>_1.nohuman.fq.gz and <sample>_2.nohuman.fq.gz
"""

import argparse
import subprocess
import os
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
    ap.add_argument("--r1", required=True, help="Input R1 FASTQ(.gz)")
    ap.add_argument("--r2", required=True, help="Input R2 FASTQ(.gz)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--trimmomatic_path", default="trimmomatic", help="trimmomatic executable (default: trimmomatic)")
    ap.add_argument("--trimmomatic_adapters", required=True, help="Adapter fasta (e.g., TruSeq3-PE.fa)")
    ap.add_argument("--human_bowtie2_index", required=True, help="Bowtie2 index prefix for human (e.g., /db/human/hg38)")
    ap.add_argument("--threads", type=int, default=8)
    # Trimmomatic knobs
    ap.add_argument("--minlen", type=int, default=50)
    ap.add_argument("--leading", type=int, default=3)
    ap.add_argument("--trailing", type=int, default=3)
    ap.add_argument("--window", type=int, default=4)
    ap.add_argument("--qual", type=int, default=15)
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / f"{args.sample}.preprocess.log"

    # 1) Trimming
    trimmed_1 = outdir / f"{args.sample}_1.trim.fq.gz"
    trimmed_2 = outdir / f"{args.sample}_2.trim.fq.gz"
    unpair_1  = outdir / f"{args.sample}_unpair1.fq.gz"
    unpair_2  = outdir / f"{args.sample}_unpair2.fq.gz"

    if not (trimmed_1.exists() and trimmed_2.exists()):
        cmd_trim = (
            f"{args.trimmomatic_path} PE -threads {args.threads} -phred33 "
            f"{args.r1} {args.r2} "
            f"{trimmed_1} {unpair_1} {trimmed_2} {unpair_2} "
            f"ILLUMINACLIP:{args.trimmomatic_adapters}:2:30:10 "
            f"LEADING:{args.leading} TRAILING:{args.trailing} "
            f"SLIDINGWINDOW:{args.window}:{args.qual} MINLEN:{args.minlen}"
        )
        run(cmd_trim, log=str(log_file))
    else:
        print(f"[SKIP] Trimming exists: {trimmed_1.name}, {trimmed_2.name}")

    # 2) Remove human reads (keep UNMAPPED pairs)
    # Bowtie2 will write unmapped pairs to <prefix>.1 and <prefix>.2
    nohuman_prefix = outdir / f"{args.sample}_nohuman"
    nohuman_1 = outdir / f"{args.sample}_1.nohuman.fq.gz"
    nohuman_2 = outdir / f"{args.sample}_2.nohuman.fq.gz"
    bam_out = outdir / f"{args.sample}.human.bam"

    if not (nohuman_1.exists() and nohuman_2.exists()):
        cmd_bt2 = (
            f"bowtie2 -p {args.threads} "
            f"-x {args.human_bowtie2_index} "
            f"-1 {trimmed_1} -2 {trimmed_2} "
            f"--very-sensitive "
            f"--un-conc-gz {nohuman_prefix} "
            f"| samtools view -b -o {bam_out}"
        )
        run(cmd_bt2, log=str(log_file))

        # Rename unmapped outputs
        cand1 = str(nohuman_prefix) + ".1"
        cand2 = str(nohuman_prefix) + ".2"
        if os.path.exists(cand1):
            os.replace(cand1, nohuman_1)
        if os.path.exists(cand2):
            os.replace(cand2, nohuman_2)
    else:
        print(f"[SKIP] No-human outputs exist: {nohuman_1.name}, {nohuman_2.name}")

    # 3) Stats
    for fq, tag in [(nohuman_1, "R1"), (nohuman_2, "R2")]:
        stat_out = outdir / f"{args.sample}_{tag}.nohuman.seqstat.txt"
        if not stat_out.exists():
            run(f"seqkit stats -a -T {fq} > {stat_out}", log=str(log_file))
        else:
            print(f"[SKIP] Stats exist: {stat_out.name}")

    print("[DONE] Short-read preprocessing complete.")

if __name__ == "__main__":
    main()
