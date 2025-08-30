#!/usr/bin/env python3
"""
Long-read preprocessing (ONT)

Steps:
  1) Adapter trimming (porechop)  [--no-porechop to skip]
  2) Length/quality filtering (filtlong)
  3) Remove human reads (two modes):
      - default: keep unmapped reads only (samtools -f 4)
      - paf_filter: use paf_parser.py with identity/coverage cutoffs
  4) Basic stats (seqkit)

Usage:
  Default (unmapped-only):
    python preprocess_longreads.py \
      --sample SAMPLE01 \
      --in_fastq SAMPLE01.fastq.gz \
      --outdir preprocess_out \
      --human_mmi /path/to/human.mmi

  With paf_filter (identity 80, coverage 30):
    python preprocess_longreads.py \
      --sample SAMPLE01 \
      --in_fastq SAMPLE01.fastq.gz \
      --outdir preprocess_out \
      --human_mmi /path/to/human.mmi \
      --paf_filter 80 30
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
    ap.add_argument("--paf_filter", nargs=2, metavar=("IDENT", "COV"), type=float,
                    help="Use paf_parser.py filtering with identity%% and coverage%% cutoffs")
    ap.add_argument("--paf_parser", default="paf_parser.py",
                    help="Path to paf_parser.py (default: in current dir)")
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / f"{args.sample}.preprocess.log"

    # 1) porechop (optional)
    porechop_out = outdir / f"{args.sample}.porechop.fastq.gz"
    if args.no_porechop:
        if not porechop_out.exists():
            run(f"ln -s {args.in_fastq} {porechop_out}", log=str(log_file))
    else:
        if not porechop_out.exists():
            run(
                f"porechop -i {args.in_fastq} -o {porechop_out} -t {args.threads} "
                f"2> {outdir}/{args.sample}.porechop.log",
                log=str(log_file),
            )

    # 2) filtlong
    filtlong_out = outdir / f"{args.sample}.filtlong.fastq.gz"
    if not filtlong_out.exists():
        run(
            f"filtlong --min_length {args.min_length} --keep_percent {args.keep_percent} "
            f"{porechop_out} | gzip > {filtlong_out} 2> {outdir}/{args.sample}.filtlong.log",
            log=str(log_file),
        )

    # 3) remove human
    bam_out = outdir / f"{args.sample}.human.bam"
    if not bam_out.exists():
        run(
            f"minimap2 -t {args.threads} -x map-ont {args.human_mmi} {filtlong_out} "
            f"| samtools view -b -o {bam_out}",
            log=str(log_file),
        )

    nohuman_fastq = outdir / f"{args.sample}.nohuman.fastq.gz"

    if args.paf_filter:
        ident, cov = args.paf_filter
        paf_out = outdir / f"{args.sample}.human.paf"
        paf_filtered = outdir / f"{args.sample}.human.filtered.paf"

        # Generate PAF file
        if not paf_out.exists():
            run(f"minimap2 -t {args.threads} -c -x map-ont {args.human_mmi} {filtlong_out} > {paf_out}",
                log=str(log_file))

        # Filter with paf_parser.py
        run(f"python {args.paf_parser} {paf_out} -i {ident} -q {cov} > {paf_filtered}",
            log=str(log_file))

        # Extract read IDs
        ids_file = outdir / f"{args.sample}.human.ids"
        run(f"cut -f1 {paf_filtered} | uniq > {ids_file}", log=str(log_file))

        # Remove human reads by ID
        run(f"seqkit grep -v -f {ids_file} {filtlong_out} | gzip > {nohuman_fastq}",
            log=str(log_file))

    else:
        # default: keep unmapped only
        if not nohuman_fastq.exists():
            run(f"samtools fastq -f 4 {bam_out} | gzip > {nohuman_fastq}", log=str(log_file))

    # 4) stats
    stat_out = outdir / f"{args.sample}.nohuman.seqstat.txt"
    if not stat_out.exists():
        run(f"seqkit stats -a -T {nohuman_fastq} > {stat_out}", log=str(log_file))

    print("[DONE] Long-read preprocessing complete.")

if __name__ == "__main__":
    main()
