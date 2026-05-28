#!/bin/bash
# =============================================================================
# run_gtdbtk.sh
#
# Taxonomic classification of genome bins using GTDB-Tk (classify_wf),
# as used in the manuscript (GTDB-Tk with the GTDB r226 database).
#
# Usage:
#   export GTDBTK_DATA_PATH=/path/to/gtdbtk/db    # GTDB release database
#   bash 05_Taxonomy/run_gtdbtk.sh <genome_dir> <out_dir> [prefix]
#
# Example:
#   export GTDBTK_DATA_PATH=/home/.../gtdbtk-2.4.1/db/
#   bash 05_Taxonomy/run_gtdbtk.sh \
#     /path/to/bins 05_Taxonomy/gtdbtk_release226_output allbins
#
# Input:
#   <genome_dir> : directory of binned genomes (FASTA, .fna)
# Output:
#   <out_dir>    : GTDB-Tk classification results
#                  (<prefix>.bac120.summary.tsv, etc.)
#
# Requires GTDBTK_DATA_PATH to point to the GTDB release database.
# =============================================================================

set -euo pipefail

GENOME_DIR="${1:?Usage: run_gtdbtk.sh <genome_dir> <out_dir> [prefix]}"
OUT_DIR="${2:?Usage: run_gtdbtk.sh <genome_dir> <out_dir> [prefix]}"
PREFIX="${3:-allbins}"

if [ -z "${GTDBTK_DATA_PATH:-}" ]; then
  echo "ERROR: set GTDBTK_DATA_PATH to the GTDB database directory before running." >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

gtdbtk classify_wf \
  --genome_dir "${GENOME_DIR}" \
  --out_dir "${OUT_DIR}" \
  -x fna \
  --cpus 40 \
  --pplacer_cpus 10 \
  --keep_intermediates \
  --write_single_copy_genes \
  --skip_ani_screen \
  --prefix "${PREFIX}"

echo "GTDB-Tk classification complete. Results in: ${OUT_DIR}"
