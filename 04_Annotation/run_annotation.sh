#!/bin/bash

if [[ $# -lt 1 ]]; then
  echo "Usage: bash $0 <sample_id>" >&2
  exit 1
fi

SAMPLE="$1"
THREADS="${THREADS:-20}"     # per-bin thread count
NPROC="${NPROC:-4}"          # number of bins to run in parallel

BINS_DIR="/home/caefs/microbiome/projects/metagenome_analysis/China/operams/${SAMPLE}.bins/refined/metawrap_50_10_bins/"
ANNOT_BASE="/home/caefs/microbiome/projects/metagenome_analysis/China/annotation"
ANNOT_DIR="${ANNOT_BASE}/${SAMPLE}"

mkdir -p "$ANNOT_DIR"
shopt -s nullglob
bins=( "${BINS_DIR}"/bin.*.fa )
shopt -u nullglob

if [[ ${#bins[@]} -eq 0 ]]; then
  echo "[WARN] No bins found in ${BINS_DIR}"
  exit 0
fi

echo "[INFO] Sample: ${SAMPLE}, Bins: ${#bins[@]}, Threads/bin: ${THREADS}, Parallel bins: ${NPROC}"

# 함수 정의: 개별 bin 처리
run_bin() {
  fasta="$1"
  SAMPLE="$2"
  THREADS="$3"
  ANNOT_DIR="$4"

  bn="$(basename "$fasta")"
  cluster_index="${bn##bin.}"
  cluster_index="${cluster_index%.fa}"
  bin_id="${SAMPLE}C${cluster_index}"

  prokka_dir="${ANNOT_DIR}/${bin_id}.prokka"
  prokka_log="${ANNOT_DIR}/${bin_id}.prokka.log"
  eggnog_dir="${ANNOT_DIR}/${bin_id}.em"
  eggnog_log="${ANNOT_DIR}/${bin_id}.em.log"

  # Prokka
  prokka_dir="${ANNOT_DIR}/${bin_id}.prokka"
  prokka_log="${ANNOT_DIR}/${bin_id}.prokka.log"

  echo "[RUN ][PROKKA] $bin_id -> $prokka_dir"
  prokka --cpus "$THREADS" \
         --outdir "$prokka_dir" \
         --prefix "$bin_id" \
         --locustag "$bin_id" \
         --metagenome \
         "$fasta" &> "$prokka_log"
  # eggNOG
  faa="${prokka_dir}/${bin_id}.faa"
  if [[ -s "$faa" ]]; then
    if [[ -d "$eggnog_dir" && -s "${eggnog_dir}/${bin_id}.annotations" ]]; then
      echo "[SKIP][EMAPPER] $bin_id"
    else
      mkdir -p "$eggnog_dir"
      emapper.py --cpu "$THREADS" --dbmem \
                 -i "$faa" \
                 --output_dir "$eggnog_dir" \
                 --output "$bin_id" &> "$eggnog_log"
    fi
  else
    echo "[WARN] $bin_id skipped (no faa)"
  fi
}

export -f run_bin
export SAMPLE THREADS ANNOT_DIR

printf "%s\n" "${bins[@]}" | xargs -n1 -P "$NPROC" -I{} bash -c 'run_bin "$@"' _ {} "$SAMPLE" "$THREADS" "$ANNOT_DIR"

echo "[INFO] Done: $SAMPLE"
