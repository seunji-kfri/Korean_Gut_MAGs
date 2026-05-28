# =============================================================================
# 01_build_linked_tables.R
#
# Build sample x KO tables for the link-based (reference-matching) baseline,
# across six ASV-MAG linking configurations that vary alignment stringency
# (PID), best-hit pruning, and KO weighting (binary vs prevalence).
# This baseline is compared against the ML framework in the manuscript
# (Fig. S2; benchmarked in 07_ML_Prediction/03_benchmark.R).
#
# Inputs (in 06_Linking/input/; produced by build_asv_features.py and the
# ASV-MAG BLAST step):
#   00_hits.rds          ASV-MAG 16S hits (ASV, BIN, pid, qcov, ...)
#   00_bin_ko.rds        per-bin KO content (BIN_ID, KO_ID)
#   00_asv_abund_long.rds  long ASV abundance (Sample, ASV, abundance)
#   00_run_manifest.rds  config manifest (run_id, pid_cut, pruning, weighting)
#
# Run from the repository root:
#   Rscript 06_Linking/01_build_linked_tables.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Paths (relative to repository root) ------------------------------------
ROOT_DIR  <- Sys.getenv("KGM_ROOT", unset = getwd())
INPUT_DIR <- file.path(ROOT_DIR, "06_Linking", "input")
OUT_DIR   <- file.path(ROOT_DIR, "06_Linking", "output", "01_linked_tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_CONFIGS <- c(
  "pid99_d02_bin",  "pid99_d02_prev",
  "pid99_tau1_bin", "pid99_tau1_prev",
  "pid100_tau1_bin", "pid100_tau1_prev"
)

qcov_cut      <- 0.90
delta_default <- 0.2
tau_default   <- 1.0

cat("Loading inputs...\n")
hits         <- readRDS(file.path(INPUT_DIR, "00_hits.rds"))
bin_ko       <- readRDS(file.path(INPUT_DIR, "00_bin_ko.rds"))
asv_abund    <- readRDS(file.path(INPUT_DIR, "00_asv_abund_long.rds"))
run_manifest <- readRDS(file.path(INPUT_DIR, "00_run_manifest.rds"))
run_manifest <- run_manifest[run_id %in% TARGET_CONFIGS]

setnames(bin_ko, c("BIN_ID", "KO_ID"))
bin_ko2 <- copy(bin_ko); setnames(bin_ko2, "BIN_ID", "BIN")

# ---- Pruning / aggregation helpers ------------------------------------------
filter_by_pid_qcov <- function(dt, pid_cut, qcov_cut) dt[pid >= pid_cut & qcov >= qcov_cut]

prune_top1 <- function(dt) {
  setorderv(dt, c("ASV", "pid", "qcov", "bits", "alnlen"), order = c(1, -1, -1, -1, -1))
  dt[, rn := seq_len(.N), by = ASV]
  out <- dt[rn == 1]; out[, rn := NULL]; out[]
}
prune_best_delta <- function(dt, delta) {
  best <- dt[, .(best_pid = max(pid, na.rm = TRUE)), by = ASV]
  dt <- merge(dt, best, by = "ASV")
  out <- dt[(best_pid - pid) <= delta]; out[, best_pid := NULL]; out[]
}
prune_soft <- function(dt, tau) {
  best <- dt[, .(best_pid = max(pid, na.rm = TRUE), best_qcov = max(qcov, na.rm = TRUE)), by = ASV]
  dt <- merge(dt, best, by = "ASV")
  dt[, `:=`(pid_gap = best_pid - pid, qcov_gap = best_qcov - qcov)]
  dt[, soft_score := pid_gap + 100 * qcov_gap]
  out <- dt[soft_score <= tau]
  out[, c("best_pid", "best_qcov", "pid_gap", "qcov_gap", "soft_score") := NULL]; out[]
}
prune_carriers <- function(dt, mode, delta, tau) {
  if (nrow(dt) == 0) return(copy(dt))
  switch(mode,
         top1 = prune_top1(copy(dt)),
         best_delta = prune_best_delta(copy(dt), delta),
         soft = prune_soft(copy(dt), tau),
         stop("Unknown pruning mode"))
}

build_asv_ko <- function(map_dt, bin_ko_dt, weighting) {
  map_unique <- unique(map_dt[, .(ASV, BIN)])
  joined <- merge(map_unique, bin_ko_dt, by = "BIN", allow.cartesian = TRUE)
  joined <- unique(joined, by = c("ASV", "BIN", "KO_ID"))
  if (weighting == "binary") {
    out <- unique(joined[, .(ASV, KO_ID)]); out[, weight := 1]; return(out[])
  }
  if (weighting == "prevalence") {
    carrier_n <- map_unique[, .N, by = ASV]; setnames(carrier_n, "N", "n_bin")
    out <- joined[, .N, by = .(ASV, KO_ID)]
    out <- merge(out, carrier_n, by = "ASV")
    out[, weight := N / n_bin]; out[, c("N", "n_bin") := NULL]; return(out[])
  }
  stop("Unknown weighting")
}
build_sample_ko <- function(asv_abund_dt, asv_ko_dt) {
  if (nrow(asv_ko_dt) == 0)
    return(data.table(Sample = character(), KO_ID = character(), value = numeric()))
  x <- merge(asv_abund_dt, asv_ko_dt, by = "ASV", allow.cartesian = TRUE)
  x[, contrib := abundance * weight]
  x[, .(value = sum(contrib)), by = .(Sample, KO_ID)][]
}
normalize_sample <- function(dt) {
  if (nrow(dt) == 0) return(copy(dt))
  totals <- dt[, .(sample_sum = sum(value)), by = Sample]
  out <- merge(dt, totals, by = "Sample")
  out[, value := fifelse(sample_sum > 0, value / sample_sum, 0)]
  out[, sample_sum := NULL]; out[]
}
to_wide <- function(dt, row_col, col_col, value_col) {
  if (nrow(dt) == 0) return(data.table())
  dcast(dt, as.formula(paste(row_col, "~", col_col)), value.var = value_col, fill = 0)
}

# ---- Main loop --------------------------------------------------------------
for (i in seq_len(nrow(run_manifest))) {
  rr <- run_manifest[i]; run_id <- rr$run_id
  message("\nRUN: ", run_id)
  hit_sub <- filter_by_pid_qcov(hits, rr$pid_cut, qcov_cut)
  pruned  <- unique(prune_carriers(hit_sub, rr$pruning, delta_default, tau_default)[, .(ASV, BIN)])
  asv_ko  <- build_asv_ko(pruned, bin_ko2, rr$weighting)
  sample_rel <- normalize_sample(build_sample_ko(asv_abund, asv_ko))
  sample_rel_wide <- to_wide(sample_rel, "Sample", "KO_ID", "value")
  out_file <- file.path(OUT_DIR, paste0("sample_x_ko.", run_id, ".rel.wide.tsv"))
  fwrite(sample_rel_wide, out_file, sep = "\t")
  cat(sprintf("  Saved %s (%d samples x %d KOs)\n", basename(out_file),
              nrow(sample_rel_wide), ncol(sample_rel_wide) - 1))
}
cat("\n=== 01 DONE ===\n")
