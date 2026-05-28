# =============================================================================
# claim4_aggregation.R
#
# Design-choice analysis (manuscript Fig. 3F): converting predicted cluster
# abundance to KO abundance through a BINARY carrier matrix outperforms
# FPKM-weighted (mean or sum) aggregation. Builds the FPKM mean/sum variants
# from the XGB CV predictions, runs MaAsLin2 on each, and compares F1 vs gold.
#
# Inputs (produced by earlier steps):
#   - XGB CV predictions: 01_train_step1_cluster/cv_xgb_prev10_log.rds
#   - bin x KO FPKM table (config FILES$bin_ko_fpkm)
#   - dRep95 cluster map  (config FILES$drep95)
#   - the binary-carrier KO tables (from 02_predict_step2_ko output) for
#     XGB and, optionally, Ridge
#   - gold MaAsLin2 results (03_benchmark/internal_gold/all_results.tsv)
#
# Run from the repository root (after 01, 02, 03):
#   Rscript 07_ML_Prediction/design_choices/claim4_aggregation.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Maaslin2)
})

source(file.path("07_ML_Prediction", "config.R"))

step1_dir <- file.path(OUTPUT_DIR, "01_train_step1_cluster")
step2_dir <- file.path(OUTPUT_DIR, "02_predict_step2_ko")
bench_dir <- file.path(OUTPUT_DIR, "03_benchmark")
out_dir   <- file.path(OUTPUT_DIR, "design_choices", "claim4_aggregation")
maaslin_dir   <- file.path(out_dir, "maaslin2")
ko_matrix_dir <- file.path(out_dir, "ko_matrices")
dir.create(maaslin_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ko_matrix_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Build XGB FPKM-mean and FPKM-sum aggregation variants ---------------
cat("[1/3] Building XGB FPKM aggregation variants...\n")
xgb_mean_file <- file.path(ko_matrix_dir, "sample_x_ko.xgb_mean.rel.wide.tsv")
xgb_sum_file  <- file.path(ko_matrix_dir, "sample_x_ko.xgb_sum.rel.wide.tsv")

if (!file.exists(xgb_mean_file) || !file.exists(xgb_sum_file)) {
  cv_xgb <- readRDS(file.path(step1_dir, "cv_xgb_prev10_log.rds"))
  Y_pred <- exp(cv_xgb$Y_pred)  # back to linear

  bin_ko <- fread(FILES$bin_ko_fpkm)
  bin_ids <- bin_ko[[1]]
  ko_mat <- as.matrix(bin_ko[, -1]); rownames(ko_mat) <- bin_ids
  colnames(ko_mat) <- sub("^ko[.:]", "", colnames(ko_mat))

  drep95 <- read.csv(FILES$drep95)
  bin_col <- intersect(names(drep95), c("genome", "user_genome", "bin_id", "BIN_ID"))[1]
  cl_col  <- intersect(names(drep95), c("secondary_cluster", "cluster"))[1]
  bin_ids_drep <- sub("\\.(fna|fa|fasta)$", "", drep95[[bin_col]])
  bin_to_cluster <- setNames(drep95[[cl_col]], bin_ids_drep)

  shared_bins <- intersect(rownames(ko_mat), names(bin_to_cluster))
  ko_mat <- ko_mat[shared_bins, ]
  clusters_vec <- bin_to_cluster[shared_bins]

  cluster_ko_mean <- rowsum(ko_mat, group = clusters_vec, na.rm = TRUE) /
    as.numeric(table(clusters_vec)[unique(clusters_vec)])
  cluster_ko_sum  <- rowsum(ko_mat, group = clusters_vec, na.rm = TRUE)

  common_cl <- intersect(rownames(cluster_ko_mean), colnames(Y_pred))
  cluster_ko_mean <- cluster_ko_mean[common_cl, , drop = FALSE]
  cluster_ko_sum  <- cluster_ko_sum[common_cl, , drop = FALSE]
  Y_pred_aligned  <- Y_pred[, common_cl, drop = FALSE]

  to_rel <- function(m) { rs <- rowSums(m, na.rm = TRUE); rs[rs == 0] <- 1; m / rs }
  sample_ko_mean <- to_rel(Y_pred_aligned %*% cluster_ko_mean)
  sample_ko_sum  <- to_rel(Y_pred_aligned %*% cluster_ko_sum)

  fwrite(data.table(sample_id = rownames(sample_ko_mean), as.data.table(sample_ko_mean)),
         xgb_mean_file, sep = "\t")
  fwrite(data.table(sample_id = rownames(sample_ko_sum), as.data.table(sample_ko_sum)),
         xgb_sum_file, sep = "\t")
} else {
  cat("  Variant files already exist, skipping build.\n")
}

# ---- 2. MaAsLin2 for aggregation variants -----------------------------------
cat("[2/3] Running MaAsLin2...\n")
meta <- fread(FILES$meta_internal)
if ("Group" %in% names(meta)) setnames(meta, "Group", "group")
meta_df <- as.data.frame(meta); rownames(meta_df) <- meta_df$sample_id

load_ko <- function(f) { df <- as.data.frame(fread(f)); rownames(df) <- df[[1]]; df[[1]] <- NULL; df }
run_one <- function(ko_file, out_subdir) {
  if (file.exists(file.path(out_subdir, "all_results.tsv"))) return(invisible())
  if (!file.exists(ko_file)) { cat(sprintf("  SKIP (not found): %s\n", ko_file)); return(invisible()) }
  ko_df <- load_ko(ko_file)
  common <- intersect(rownames(ko_df), rownames(meta_df))
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    Maaslin2(input_data = ko_df[common, , drop = FALSE],
             input_metadata = meta_df[common, , drop = FALSE], output = out_subdir,
             fixed_effects = c("group"), reference = c("group,Healthy"),
             normalization = "NONE", transform = "LOG", analysis_method = "LM",
             min_prevalence = 0.1, max_significance = 1.0,
             plot_heatmap = FALSE, plot_scatter = FALSE),
    error = function(e) cat(sprintf("  ERROR: %s\n", e$message)))
}

# Binary carrier (winner) comes from 02_predict_step2_ko output
entries_agg <- list(
  xgb_carrier = file.path(step2_dir, "sample_x_ko.ours.rel.wide.tsv"),
  xgb_mean    = xgb_mean_file,
  xgb_sum     = xgb_sum_file)
for (e in names(entries_agg)) run_one(entries_agg[[e]], file.path(maaslin_dir, e))

# ---- 3. F1 comparison -------------------------------------------------------
cat("[3/3] F1 by aggregation strategy...\n")
gold <- fread(file.path(bench_dir, "internal_gold", "all_results.tsv"))
gold <- gold[metadata == "group" & value == "MASLD"]
gold[, feature := sub("^ko[.:]", "", feature)]
gold_sig <- gold[qval < 0.05, feature]

compute_f1 <- function(sig_kos, gold_sig_kos, all_kos) {
  truth <- all_kos %in% gold_sig_kos; pred <- all_kos %in% sig_kos
  TP <- sum(truth & pred); FP <- sum(!truth & pred); FN <- sum(truth & !pred)
  P <- if (TP + FP > 0) TP / (TP + FP) else 0
  R <- if (TP + FN > 0) TP / (TP + FN) else 0
  data.table(precision = P, recall = R,
             F1 = if (P + R > 0) 2 * P * R / (P + R) else 0)
}

variant_info <- data.table(
  entry = names(entries_agg),
  aggregation = c("Binary carrier", "FPKM mean", "FPKM sum"))
f1_list <- list()
for (e in names(entries_agg)) {
  f <- file.path(maaslin_dir, e, "all_results.tsv")
  if (!file.exists(f)) { cat(sprintf("  SKIP %s\n", e)); next }
  r <- fread(f)[metadata == "group" & value == "MASLD"]
  r[, feature := sub("^ko[.:]", "", feature)]
  m <- compute_f1(r[qval < 0.05, feature], gold_sig, r$feature)
  m[, `:=`(entry = e, n_sig = r[qval < 0.05, .N])]
  f1_list[[e]] <- m
}
f1_agg <- merge(rbindlist(f1_list), variant_info, by = "entry")
print(f1_agg[, .(aggregation, F1 = round(F1, 3),
                 precision = round(precision, 3), recall = round(recall, 3), n_sig)])
fwrite(f1_agg, file.path(out_dir, "claim4_f1_aggregation.tsv"), sep = "\t")
cat("\nDone.\n")
