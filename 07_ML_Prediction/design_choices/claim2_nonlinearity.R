# =============================================================================
# claim2_nonlinearity.R
#
# Design-choice analysis (manuscript Fig. 3B-E): the non-linear cluster
# predictor (XGBoost) outperforms the linear one (Ridge) at the per-cluster
# level, while the two converge at the per-sample level. Compares per-cluster
# and per-sample Spearman between the two models with paired Wilcoxon tests.
#
# Inputs (produced by earlier steps):
#   - per-cluster Spearman table with XGB and Ridge columns
#     (config FILES$xgb_vs_ridge_spearman; derived from the CV outputs of
#      01_train_step1_cluster.R for XGBoost and the Ridge baseline)
#   - per-sample similarity table with a `method` column and `spearman`
#     (config FILES$per_sample_similarity)
#
# Run from the repository root:
#   Rscript 07_ML_Prediction/design_choices/claim2_nonlinearity.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("07_ML_Prediction", "config.R"))

out_dir <- file.path(OUTPUT_DIR, "design_choices", "claim2_nonlinearity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Per-cluster Spearman (XGB vs Ridge) ---------------------------------
cat("[1/2] Per-cluster XGB vs Ridge...\n")
cluster_dt <- fread(FILES$xgb_vs_ridge_spearman)
ridge_col <- grep("ridge|Ridge", names(cluster_dt), value = TRUE)[1]
xgb_col   <- grep("xgb|XGB|XGBoost", names(cluster_dt), value = TRUE)[1]
cat(sprintf("  Ridge col: %s ; XGB col: %s\n", ridge_col, xgb_col))

cluster_dt[, delta := get(xgb_col) - get(ridge_col)]
cluster_dt[, xgb_wins := delta > 0]
cluster_dt[, strong_advantage := delta > 0.1]
wt_cluster <- wilcox.test(cluster_dt[[xgb_col]], cluster_dt[[ridge_col]],
                          paired = TRUE, exact = FALSE)
cluster_summary <- data.table(
  level = "per-cluster", n_total = nrow(cluster_dt),
  xgb_median   = median(cluster_dt[[xgb_col]], na.rm = TRUE),
  ridge_median = median(cluster_dt[[ridge_col]], na.rm = TRUE),
  median_delta = median(cluster_dt$delta, na.rm = TRUE),
  xgb_wins_n   = sum(cluster_dt$xgb_wins, na.rm = TRUE),
  xgb_wins_pct = 100 * mean(cluster_dt$xgb_wins, na.rm = TRUE),
  strong_advantage_n = sum(cluster_dt$strong_advantage, na.rm = TRUE),
  wilcoxon_p   = wt_cluster$p.value)
print(cluster_summary)
fwrite(cluster_dt, file.path(out_dir, "per_cluster_xgb_vs_ridge.tsv"), sep = "\t")

# ---- 2. Per-sample Spearman (XGB vs Ridge) ----------------------------------
cat("[2/2] Per-sample XGB vs Ridge...\n")
per_sample <- fread(FILES$per_sample_similarity)
xgb_spr   <- per_sample[grepl("XGB", method), .(sample_id, xgb = spearman)]
ridge_spr <- per_sample[grepl("Ridge", method), .(sample_id, ridge = spearman)]
paired <- merge(xgb_spr, ridge_spr, by = "sample_id")
paired[, delta := xgb - ridge]
paired[, xgb_wins := delta > 0]
wt_sample <- wilcox.test(paired$xgb, paired$ridge, paired = TRUE, exact = FALSE)
sample_summary <- data.table(
  level = "per-sample", n_total = nrow(paired),
  xgb_median = median(paired$xgb), ridge_median = median(paired$ridge),
  median_delta = median(paired$delta),
  xgb_wins_n = sum(paired$xgb_wins),
  xgb_wins_pct = 100 * mean(paired$xgb_wins),
  strong_advantage_n = sum(paired$delta > 0.05),
  wilcoxon_p = wt_sample$p.value)
print(sample_summary)
fwrite(paired, file.path(out_dir, "per_sample_xgb_vs_ridge.tsv"), sep = "\t")

fwrite(rbind(cluster_summary, sample_summary),
       file.path(out_dir, "claim2_summary.tsv"), sep = "\t")
cat("\nDone.\n")
