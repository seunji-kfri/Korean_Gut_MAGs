# =============================================================================
# 01_train_step1_cluster.R
#
# Step 1 of the two-step framework: predict dRep95 species-cluster abundance
# from CLR-transformed ASV abundance, using per-cluster XGBoost regressors.
#
# Performance is estimated by 5-fold group-stratified cross-validation, and an
# optional Ridge baseline is compared (manuscript Fig. 3B-E). The CV predictions
# are saved for the downstream KO step (02_predict_step2_ko.R).
#
# NOTE on models:
#   - This script produces CROSS-VALIDATED predictions for performance
#     evaluation (Fig. 2/3).
#   - The reusable FULL-cohort model applied to external data is produced by
#     04_external_apply.R and stored in models/full_xgb_models.rds.
#
# Run from the repository root:
#   Rscript 07_ML_Prediction/01_train_step1_cluster.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(xgboost)
})

source(file.path("07_ML_Prediction", "config.R"))

outdir <- file.path(OUTPUT_DIR, "01_train_step1_cluster")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

N_FOLDS <- 5

# ---- 1. Load ----------------------------------------------------------------
cat("[1/6] Loading inputs...\n")
asv_raw   <- read_csv(FILES$asv_internal, show_col_types = FALSE)
drep95_df <- read_csv(FILES$drep95, show_col_types = FALSE) %>%
  mutate(bin_id = str_remove(genome, "\\.fna$"))
cpm_ra    <- fread(FILES$cpm_taxonomy)

# Sample metadata used for group-stratified folds
meta_dt <- fread(FILES$meta_internal)
if ("Group" %in% names(meta_dt)) setnames(meta_dt, "Group", "group")
meta_dt <- meta_dt[group %in% c("Healthy", "MASLD"), .(sample_id, group)]

# ---- 2. Build cluster abundance target --------------------------------------
cat("[2/6] Building cluster abundance target...\n")
bin2cluster <- drep95_df %>% select(bin_id, drep95_cluster = secondary_cluster)

cpm_cluster <- cpm_ra %>%
  select(sample_id, bin_id, CPM_RA) %>%
  inner_join(bin2cluster, by = "bin_id") %>%
  group_by(sample_id, drep95_cluster) %>%
  summarise(abund = sum(CPM_RA), .groups = "drop")

prev95 <- cpm_cluster %>%
  filter(abund > 0) %>%
  group_by(drep95_cluster) %>%
  summarise(n_samples = n_distinct(sample_id), .groups = "drop")

# ---- 3. ASV CLR input -------------------------------------------------------
cat("[3/6] Building CLR-transformed ASV input...\n")
asv_mat <- asv_raw %>% select(-ASV) %>% as.matrix()
rownames(asv_mat) <- asv_raw$ASV
asv_keep <- rowSums(asv_mat > 0) >= PARAMS$asv_prev_cut
asv_mat <- asv_mat[asv_keep, ]
cat(sprintf("  ASVs after prev>=%d filter: %d\n", PARAMS$asv_prev_cut, nrow(asv_mat)))

pseudo_asv <- min(asv_mat[asv_mat > 0]) / 2
asv_clr <- apply(asv_mat + pseudo_asv, 2, function(x) log(x) - mean(log(x)))
asv_clr_t <- t(asv_clr)   # samples x ASVs

# ---- 4. Sample alignment & folds -------------------------------------------
cat("[4/6] Aligning samples & assigning folds...\n")
all_samples <- sort(intersect(rownames(asv_clr_t), unique(cpm_cluster$sample_id)))
X <- asv_clr_t[all_samples, ]

meta <- as.data.table(meta_dt)[sample_id %in% all_samples]
meta <- meta[match(all_samples, sample_id)]
stopifnot(all(meta$sample_id == all_samples))
cat(sprintf("  samples: %d\n", length(all_samples)))

make_folds <- function(group_vec, k = 5, seed = 42) {
  set.seed(seed)
  folds <- integer(length(group_vec))
  for (g in unique(group_vec)) {
    idx <- which(group_vec == g)
    folds[idx] <- sample(rep(1:k, length.out = length(idx)))
  }
  folds
}
folds <- make_folds(meta$group, k = N_FOLDS)

# ---- 5. Helpers -------------------------------------------------------------
build_Y_log <- function(clusters_keep, samples) {
  Y <- cpm_cluster %>%
    filter(drep95_cluster %in% clusters_keep, sample_id %in% samples) %>%
    complete(sample_id = samples, drep95_cluster = clusters_keep,
             fill = list(abund = 0)) %>%
    pivot_wider(names_from = drep95_cluster, values_from = abund,
                values_fill = 0) %>%
    column_to_rownames("sample_id") %>% as.matrix()
  Y <- Y[samples, , drop = FALSE]
  nz <- Y[Y > 0]
  pseudo <- if (length(nz)) min(nz) / 2 else 1e-6
  log(Y + pseudo)
}

predict_xgb <- function(X_tr, Y_tr, X_te,
                        nrounds = 200, early_stopping_rounds = 20) {
  n_mag  <- ncol(Y_tr)
  Y_pred <- matrix(0, nrow = nrow(X_te), ncol = n_mag)
  colnames(Y_pred) <- colnames(Y_tr)
  for (j in seq_len(n_mag)) {
    y <- Y_tr[, j]
    if (sd(y, na.rm = TRUE) < 1e-10) { Y_pred[, j] <- mean(y, na.rm = TRUE); next }
    n_tr <- length(y)
    inner_val <- sample(seq_len(n_tr), size = max(8, round(n_tr * 0.2)))
    inner_tr  <- setdiff(seq_len(n_tr), inner_val)
    dtrain <- xgb.DMatrix(data = X_tr[inner_tr, , drop = FALSE], label = y[inner_tr])
    dval   <- xgb.DMatrix(data = X_tr[inner_val, , drop = FALSE], label = y[inner_val])
    fit <- tryCatch(
      xgb.train(params = PARAMS$xgb_params, data = dtrain, nrounds = nrounds,
                watchlist = list(val = dval),
                early_stopping_rounds = early_stopping_rounds, verbose = 0),
      error = function(e) NULL)
    if (is.null(fit)) { Y_pred[, j] <- mean(y); next }
    Y_pred[, j] <- predict(fit, X_te)
  }
  Y_pred
}

run_cv <- function(X, Y, folds) {
  Y_pred <- matrix(NA_real_, nrow(Y), ncol(Y),
                   dimnames = list(rownames(Y), colnames(Y)))
  for (k in seq_len(N_FOLDS)) {
    tk <- Sys.time()
    te <- which(folds == k)
    tr <- setdiff(seq_len(nrow(X)), te)
    Y_pred[te, ] <- predict_xgb(X[tr, ], Y[tr, ], X[te, ])
    cat(sprintf("    fold %d done (%.1f min)\n", k,
                as.numeric(difftime(Sys.time(), tk, units = "mins"))))
  }
  Y_pred
}

evaluate_cv <- function(Y_true, Y_pred) {
  spr <- sapply(seq_len(ncol(Y_true)), function(j)
    suppressWarnings(cor(Y_true[, j], Y_pred[, j], method = "spearman")))
  r2  <- sapply(seq_len(ncol(Y_true)), function(j) {
    y <- Y_true[, j]; p <- Y_pred[, j]
    ss_res <- sum((y - p)^2); ss_tot <- sum((y - mean(y))^2)
    if (ss_tot == 0) return(NA_real_)
    1 - ss_res / ss_tot
  })
  spr_samp <- sapply(seq_len(nrow(Y_true)), function(i)
    suppressWarnings(cor(Y_true[i, ], Y_pred[i, ], method = "spearman")))
  list(spr_mag_median = median(spr, na.rm = TRUE),
       spr_mag_mean   = mean(spr, na.rm = TRUE),
       r2_mag_median  = median(r2, na.rm = TRUE),
       spr_samp_median = median(spr_samp, na.rm = TRUE),
       rmse_global    = sqrt(mean((Y_true - Y_pred)^2)),
       per_mag_spr = spr, per_mag_r2 = r2, per_samp_spr = spr_samp)
}

# ---- 6. Run CV for the prevalence cutoffs -----------------------------------
cat("[5/6] Running XGBoost CV...\n")
# Main analysis uses prev10; prev20 is reported as a sensitivity comparison.
prev_cutoffs <- c(PARAMS$cluster_prev, 20)
prev_cutoffs <- unique(prev_cutoffs)

summary_list <- list()
for (prev_cut in prev_cutoffs) {
  clusters_keep <- prev95 %>% filter(n_samples >= prev_cut) %>% pull(drep95_cluster)
  cat(sprintf("\n  === prev>=%d : %d clusters ===\n", prev_cut, length(clusters_keep)))

  Y <- build_Y_log(clusters_keep, all_samples)
  t0 <- Sys.time()
  Y_pred <- run_cv(X, Y, folds)
  cat(sprintf("  CV total: %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  ev <- evaluate_cv(Y, Y_pred)
  cat(sprintf("  per-cluster rho median=%.3f  per-sample rho median=%.3f  R2=%.3f\n",
              ev$spr_mag_median, ev$spr_samp_median, ev$r2_mag_median))

  lbl <- sprintf("xgb_prev%d_log", prev_cut)
  saveRDS(list(Y = Y, Y_pred = Y_pred, eval = ev),
          file.path(outdir, sprintf("cv_%s.rds", lbl)))

  summary_list[[lbl]] <- tibble(
    label = lbl, prev_cut = prev_cut, n_clusters = ncol(Y),
    spearman_per_mag_median    = ev$spr_mag_median,
    spearman_per_mag_mean      = ev$spr_mag_mean,
    spearman_per_sample_median = ev$spr_samp_median,
    r2_per_mag_median          = ev$r2_mag_median,
    rmse_global                = ev$rmse_global)
}

summary_tbl <- bind_rows(summary_list)
write_tsv(summary_tbl, file.path(outdir, "summary.tsv"))
cat("\n[6/6] === CV summary ===\n")
print(summary_tbl)
cat(sprintf("\nSaved CV predictions and summary to: %s\n", outdir))
cat("Next: 02_predict_step2_ko.R (cluster -> KO via carrier matrix)\n")
