# =============================================================================
# 02_predict_step2_ko.R
#
# Step 2 of the two-step framework: convert predicted dRep95 cluster abundances
# (from 01_train_step1_cluster.R) into sample x KO functional profiles via a
# binary carrier matrix, and evaluate against the shotgun-derived gold standard
# (F1 / precision / recall / effect-size Spearman).
#
# The carrier matrix encodes, for each (cluster, KO), the fraction of cluster
# member MAGs that carry the KO. A carrier threshold (default 0.25) binarises
# this matrix; a small threshold scan is run to confirm the chosen value.
#
# Run from the repository root (after 01_train_step1_cluster.R):
#   Rscript 07_ML_Prediction/02_predict_step2_ko.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

source(file.path("07_ML_Prediction", "config.R"))

step1_dir <- file.path(OUTPUT_DIR, "01_train_step1_cluster")
outdir    <- file.path(OUTPUT_DIR, "02_predict_step2_ko")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

CARRIER_THRS <- c(0, 0.25, 0.5, 0.67, 0.75)  # scan; manuscript winner = 0.25

# ---- 1. Load cluster map and bin-KO content ---------------------------------
cat("[1/5] Loading cluster map and bin-KO content...\n")
drep95_df <- read_csv(FILES$drep95, show_col_types = FALSE) %>%
  mutate(bin_id = str_remove(genome, "\\.fna$"))
bin2cluster <- as.data.table(drep95_df)[, .(bin_id, drep95_cluster = secondary_cluster)]

bin_ko_wide <- fread(FILES$bin_ko_fpkm)
if (names(bin_ko_wide)[1] != "bin_id") setnames(bin_ko_wide, 1, "bin_id")
ko_cols  <- setdiff(names(bin_ko_wide), "bin_id")
ko_clean <- sub("^ko:", "", ko_cols)
setnames(bin_ko_wide, ko_cols, ko_clean)

# ---- 2. Gold standard + metadata --------------------------------------------
cat("[2/5] Loading gold standard and metadata...\n")
fpkm_gold <- fread(FILES$gold_fpkm)
setnames(fpkm_gold, c("KO", "fpkm_sum"), c("KO_ID", "abundance"), skip_absent = TRUE)
fpkm_gold[, abundance := abundance / sum(abundance, na.rm = TRUE), by = sample_id]

meta_dt <- fread(FILES$meta_internal)
if ("Group" %in% names(meta_dt)) setnames(meta_dt, "Group", "group")
meta_dt <- meta_dt[group %in% c("Healthy", "MASLD"), .(sample_id, group)]

# ---- 3. Helpers -------------------------------------------------------------
make_carrier <- function(clusters_keep) {
  bin2c <- bin2cluster[drep95_cluster %in% clusters_keep]
  bk <- merge(bin_ko_wide[bin_id %in% unique(bin2c$bin_id)], bin2c, by = "bin_id")
  ko_mat <- as.matrix(bk[, ..ko_clean])
  cvec <- bk$drep95_cluster
  presence <- (ko_mat > 0) * 1L
  carrier <- rowsum(presence, group = cvec, reorder = TRUE) /
    as.vector(table(cvec)[sort(unique(cvec))])
  carrier[clusters_keep, , drop = FALSE]
}

log_to_linear <- function(Y_pred_log, pseudo) {
  M <- exp(Y_pred_log) - pseudo
  M[M < 0] <- 0
  rs <- rowSums(M)
  sweep(M, 1, ifelse(rs > 0, rs / 100, 1), "/")
}

run_feature_test <- function(wide_dt, sample_col, group_dt) {
  dt <- merge(copy(wide_dt), group_dt[, .(sample_id, group)],
              by.x = sample_col, by.y = "sample_id")
  feature_cols <- setdiff(names(dt), c(sample_col, "group"))
  res <- rbindlist(lapply(feature_cols, function(feat) {
    x_h <- dt[group == "Healthy"][[feat]]
    x_m <- dt[group == "MASLD"][[feat]]
    beta <- mean(x_m, na.rm = TRUE) - mean(x_h, na.rm = TRUE)
    pval <- tryCatch(wilcox.test(x_m, x_h)$p.value, error = function(e) NA_real_)
    data.table(feature_id = feat, effect_size = beta, pvalue = pval)
  }))
  res[, padj := p.adjust(pvalue, method = "BH")]
  res[, significant := !is.na(padj) & padj < 0.05]
  res[]
}

compare_to_gold <- function(gold_stats, model_stats, label) {
  x <- merge(
    gold_stats[, .(feature_id, gs = significant, ge = effect_size)],
    model_stats[, .(feature_id, ms = significant, me = effect_size)],
    by = "feature_id", all = TRUE)
  x[is.na(gs), gs := FALSE]; x[is.na(ms), ms := FALSE]
  x[is.na(ge), ge := 0];     x[is.na(me), me := 0]
  TP <- x[gs & ms, .N]; FP <- x[!gs & ms, .N]
  FN <- x[gs & !ms, .N]; TN <- x[!gs & !ms, .N]
  P  <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
  R  <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  F1 <- if (!is.na(P) && !is.na(R) && P + R > 0) 2 * P * R / (P + R) else NA_real_
  sp <- suppressWarnings(cor(x$ge, x$me, method = "spearman"))
  data.table(method = label, TP = TP, FP = FP, FN = FN, TN = TN,
             precision = round(P, 4), recall = round(R, 4),
             F1 = round(F1, 4), spearman_effect = round(sp, 4),
             n_model_sig = sum(x$ms))
}

# Recover the cluster pseudo-count used at training time (from raw CPM data)
cpm_ra <- fread(FILES$cpm_taxonomy)
cpm_cluster <- cpm_ra %>%
  select(sample_id, bin_id, CPM_RA) %>%
  inner_join(as_tibble(bin2cluster), by = "bin_id") %>%
  group_by(sample_id, drep95_cluster) %>%
  summarise(abund = sum(CPM_RA), .groups = "drop")

recover_pseudo <- function(clusters_keep, samples) {
  Y <- cpm_cluster %>%
    filter(drep95_cluster %in% clusters_keep, sample_id %in% samples) %>%
    complete(sample_id = samples, drep95_cluster = clusters_keep,
             fill = list(abund = 0)) %>%
    pivot_wider(names_from = drep95_cluster, values_from = abund,
                values_fill = 0) %>%
    column_to_rownames("sample_id") %>% as.matrix()
  nz <- Y[Y > 0]; min(nz) / 2
}

# Gold differential abundance (computed once)
gold_wide  <- dcast(fpkm_gold, sample_id ~ KO_ID, value.var = "abundance", fill = 0)
gold_stats <- run_feature_test(gold_wide, "sample_id", meta_dt)
cat(sprintf("  Gold significant KOs: %d / %d\n",
            sum(gold_stats$significant), nrow(gold_stats)))

# ---- 4. cluster -> KO for each CV prediction + threshold scan ----------------
cat("[3/5] Converting cluster predictions to KO profiles...\n")
cv_files <- list.files(step1_dir, pattern = "^cv_xgb_prev.*\\.rds$", full.names = TRUE)
if (length(cv_files) == 0) {
  stop("No CV prediction files found in ", step1_dir,
       ". Run 01_train_step1_cluster.R first.", call. = FALSE)
}

ko_results <- list()
for (cvf in cv_files) {
  lbl_base <- sub("^cv_", "", sub("\\.rds$", "", basename(cvf)))
  cv <- readRDS(cvf)
  clusters_keep <- colnames(cv$Y)
  samples <- rownames(cv$Y)

  pseudo  <- recover_pseudo(clusters_keep, samples)
  Y_pred_linear <- log_to_linear(cv$Y_pred, pseudo)
  carrier <- make_carrier(clusters_keep)

  for (thr in CARRIER_THRS) {
    lbl <- sprintf("%s_thr%.2f", lbl_base, thr)
    pm <- carrier; if (thr > 0) pm[pm < thr] <- 0
    sko <- Y_pred_linear %*% pm
    rs <- rowSums(sko)
    sko <- sweep(sko, 1, ifelse(rs > 0, rs, 1), "/")
    wide_dt <- as.data.table(sko, keep.rownames = "sample_id")
    stats <- run_feature_test(wide_dt, "sample_id", meta_dt)
    res <- compare_to_gold(gold_stats, stats, lbl)
    res[, `:=`(model = lbl_base, thr = thr, n_clusters = ncol(cv$Y))]
    cat(sprintf("  %-28s F1=%.3f P=%.3f R=%.3f rho=%.3f\n",
                lbl, res$F1, res$precision, res$recall, res$spearman_effect))
    ko_results[[lbl]] <- res

    # Save the winner-config profile (prev10, thr 0.25) for reuse downstream
    if (grepl("prev10", lbl_base) && abs(thr - PARAMS$carrier_thr) < 1e-9) {
      fwrite(wide_dt, file.path(outdir, "sample_x_ko.ours.rel.wide.tsv"), sep = "\t")
    }
  }
}

# ---- 5. Summary -------------------------------------------------------------
final <- rbindlist(ko_results, fill = TRUE)
cat("\n[4/5] === KO-level results (sorted by F1) ===\n")
print(final[order(-F1), .(method, thr, precision, recall, F1, spearman_effect)])
fwrite(final, file.path(outdir, "step2_ko_results.tsv"), sep = "\t")

cat(sprintf("\n[5/5] Saved results to: %s\n", outdir))
cat("  step2_ko_results.tsv        (threshold scan)\n")
cat("  sample_x_ko.ours.rel.wide.tsv  (winner config: prev10, thr 0.25)\n")
