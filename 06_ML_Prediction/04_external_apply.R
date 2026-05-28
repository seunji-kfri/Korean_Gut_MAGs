# =============================================================================
# 04_external_apply.R
#
# Train the full-cohort XGBoost models on all 224 internal samples, save them
# for reuse, then apply them to an independent amplicon-only external cohort.
# Produces external sample x KO profiles and a cross-cohort consistency table
# (manuscript Fig. 4).
#
# Pipeline:
#   PART 0  Load internal + external data; normalise PICRUSt2 external
#   PART 1  Train full-224 XGBoost models (prev>=10); save to MODEL_DIR
#   PART 2  Align external ASVs to training features; predict cluster abundance;
#           build binary carrier profile; convert to sample x KO
#   PART 3  Wilcoxon Healthy vs MASLD on external + internal
#   PART 4  Cross-cohort consistency (recovery, direction agreement, effect rho)
#
# Run from the repository root:
#   Rscript 07_ML_Prediction/04_external_apply.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(xgboost)
})

source(file.path("07_ML_Prediction", "config.R"))

outdir <- file.path(OUTPUT_DIR, "04_external_apply")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# PART 0: Setup
# =============================================================================
cat("=== PART 0: Setup ===\n")

# Internal ASV relative abundance (ASV x sample)
asv_int_raw <- read_csv(FILES$asv_internal, show_col_types = FALSE)

# genome -> dRep95 cluster map
drep95_df <- read_csv(FILES$drep95, show_col_types = FALSE) %>%
  mutate(bin_id = str_remove(genome, "\\.fna$"))
bin2cluster <- drep95_df %>% select(bin_id, drep95_cluster = secondary_cluster)

# Per-sample MAG CPM -> cluster-level abundance
cpm_ra <- fread(FILES$cpm_taxonomy)
cpm_cluster <- cpm_ra %>%
  select(sample_id, bin_id, CPM_RA) %>%
  inner_join(bin2cluster, by = "bin_id") %>%
  group_by(sample_id, drep95_cluster) %>%
  summarise(abund = sum(CPM_RA), .groups = "drop")

prev95 <- cpm_cluster %>%
  filter(abund > 0) %>%
  group_by(drep95_cluster) %>%
  summarise(n_samples = n_distinct(sample_id), .groups = "drop")

# Internal metadata
meta_int <- fread(FILES$meta_internal)
if ("Group" %in% names(meta_int)) setnames(meta_int, "Group", "group")
meta_int <- meta_int[group %in% c("Healthy", "MASLD"), .(sample_id, group)]

# External data
asv_ext_raw <- fread(FILES$asv_external)
setnames(asv_ext_raw, names(asv_ext_raw)[1], "ASV")
meta_ext <- fread(FILES$meta_external)
setnames(meta_ext, c("SampleID", "Group"), c("sample_id", "group"))
meta_ext <- meta_ext[group %in% c("Healthy", "MASLD"), .(sample_id, group)]

# PICRUSt2 external (raw count -> relative abundance)
pic_ext_raw <- fread(FILES$picrust_external)
setnames(pic_ext_raw, names(pic_ext_raw)[1], "sample_id")
pic_ext_mat <- as.matrix(pic_ext_raw[, -1])
rownames(pic_ext_mat) <- pic_ext_raw$sample_id
rs <- rowSums(pic_ext_mat)
pic_ext_mat <- sweep(pic_ext_mat, 1, ifelse(rs > 0, rs, 1), "/")
pic_ext <- as.data.table(pic_ext_mat, keep.rownames = "sample_id")
cat(sprintf("  PICRUSt2 external normalised: %d samples x %d KOs\n",
            nrow(pic_ext), ncol(pic_ext) - 1))

# PICRUSt2 internal
pic_int <- fread(FILES$picrust_internal)
setnames(pic_int, names(pic_int)[1], "sample_id")

# Gold (shotgun-derived) KO relative abundance, internal
fpkm_gold <- fread(FILES$gold_fpkm)
setnames(fpkm_gold, c("KO", "fpkm_sum"), c("KO_ID", "abundance"), skip_absent = TRUE)
fpkm_gold[, abundance := abundance / sum(abundance, na.rm = TRUE), by = sample_id]
gold_wide <- dcast(fpkm_gold, sample_id ~ KO_ID, value.var = "abundance", fill = 0)

# =============================================================================
# PART 1: Train full-224 XGBoost models
# =============================================================================
cat("\n=== PART 1: Training full-224 XGBoost models ===\n")

# ASV CLR transform (internal)
asv_int_ids <- asv_int_raw$ASV
asv_int_mat <- asv_int_raw %>% select(-ASV) %>% as.matrix()
rownames(asv_int_mat) <- asv_int_ids
asv_keep <- rowSums(asv_int_mat > 0) >= PARAMS$asv_prev_cut
asv_int_filt <- asv_int_mat[asv_keep, ]
cat(sprintf("  Internal ASV after prev>=%d filter: %d\n",
            PARAMS$asv_prev_cut, nrow(asv_int_filt)))

pseudo_asv <- min(asv_int_filt[asv_int_filt > 0]) / 2
asv_int_clr <- apply(asv_int_filt + pseudo_asv, 2,
                     function(x) log(x) - mean(log(x)))
X_int_full <- t(asv_int_clr)  # samples x ASVs
all_samples <- rownames(X_int_full)

# Build log-space target matrix for a set of clusters
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
  list(Y_log = log(Y + pseudo), pseudo = pseudo, Y_linear = Y)
}

train_full_xgb <- function(X_tr, Y_tr, nrounds) {
  n_mag <- ncol(Y_tr)
  models <- vector("list", n_mag)
  names(models) <- colnames(Y_tr)
  for (j in seq_len(n_mag)) {
    y <- Y_tr[, j]
    if (sd(y, na.rm = TRUE) < 1e-10) {
      models[[j]] <- list(constant = mean(y, na.rm = TRUE))
      next
    }
    dtrain <- xgb.DMatrix(data = X_tr, label = y)
    fit <- tryCatch(
      xgb.train(params = PARAMS$xgb_params, data = dtrain,
                nrounds = nrounds, verbose = 0),
      error = function(e) NULL
    )
    models[[j]] <- if (is.null(fit)) list(constant = mean(y)) else fit
  }
  models
}

predict_with_models <- function(models, X_te) {
  out <- matrix(0, nrow = nrow(X_te), ncol = length(models))
  colnames(out) <- names(models)
  for (j in seq_along(models)) {
    m <- models[[j]]
    if (!is.null(m$constant)) out[, j] <- m$constant
    else out[, j] <- predict(m, X_te)
  }
  out
}

# Train on the main prevalence cutoff (prev10) described in the manuscript.
prev_cut <- PARAMS$cluster_prev
cat(sprintf("  Training for prev>=%d... ", prev_cut))
t0 <- Sys.time()
clusters_keep <- prev95 %>% filter(n_samples >= prev_cut) %>% pull(drep95_cluster)
yi <- build_Y_log(clusters_keep, all_samples)
full_models <- train_full_xgb(X_int_full, yi$Y_log, nrounds = PARAMS$xgb_nrounds)
cat(sprintf("%d clusters, %.1f min\n", length(clusters_keep),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Save trained models for external reuse
saveRDS(list(models = full_models, pseudo = yi$pseudo,
             feature_names = colnames(X_int_full)),
        MODEL_FILES$xgb_full)
cat(sprintf("  Saved trained models -> %s\n", MODEL_FILES$xgb_full))

# =============================================================================
# PART 2: Apply to external cohort
# =============================================================================
cat("\n=== PART 2: Applying to external cohort ===\n")

training_asv <- colnames(X_int_full)
ext_asv_ids <- asv_ext_raw$ASV
ext_in_train <- intersect(ext_asv_ids, training_asv)
missing_in_ext <- setdiff(training_asv, ext_asv_ids)
cat(sprintf("  External has %d / %d training ASVs (missing: %d -> 0-pad)\n",
            length(ext_in_train), length(training_asv), length(missing_in_ext)))

asv_ext_mat <- as.matrix(asv_ext_raw[, -1])
rownames(asv_ext_mat) <- asv_ext_raw$ASV
ext_samples <- intersect(colnames(asv_ext_mat), meta_ext$sample_id)
cat(sprintf("  External samples with metadata: %d\n", length(ext_samples)))

# Align external ASVs to training feature space (0-pad missing)
X_ext_aligned <- matrix(0, nrow = length(ext_samples), ncol = length(training_asv))
rownames(X_ext_aligned) <- ext_samples
colnames(X_ext_aligned) <- training_asv
for (asv in ext_in_train) {
  X_ext_aligned[, asv] <- asv_ext_mat[asv, ext_samples]
}

# CLR transform on external (per-sample), using internal pseudo-count
X_ext_clr <- t(apply(X_ext_aligned + pseudo_asv, 1,
                     function(x) log(x) - mean(log(x))))
cat(sprintf("  External CLR matrix: %d x %d\n", nrow(X_ext_clr), ncol(X_ext_clr)))

# log -> linear (relative abundance summing to 100)
log_to_linear <- function(Y_pred_log, pseudo) {
  M <- exp(Y_pred_log) - pseudo
  M[M < 0] <- 0
  rs <- rowSums(M)
  sweep(M, 1, ifelse(rs > 0, rs / 100, 1), "/")
}

# Build cluster x KO binary carrier profile
cat("  Building carrier profile...\n")
bin_ko_wide <- fread(FILES$bin_ko_fpkm)
if (names(bin_ko_wide)[1] != "bin_id") setnames(bin_ko_wide, 1, "bin_id")
ko_cols <- setdiff(names(bin_ko_wide), "bin_id")
ko_clean <- sub("^ko:", "", ko_cols)
setnames(bin_ko_wide, ko_cols, ko_clean)

make_carrier <- function(clusters_keep) {
  bin2c <- as.data.table(bin2cluster)[drep95_cluster %in% clusters_keep]
  bins_in <- unique(bin2c$bin_id)
  bk <- bin_ko_wide[bin_id %in% bins_in]
  bk <- merge(bk, bin2c, by = "bin_id")
  ko_mat <- as.matrix(bk[, ..ko_clean])
  cvec <- bk$drep95_cluster
  presence <- (ko_mat > 0) * 1L
  carrier <- rowsum(presence, group = cvec, reorder = TRUE) /
    as.vector(table(cvec)[sort(unique(cvec))])
  carrier[clusters_keep, , drop = FALSE]
}

# Predict external cluster abundance, then convert to sample x KO
Y_pred_log_ext <- predict_with_models(full_models, X_ext_clr)
rownames(Y_pred_log_ext) <- rownames(X_ext_clr)
Y_pred_linear_ext <- log_to_linear(Y_pred_log_ext, yi$pseudo)

clusters_pred <- colnames(Y_pred_linear_ext)
carrier <- make_carrier(clusters_pred)
carrier_thr <- carrier
carrier_thr[carrier_thr < PARAMS$carrier_thr] <- 0

# Save carrier profile (cluster x KO) for reuse / sharing
saveRDS(carrier_thr, MODEL_FILES$carrier)
cat(sprintf("  Saved carrier profile -> %s\n", MODEL_FILES$carrier))

sko_ext <- Y_pred_linear_ext %*% carrier_thr
rs <- rowSums(sko_ext)
sko_ext <- sweep(sko_ext, 1, ifelse(rs > 0, rs, 1), "/")

ext_ko_wide <- as.data.table(sko_ext, keep.rownames = "sample_id")
fwrite(ext_ko_wide, file.path(outdir, "ext_sample_x_ko.ours.rel.wide.tsv"),
       sep = "\t")
cat(sprintf("  Saved external KO profile (%d x %d)\n",
            nrow(ext_ko_wide), ncol(ext_ko_wide) - 1))

# =============================================================================
# PART 3: External + internal Wilcoxon (Healthy vs MASLD)
# =============================================================================
cat("\n=== PART 3: Wilcoxon Healthy vs MASLD ===\n")

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

ext_stats <- list()
ext_stats[["ours"]]     <- run_feature_test(ext_ko_wide, "sample_id", meta_ext)
ext_stats[["picrust2"]] <- run_feature_test(pic_ext,     "sample_id", meta_ext)

int_stats <- list()
int_stats[["picrust2"]] <- run_feature_test(pic_int,   "sample_id", meta_int)
int_stats[["gold"]]     <- run_feature_test(gold_wide, "sample_id", meta_int)

# Internal "ours" profile (produced by an earlier step in the pipeline).
# Falls back gracefully if the file is not present.
int_ours_path <- file.path(OUTPUT_DIR, "02_predict_step2_ko",
                           "sample_x_ko.ours.rel.wide.tsv")
if (file.exists(int_ours_path)) {
  int_ours <- fread(int_ours_path)
  setnames(int_ours, names(int_ours)[1], "sample_id")
  int_stats[["ours"]] <- run_feature_test(int_ours, "sample_id", meta_int)
} else {
  cat(sprintf("  [note] internal 'ours' profile not found at %s; ",
              int_ours_path))
  cat("skipping internal-ours consistency.\n")
}

# =============================================================================
# PART 4: Cross-cohort consistency (Fig. 4)
# =============================================================================
cat("\n=== PART 4: Cross-cohort consistency ===\n")

methods_for_consistency <- intersect(names(ext_stats), names(int_stats))
consistency_results <- list()
for (nm in methods_for_consistency) {
  merged <- merge(
    int_stats[[nm]][, .(feature_id, int_sig = significant, int_eff = effect_size)],
    ext_stats[[nm]][, .(feature_id, ext_sig = significant, ext_eff = effect_size)],
    by = "feature_id", all = FALSE)

  int_sig_kos <- merged[int_sig == TRUE, .N]
  both_sig <- merged[int_sig & ext_sig, .N]
  recovery <- if (int_sig_kos > 0) both_sig / int_sig_kos else NA

  sig_int <- merged[int_sig == TRUE]
  dir_agree <- if (nrow(sig_int) > 0) {
    sum(sign(sig_int$int_eff) == sign(sig_int$ext_eff), na.rm = TRUE) / nrow(sig_int)
  } else NA

  ef_cor_all <- suppressWarnings(
    cor(merged$int_eff, merged$ext_eff, method = "spearman"))
  ef_cor_intsig <- if (nrow(sig_int) > 1) {
    suppressWarnings(cor(sig_int$int_eff, sig_int$ext_eff, method = "spearman"))
  } else NA

  consistency_results[[nm]] <- data.table(
    method = nm,
    n_common_KO = nrow(merged),
    int_sig = int_sig_kos,
    ext_sig = merged[ext_sig == TRUE, .N],
    recovered_in_ext = both_sig,
    recovery_rate = round(recovery, 3),
    direction_agreement = round(dir_agree, 3),
    effect_cor_all = round(ef_cor_all, 3),
    effect_cor_intsig = round(ef_cor_intsig, 3)
  )
}
consistency_tbl <- rbindlist(consistency_results)
print(consistency_tbl)
fwrite(consistency_tbl, file.path(outdir, "cross_cohort_consistency.tsv"),
       sep = "\t")

# Save all Wilcoxon results
for (nm in names(ext_stats)) ext_stats[[nm]][, source := "external"]
for (nm in names(int_stats)) int_stats[[nm]][, source := "internal"]
all_sig_out <- rbindlist(c(ext_stats, int_stats), fill = TRUE, idcol = "method_src")
fwrite(all_sig_out, file.path(outdir, "all_wilcoxon_results.tsv"), sep = "\t")

cat("\n=== DONE ===\n")
cat(sprintf("Outputs in: %s\n", outdir))
cat("Key files:\n")
cat("  ext_sample_x_ko.ours.rel.wide.tsv\n")
cat("  cross_cohort_consistency.tsv\n")
cat("  all_wilcoxon_results.tsv\n")
cat(sprintf("Reusable model: %s\n", MODEL_FILES$xgb_full))
cat(sprintf("Carrier profile: %s\n", MODEL_FILES$carrier))
