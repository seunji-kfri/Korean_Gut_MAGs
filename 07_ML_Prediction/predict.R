#!/usr/bin/env Rscript
# =============================================================================
# predict.R
#
# Standalone prediction entry point.
#
# Given a user-supplied 16S rRNA ASV relative-abundance table, apply the
# pre-trained XGBoost models (full_xgb_models.rds, trained on the 224 paired
# Korean internal cohort) to infer dRep95 cluster abundances, then convert them
# to a sample x KO functional profile via a binary carrier matrix.
#
# This is the external-facing tool described in the manuscript: a user with
# amplicon data alone can reproduce the functional inference without retraining.
#
# -----------------------------------------------------------------------------
# USAGE
#   Rscript predict.R \
#       --asv      path/to/your_asv_relabund.csv \
#       --models   07_ML_Prediction/models/full_xgb_models.rds \
#       --carrier  07_ML_Prediction/models/ko_cluster_table.tsv \
#       --out      predicted_sample_x_ko.tsv \
#       [--prev    prev10]            # which trained model set (prev10 | prev20)
#       [--thr     0.25]              # carrier-fraction threshold
#
# If --carrier is omitted, the script will look for either
#   07_ML_Prediction/models/ko_cluster_table.tsv   (cluster, KO, carrier_fraction)
#   07_ML_Prediction/models/carrier_profile.rds    (cluster x KO matrix)
#
# INPUT ASV TABLE FORMAT
#   - CSV (or TSV) with the FIRST column = ASV id, remaining columns = samples.
#   - Values = relative abundance (or counts; they are CLR-transformed per
#     sample, so either works as long as it is consistent).
#   - ASV ids must match those used in training. ASVs missing from the user
#     table are zero-padded; ASVs not seen in training are ignored.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(xgboost)
})

# ---- Minimal argument parser (no external dependency) -----------------------
parse_args <- function(args) {
  out <- list(prev = "prev10", thr = 0.25,
              models = "07_ML_Prediction/models/full_xgb_models.rds",
              carrier = NULL, asv = NULL, out = "predicted_sample_x_ko.tsv")
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    val <- args[i + 1]
    out[[key]] <- val
    i <- i + 2
  }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))

if (is.null(args$asv)) {
  stop("Missing --asv. See header of predict.R for usage.", call. = FALSE)
}
# Normalise prev key (accept "10" or "prev10")
if (!grepl("^prev", args$prev)) args$prev <- paste0("prev", args$prev)
thr <- as.numeric(args$thr)

cat("predict.R\n")
cat(sprintf("  ASV table : %s\n", args$asv))
cat(sprintf("  Models    : %s (%s)\n", args$models, args$prev))
cat(sprintf("  Carrier   : %s\n", ifelse(is.null(args$carrier), "(auto-detect)", args$carrier)))
cat(sprintf("  Threshold : %.2f\n", thr))
cat(sprintf("  Output    : %s\n\n", args$out))

# ---- 1. Load trained models -------------------------------------------------
bundle <- readRDS(args$models)
if (!all(c("models", "Y_info", "feature_names") %in% names(bundle))) {
  stop("Model file does not have expected structure ",
       "(models / Y_info / feature_names).", call. = FALSE)
}
if (!args$prev %in% names(bundle$models)) {
  stop(sprintf("Model set '%s' not found. Available: %s",
               args$prev, paste(names(bundle$models), collapse = ", ")),
       call. = FALSE)
}
models <- bundle$models[[args$prev]]
pseudo_mag <- bundle$Y_info[[args$prev]]$pseudo
training_asv <- bundle$feature_names
cat(sprintf("  Loaded %d cluster models; %d training ASV features.\n",
            length(models), length(training_asv)))

# ---- 2. Load and align the user ASV table -----------------------------------
read_table_any <- function(path) {
  sep <- if (grepl("\\.tsv$|\\.txt$", path)) "\t" else ","
  fread(path, sep = sep)
}
asv_raw <- read_table_any(args$asv)
setnames(asv_raw, names(asv_raw)[1], "ASV")

asv_mat <- as.matrix(asv_raw[, -1])
rownames(asv_mat) <- asv_raw$ASV
samples <- colnames(asv_mat)

in_train <- intersect(rownames(asv_mat), training_asv)
missing  <- setdiff(training_asv, rownames(asv_mat))
cat(sprintf("  User table: %d ASVs x %d samples\n",
            nrow(asv_mat), length(samples)))
cat(sprintf("  Overlap with training features: %d / %d (%d zero-padded)\n",
            length(in_train), length(training_asv), length(missing)))

# Align to training feature space (samples x training_asv), zero-pad missing
X_aligned <- matrix(0, nrow = length(samples), ncol = length(training_asv),
                    dimnames = list(samples, training_asv))
X_aligned[, in_train] <- t(asv_mat[in_train, , drop = FALSE])

# ---- 3. CLR transform (per sample) ------------------------------------------
# Use a pseudo-count of half the global minimum non-zero value, matching the
# training-time transform.
nz <- X_aligned[X_aligned > 0]
pseudo_asv <- if (length(nz)) min(nz) / 2 else 1e-6
X_clr <- t(apply(X_aligned + pseudo_asv, 1,
                 function(x) log(x) - mean(log(x))))

# ---- 4. Predict cluster abundance -------------------------------------------
predict_with_models <- function(models, X_te) {
  out <- matrix(0, nrow = nrow(X_te), ncol = length(models),
                dimnames = list(rownames(X_te), names(models)))
  for (j in seq_along(models)) {
    m <- models[[j]]
    if (!is.null(m$constant)) out[, j] <- m$constant
    else out[, j] <- predict(m, X_te)
  }
  out
}
Y_log <- predict_with_models(models, X_clr)

# log -> linear, clip negatives, renormalise per sample to sum 100
M <- exp(Y_log) - pseudo_mag
M[M < 0] <- 0
rs <- rowSums(M)
cluster_abund <- sweep(M, 1, ifelse(rs > 0, rs / 100, 1), "/")
cat(sprintf("  Predicted cluster abundance: %d samples x %d clusters\n",
            nrow(cluster_abund), ncol(cluster_abund)))

# ---- 5. Load / build carrier matrix (cluster x KO) --------------------------
clusters_pred <- colnames(cluster_abund)

load_carrier <- function(path, clusters_pred) {
  if (grepl("\\.rds$", path)) {
    carrier <- readRDS(path)                      # expect cluster x KO matrix
    return(as.matrix(carrier))
  }
  # Long format: columns cluster, KO, carrier_fraction (names auto-detected)
  dt <- fread(path)
  nm <- tolower(names(dt))
  ccol <- names(dt)[which(nm %in% c("cluster", "drep95_cluster", "secondary_cluster"))][1]
  kcol <- names(dt)[which(nm %in% c("ko", "ko_id", "kegg", "feature_id"))][1]
  vcol <- names(dt)[which(nm %in% c("carrier_frac", "carrier_fraction", "carrier",
                                    "fraction", "value"))][1]
  if (is.na(ccol) || is.na(kcol)) {
    stop("Could not detect cluster/KO columns in carrier table: ",
         paste(names(dt), collapse = ", "), call. = FALSE)
  }
  if (is.na(vcol)) { dt[, .__v := 1]; vcol <- ".__v" }  # presence-only table
  wide <- dcast(dt, get(ccol) ~ get(kcol), value.var = vcol, fill = 0)
  m <- as.matrix(wide[, -1]); rownames(m) <- wide[[1]]
  m
}

carrier_path <- args$carrier
if (is.null(carrier_path)) {
  cand <- c("07_ML_Prediction/models/ko_cluster_table.tsv",
            "07_ML_Prediction/models/carrier_profile.rds")
  carrier_path <- cand[file.exists(cand)][1]
  if (is.na(carrier_path)) {
    stop("No carrier file given and none found. Provide --carrier ",
         "(ko_cluster_table.tsv or carrier_profile.rds).", call. = FALSE)
  }
  cat(sprintf("  Auto-detected carrier: %s\n", carrier_path))
}
carrier <- load_carrier(carrier_path, clusters_pred)

# Keep only clusters present in the prediction; align order
common_clusters <- intersect(clusters_pred, rownames(carrier))
if (length(common_clusters) == 0) {
  stop("No overlap between predicted clusters and carrier rows. ",
       "Check that the carrier table uses the same cluster IDs.", call. = FALSE)
}
cat(sprintf("  Carrier matrix: %d clusters x %d KOs (%d clusters matched)\n",
            nrow(carrier), ncol(carrier), length(common_clusters)))

carrier <- carrier[common_clusters, , drop = FALSE]
cluster_abund <- cluster_abund[, common_clusters, drop = FALSE]

# Apply carrier threshold (winner config = 0.25)
carrier[carrier < thr] <- 0

# ---- 6. cluster abundance x carrier -> sample x KO --------------------------
sko <- cluster_abund %*% carrier
rs <- rowSums(sko)
sko <- sweep(sko, 1, ifelse(rs > 0, rs, 1), "/")   # relative abundance per sample

out_dt <- as.data.table(sko, keep.rownames = "sample_id")
sep_out <- if (grepl("\\.csv$", args$out)) "," else "\t"
fwrite(out_dt, args$out, sep = sep_out)

cat(sprintf("\nDone. Predicted KO profile: %d samples x %d KOs\n",
            nrow(out_dt), ncol(out_dt) - 1))
cat(sprintf("Saved -> %s\n", args$out))
