# =============================================================================
# claim3_cluster_resolution.R
#
# Design-choice analysis (manuscript Fig. 3A-B): species-level dRep95 clustering
# matches per-sample amplicon richness better than the finer dRep99 resolution.
# Compares F1 vs gold across cluster catalogs:
#   - dRep95 prev>=10 (350 clusters; the main/winner catalog)
#   - dRep99 prev>=10 (dRep99A, ~130 clusters)
#   - dRep99 prev>=5  (dRep99B, ~538 clusters)
# for both XGBoost and Ridge.
#
# Inputs:
#   - MaAsLin2 results for dRep95 XGB/Ridge (from 03_benchmark output)
#   - MaAsLin2 results for dRep99 variants. These require building the dRep99
#     KO tables first (re-cluster at 99% ANI, retrain, aggregate to KO). The
#     pre-built dRep99 KO tables are expected under design_choices input as:
#       sample_x_ko.drep99A_xgb.rel.wide.tsv,  sample_x_ko.drep99A_ridge...
#       sample_x_ko.drep99B_xgb.rel.wide.tsv,  sample_x_ko.drep99B_ridge...
#   - gold MaAsLin2 results (03_benchmark/internal_gold/all_results.tsv)
#
# Run from the repository root (after 03_benchmark.R):
#   Rscript 07_ML_Prediction/design_choices/claim3_cluster_resolution.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Maaslin2)
})

source(file.path("07_ML_Prediction", "config.R"))

bench_dir <- file.path(OUTPUT_DIR, "03_benchmark")
dc_input  <- file.path(INPUT_DIR, "design_choices")   # pre-built dRep99 KO tables
out_dir   <- file.path(OUTPUT_DIR, "design_choices", "claim3_cluster_resolution")
maaslin_dir <- file.path(out_dir, "maaslin2")
dir.create(maaslin_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. MaAsLin2 for dRep99 variants ---------------------------------------
cat("[1/3] Running MaAsLin2 for dRep99 variants...\n")
meta <- fread(FILES$meta_internal)
if ("Group" %in% names(meta)) setnames(meta, "Group", "group")
meta_df <- as.data.frame(meta); rownames(meta_df) <- meta_df$sample_id

load_ko_wide <- function(f) {
  if (!file.exists(f)) return(NULL)
  df <- as.data.frame(fread(f)); rownames(df) <- df[[1]]; df[[1]] <- NULL; df
}

run_maaslin <- function(ko_file, out_subdir) {
  if (file.exists(file.path(out_subdir, "all_results.tsv"))) return(invisible())
  ko_df <- load_ko_wide(ko_file)
  if (is.null(ko_df)) { cat(sprintf("  SKIP (not found): %s\n", ko_file)); return(invisible()) }
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

drep99_entries <- c(
  drep99A_xgb = "sample_x_ko.drep99A_xgb.rel.wide.tsv",
  drep99A_ridge = "sample_x_ko.drep99A_ridge.rel.wide.tsv",
  drep99B_xgb = "sample_x_ko.drep99B_xgb.rel.wide.tsv",
  drep99B_ridge = "sample_x_ko.drep99B_ridge.rel.wide.tsv")
for (e in names(drep99_entries)) {
  run_maaslin(file.path(dc_input, drep99_entries[e]), file.path(maaslin_dir, e))
}

# ---- 2. F1 for all variants -------------------------------------------------
cat("[2/3] Computing F1 vs gold...\n")
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
             F1 = if (P + R > 0) 2 * P * R / (P + R) else 0, TP = TP, FP = FP, FN = FN)
}

entries_drep <- list(
  drep95_xgb   = file.path(bench_dir, "internal_xgb"),
  drep95_ridge = file.path(bench_dir, "internal_ridge"),       # if Ridge benchmark was run
  drep99A_xgb  = file.path(maaslin_dir, "drep99A_xgb"),
  drep99A_ridge = file.path(maaslin_dir, "drep99A_ridge"),
  drep99B_xgb  = file.path(maaslin_dir, "drep99B_xgb"),
  drep99B_ridge = file.path(maaslin_dir, "drep99B_ridge"))
variant_info <- data.table(
  entry = names(entries_drep),
  catalog = c("dRep95", "dRep95", "dRep99A", "dRep99A", "dRep99B", "dRep99B"),
  n_clusters = c(350, 350, 130, 130, 538, 538),
  model = c("XGBoost", "Ridge", "XGBoost", "Ridge", "XGBoost", "Ridge"))

f1_list <- list()
for (e in names(entries_drep)) {
  f <- file.path(entries_drep[[e]], "all_results.tsv")
  if (!file.exists(f)) { cat(sprintf("  SKIP %s\n", e)); next }
  r <- fread(f)[metadata == "group" & value == "MASLD"]
  r[, feature := sub("^ko[.:]", "", feature)]
  m <- compute_f1(r[qval < 0.05, feature], gold_sig, r$feature)
  m[, `:=`(entry = e, n_sig = r[qval < 0.05, .N])]
  f1_list[[e]] <- m
}
f1_drep <- merge(rbindlist(f1_list), variant_info, by = "entry")
print(f1_drep[, .(catalog, n_clusters, model, F1 = round(F1, 3),
                  precision = round(precision, 3), recall = round(recall, 3), n_sig)])
fwrite(f1_drep, file.path(out_dir, "claim3_f1_drep_comparison.tsv"), sep = "\t")

# ---- 3. Summary -------------------------------------------------------------
cat("[3/3] XGBoost across dRep variants:\n")
print(f1_drep[model == "XGBoost", .(catalog, n_clusters, F1 = round(F1, 3))])
cat("\nDone.\n")
