# =============================================================================
# 03_benchmark.R
#
# Benchmark the four amplicon-based functional inference methods against the
# shotgun-derived gold standard, using MaAsLin2 differential abundance with
# Age + Sex + BMI adjustment (manuscript Fig. 2, Fig. 4, Fig. S5).
#
# Methods compared:
#   - gold        shotgun-derived KO reference (internal; the target)
#   - xgb         our two-step XGBoost framework (winner config: prev10, thr0.25)
#   - baseline    link-based inference (06_Linking output)
#   - ridge_avg   Ridge prediction averaged with the link baseline (optional)
#   - picrust2    PICRUSt2 prediction
#
# Outputs:
#   PART 1-2  MaAsLin2 for each method (internal: 4-5 methods; external: xgb, picrust2)
#   PART 3    parse results, significant-feature counts
#   PART 4    internal F1 / precision / recall / effect Spearman vs gold (Fig. 2)
#   PART 5    cross-cohort consistency, internal vs external (Fig. 4)
#
# Inputs expected (produced by earlier steps; see README):
#   02_predict_step2_ko/sample_x_ko.ours.rel.wide.tsv   (our XGB profile)
#   04_external_apply/ext_sample_x_ko.ours.rel.wide.tsv (our external profile)
#   a link-based baseline KO table (06_Linking output)  -> FILES$baseline_internal
#   PICRUSt2 internal/external                           -> FILES$picrust_*
#   gold (FILES$gold_fpkm)
#
# Metadata (FILES$meta_internal / meta_external) must contain columns
#   sample_id, group (Healthy/MASLD), Age, Sex, BMI
# (an AgeGroup column is derived from Age using AGE_CUTOFF).
#
# Run from the repository root (after 02 and 04):
#   Rscript 07_ML_Prediction/03_benchmark.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Maaslin2)
})

source(file.path("07_ML_Prediction", "config.R"))

outdir <- file.path(OUTPUT_DIR, "03_benchmark")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

N_CORES    <- PARAMS$xgb_params$nthread
AGE_CUTOFF <- 50

# ---- Metadata prep ----------------------------------------------------------
cat("[1/6] Preparing metadata (Age/Sex/BMI adjustment)...\n")

prep_meta_dt <- function(path, is_external = FALSE) {
  mt <- fread(path)
  if (is_external) setnames(mt, c("SampleID", "Group"),
                            c("sample_id", "group"), skip_absent = TRUE)
  if ("Group" %in% names(mt)) setnames(mt, "Group", "group")
  mt <- mt[group %in% c("Healthy", "MASLD")]
  if (!"AgeGroup" %in% names(mt) && "Age" %in% names(mt)) {
    mt[, AgeGroup := fifelse(Age < AGE_CUTOFF, "Younger", "Older")]
  }
  mt[!is.na(BMI) & !is.na(AgeGroup) & !is.na(Sex)]
}

meta_int <- prep_meta_dt(FILES$meta_internal, is_external = FALSE)
meta_ext <- prep_meta_dt(FILES$meta_external, is_external = TRUE)
cat(sprintf("  Internal: %d samples; External: %d samples\n",
            nrow(meta_int), nrow(meta_ext)))

prep_md <- function(mt) {
  md <- as.data.frame(mt[, .(group, AgeGroup, Sex, BMI)])
  rownames(md) <- mt$sample_id
  md
}
md_int <- prep_md(meta_int)
md_ext <- prep_md(meta_ext)

# ---- MaAsLin2 runner --------------------------------------------------------
run_maaslin <- function(wide_dt, md, cohort_tag, method_tag) {
  out_path <- file.path(outdir, sprintf("%s_%s", cohort_tag, method_tag))
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  common <- intersect(wide_dt$sample_id, rownames(md))
  df <- as.data.frame(wide_dt[sample_id %in% common])
  rownames(df) <- df$sample_id; df$sample_id <- NULL
  df <- df[, apply(df, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
  md_sub <- md[common, , drop = FALSE]
  cat(sprintf("  -- %s / %s: %d samples x %d features\n",
              cohort_tag, method_tag, nrow(df), ncol(df)))
  tryCatch(
    Maaslin2(input_data = df, input_metadata = md_sub, output = out_path,
             analysis_method = "LM", normalization = "NONE", transform = "LOG",
             fixed_effects = c("group", "AgeGroup", "Sex", "BMI"),
             reference = c("group,Healthy"),
             min_abundance = 0.0, min_prevalence = 0.1, cores = N_CORES,
             plot_heatmap = FALSE, plot_scatter = FALSE, standardize = FALSE),
    error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL })
  out_path
}

parse_maaslin <- function(out_path, cohort, method) {
  f <- file.path(out_path, "all_results.tsv")
  if (!file.exists(f)) return(NULL)
  res <- fread(f)[metadata == "group" & value == "MASLD",
                  .(feature, coef, stderr, pval, qval)]
  setnames(res, c("feature", "coef", "pval", "qval"),
           c("feature_id", "effect_adj", "pval_adj", "qval_adj"))
  res[, `:=`(cohort = cohort, method = method)][]
}

load_wide <- function(path) {
  dt <- fread(path); setnames(dt, names(dt)[1], "sample_id"); dt
}
to_rel <- function(dt) {  # row-normalise a count table to relative abundance
  m <- as.matrix(dt[, -1]); rownames(m) <- dt$sample_id
  rs <- rowSums(m); m <- sweep(m, 1, ifelse(rs > 0, rs, 1), "/")
  out <- as.data.table(m, keep.rownames = "sample_id"); out
}

# ---- PART 1: Internal methods ----------------------------------------------
cat("[2/6] === Internal MaAsLin2 ===\n")

# Gold
fpkm_gold <- fread(FILES$gold_fpkm)
setnames(fpkm_gold, c("KO", "fpkm_sum"), c("KO_ID", "abundance"), skip_absent = TRUE)
fpkm_gold[, abundance := abundance / sum(abundance, na.rm = TRUE), by = sample_id]
gold_wide <- dcast(fpkm_gold, sample_id ~ KO_ID, value.var = "abundance", fill = 0)

# Our XGBoost framework (from 02_predict_step2_ko.R)
xgb_path <- file.path(OUTPUT_DIR, "02_predict_step2_ko", "sample_x_ko.ours.rel.wide.tsv")
xgb_winner <- load_wide(xgb_path)

# Link-based baseline (from 06_Linking); path configured in config.R
baseline <- if (!is.null(FILES$baseline_internal) && file.exists(FILES$baseline_internal)) {
  load_wide(FILES$baseline_internal)
} else { cat("  [note] baseline KO table not found; skipping link-based method.\n"); NULL }

# PICRUSt2 internal
pic_int <- load_wide(FILES$picrust_internal)

out_paths <- list()
out_paths$internal_gold    <- run_maaslin(gold_wide,  md_int, "internal", "gold")
out_paths$internal_xgb     <- run_maaslin(xgb_winner, md_int, "internal", "xgb")
if (!is.null(baseline))
  out_paths$internal_baseline <- run_maaslin(baseline, md_int, "internal", "baseline")
out_paths$internal_picrust2 <- run_maaslin(pic_int,   md_int, "internal", "picrust2")

# ---- PART 2: External methods ----------------------------------------------
cat("[3/6] === External MaAsLin2 ===\n")
xgb_ext_path <- file.path(OUTPUT_DIR, "04_external_apply",
                          "ext_sample_x_ko.ours.rel.wide.tsv")
xgb_ext <- load_wide(xgb_ext_path)
pic_ext <- to_rel(load_wide(FILES$picrust_external))

out_paths$external_xgb      <- run_maaslin(xgb_ext, md_ext, "external", "xgb")
out_paths$external_picrust2 <- run_maaslin(pic_ext, md_ext, "external", "picrust2")

# ---- PART 3: Parse ----------------------------------------------------------
cat("[4/6] === Parsing results ===\n")
all_res <- rbindlist(lapply(names(out_paths), function(key) {
  parts <- strsplit(key, "_")[[1]]
  parse_maaslin(out_paths[[key]], parts[1], paste(parts[-1], collapse = "_"))
}), fill = TRUE)
fwrite(all_res, file.path(outdir, "all_maaslin_results.tsv"), sep = "\t")

# ---- PART 4: Internal F1 vs gold (Fig. 2) -----------------------------------
cat("[5/6] === Internal F1 vs gold ===\n")
gold_int <- all_res[cohort == "internal" & method == "gold"]
for (thresh_name in c("q0.25", "q0.05")) {
  thresh <- if (thresh_name == "q0.25") 0.25 else 0.05
  gold_sig <- gold_int[qval_adj < thresh, feature_id]
  f1_list <- list()
  for (m in c("xgb", "baseline", "picrust2")) {
    mr <- all_res[cohort == "internal" & method == m]
    if (nrow(mr) == 0) next
    m_sig <- mr[qval_adj < thresh, feature_id]
    all_feat <- union(gold_int$feature_id, mr$feature_id)
    TP <- length(intersect(gold_sig, m_sig)); FP <- length(setdiff(m_sig, gold_sig))
    FN <- length(setdiff(gold_sig, m_sig));   TN <- length(all_feat) - TP - FP - FN
    P <- if (TP + FP > 0) TP / (TP + FP) else NA
    R <- if (TP + FN > 0) TP / (TP + FN) else NA
    F1 <- if (!is.na(P) && !is.na(R) && P + R > 0) 2 * P * R / (P + R) else NA
    merged <- merge(gold_int[, .(feature_id, ge = effect_adj)],
                    mr[, .(feature_id, me = effect_adj)], by = "feature_id")
    spr_all <- suppressWarnings(cor(merged$ge, merged$me, method = "spearman"))
    f1_list[[m]] <- data.table(threshold = thresh_name, method = m,
                               n_gold_sig = length(gold_sig), n_method_sig = length(m_sig),
                               TP = TP, FP = FP, FN = FN,
                               precision = round(P, 3), recall = round(R, 3),
                               F1 = round(F1, 3), spearman_all = round(spr_all, 3))
  }
  f1_tbl <- rbindlist(f1_list)
  print(f1_tbl)
  fwrite(f1_tbl, file.path(outdir, sprintf("internal_f1_%s.tsv", thresh_name)), sep = "\t")
}

# ---- PART 5: Cross-cohort consistency (Fig. 4) ------------------------------
cat("[6/6] === Cross-cohort consistency ===\n")
for (thresh_name in c("q0.25", "q0.05")) {
  thresh <- if (thresh_name == "q0.25") 0.25 else 0.05
  cc_list <- list()
  for (m in c("xgb", "picrust2")) {
    ir <- all_res[cohort == "internal" & method == m]
    er <- all_res[cohort == "external" & method == m]
    if (nrow(ir) == 0 || nrow(er) == 0) next
    merged <- merge(ir[, .(feature_id, int_eff = effect_adj, int_q = qval_adj)],
                    er[, .(feature_id, ext_eff = effect_adj, ext_q = qval_adj)],
                    by = "feature_id")
    merged[, `:=`(int_sig = int_q < thresh, ext_sig = ext_q < thresh)]
    n_int <- sum(merged$int_sig, na.rm = TRUE)
    both  <- sum(merged$int_sig & merged$ext_sig, na.rm = TRUE)
    is_sub <- merged[int_sig == TRUE]
    dir_agree <- if (nrow(is_sub) > 0)
      sum(sign(is_sub$int_eff) == sign(is_sub$ext_eff), na.rm = TRUE) / nrow(is_sub) else NA
    both_sub <- merged[int_sig & ext_sig]
    spr_all  <- suppressWarnings(cor(merged$int_eff, merged$ext_eff, method = "spearman"))
    spr_both <- if (nrow(both_sub) > 1)
      suppressWarnings(cor(both_sub$int_eff, both_sub$ext_eff, method = "spearman")) else NA
    cc_list[[m]] <- data.table(threshold = thresh_name, method = m,
                               n_common = nrow(merged), int_sig = n_int,
                               recovered = both,
                               recovery_rate = round(if (n_int > 0) both / n_int else NA, 3),
                               dir_agree_intsig = round(dir_agree, 3),
                               spearman_all = round(spr_all, 3),
                               spearman_both = round(spr_both, 3))
  }
  cc_tbl <- rbindlist(cc_list)
  print(cc_tbl)
  fwrite(cc_tbl, file.path(outdir, sprintf("cross_cohort_%s.tsv", thresh_name)), sep = "\t")
}

cat(sprintf("\nDone. Outputs in: %s\n", outdir))
cat("  all_maaslin_results.tsv, internal_f1_q*.tsv, cross_cohort_q*.tsv\n")
