# =============================================================================
# 05_biology.R
#
# Biological interpretation of the MASLD functional signature recovered by the
# XGBoost framework (manuscript Fig. 5). Combines:
#   - top differentially abundant KOs (MASLD-up / MASLD-down)
#   - cross-cohort consistent KOs (internal & external, same direction)
#   - KEGG pathway over-representation (hypergeometric test, Levels 1-3)
#   - targeted check of pathways with established roles in MASLD
#     (LPS/endotoxin, bile acid, SCFA, tryptophan, bacterial secretion)
#   - the same enrichment applied to the external cohort
#
# Inputs (produced by 03_benchmark.R MaAsLin2 runs):
#   03_benchmark/internal_gold/all_results.tsv
#   03_benchmark/internal_xgb/all_results.tsv
#   03_benchmark/external_xgb/all_results.tsv
#   KEGG hierarchy table (config FILES$ko_pathway_hierarchy)
#
# Run from the repository root (after 03_benchmark.R):
#   Rscript 07_ML_Prediction/05_biology.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

source(file.path("07_ML_Prediction", "config.R"))

bench_dir <- file.path(OUTPUT_DIR, "03_benchmark")
out_dir   <- file.path(OUTPUT_DIR, "05_biology")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Q_SIG <- 0.05

# ---- 1. Load MaAsLin2 results -----------------------------------------------
cat("[1/5] Loading MaAsLin2 results...\n")
parse_res <- function(folder) {
  x <- fread(file.path(folder, "all_results.tsv"))
  x <- x[metadata == "group" & value == "MASLD"]
  setnames(x, c("feature", "coef", "qval"),
           c("feature_id", "effect", "qval"), skip_absent = TRUE)
  x[, feature_id := sub("^ko[.:]", "", feature_id)]
  x
}
gold_int <- parse_res(file.path(bench_dir, "internal_gold"))
xgb_int  <- parse_res(file.path(bench_dir, "internal_xgb"))
xgb_ext  <- parse_res(file.path(bench_dir, "external_xgb"))
cat(sprintf("  gold int=%d sig; xgb int=%d sig; xgb ext=%d sig\n",
            sum(gold_int$qval < Q_SIG), sum(xgb_int$qval < Q_SIG),
            sum(xgb_ext$qval < Q_SIG)))

# ---- 2. KEGG hierarchy ------------------------------------------------------
cat("[2/5] Loading KEGG pathway hierarchy...\n")
ko_path <- fread(FILES$ko_pathway_hierarchy)
# Expected columns: feature (KO), Pathway, Level1, Level2, Level3

# ---- 3. Top KOs + cross-cohort consistency ----------------------------------
cat("[3/5] Top KOs and cross-cohort consistent KOs...\n")
gold_sig <- gold_int[qval < Q_SIG][order(-abs(effect))]
gold_sig_path <- merge(gold_sig,
  ko_path[, .(feature_id = feature, Level1, Level2, Level3, Pathway)],
  by = "feature_id", all.x = TRUE)

top_up   <- gold_sig_path[effect > 0][order(-effect)][!duplicated(feature_id)]
top_down <- gold_sig_path[effect < 0][order(effect)][!duplicated(feature_id)]
fwrite(top_up[, .(feature_id, effect, qval, Pathway, Level3, Level2, Level1)],
       file.path(out_dir, "top_ko_masld_up.tsv"), sep = "\t")
fwrite(top_down[, .(feature_id, effect, qval, Pathway, Level3, Level2, Level1)],
       file.path(out_dir, "top_ko_masld_down.tsv"), sep = "\t")
cat(sprintf("  MASLD-up: %d, MASLD-down: %d\n", nrow(top_up), nrow(top_down)))

cc <- merge(xgb_int[, .(feature_id, int_effect = effect, int_qval = qval)],
            xgb_ext[, .(feature_id, ext_effect = effect, ext_qval = qval)],
            by = "feature_id")
cc[, `:=`(int_sig = int_qval < Q_SIG, ext_sig = ext_qval < Q_SIG,
          same_direction = sign(int_effect) == sign(ext_effect))]
both_consistent <- cc[int_sig & ext_sig & same_direction]
both_consistent <- merge(both_consistent,
  ko_path[, .(feature_id = feature, Level1, Level2, Level3, Pathway)],
  by = "feature_id", all.x = TRUE)[!duplicated(feature_id)]
both_consistent[, mean_effect := (int_effect + ext_effect) / 2]
setorder(both_consistent, -abs(mean_effect))
fwrite(both_consistent[, .(feature_id,
         int_effect = round(int_effect, 3), ext_effect = round(ext_effect, 3),
         int_qval = signif(int_qval, 3), ext_qval = signif(ext_qval, 3),
         direction = ifelse(mean_effect > 0, "up (MASLD)", "down (MASLD)"),
         Pathway, Level3, Level2, Level1)],
       file.path(out_dir, "cross_cohort_consistent_ko.tsv"), sep = "\t")
cat(sprintf("  cross-cohort consistent KOs: %d\n", nrow(both_consistent)))

# ---- 4. KEGG over-representation (hypergeometric) ---------------------------
cat("[4/5] Pathway over-representation analysis...\n")
run_ora <- function(stats_dt, level_col, tag) {
  sig_ids <- unique(stats_dt[qval < Q_SIG, feature_id])
  universe <- unique(stats_dt$feature_id)
  kp <- unique(ko_path[!is.na(get(level_col)),
                       .(feature_id = feature, pathway = get(level_col))])
  kp <- kp[feature_id %in% universe]
  res <- rbindlist(lapply(unique(kp$pathway), function(pth) {
    pth_kos <- kp[pathway == pth, feature_id]
    n_total <- length(pth_kos); n_sig <- sum(pth_kos %in% sig_ids)
    if (n_total < 5 || n_sig < 2) return(NULL)
    pval <- phyper(n_sig - 1, n_total, length(universe) - n_total,
                   length(sig_ids), lower.tail = FALSE)
    er <- (n_sig / length(sig_ids)) / (n_total / length(universe))
    avg_eff <- mean(stats_dt[feature_id %in% intersect(pth_kos, sig_ids), effect],
                    na.rm = TRUE)
    data.table(pathway = pth, n_total = n_total, n_sig = n_sig,
               enrichment = round(er, 2), avg_effect = round(avg_eff, 3),
               direction = ifelse(avg_eff > 0, "up (MASLD)", "down (MASLD)"),
               pvalue = pval)
  }), fill = TRUE)
  if (nrow(res)) { res[, qvalue := p.adjust(pvalue, method = "BH")]; setorder(res, pvalue) }
  fwrite(res, file.path(out_dir, sprintf("pathway_enrichment_%s_%s.tsv",
                                         tag, tolower(level_col))), sep = "\t")
  res
}
for (lv in c("Level1", "Level2", "Level3")) {
  run_ora(gold_int, lv, "internal_gold")
  run_ora(xgb_ext,  lv, "external_xgb")
}

# ---- 5. Known MASLD pathways ------------------------------------------------
cat("[5/5] Known MASLD-related pathways...\n")
known_masld <- data.table(
  category = c("LPS / endotoxin", "LPS / endotoxin", "Bile acid", "Bile acid",
               "SCFA (butyrate)", "SCFA (propionate)", "Amino acid / BCAA",
               "Tryptophan metabolism", "Bacterial secretion"),
  pathway_id = c("ko00540", "ko00541", "ko00121", "ko00120", "ko00650",
                 "ko00640", "ko00280", "ko00380", "ko03070"))

check_known <- function(stats_dt, tag) {
  sig_ids <- unique(stats_dt[qval < Q_SIG, feature_id])
  universe <- unique(stats_dt$feature_id)
  res <- rbindlist(lapply(seq_len(nrow(known_masld)), function(i) {
    pth_id <- known_masld$pathway_id[i]
    pth_kos <- intersect(ko_path[Pathway == pth_id, unique(feature)], universe)
    if (length(pth_kos) == 0)
      return(data.table(category = known_masld$category[i], pathway_id = pth_id,
                        pathway_name = "not in data", n_total = 0, n_sig = 0,
                        frac_sig = NA, avg_effect = NA))
    sig_kos <- intersect(pth_kos, sig_ids)
    data.table(category = known_masld$category[i], pathway_id = pth_id,
               pathway_name = ko_path[Pathway == pth_id, Level3][1],
               n_total = length(pth_kos), n_sig = length(sig_kos),
               frac_sig = round(length(sig_kos) / length(pth_kos), 3),
               avg_effect = if (length(sig_kos))
                 round(mean(stats_dt[feature_id %in% sig_kos, effect], na.rm = TRUE), 3)
                 else NA_real_)
  }), fill = TRUE)
  fwrite(res, file.path(out_dir, sprintf("masld_known_pathways_%s.tsv", tag)), sep = "\t")
  res
}
cat("\n  Internal (gold):\n"); print(check_known(gold_int, "internal_gold"))
cat("\n  External (xgb):\n");  print(check_known(xgb_ext, "external_xgb"))

cat(sprintf("\nDone. Outputs in: %s\n", out_dir))
