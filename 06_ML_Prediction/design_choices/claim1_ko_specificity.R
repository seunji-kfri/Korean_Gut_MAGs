# =============================================================================
# claim1_ko_specificity.R
#
# Design-choice analysis (manuscript Fig. 3G): the machine-learning advantage
# over link-based inference and PICRUSt2 is largest for moderate-specificity
# KOs. KOs are grouped by the number of MAG bins carrying them (Specialist
# <=10, Moderate 11-100, Generalist >100), and F1 vs the gold reference is
# computed per group for each method.
#
# Inputs (produced by earlier steps):
#   - bin x KO FPKM table (config FILES$bin_ko_fpkm)
#   - MaAsLin2 results per method (from 03_benchmark.R output:
#     03_benchmark/internal_<method>/all_results.tsv)
#
# Run from the repository root (after 03_benchmark.R):
#   Rscript 07_ML_Prediction/design_choices/claim1_ko_specificity.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("07_ML_Prediction", "config.R"))

bench_dir <- file.path(OUTPUT_DIR, "03_benchmark")
out_dir   <- file.path(OUTPUT_DIR, "design_choices", "claim1_ko_specificity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Count bin carriers per KO -> specificity group ----------------------
cat("[1/4] Loading bin x KO FPKM matrix...\n")
bin_ko <- fread(FILES$bin_ko_fpkm)
setnames(bin_ko, names(bin_ko)[1], "BIN_ID")
bin_ko_long <- melt(bin_ko, id.vars = "BIN_ID", variable.name = "KO_ID",
                    value.name = "FPKM")
bin_ko_long <- bin_ko_long[FPKM > 0]
bin_ko_long[, KO_ID := sub("^ko[.:]", "", as.character(KO_ID))]

ko_carrier <- bin_ko_long[, .(n_bin_carriers = uniqueN(BIN_ID)), by = KO_ID]
ko_carrier[, spec_group := fcase(
  n_bin_carriers <= 10,  "Specialist",
  n_bin_carriers <= 100, "Moderate",
  default = "Generalist")]
ko_carrier[, spec_group := factor(spec_group,
           levels = c("Specialist", "Moderate", "Generalist"))]
cat("  KO specificity distribution:\n")
print(ko_carrier[, .N, by = spec_group][order(spec_group)])

# ---- 2. Load MaAsLin2 results for each method -------------------------------
cat("[2/4] Loading MaAsLin2 results (from 03_benchmark)...\n")
ENTRIES <- c(xgb = "XGBoost", baseline = "Link-based",
             picrust2 = "PICRUSt2", gold = "Gold")

maaslin_list <- list()
for (e in names(ENTRIES)) {
  f <- file.path(bench_dir, sprintf("internal_%s", e), "all_results.tsv")
  if (!file.exists(f)) { cat(sprintf("  SKIP %s: not found\n", e)); next }
  res <- fread(f)[metadata == "group" & value == "MASLD"]
  res[, KO_ID := sub("^ko[.:]", "", feature)]
  res[, method := ENTRIES[e]]
  maaslin_list[[e]] <- res[, .(method, KO_ID, coef, qval)]
}
maaslin_all <- rbindlist(maaslin_list)
maaslin_all <- merge(maaslin_all, ko_carrier[, .(KO_ID, spec_group)],
                     by = "KO_ID", all.x = TRUE)
maaslin_all <- maaslin_all[!is.na(spec_group)]

# ---- 3. F1 per (method, specificity group) ----------------------------------
cat("[3/4] Computing F1 per specificity group...\n")
gold_sig_by_grp <- maaslin_all[method == "Gold" & qval < 0.05, .(KO_ID, spec_group)]

compute_metrics <- function(sig_kos, gold_sig_kos, all_kos) {
  truth <- all_kos %in% gold_sig_kos; pred <- all_kos %in% sig_kos
  TP <- sum(truth & pred); FP <- sum(!truth & pred); FN <- sum(truth & !pred)
  P <- if (TP + FP > 0) TP / (TP + FP) else 0
  R <- if (TP + FN > 0) TP / (TP + FN) else 0
  F1 <- if (P + R > 0) 2 * P * R / (P + R) else 0
  list(precision = P, recall = R, F1 = F1, TP = TP, FP = FP, FN = FN)
}

f1_list <- list()
for (m in c("XGBoost", "Link-based", "PICRUSt2")) {
  for (grp in c("Specialist", "Moderate", "Generalist")) {
    grp_kos <- maaslin_all[method == m & spec_group == grp, KO_ID]
    sig_kos <- maaslin_all[method == m & spec_group == grp & qval < 0.05, KO_ID]
    grp_gold <- gold_sig_by_grp[spec_group == grp, KO_ID]
    if (length(grp_kos) == 0) {
      f1_list[[paste(m, grp)]] <- data.table(method = m, spec_group = grp,
        n_KOs = 0, n_gold_sig = length(grp_gold),
        precision = NA, recall = NA, F1 = NA); next
    }
    mt <- compute_metrics(sig_kos, grp_gold, grp_kos)
    f1_list[[paste(m, grp)]] <- data.table(method = m, spec_group = grp,
      n_KOs = length(grp_kos), n_gold_sig = length(grp_gold),
      precision = mt$precision, recall = mt$recall, F1 = mt$F1)
  }
}
f1_by_grp <- rbindlist(f1_list)
print(f1_by_grp[, .(method, spec_group, n_KOs, n_gold_sig,
                    F1 = round(F1, 3), precision = round(precision, 3),
                    recall = round(recall, 3))])
fwrite(f1_by_grp, file.path(out_dir, "f1_by_specificity.tsv"), sep = "\t")
fwrite(ko_carrier, file.path(out_dir, "ko_specificity_groups.tsv"), sep = "\t")

# ---- 4. ML advantage summary ------------------------------------------------
cat("[4/4] ML advantage by group...\n")
f1_wide <- dcast(f1_by_grp[, .(method, spec_group, F1)],
                 spec_group ~ method, value.var = "F1")
print(f1_wide)
fwrite(f1_wide, file.path(out_dir, "ml_advantage_summary.tsv"), sep = "\t")
cat("\nDone.\n")
