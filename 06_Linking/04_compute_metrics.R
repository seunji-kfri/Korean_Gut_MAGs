# =============================================================================
# 04_compute_metrics.R
#
# For the link-based baseline, compute (1) F1 / precision / recall vs the gold
# standard for each method, under both MaAsLin2 and Wilcoxon definitions of
# significance, and (2) per-sample Spearman similarity to the gold profile.
#
# The per-sample similarity table is also consumed by the design-choice
# analysis 07_ML_Prediction/design_choices/claim2_nonlinearity.R.
#
# Run from the repository root (after 02 and 03):
#   Rscript 06_Linking/04_compute_metrics.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

ROOT_DIR     <- Sys.getenv("KGM_ROOT", unset = getwd())
INPUT_DIR    <- file.path(ROOT_DIR, "06_Linking", "input")
linked_dir   <- file.path(ROOT_DIR, "06_Linking", "output", "01_linked_tables")
maaslin_dir  <- file.path(ROOT_DIR, "06_Linking", "output", "02_maaslin2")
wilcoxon_file<- file.path(ROOT_DIR, "06_Linking", "output", "03_wilcoxon", "wilcoxon_all_results.tsv")
out_dir      <- file.path(ROOT_DIR, "06_Linking", "output", "04_metrics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

CONFIGS <- c("pid99_d02_bin", "pid99_d02_prev", "pid99_tau1_bin",
             "pid99_tau1_prev", "pid100_tau1_bin", "pid100_tau1_prev")
ENTRIES_FOR_F1 <- c(CONFIGS, "picrust2")
DISPLAY_LABELS <- c(pid99_d02_bin = "pid99_d02_bin", pid99_d02_prev = "pid99_d02_prev",
                    pid99_tau1_bin = "pid99_bin", pid99_tau1_prev = "pid99_prev",
                    pid100_tau1_bin = "pid100_bin", pid100_tau1_prev = "pid100_prev",
                    picrust2 = "PICRUSt2")
entry_files <- setNames(
  c(file.path(linked_dir, sprintf("sample_x_ko.%s.rel.wide.tsv", CONFIGS)),
    file.path(INPUT_DIR, "pic_rel.csv")),
  ENTRIES_FOR_F1)

compute_metrics <- function(sig_kos, gold_sig_kos, all_kos) {
  truth <- all_kos %in% gold_sig_kos; pred <- all_kos %in% sig_kos
  TP <- sum(truth & pred); FP <- sum(!truth & pred)
  FN <- sum(truth & !pred); TN <- sum(!truth & !pred)
  P <- if (TP + FP > 0) TP / (TP + FP) else 0
  R <- if (TP + FN > 0) TP / (TP + FN) else 0
  data.table(TP = TP, FP = FP, FN = FN, TN = TN, precision = P, recall = R,
             F1 = if (P + R > 0) 2 * P * R / (P + R) else 0)
}

# ---- MaAsLin2 F1 ------------------------------------------------------------
cat("=== MaAsLin2 F1 ===\n")
gold_maaslin <- fread(file.path(maaslin_dir, "gold", "all_results.tsv"))
gold_maaslin <- gold_maaslin[metadata == "group" & value == "MASLD"]
gold_sig_maaslin <- gold_maaslin[qval < 0.05, feature]

maaslin_f1 <- rbindlist(lapply(ENTRIES_FOR_F1, function(entry) {
  f <- file.path(maaslin_dir, entry, "all_results.tsv")
  if (!file.exists(f)) return(NULL)
  res <- fread(f)[metadata == "group" & value == "MASLD"]
  m <- compute_metrics(res[qval < 0.05, feature], gold_sig_maaslin, res$feature)
  m[, c("entry", "test", "label") := .(entry, "MaAsLin2", DISPLAY_LABELS[entry])][]
}))

# ---- Wilcoxon F1 ------------------------------------------------------------
cat("=== Wilcoxon F1 ===\n")
wilcoxon_all <- fread(wilcoxon_file)
gold_wilcoxon <- wilcoxon_all[entry == "gold"]
gold_sig_wilcoxon <- gold_wilcoxon[!is.na(qval) & qval < 0.05, feature]
wilcox_f1 <- rbindlist(lapply(ENTRIES_FOR_F1, function(e) {
  res <- wilcoxon_all[entry == e]
  if (nrow(res) == 0) return(NULL)
  m <- compute_metrics(res[!is.na(qval) & qval < 0.05, feature], gold_sig_wilcoxon, res$feature)
  m[, c("entry", "test", "label") := .(e, "Wilcoxon", DISPLAY_LABELS[e])][]
}))

F1_summary <- rbind(maaslin_f1, wilcox_f1)
fwrite(F1_summary, file.path(out_dir, "F1_summary.tsv"), sep = "\t")
print(F1_summary[, .(label, test, F1 = round(F1, 3),
                     precision = round(precision, 3), recall = round(recall, 3))])

# ---- Per-sample Spearman similarity vs gold ---------------------------------
cat("\n=== Per-sample Spearman similarity ===\n")
gold_wide <- fread(file.path(linked_dir, "sample_x_ko.gold.rel.wide.tsv"))
gold_mat <- as.matrix(gold_wide[, -1]); rownames(gold_mat) <- gold_wide$Sample

compute_sample_spr <- function(ko_file, label) {
  if (!file.exists(ko_file)) return(NULL)
  wide <- fread(ko_file); setnames(wide, names(wide)[1], "Sample")
  common <- intersect(wide$Sample, rownames(gold_mat))
  if (length(common) < 10) return(NULL)
  wide <- wide[Sample %in% common]; setorder(wide, Sample)
  method_mat <- as.matrix(wide[, -1]); rownames(method_mat) <- wide$Sample
  common_kos <- intersect(colnames(method_mat), colnames(gold_mat))
  if (length(common_kos) < 5) return(NULL)
  m <- method_mat[common, common_kos]; g <- gold_mat[common, common_kos]
  spr <- sapply(seq_len(nrow(g)), function(i)
    suppressWarnings(cor(g[i, ], m[i, ], method = "spearman")))
  data.table(sample_id = common, method = label, spearman = spr)
}
per_sample_spr <- rbindlist(lapply(ENTRIES_FOR_F1, function(e)
  compute_sample_spr(entry_files[[e]], DISPLAY_LABELS[e])), fill = TRUE)
fwrite(per_sample_spr, file.path(out_dir, "per_sample_similarity.tsv"), sep = "\t")
print(per_sample_spr[, .(median_spr = round(median(spearman, na.rm = TRUE), 3)), by = method])

cat("\n=== 04 DONE ===\n")
