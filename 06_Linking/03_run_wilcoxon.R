# =============================================================================
# 03_run_wilcoxon.R
#
# Per-KO Wilcoxon differential abundance (Healthy vs MASLD) with BH correction,
# for the link-based baseline: 6 linking configs + PICRUSt2 + Gold.
# Provides the univariate counterpart to the MaAsLin2 analysis (02).
#
# Run from the repository root (after 02, which writes the gold rel table):
#   Rscript 06_Linking/03_run_wilcoxon.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

ROOT_DIR   <- Sys.getenv("KGM_ROOT", unset = getwd())
INPUT_DIR  <- file.path(ROOT_DIR, "06_Linking", "input")
linked_dir <- file.path(ROOT_DIR, "06_Linking", "output", "01_linked_tables")
out_dir    <- file.path(ROOT_DIR, "06_Linking", "output", "03_wilcoxon")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

CONFIGS <- c("pid99_d02_bin", "pid99_d02_prev", "pid99_tau1_bin",
             "pid99_tau1_prev", "pid100_tau1_bin", "pid100_tau1_prev")
ALL_ENTRIES <- c(CONFIGS, "picrust2", "gold")

meta <- fread(file.path(INPUT_DIR, "sample_meta.csv"))
setnames(meta, "sample_id", "Sample")
meta_df <- as.data.frame(meta); rownames(meta_df) <- meta_df$Sample

entry_files <- setNames(
  c(file.path(linked_dir, sprintf("sample_x_ko.%s.rel.wide.tsv", CONFIGS)),
    file.path(INPUT_DIR, "pic_rel.csv"),
    file.path(linked_dir, "sample_x_ko.gold.rel.wide.tsv")),
  ALL_ENTRIES)

wilcoxon_one <- function(ko_df, meta_use) {
  group_factor <- factor(meta_use$group, levels = c("Healthy", "MASLD"))
  res_dt <- rbindlist(lapply(colnames(ko_df), function(ko) {
    x_h <- ko_df[group_factor == "Healthy", ko]
    x_m <- ko_df[group_factor == "MASLD", ko]
    if (sum(x_h > 0) + sum(x_m > 0) < 5) return(data.table(feature = ko, pval = NA_real_))
    wt <- tryCatch(wilcox.test(x_h, x_m, exact = FALSE), error = function(e) NULL)
    data.table(feature = ko, pval = if (is.null(wt)) NA_real_ else wt$p.value)
  }))
  res_dt[, qval := p.adjust(pval, method = "BH")]
  res_dt
}

wilcoxon_all_list <- list()
for (entry in ALL_ENTRIES) {
  ko_file <- entry_files[[entry]]
  if (!file.exists(ko_file)) { cat(sprintf("SKIP %s: not found\n", entry)); next }
  cat(sprintf("\nWilcoxon for %s...\n", entry))
  ko_df <- as.data.frame(fread(ko_file)); rownames(ko_df) <- ko_df[[1]]; ko_df[[1]] <- NULL
  common <- intersect(rownames(ko_df), rownames(meta_df))
  res <- wilcoxon_one(ko_df[common, , drop = FALSE], meta_df[common, , drop = FALSE])
  res[, entry := entry]
  wilcoxon_all_list[[entry]] <- res
  cat(sprintf("  Significant KOs (q<0.05): %d / %d\n",
              sum(res$qval < 0.05, na.rm = TRUE), nrow(res)))
}
wilcoxon_all <- rbindlist(wilcoxon_all_list)
fwrite(wilcoxon_all, file.path(out_dir, "wilcoxon_all_results.tsv"), sep = "\t")
cat("\n=== 03 DONE ===\n")
