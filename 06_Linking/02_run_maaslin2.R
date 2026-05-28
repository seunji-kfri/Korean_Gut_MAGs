# =============================================================================
# 02_run_maaslin2.R
#
# MaAsLin2 differential abundance (Healthy vs MASLD) for the link-based
# baseline: 6 linking configs + PICRUSt2 + Gold (shotgun truth).
# Also writes the gold relative-abundance table used by 03/04.
#
# Inputs:
#   06_Linking/output/01_linked_tables/sample_x_ko.<config>.rel.wide.tsv
#   06_Linking/input/sample_meta.csv     (sample_id, group)
#   06_Linking/input/fpkm_gold_long.csv  (shotgun gold; KO, fpkm_sum)
#   06_Linking/input/pic_rel.csv         (PICRUSt2 relative abundance)
#
# Run from the repository root (after 01):
#   Rscript 06_Linking/02_run_maaslin2.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Maaslin2)
})

ROOT_DIR   <- Sys.getenv("KGM_ROOT", unset = getwd())
INPUT_DIR  <- file.path(ROOT_DIR, "06_Linking", "input")
linked_dir <- file.path(ROOT_DIR, "06_Linking", "output", "01_linked_tables")
out_base   <- file.path(ROOT_DIR, "06_Linking", "output", "02_maaslin2")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

CONFIGS <- c("pid99_d02_bin", "pid99_d02_prev", "pid99_tau1_bin",
             "pid99_tau1_prev", "pid100_tau1_bin", "pid100_tau1_prev")

meta <- fread(file.path(INPUT_DIR, "sample_meta.csv"))
setnames(meta, "sample_id", "Sample")
meta_df <- as.data.frame(meta); rownames(meta_df) <- meta_df$Sample

load_ko_wide <- function(f) {
  if (!file.exists(f)) return(NULL)
  df <- as.data.frame(fread(f)); rownames(df) <- df[[1]]; df[[1]] <- NULL; df
}

# Build & save gold relative-abundance matrix
gold_long <- fread(file.path(INPUT_DIR, "fpkm_gold_long.csv"))
setnames(gold_long, c("KO", "fpkm_sum"), c("KO_ID", "abundance"), skip_absent = TRUE)
gold_wide <- dcast(gold_long, sample_id ~ KO_ID, value.var = "abundance", fill = 0)
setnames(gold_wide, "sample_id", "Sample")
gold_df <- as.data.frame(gold_wide); rownames(gold_df) <- gold_df$Sample; gold_df$Sample <- NULL
gold_df <- as.data.frame(t(apply(gold_df, 1, function(x) if (sum(x) > 0) x / sum(x) else x)))
fwrite(data.table(Sample = rownames(gold_df), gold_df),
       file.path(linked_dir, "sample_x_ko.gold.rel.wide.tsv"), sep = "\t")

run_maaslin <- function(entry_name, ko_df, meta_df, out_dir) {
  cat(sprintf("\n--- MaAsLin2: %s ---\n", entry_name))
  common <- intersect(rownames(ko_df), rownames(meta_df))
  if (length(common) < 10) { cat("  SKIP: too few samples\n"); return(NULL) }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    Maaslin2(input_data = ko_df[common, , drop = FALSE],
             input_metadata = meta_df[common, , drop = FALSE], output = out_dir,
             fixed_effects = c("group"), reference = c("group,Healthy"),
             normalization = "NONE", transform = "LOG", analysis_method = "LM",
             min_prevalence = 0.1, max_significance = 1.0,
             plot_heatmap = FALSE, plot_scatter = FALSE),
    error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL })
}

for (cfg in CONFIGS) {
  ko_df <- load_ko_wide(file.path(linked_dir, sprintf("sample_x_ko.%s.rel.wide.tsv", cfg)))
  if (is.null(ko_df)) { cat(sprintf("SKIP %s: not found\n", cfg)); next }
  run_maaslin(cfg, ko_df, meta_df, file.path(out_base, cfg))
}
pic_df <- load_ko_wide(file.path(INPUT_DIR, "pic_rel.csv"))
if (!is.null(pic_df)) run_maaslin("picrust2", pic_df, meta_df, file.path(out_base, "picrust2"))
run_maaslin("gold", gold_df, meta_df, file.path(out_base, "gold"))

cat("\n=== 02 DONE ===\n")
