# =============================================================================
# config.R  --  Shared configuration for the ML prediction pipeline
#
# Edit ONLY this file to run the pipeline in a new environment.
# All scripts in 07_ML_Prediction/ source this file at the top, so paths and
# parameters live in one place instead of being hard-coded in each script.
#
# Usage (from the repository root):
#   Rscript 07_ML_Prediction/01_train_step1_cluster.R
#   Rscript 07_ML_Prediction/02_predict_step2_ko.R
#   ...
# Each script changes its working directory assumption to the repo root, then
# does:  source("07_ML_Prediction/config.R")
# =============================================================================

# ---- 1. Root directory ------------------------------------------------------
# By default the repository root is taken from the environment variable
# KGM_ROOT if it is set; otherwise it falls back to the current working
# directory. This lets a user run the pipeline without editing any path:
#   export KGM_ROOT=/path/to/Korean_Gut_MAGs   (Linux/macOS)
# or simply run R from inside the repository root.
ROOT_DIR <- Sys.getenv("KGM_ROOT", unset = getwd())

# ---- 2. Standard sub-directories --------------------------------------------
# Inputs the user must supply (see README "Input data" section), intermediate
# files written by the pipeline, and trained models for external application.
INPUT_DIR  <- file.path(ROOT_DIR, "07_ML_Prediction", "input")   # user-supplied data
OUTPUT_DIR <- file.path(ROOT_DIR, "07_ML_Prediction", "output")  # generated results
MODEL_DIR  <- file.path(ROOT_DIR, "07_ML_Prediction", "models")  # trained models / carrier profile

# Create output/model dirs if missing (inputs must already exist).
for (d in c(OUTPUT_DIR, MODEL_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- 3. Input file names ----------------------------------------------------
# Centralised so a user with differently named files edits here, not in code.
FILES <- list(
  # Internal cohort (224 paired samples)
  asv_internal      = file.path(INPUT_DIR, "ASV_RA_224.csv"),               # ASV x sample relative abundance
  drep95            = file.path(INPUT_DIR, "drep_95_bigliver_hybrid.csv"),  # genome -> dRep95 cluster map
  cpm_taxonomy      = file.path(INPUT_DIR, "taxonomy_cpm_pivot_gtdb.csv"),  # per-sample MAG CPM abundance
  meta_internal     = file.path(INPUT_DIR, "sample_meta.csv"),             # sample_id, group (Healthy/MASLD)
  bin_ko_fpkm       = file.path(INPUT_DIR, "all_bins_KEGG_ko_FPKM_pivot.csv"), # bin x KO FPKM (for carrier profile)
  gold_fpkm         = file.path(INPUT_DIR, "fpkm_gold_long.csv"),           # shotgun gold KO abundance (long)
  picrust_internal  = file.path(INPUT_DIR, "pic_rel.csv"),                  # PICRUSt2 internal (for benchmark)

  # External cohort (amplicon-only, e.g. 394 samples)
  asv_external      = file.path(INPUT_DIR, "ASV_RA_external.csv"),
  meta_external     = file.path(INPUT_DIR, "meta_external.csv"),
  picrust_external  = file.path(INPUT_DIR, "picrust2_ext_ko.tsv")
)

# ---- 4. Trained-model / carrier output names --------------------------------
MODEL_FILES <- list(
  xgb_full   = file.path(MODEL_DIR, "full_xgb_models.rds"),   # full-224 XGBoost models (reusable)
  carrier    = file.path(MODEL_DIR, "carrier_profile.rds")    # cluster x KO binary carrier matrix
)

# ---- 5. Pipeline parameters -------------------------------------------------
# These reproduce the values described in the manuscript Methods.
PARAMS <- list(
  seed          = 42,
  asv_prev_cut  = 5,      # ASV kept if present in >= 5 samples
  cluster_prev  = 10,     # dRep95 cluster kept if present in >= 10 samples (prev10, main)
  carrier_thr   = 0.25,   # cluster-KO carrier fraction threshold (winner config)
  xgb_nrounds   = 100,
  xgb_params    = list(
    objective        = "reg:squarederror",
    eta              = 0.1,
    max_depth        = 4,
    subsample        = 0.8,
    colsample_bytree = 0.5,
    min_child_weight = 3,
    gamma            = 0,
    lambda           = 1,
    nthread          = 4
  )
)

set.seed(PARAMS$seed)

# ---- 6. Small startup report ------------------------------------------------
message("config.R loaded")
message("  ROOT_DIR  : ", ROOT_DIR)
message("  INPUT_DIR : ", INPUT_DIR)
message("  OUTPUT_DIR: ", OUTPUT_DIR)
message("  MODEL_DIR : ", MODEL_DIR)
