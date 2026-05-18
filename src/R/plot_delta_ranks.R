# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
# load required R packages
library(phiper)
library(rlang)
library(ggplot2)
library(Cairo)
library(openxlsx)
library(dplyr)
library(purrr)
library(locfdr)
library(ggpubr)
library(tidyr)
library(stringr)
library(mgcv)

# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------
# parse CLI inputs and optional parameters (e.g., paths, filters, flags),
# validate values, and populate defaults used throughout the script
# ------------------------------------------------------------------------------
N_CORES <- 5
LOG <- TRUE
LOG_FILE <- NULL
MAX_GB <- 10
ACTIVE_GROUP <- "group_test"
DEFAULT_LONGITUDINAL <- FALSE
ALL <- TRUE
PROJECT_DIR <- "IBD-Berlin"
PARQUET_NAME <- paste0(PROJECT_DIR, ".parquet")
FORCE <- FALSE

args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  if (grepl("=", arg)) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    key <- parts[1]
    value <- parts[2]
    
    if (key == "N_CORES") {
      val_num <- suppressWarnings(as.numeric(value))
      if (!is.na(val_num)) N_CORES <- val_num
      
    } else if (key == "MAX_GB") {
      val_num <- suppressWarnings(as.numeric(value))
      if (!is.na(val_num)) MAX_GB <- val_num
      
    } else if (key == "LOG") {
      val_log <- tolower(value)
      if (val_log %in% c("true", "t", "1")) {
        LOG <- TRUE
      } else if (val_log %in% c("false", "f", "0")) {
        LOG <- FALSE
      }
      
    } else if (key == "LOG_FILE") {
      LOG_FILE <- sub("^['\\\"]|['\\\"]$", "", value)
      
    } else if (key == "FORCE") {
      val_log <- tolower(value)
      if (val_log %in% c("true", "t", "1")) {
        FORCE <- TRUE
      } else if (val_log %in% c("false", "f", "0")) {
        FORCE <- FALSE
      }
      
    } else if (key == "ALL") {
      val_log <- tolower(value)
      if (val_log %in% c("true", "t", "1")) {
        ALL <- TRUE
      } else if (val_log %in% c("false", "f", "0")) {
        ALL <- FALSE
      }
      
    } else if (key == "ACTIVE_GROUP") {
      ACTIVE_GROUP <- sub("^['\\\"]|['\\\"]$", "", value)
      
    } else if (key == "DEFAULT_LONGITUDINAL") {
      val_log <- tolower(value)
      if (val_log %in% c("true", "t", "1")) {
        DEFAULT_LONGITUDINAL <- TRUE
      } else if (val_log %in% c("false", "f", "0")) {
        DEFAULT_LONGITUDINAL <- FALSE
      } else {
        stop("DEFAULT_LONGITUDINAL must be one of TRUE/FALSE/T/F/1/0")
      }
    } else if (key == "PROJECT_DIR") {
        PROJECT_DIR <- sub("^['\\\"]|['\\\"]$", "", value)
    } else if (key == "PARQUET_NAME") {
        PARQUET_NAME <- sub("^['\\\"]|['\\\"]$", "", value)
    }
  }
}


# ------------------------------------------------------------------------------
# Define main directory paths and load group configuration and helper functions
# ------------------------------------------------------------------------------
project_dir <- PROJECT_DIR
data_dir <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "results")
r_dir <- file.path(project_dir, "R")

group_config_file <- file.path(r_dir, "group_config.R")
data_long_path <- file.path(data_dir, PARQUET_NAME)

message("PROJECT_DIR: ", project_dir)
message("data_long_path: ", data_long_path)
message("results_dir: ", results_dir)
message("group_config_file: ", group_config_file)

source(file.path("src/R", "helper_functions.R"))
source(group_config_file)

# ------------------------------------------------------------------------------
# Build PHIP data object (reproducible)
# ------------------------------------------------------------------------------
# create the PHIP data object from the input files and set a fixed seed to
# ensure reproducible sampling, random splits, and any stochastic steps
# ------------------------------------------------------------------------------
withr::with_preserve_seed({
  ps <- phip_convert(
    data_long_path    = data_long_path,
    sample_id         = "sample_id",
    peptide_id        = "peptide_id",
    exist             = "exist",
    fold_change       = "fold_change",
    counts_input      = NULL,
    counts_hit        = NULL,
    peptide_library   = TRUE,
    materialise_table = TRUE,
    auto_expand       = TRUE,
    n_cores           = 10
  )
})

# ------------------------------------------------------------------------------
# Results directory + peptide library snapshot
# ------------------------------------------------------------------------------
# create a base results folder and persist the peptide library used in this run
# for provenance and reproducibility
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

peptide_library_path <- file.path(results_dir, "peptide_library.rds")
if (!file.exists(peptide_library_path)) {
  message("Creating peptide library RDS: ", peptide_library_path)
  
  peptide_library <- ps %>%
    get_peptide_library() %>%
    collect() %>%
    as.data.frame()
  
  saveRDS(peptide_library, file = peptide_library_path)
} else {
  message("Using existing peptide library RDS: ", peptide_library_path)
}

# ------------------------------------------------------------------------------
# Analysis setup
# ------------------------------------------------------------------------------
group_cfg <- build_group_config(
  active_group_name = ACTIVE_GROUP,
  default_longitudinal = DEFAULT_LONGITUDINAL,
  group_definitions = group_definitions,
  fallback_palette = phip_palette,
  manual_comparisons = manual_comparisons,
  manual_longitudinal = manual_longitudinal
)

group_col <- group_cfg$group_col
groups <- group_cfg$groups
comparisons <- group_cfg$comparisons
longitudinal <- group_cfg$longitudinal
group_palette <- group_cfg$group_palette
#group_palette_full <- group_cfg$group_palette_full

message("group_col: ", group_col)
message("groups: ", paste(groups, collapse = ", "))
message("n comparisons: ", length(comparisons))
# ------------------------------------------------------------------------------
# Parallel backend (DELTA)
# ------------------------------------------------------------------------------
# configure BLAS/OpenMP to single-threaded mode to avoid CPU oversubscription
# when running multiple parallel R workers. Then set up a `future` plan that
# works across platforms (multisession on Windows; multicore/sequential on Unix)
Sys.setenv(
  OMP_NUM_THREADS     = "1",
  MKL_NUM_THREADS     = "1",
  OPENBLAS_NUM_THREADS = "1"
)

options(
  future.globals.maxSize = MAX_GB * 1024^3,
  future.scheduling      = 1
)

# store the current plan so it can be restored later.
original_plan <- future::plan()

if (.Platform$OS.type == "windows") {
  future::plan(future::multisession, workers = N_CORES)
} else if (N_CORES > 1L) {
  future::plan(future::multicore, workers = N_CORES)
} else {
  future::plan(future::sequential)
}

# ------------------------------------------------------------------------------
# Comparisons loop
# ------------------------------------------------------------------------------
# Run the same pipeline for each (var1 vs var2) comparison:
#   1) filter + export comparison subset
#   2) enrichment counts
#   3) alpha diversity
#   4) beta diversity (distances, ordinations, PERMANOVA, dispersion, t-SNE)
#   5) POP framework (prevalence comparison + plots)
#   6) DELTA framework (permutation-based differential prevalence)
#cmp <- comparisons[[1]]
#for (cmp in comparisons) {
for (cmp_idx in seq_along(comparisons)) {
  cmp <- comparisons[[cmp_idx]]
  var1 <- cmp[[1]]
  var2 <- cmp[[2]]
  message("Running comparison: ", var1, " vs. ", var2)
  
  is_longitudinal <- isTRUE(longitudinal[cmp_idx])
  paired_col <- if (is_longitudinal) "subject_id" else NULL
  
  # create output directory for this comparison
  label_dir <- paste(var1, "vs.", var2, sep = "_")
  out_dir   <- file.path(results_dir, group_col, label_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ----------------------------------------------------------------------------
  # filter to the two groups and add analysis-friendly grouping columns
  # ----------------------------------------------------------------------------
  ps_cmp <- ps %>%
  dplyr::filter(.data[[group_col]] %in% c(var1, var2)) %>%
  dplyr::mutate(
    group_char  = .data[[group_col]]#,
    #group_dummy = dplyr::if_else(.data[[group_col]] == var1, 1L, 0L)
  )
  
  extra_vars_to_save <- "fold_change"
  if (!is.null(paired_col)) {
  ps_cmp <- ps_cmp %>%
    dplyr::group_by(.data[[paired_col]]) %>%
    dplyr::filter(
      any(.data[[group_col]] == var1) & any(.data[[group_col]] == var2)
    ) %>%
    dplyr::ungroup()

  extra_vars_to_save <- c(extra_vars_to_save, paired_col)
  } 
  # persist the filtered dataset used for downstream frameworks.
  make_and_save(
    data       = ps_cmp,
    out_path   = file.path(out_dir, paste0(label_dir, "_data.rds")),
    extra_vars = extra_vars_to_save
    )
  
  # ----------------------------------------------------------------------------
  # 
  # ----------------------------------------------------------------------------
  data_frameworks <- readRDS(file.path(out_dir, paste0(label_dir, "_data.rds")))
  data_frameworks$group_char <- factor(data_frameworks$group_char, levels = rev(cmp))
  
  peplib <- readRDS(file.path("BC-Engl/results", "peptide_library.rds"))

  # paired_prev <- if (is_longitudinal) "subject_id" else FALSE
  # prev_res_pep <- phiper::ph_prevalence_compare(
  #   x                 = data_frameworks,
  #   group_cols        = "group_char",
  #   rank_cols         = "peptide_id",
  #   compute_ratios_db = TRUE,
  #   paired            = paired_prev,
  #   parallel          = TRUE,
  #   collect           = TRUE
  # )
  # pep_tbl <- extract_tbl(prev_res_pep)
  
  #ranks_tax <- c("phylum", "class", "order", "family", "genus", "species")
  # ranks_tax <- c( "species")

  # prev_res_rank <- phiper::ph_prevalence_compare(
  #   x                 = data_frameworks,
  #   group_cols        = "group_char",
  #   rank_cols         = ranks_tax,
  #   compute_ratios_db = FALSE,
  #   paired            = paired_prev,
  #   parallel          = TRUE,
  #   peptide_library   = peplib,
  #   collect           = TRUE
  # )

  # rank_tbl <- extract_tbl(prev_res_rank)  
  # ranks_combined <- c(ranks_tax, "peptide_id")  
  # for (rank_name in ranks_combined) {
  #   #rank_name <- "species"
  #   rank_chr <- as.character(rank_name)
  #   df_rank <- if (rank_chr == "peptide_id") {
  #     pep_tbl
  #   } else {
  #     rank_tbl %>% filter(rank == rank_chr)
  #   }
  # }
  # ----------------------------------------------------------------------------
  # DELTA framework
  # ----------------------------------------------------------------------------
  dir.create(file.path(out_dir, "DELTA_framework"), recursive = TRUE,
             showWarnings = FALSE)
  
  dir.create(file.path(out_dir, "DELTA_framework", "species"), recursive = TRUE,
             showWarnings = FALSE)

  data_frameworks$subject_id <- data_frameworks$sample_id
  peplib[] <- lapply(peplib, as.character)
  log_file_current <- if (is.null(LOG_FILE)) {
    file.path(out_dir, "DELTA_framework",  "species", "log.txt")
  } else {
    LOG_FILE
  }
  
  stat_mode <- if (is_longitudinal) {
    "srlr_paired"
  } else {
    "srlr"
  }
  paired_by <- if (is_longitudinal) paired_col else NULL


  res <- read.csv(file.path(out_dir, "DELTA_framework",
                                "delta_table.csv"))

  for (taxon in c("order", "family", "genus", "species", "all")){
    if(taxon %in% "all"){
        std_idx <- res$rank %in% c("domain", "kingdom", "phylum", "class", "order", "family",
                 "genus", "species")
        res$feature[!std_idx] <- as.character(res$rank[!std_idx])
        res$rank <- taxon
    }
    
    svglite::svglite(
      file.path(out_dir, "DELTA_framework",
                paste0("significant_static_", taxon, "_forestplot.svg")),
      width = 20/2.54,
      height = 20/2.54,
      bg = "white"
    )

    message("Creating forest plot for taxon: ", taxon, ", group1 = ", res$group1[1], ", group2 = ", res$group2[1])
    
    p_forest_unc <- phiper::forestplot(
      results_tbl          = res,
      rank_of_interest     = taxon,
      use_diverging_colors = TRUE,
      filter_significant   = "p_perm",
      left_label           = paste0("More in ", as.character(res$group1[1])),#df_rank$group1[1]),
      right_label          = paste0("More in ", as.character(res$group2[1])),#df_rank$group2[1]),
      label_vjust           = -0.9,
      y_pad                 = 0.3,
      label_x_gap_frac      = -0.3,
      statistic_to_plot     = "T_stand"
    )
    print(p_forest_unc)
    dev.off()
    
    p_inter <- phiper::forestplot_interactive(
      results_tbl            = res,
      rank_of_interest       = taxon,
      statistic_to_plot      = "T_stand",
      use_diverging_colors   = TRUE,
      filter_significant     = "p_perm",
      left_label           = paste0("More in ", as.character(res$group1[1])),#df_rank$group1[1]),
      right_label          = paste0("More in ", as.character(res$group2[1])),#df_rank$group2[1]),
      arrow_length_frac      = 0.35,
      label_x_gap_frac       = -0.3,
      label_y_offset        = -0.9
    )$plot
    htmlwidgets::saveWidget(
      p_inter,
      file = file.path(out_dir, "DELTA_framework",
                      paste0("significant_interactive_", taxon, "_forestplot.html")),
      selfcontained = TRUE
    )
  }
}



# ------------------------------------------------------------------------------
# restore original future plan
# ------------------------------------------------------------------------------
future::plan(original_plan)
