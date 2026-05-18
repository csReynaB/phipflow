#!/usr/bin/env Rscript

message("R version: ", R.version.string)
message("R executable: ", Sys.which("R"))
message("Rscript executable: ", Sys.which("Rscript"))
message("Library paths: ", paste(.libPaths(), collapse = " | "))

## ---------------------------- ARGUMENTS --------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL, required = FALSE) {
  idx <- match(flag, args)
  if (is.na(idx)) {
    if (required) stop("Missing required argument: ", flag, call. = FALSE)
    return(default)
  }
  if (idx == length(args)) stop("Missing value for argument: ", flag, call. = FALSE)
  args[[idx + 1]]
}

base_dir       <- get_arg("--base_dir", required = TRUE)
project_name   <- get_arg("--project_name", required = TRUE)

exist_file     <- get_arg("--exist_file", "exist.csv")
fold_file      <- get_arg("--fold_file", "fold.csv")
metadata_file  <- get_arg("--metadata_file", required = TRUE)

out_parquet    <- get_arg("--out_parquet", required = TRUE)

sample_prefix  <- get_arg("--sample_prefix", "R")
peptide_col    <- get_arg("--peptide_col", "peptide_name")

replace_inf    <- get_arg("--replace_inf", "TRUE")
replace_inf    <- toupper(replace_inf) %in% c("TRUE", "T", "1", "YES", "Y")

## ---------------------------- SETUP ------------------------------------------

set.seed(16748991)

packages <- c(
  "dplyr",
  "tidyr",
  "data.table",
  "DBI",
  "duckdb"
)

missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them in renv/container before running the workflow.",
    call. = FALSE
  )
}

invisible(lapply(packages, library, character.only = TRUE))
## ---------------------------- PATHS ------------------------------------------

project_dir <- file.path(base_dir, project_name)
data_dir    <- file.path(project_dir, "Data")
meta_dir    <- file.path(project_dir, "Metadata")

exist_path  <- file.path(data_dir, exist_file)
fold_path   <- file.path(data_dir, fold_file)
meta_path   <- file.path(meta_dir, metadata_file)

out_path    <- file.path(data_dir, out_parquet)

message("Project directory: ", project_dir)
message("Reading exist matrix: ", exist_path)
message("Reading fold matrix:  ", fold_path)
message("Reading metadata:     ", meta_path)
message("Writing parquet:      ", out_path)

if (!file.exists(exist_path)) stop("exist file not found: ", exist_path, call. = FALSE)
if (!file.exists(fold_path)) stop("fold file not found: ", fold_path, call. = FALSE)
if (!file.exists(meta_path)) stop("metadata file not found: ", meta_path, call. = FALSE)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

## ---------------------------- HELPERS ----------------------------------------

detect_sample_columns <- function(df, sample_prefix = "R") {
  start_idx <- which(startsWith(names(df), sample_prefix))[1]

  if (is.na(start_idx)) {
    stop(
      "Could not detect sample columns using prefix: ", sample_prefix,
      "\nAvailable columns are:\n",
      paste(names(df), collapse = ", "),
      call. = FALSE
    )
  }

  names(df)[start_idx:ncol(df)]
}

matrix_to_long <- function(df, value_col, peptide_col, sample_prefix = "R") {
  if (!peptide_col %in% names(df)) {
    stop(
      "Peptide column '", peptide_col, "' not found.\n",
      "Available columns are:\n",
      paste(names(df), collapse = ", "),
      call. = FALSE
    )
  }

  sample_cols <- detect_sample_columns(df, sample_prefix = sample_prefix)

  df |>
    dplyr::select(
      dplyr::all_of(peptide_col),
      dplyr::all_of(sample_cols)
    ) |>
    tidyr::pivot_longer(
      cols      = -dplyr::all_of(peptide_col),
      names_to  = "sample_id",
      values_to = value_col
    ) |>
    dplyr::rename(peptide_id = dplyr::all_of(peptide_col))
}

## ---------------------------- READ DATA --------------------------------------

data_exist <- data.table::fread(exist_path)
data_fold  <- data.table::fread(fold_path)

long_exist <- matrix_to_long(
  df            = data_exist,
  value_col     = "exist",
  peptide_col   = peptide_col,
  sample_prefix = sample_prefix
)

long_fc <- matrix_to_long(
  df            = data_fold,
  value_col     = "fold_change",
  peptide_col   = peptide_col,
  sample_prefix = sample_prefix
)

long_df <- long_exist |>
  dplyr::left_join(
    long_fc,
    by = c("sample_id", "peptide_id")
  )

## ---------------------------- METADATA ---------------------------------------

metadata <- data.table::fread(meta_path)

colnames(metadata)[1] <- "sample_id"
names(metadata) <- make.unique(names(metadata))

long_df <- long_df |>
  dplyr::left_join(metadata, by = "sample_id")

missing_meta <- long_df |>
  dplyr::filter(is.na(.data[[names(metadata)[2]]])) |>
  dplyr::distinct(sample_id)

if (nrow(missing_meta) > 0) {
  warning(
    "Some samples in the matrices were not found in metadata. Example missing samples: ",
    paste(head(missing_meta$sample_id, 10), collapse = ", ")
  )
}

## ---------------------------- CLEAN FOLD CHANGE ------------------------------

if (replace_inf) {
  max_fc <- max(long_df$fold_change[is.finite(long_df$fold_change)], na.rm = TRUE)

  if (is.finite(max_fc)) {
    idx_inf <- is.infinite(long_df$fold_change)
    n_inf <- sum(idx_inf, na.rm = TRUE)

    if (n_inf > 0) {
      message("Replacing ", n_inf, " infinite fold_change values with max finite value: ", max_fc)
      long_df$fold_change[idx_inf] <- max_fc
    }
  } else {
    warning("No finite fold_change values found. Infinite values were not replaced.")
  }
}

## ---------------------------- SAFETY CHECKS ----------------------------------

data_export <- long_df |>
  dplyr::distinct(sample_id, peptide_id, .keep_all = TRUE)

n_before <- nrow(long_df)
n_after  <- nrow(data_export)

if (n_before != n_after) {
  warning(
    "Removed duplicated sample_id/peptide_id rows: ",
    n_before - n_after
  )
}

message("Number of samples:  ", dplyr::n_distinct(data_export$sample_id))
message("Number of peptides: ", dplyr::n_distinct(data_export$peptide_id))
message("Number of rows:     ", nrow(data_export))

## ---------------------------- EXPORT PARQUET ---------------------------------

con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:", read_only = FALSE)

on.exit({
  try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
}, add = TRUE)

duckdb::duckdb_register(con, "df_tbl", data_export)

DBI::dbExecute(
  con,
  paste0(
    "COPY df_tbl TO ",
    DBI::dbQuoteString(con, out_path),
    " (FORMAT PARQUET);"
  )
)

message("Done. Parquet written to: ", out_path)

rm(
  list = c(
    "long_df",
    "data_export",
    "con",
    "data_exist",
    "data_fold",
    "metadata",
    "long_exist",
    "long_fc"
  )
)

gc()