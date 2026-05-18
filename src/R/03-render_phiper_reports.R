#!/usr/bin/env Rscript

message("R version: ", R.version.string)
message("R executable: ", Sys.which("R"))
message("Rscript executable: ", Sys.which("Rscript"))
message("Library paths: ", paste(.libPaths(), collapse = " | "))
message("Quarto: ", Sys.which("quarto"))
message("Pandoc: ", Sys.which("pandoc"))

suppressPackageStartupMessages({
  library(fs)
  library(quarto)
})

args <- commandArgs(trailingOnly = TRUE)

get_kv_arg <- function(key, default = NULL, required = FALSE) {
  pattern <- paste0("^", key, "=")
  hit <- grep(pattern, args, value = TRUE)

  if (length(hit) == 0) {
    if (required) {
      stop("Missing required argument: ", key, "=...", call. = FALSE)
    }
    return(default)
  }

  value <- sub(pattern, "", hit[[1]])
  value <- sub("^['\\\"]|['\\\"]$", "", value)
  trimws(value)
}

parse_csv_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  x <- strsplit(x, ",", fixed = TRUE)[[1]]
  x <- trimws(x)
  x[nzchar(x)]
}

base_dir <- get_kv_arg("BASE_DIR", required = TRUE)
group_cols <- parse_csv_arg(get_kv_arg("GROUP_COLS", required = TRUE))
PHIPFLOW_SRC <- get_kv_arg("PHIPFLOW_SRC", required = TRUE)

template <- get_kv_arg(
  "TEMPLATE",
  default = file.path(dirname(PHIPFLOW_SRC), "template", "phiper_summary_report.qmd")
)

if (length(group_cols) == 0) {
  stop("Please provide GROUP_COLS=group1,group2,...", call. = FALSE)
}

base_dir <- fs::path_abs(base_dir)
PHIPFLOW_SRC <- fs::path_abs(PHIPFLOW_SRC)
template <- fs::path_abs(template)

if (!fs::dir_exists(base_dir)) {
  stop("BASE_DIR does not exist: ", base_dir, call. = FALSE)
}

if (!fs::dir_exists(PHIPFLOW_SRC)) {
  stop("PHIPFLOW_SRC does not exist: ", PHIPFLOW_SRC, call. = FALSE)
}

if (!fs::file_exists(template)) {
  stop("Template not found: ", template, call. = FALSE)
}

message("Base dir: ", base_dir)
message("PHIPFLOW_SRC: ", PHIPFLOW_SRC)
message("Template: ", template)
message("Group columns: ", paste(group_cols, collapse = ", "))

for (gc in group_cols) {
  out_dir <- fs::path(base_dir, gc)

  if (!fs::dir_exists(out_dir)) {
    message("Skipping missing group column folder: ", out_dir)
    next
  }

  qmd_copy <- fs::path(out_dir, "phiper_summary_report.qmd")
  fs::file_copy(template, qmd_copy, overwrite = TRUE)

  old_wd <- getwd()
  setwd(out_dir)

  tryCatch({
    quarto::quarto_render(
      input = qmd_copy,
      output_file = paste0("summary_report_", gc, ".html"),
      execute_dir = out_dir,
      execute_params = list(
        base_dir = base_dir,
        group_col = gc,
        top_n = 500,
        tables_open = FALSE,
        include_beta_tables = FALSE,
        include_alpha_tables = FALSE,
        include_tsne3d = TRUE,
        include_longitudinal_stability = TRUE,
        include_pop_tables = TRUE,
        include_pop_plots = TRUE,
        include_delta_tables = TRUE,
        include_delta_feature_plots = TRUE,
        pop_interactive_mode = "link",
        delta_interactive_mode = "embed",
        max_delta_feature_plots = 5
      )
    )

    message("Rendered report for: ", gc)

  }, error = function(e) {
    message("Failed for ", gc, ": ", conditionMessage(e))

  }, finally = {
    setwd(old_wd)
    if (fs::file_exists(qmd_copy)) {
      fs::file_delete(qmd_copy)
    }
  })
}