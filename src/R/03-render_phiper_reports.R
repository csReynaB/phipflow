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


resolve_output_group_name <- function(active_group, output_group_mode, group_definitions = NULL) {
  if (!output_group_mode %in% c("group_name", "group_col")) {
    stop(
      "OUTPUT_GROUP_MODE must be either 'group_name' or 'group_col'. Got: ",
      output_group_mode,
      call. = FALSE
    )
  }

  if (output_group_mode == "group_name") {
    return(active_group)
  }

  # Legacy/current behavior:
  # use group_col from group_config.R if active_group is a group definition.
  if (!is.null(group_definitions) && active_group %in% names(group_definitions)) {
    return(group_definitions[[active_group]]$group_col)
  }

  # If active_group is already a simple metadata column, use it directly.
  active_group
}
  
base_dir <- get_kv_arg("BASE_DIR", required = TRUE)
group_cols <- parse_csv_arg(get_kv_arg("GROUP_COLS", required = TRUE))
OUTPUT_GROUP_MODE <- get_kv_arg("OUTPUT_GROUP_MODE", default = "group_name")
PHIPFLOW_SRC <- get_kv_arg("PHIPFLOW_SRC", required = TRUE)

template <- get_kv_arg(
  "TEMPLATE",
  default = file.path(dirname(PHIPFLOW_SRC), "template", "phiper_summary_report.qmd")
)

if (length(group_cols) == 0) {
  stop("Please provide GROUP_COLS=group1,group2,...", call. = FALSE)
}

if (!OUTPUT_GROUP_MODE %in% c("group_name", "group_col")) {
  stop(
    "OUTPUT_GROUP_MODE must be either 'group_name' or 'group_col'. Got: ",
    OUTPUT_GROUP_MODE,
    call. = FALSE
  )
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

# BASE_DIR is expected to be:
#   <project_dir>/results
# Therefore group_config.R should be:
#   <project_dir>/R/group_config.R
project_dir <- fs::path_dir(base_dir)
group_config_file <- fs::path(project_dir, "R", "group_config.R")

group_definitions <- NULL

if (fs::file_exists(group_config_file)) {
  message("Loading group_config.R: ", group_config_file)
  group_env <- new.env(parent = globalenv())
  sys.source(group_config_file, envir = group_env)

  if (exists("group_definitions", envir = group_env, inherits = FALSE)) {
    group_definitions <- get("group_definitions", envir = group_env)
  }
} else {
  message("group_config.R not found: ", group_config_file)
  message("Proceeding without group_config.R; active group names will be used directly.")
}

message("Base dir: ", base_dir)
message("Project dir: ", project_dir)
message("PHIPFLOW_SRC: ", PHIPFLOW_SRC)
message("Template: ", template)
message("Group columns / active groups: ", paste(group_cols, collapse = ", "))
message("OUTPUT_GROUP_MODE: ", OUTPUT_GROUP_MODE)

for (gc in group_cols) {
  output_group_name <- resolve_output_group_name(
    active_group = gc,
    output_group_mode = OUTPUT_GROUP_MODE,
    group_definitions = group_definitions
  )

  out_dir <- fs::path(base_dir, output_group_name)

  message("------------------------------------------------------------")
  message("Active group      : ", gc)
  message("Output group name : ", output_group_name)
  message("Report folder     : ", out_dir)

  if (!fs::dir_exists(out_dir)) {
    message("Skipping missing report folder: ", out_dir)
    next
  }

  qmd_copy <- fs::path(out_dir, "phiper_summary_report.qmd")
  fs::file_copy(template, qmd_copy, overwrite = TRUE)

  old_wd <- getwd()
  setwd(out_dir)

  tryCatch({
    quarto::quarto_render(
      input = qmd_copy,
      output_file = paste0("summary_report_", output_group_name, ".html"),
      execute_dir = out_dir,
      execute_params = list(
        base_dir = base_dir,
        group_col = output_group_name,
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

    message("Rendered report for active group: ", gc)
    message("Report output folder: ", out_dir)

  }, error = function(e) {
    message("Failed for ", gc, ": ", conditionMessage(e))

  }, finally = {
    setwd(old_wd)
    if (fs::file_exists(qmd_copy)) {
      fs::file_delete(qmd_copy)
    }
  })
}