options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  Ncpus = max(1, parallel::detectCores() - 1)
)

cat("Installing pak...\n")
install.packages("pak")

cran_pkgs <- c(
  # core IO / data
  "data.table",
  "DBI",
  "dbplyr",
  "dplyr",
  "duckdb",
  "fs",
  "readr",
  "readxl",
  "tidyr",
  "tibble",
  "purrr",
  "stringr",
  "forcats",
  "rlang",
  "withr",
  "yaml",
  "glue",

  # statistics / ecology / modeling helpers
  "vegan",
  "permute",
  "mgcv",
  "locfdr",
  "future",
  "future.apply",

  # plotting / reports
  "ggplot2",
  "ggpubr",
  "ggsignif",
  "ggtext",
  "patchwork",
  "plotly",
  "Cairo",
  "svglite",
  "showtext",
  "sysfonts",
  "scales",
  "viridisLite",
  "htmlwidgets",
  "htmltools",
  "shiny",
	
  # tables / office / rendering
  "knitr",
  "DT",
  "jsonlite",
  "kableExtra",
  "rmarkdown",
  "quarto",
  "openxlsx",
	
  # package/helper dependencies
  "chk",
  "cli",
  "tidyselect",
  "Rcpp",
  "RcppParallel",
  "Rtsne"
)

cat("Installing CRAN packages...\n")
install.packages(setdiff(cran_pkgs, rownames(installed.packages())))

cat("Installing phiperio and phiper from GitHub...\n")
pak::pak(c(
  "Polymerase3/phiperio",
  "Polymerase3/phiper"
))

cat("Running package checks...\n")
required <- c(
  "phiper",
  "phiperio",
  "duckdb",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "ggtext",
  "plotly",
  "Cairo",
  "svglite",
  "quarto",
  "rmarkdown",
  "openxlsx",
  "vegan",
  "locfdr",
  "future",
  "future.apply",
  "showtext",
  "sysfonts",
  "viridisLite",
  "permute",
  "dbplyr",
  "Rtsne",
  "htmlwidgets",
  "htmltools",
  "knitr",
  "shiny",
  "DT",
  "jsonlite"
)

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  stop("Missing required packages: ", paste(missing, collapse = ", "))
}

cat("Installed versions:\n")
cat("R:", R.version.string, "\n")
cat("phiper:", as.character(packageVersion("phiper")), "\n")
cat("phiperio:", as.character(packageVersion("phiperio")), "\n")
cat("duckdb:", as.character(packageVersion("duckdb")), "\n")
cat("quarto R package:", as.character(packageVersion("quarto")), "\n")
