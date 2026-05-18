#!/usr/bin/env Rscript

message("R version: ", R.version.string)
message("R executable: ", Sys.which("R"))
message("Rscript executable: ", Sys.which("Rscript"))
message("Library paths: ", paste(.libPaths(), collapse = " | "))
message("Quarto: ", Sys.which("quarto"))
message("Pandoc: ", Sys.which("pandoc"))

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
required_packages <- c(
  "phiper", "phiperio", "rlang", "ggplot2", "openxlsx",
  "dplyr", "purrr", "ggpubr", "tidyr", "ggtext", "patchwork",
  "withr", "future", "svglite", "htmlwidgets", "plotly", "viridisLite",
  "forcats", "tibble", "vegan", "permute", "dbplyr"
)


missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them in renv/container before running the workflow.",
    call. = FALSE
  )
}

attach_packages <- c(
  "phiper",
  "rlang",
  "ggplot2",
  "openxlsx",
  "dplyr",
  "purrr",
  "ggpubr",
  "tidyr",
  "ggtext",
  "patchwork"
)

invisible(lapply(attach_packages, library, character.only = TRUE))
# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

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
  value <- trimws(value)

  value
}

parse_bool <- function(x, default = FALSE, name = "value") {
  if (is.null(x) || !nzchar(x)) return(default)

  x <- tolower(trimws(as.character(x)))

  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)

  stop(name, " must be one of TRUE/FALSE/T/F/1/0/YES/NO", call. = FALSE)
}

parse_nullable <- function(x) {
  if (is.null(x)) return(NULL)

  x <- trimws(as.character(x))

  if (toupper(x) %in% c("", "NULL", "NA", "NONE")) return(NULL)

  x
}

N_CORES <- as.integer(get_kv_arg("N_CORES", default = "1"))
MAX_GB  <- as.numeric(get_kv_arg("MAX_GB", default = "40"))

if (is.na(N_CORES) || N_CORES < 1L) {
  stop("N_CORES must be a positive integer.", call. = FALSE)
}

if (is.na(MAX_GB) || MAX_GB <= 0) {
  stop("MAX_GB must be positive.", call. = FALSE)
}

IO_CORES   <- max(1L, min(N_CORES, 10L))
DIST_CORES <- max(1L, min(N_CORES, 8L))

LOG   <- parse_bool(get_kv_arg("LOG", default = "TRUE"), default = TRUE, name = "LOG")
FORCE <- parse_bool(get_kv_arg("FORCE", default = "FALSE"), default = FALSE, name = "FORCE")
ALL   <- parse_bool(get_kv_arg("ALL", default = "FALSE"), default = FALSE, name = "ALL")

DEFAULT_LONGITUDINAL <- parse_bool(
  get_kv_arg("DEFAULT_LONGITUDINAL", default = "FALSE"),
  default = FALSE,
  name = "DEFAULT_LONGITUDINAL"
)

LOG_FILE <- parse_nullable(get_kv_arg("LOG_FILE", default = NULL))

ACTIVE_GROUP <- get_kv_arg("ACTIVE_GROUP", required = TRUE)
PROJECT_DIR  <- get_kv_arg("PROJECT_DIR", required = TRUE)
PARQUET_NAME <- get_kv_arg("PARQUET_NAME", default = paste0(PROJECT_DIR, ".parquet"))

MANUAL_COMPARISON_FILE <- parse_nullable(
  get_kv_arg("MANUAL_COMPARISON_FILE", default = NULL)
)

PHIPFLOW_SRC <- get_kv_arg("PHIPFLOW_SRC", required = TRUE)
PEPTIDE_LIBRARY <- get_kv_arg("PEPTIDE_LIBRARY", required = TRUE)
# ------------------------------------------------------------------------------
# Define main directory paths and load group configuration/helper functions
# ------------------------------------------------------------------------------

project_dir <- PROJECT_DIR
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "results")
r_dir       <- file.path(project_dir, "R")

group_config_file <- file.path(r_dir, "group_config.R")
data_long_path    <- file.path(data_dir, PARQUET_NAME)

manual_comparison_file <- NULL

if (!is.null(MANUAL_COMPARISON_FILE) && nzchar(MANUAL_COMPARISON_FILE)) {
  manual_comparison_file <- file.path(r_dir, MANUAL_COMPARISON_FILE)
}

helper_file <- file.path(PHIPFLOW_SRC, "helper_functions.R")
peptide_library_path <- normalizePath(PEPTIDE_LIBRARY, mustWork = FALSE)

message("PROJECT_DIR: ", project_dir)
message("ACTIVE_GROUP: ", ACTIVE_GROUP)
message("PARQUET_NAME: ", PARQUET_NAME)
message("data_long_path: ", data_long_path)
message("results_dir: ", results_dir)
message("group_config_file: ", group_config_file)
message("PHIPFLOW_SRC: ", PHIPFLOW_SRC)
message("PEPTIDE_LIBRARY: ", peptide_library_path)
message("helper_file: ", helper_file)
message("N_CORES: ", N_CORES)
message("IO_CORES: ", IO_CORES)
message("DIST_CORES: ", DIST_CORES)
message("MAX_GB: ", MAX_GB)
message("ALL: ", ALL)
message("DEFAULT_LONGITUDINAL: ", DEFAULT_LONGITUDINAL)
message("FORCE: ", FORCE)


if (!file.exists(data_long_path)) {
  stop("Input parquet does not exist: ", data_long_path, call. = FALSE)
}

if (!file.exists(group_config_file)) {
  stop("group_config.R does not exist: ", group_config_file, call. = FALSE)
}

if (!file.exists(helper_file)) {
  stop("helper_functions.R does not exist: ", helper_file, call. = FALSE)
}

source(helper_file)
source(group_config_file)
# ------------------------------------------------------------------------------
# Optional manual comparisons
# ------------------------------------------------------------------------------

comparison_cfg <- load_manual_comparison_config(
  manual_comparison_file = manual_comparison_file
)

manual_comparisons <- comparison_cfg$manual_comparisons
manual_longitudinal <- comparison_cfg$manual_longitudinal

if (is.null(manual_comparisons)) {
  message(
    "No external manual comparison file loaded. ",
    "Using comparisons from group_config.R if available; ",
    "otherwise generating automatic pairwise comparisons."
  )
} else {
  message("External manual comparisons loaded: ", length(manual_comparisons))
}

# ------------------------------------------------------------------------------
# Build PHIP data object (reproducible)
# ------------------------------------------------------------------------------
# create the PHIP data object from the input files and set a fixed seed to
# ensure reproducible sampling, random splits, and any stochastic steps
# ------------------------------------------------------------------------------
withr::with_preserve_seed({
  ps <- phiperio::convert_standard(
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
    n_cores           = IO_CORES
  )
})

# ------------------------------------------------------------------------------
# Results directory + peptide library
# ------------------------------------------------------------------------------
# The peptide library is now a workflow-level resource stored in:
#   phipflow/peplib/peptide_library.rds
#
# It is required and is not created by this analysis script.
# ------------------------------------------------------------------------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

message("Using workflow peptide library RDS: ", peptide_library_path)
peplib <- as.data.frame(readRDS(peptide_library_path))

# ------------------------------------------------------------------------------
# Analysis setup
# ------------------------------------------------------------------------------
phip_palette <- extend_palette_distinct(phip_palette, n = 50)
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

for (i in seq_along(comparisons)) {
  message(
    "comparison ", i, ": ",
    comparisons[[i]][1], " vs ", comparisons[[i]][2],
    " | longitudinal = ", longitudinal[[i]]
  )
}
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

# store the current plan so it can be restored even if the script fails
original_plan <- future::plan()
on.exit(future::plan(original_plan), add = TRUE)

if (.Platform$OS.type == "windows") {
  future::plan(future::multisession, workers = N_CORES)
} else if (N_CORES > 1L) {
  future::plan(future::multicore, workers = N_CORES)
} else {
  future::plan(future::sequential)
}


# ------------------------------------------------------------------------------
# Using all group test for alpha and beta
# ------------------------------------------------------------------------------
if(ALL){
    message("Running analysis on all ", group_col)

    out_dir   <- file.path(results_dir, group_col, "all")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # ----------------------------------------------------------------------------
    # enrichment counts
    # ----------------------------------------------------------------------------
    svglite::svglite(
    file.path(out_dir, "enrichment_counts.svg"),
    width = 30/2.54,
    height = 20/2.54,
    bg = "white"
    )

    p_enrich <- plot_enrichment_counts(
    ps,
    group_cols = group_col,
    custom_colors = group_palette,
    annotation_size = 4
    )
    p_enrich$data <- p_enrich$data %>%
    dplyr::filter(Cohort %in% groups) %>%
    dplyr::mutate(Cohort = factor(Cohort, levels = groups))
    for (i in seq_along(p_enrich$layers)) {
    layer_data <- p_enrich$layers[[i]]$data
    if (!is.null(layer_data) && "Cohort" %in% colnames(layer_data)) {
        p_enrich$layers[[i]]$data <- layer_data %>%
        dplyr::filter(Cohort %in% groups) %>%
        dplyr::mutate(Cohort = factor(Cohort, levels = groups))
    }
    }

    p_enrich <- p_enrich +
    labs(title = "Enrichment counts") +
    theme(text = element_text(family = "Montserrat")) +
    facet_wrap(
        vars(Cohort),
        ncol = 2,
        scales = "free_x"
    )
    print(p_enrich)
    dev.off()


    # ----------------------------------------------------------------------------
    # alpha diversity
    # ----------------------------------------------------------------------------
    alpha_div <- compute_alpha(ps,
                                        group_cols = group_col,
                                        carry_cols = c("Sex", "Age"))

    dir.create(file.path(out_dir, "alpha_diversity"), 
            recursive = T, showWarnings = F)
    write.xlsx(alpha_div[[group_col]], file.path(out_dir, 
                                                "alpha_diversity", "table.xlsx"))
    pv <- kruskal.test(reformulate(group_col, response = "richness"), 
                    data = alpha_div[[group_col]])$p.value

    svglite::svglite(
    file.path(out_dir, "alpha_diversity", "plot.svg"),
    width = 20/2.54,
    height = 13/2.54,
    bg = "white"
    )
    p_alpha <- plot_alpha(alpha_div,
                                    metric = "richness",
                                    group_col = group_col,
                                    x_order = groups,
                                    custom_colors = group_palette,
                                    text_size = 14) +
    labs(
        x = NULL#,
        #y = "New Y Label"
    ) +
    annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste("Kruskal-Wallis p =", format_pval(pv)),
        hjust = 1,
        vjust = 1,
        size = 4.5
    ) +
    coord_cartesian(clip = "off") +
    theme(
        axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
        )
    )
    p_alpha$layers$geom_boxplot$stat_params$width <- 0.4
    p_alpha$layers$geom_point$position$width <- 0.1
    print(p_alpha)
    dev.off()

    # species alpha diversity
    alpha_div <- compute_alpha(ps,
                                        group_cols = group_col,
                                        carry_cols = c("Sex", "Age"),
                                        ranks = "species")
    write.xlsx(alpha_div[[group_col]], file.path(out_dir, 
                                                "alpha_diversity", "table_species.xlsx"))
    # Calculate and format p-values first
    pv <- kruskal.test(reformulate(group_col, response = "richness"), 
                    data = alpha_div[[group_col]])$p.value
    svglite::svglite(
    file.path(out_dir, "alpha_diversity", "plot_species.svg"),
    width = 20/2.54,
    height = 13/2.54,
    bg = "white"
    )
    p_alpha_sh <- plot_alpha(
    alpha_div,
    #metric = "shannon_diversity",
    group_col = group_col,
    custom_colors = group_palette,
    x_order = groups,
    text_size = 14
    ) +
    labs(
        x = NULL#,
        #y = "New Y Label"
    ) +
    annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste("Kruskal-Wallis p =", format_pval(pv)),
        hjust = 1,
        vjust = 1,
        size = 4.5
    ) +
    coord_cartesian(clip = "off") +
    #scale_fill_manual(values = group_palette) +
    theme(
        axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
        )
    ) 
    p_alpha_sh$layers$geom_boxplot$stat_params$width <- 0.4
    p_alpha_sh$layers$geom_point$position$width <- 0.1
    print(p_alpha_sh)
    dev.off()


    # ----------------------------------------------------------------------------
    # beta diversity
    # ----------------------------------------------------------------------------
    dist_bc <- phiper:::compute_distance(ps, value_col = "exist",
                                        method_normalization = "auto",
                                        distance = "kulczynski", n_threads = DIST_CORES)                                      

    dir.create(file.path(out_dir, "beta_diversity"), 
            recursive = T, showWarnings = F)


    dist_mat <- as.matrix(dist_bc)
    openxlsx::write.xlsx(
    dist_mat,
    file = file.path(out_dir, "beta_diversity", "distance_matrix.xlsx"),
    rowNames = TRUE
    )

    cap_res <- phiper:::compute_capscale(dist_bc,
                                        ps = ps,
                                        formula = reformulate(group_col))
    saveRDS(cap_res, file.path(out_dir, "beta_diversity", "capscale_results.rds"))

    permanova_res <- phiper:::compute_permanova(dist_bc,
                                                ps = ps,
                                                group_col = group_col)
    saveRDS(permanova_res, file.path(out_dir, 
                                    "beta_diversity", "permanova_results.rds"))


    disp_res <- phiper:::compute_dispersion(dist_bc,
                                            ps = ps,
                                            group_col = group_col)
    saveRDS(disp_res, file.path(out_dir, 
                                "beta_diversity", "dispersion_results.rds"))

    #print(disp_res)
    n_axes_pcoa <- max(2L, min(109L, length(disp_res$distances$sample_id) - 1L))
    pcoa_res <- phiper:::compute_pcoa(dist_bc,
                                    neg_correction = "none",
                                    n_axes = n_axes_pcoa)
    saveRDS(pcoa_res, file.path(out_dir, "beta_diversity", "pcoa_results.rds"))


    tsne_perplexity_all <- max(2L, min(15L, length(disp_res$distances$sample_id) - 1L))
    tsne_res <- phiper:::compute_tsne(ps = ps,
                                    dist_obj = dist_bc,
                                    dims = 2L,
                                    perplexity = tsne_perplexity_all,
                                    meta_cols = group_col)
    openxlsx::write.xlsx(tsne_res, 
                        file = file.path(out_dir, "beta_diversity", "tsne2d_results.xlsx"),
                        rowNames = TRUE)


    svglite::svglite(
    file.path(out_dir, "beta_diversity", "tsne2d_plot.svg"),
    width = 20/2.54,
    height = 20/2.54,
    bg = "white"
    )

    p_tsne2d <- phiper:::plot_tsne(tsne_res,
                                view = "2d",
                                colour = group_col,
                                palette = group_palette) +
    scale_color_manual(
        values = group_palette,
        breaks = groups,
        name = "Groups"
    ) +
        theme(
        axis.text  = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)
        )
    p_tsne2d$data <- p_tsne2d$data %>% 
    dplyr::filter(.data[[group_col]] %in% groups) %>% 
    dplyr::mutate(
        !!group_col := factor(.data[[group_col]], levels = groups)
    )
    print(p_tsne2d)
    dev.off()

    tsne_perplexity_all <- max(2L, min(20L, length(disp_res$distances$sample_id) - 1L))
    tsne_res <- phiper:::compute_tsne(ps = ps,
                                    dist_obj = dist_bc,
                                    dims = 3L,
                                    perplexity = tsne_perplexity_all,
                                    meta_cols = group_col)
    tsne_res <- tsne_res %>% 
    dplyr::filter(.data[[group_col]] %in% groups) %>% 
    dplyr::mutate(
        !!group_col := factor(.data[[group_col]], levels = groups)
    )
    openxlsx::write.xlsx(tsne_res, file = file.path(out_dir,
                                                    "beta_diversity",
                                                    "tsne3d_results.xlsx"),
                        rowNames = TRUE)

    p3d <- phiper:::plot_tsne(tsne_res,
                            view = "3d",
                            colour = group_col,
                            palette = group_palette)
    htmlwidgets::saveWidget(p3d, file = file.path(out_dir, "beta_diversity",
                                                "tsne3d_plot.html"),
                            selfcontained = TRUE)

    # add group information to PCoA sample coordinates
    pcoa_res$sample_coords <- pcoa_res$sample_coords %>%
    dplyr::left_join(
        ps$data_long %>%
        dplyr::select(sample_id, all_of(group_col) ) %>%
        dplyr::distinct(),
        by = "sample_id",
        copy = TRUE
    )

    lab_perm <- paste0("PERMANOVA p = ", format_pval(permanova_res$p_adjust))
    lab_disp <- paste0("Dispersion p = ", format_pval(disp_res$tests$p_adjust))


    # PCoA plot with group centroids and ellipses
    svglite::svglite(
    file.path(out_dir, "beta_diversity", "pcoa_plot.svg"),
    width = 20/2.54,
    height = 20/2.54,
    bg = "white"
    )
    p_pcoa <- phiper:::plot_pcoa(pcoa_res,
                                axes = c(1, 2),
                                group_col = group_col,
                                ellipse_by = "group",
                                show_centroids = TRUE,
                                point_size = 2,
                                centroid_size=4) +
    theme(
        axis.title = element_text(size =14),
        axis.text = element_text(size =14),
        legend.title = element_blank(),
        legend.text = element_text(size = 13)
    ) +
    annotate(
        "text",
        x = Inf,
        y = -Inf,
        label = paste(lab_perm, lab_disp, sep = "\n"),
        hjust = 1.1,  # nudge left a bit
        vjust = -0.1,  # nudge up a bit
        size = 4
    ) +
    scale_color_manual(
        values = group_palette,
        breaks = groups,
        name = "Groups"
    )

    ellipse_idx <- which(vapply(p_pcoa$layers, function(x) inherits(x$stat, "StatEllipse"), logical(1)))
    ellipse_groups <- sapply(ellipse_idx, function(i) {
    dat_i <- p_pcoa$layers[[i]]$data
    if (is.null(dat_i) || !group_col %in% names(dat_i)) return(NA_character_)
    unique(as.character(dat_i[[group_col]]))[1]
    })
    for (j in seq_along(ellipse_idx)) {
    g <- ellipse_groups[j]
    if (!is.na(g) && g %in% names(group_palette)) {
        p_pcoa$layers[[ellipse_idx[j]]]$aes_params$colour <- unname(group_palette[g])
      }
    }

    print(p_pcoa)
    dev.off()

    # scree plot for first 15 axes of PCoA
    svglite::svglite(
    file.path(out_dir, "beta_diversity", "scree_plot.svg"),
    width = 20/2.54,
    height = 20/2.54,
    bg = "white"
    )

    p_scree <- phiper:::plot_scree(pcoa_res,
                                n_axes = tsne_perplexity_all,
                                type = "line") +
                                    theme(
        axis.text  = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)
        )
    print(p_scree)
    dev.off()

    # determine which contrast label actually exists in the dispersion object
    available_contrasts <- unique(disp_res$distances$contrast)

    # preferred pairwise label: var1 vs var2 (e.g. "control vs MCI")
    pair_contrast <- NA

    if (pair_contrast %in% available_contrasts) {
    contrast_to_use <- pair_contrast
    } else if ("<global>" %in% available_contrasts) {
    # fallback: use global dispersion if pairwise distances are not stored
    contrast_to_use <- "<global>"
    } else {
    # last resort: just take the first available contrast and warn
    contrast_to_use <- available_contrasts[1]
    message(
        "Warning: requested contrast '", pair_contrast,
        "' not found in disp_res$distances$contrast. Using '",
        contrast_to_use, "' instead."
    )
    }

    svglite::svglite(
    file.path(out_dir, "beta_diversity", "dispersion_plot.svg"),
    width = 20/2.54,
    height = 20/2.54,
    bg = "white"
    )

    p_disp <- phiper:::plot_dispersion(
    disp_res,
    scope        = "group",
    contrast     = contrast_to_use,
    show_violin  = TRUE,
    show_box     = TRUE,
    show_points  = TRUE
    )
    p_disp$data <- p_disp$data %>%
    dplyr::filter(level %in% groups) %>%
    dplyr::mutate(level = factor(level, levels = groups))

    #cpair_group  <- comparisons
    # pv <- p_disp$data %>%
    # rstatix::pairwise_wilcox_test(
    #     distance ~ level,
    #     comparisons = pair_group,
    #     p.adjust.method = "BH"
    # ) %>%
    # rstatix::add_xy_position(x = "level") %>%
    # mutate(label = sapply(p.adj, format_pval))

    p_disp <- p_disp +
    labs(x = NULL) +
    scale_colour_manual(values = group_palette,
                        breaks = groups,
                        name = "Groups") +
    scale_fill_manual(values = group_palette,
                        breaks = groups,
                        name = "Groups") +
        theme(
        axis.text  = element_text(size = 14),
        axis.text.x = element_text(
            angle = 45,
            vjust = 1,
            hjust = 1
        ),
        axis.title = element_text(size = 14),
        legend.position = "none",
        legend.text = element_text(size = 13)
        ) #+ 
    # stat_pvalue_manual(
    #     pv,
    #     label = "label",
    #     tip.length = 0.01,
    #     bracket.size = 0.5
    # )
    p_disp$layers$geom_boxplot$stat_params$width <- 0.4
    print(p_disp)
    dev.off()
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
  #cmp_idx <- 1
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
  # extra_vars_to_save <- ps_cmp$data_long %>%
  #   colnames() %>%
  #   setdiff(c("peptide_id", "sample_id", "exist", "group_char"))
  
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
  # enrichment counts
  # ----------------------------------------------------------------------------
  svglite::svglite(
    file.path(out_dir, "enrichment_counts.svg"),
    width = 20/2.54,
    height = 13/2.54,
    bg = "white"
  )
  p_enrich <- plot_enrichment_counts(ps_cmp, 
                                     group_cols = "group_char",
                                     custom_colors = group_palette, annotation_size = 4) +
    labs(title = "Enrichment counts") +
    theme(text = element_text(family = "Montserrat")) +
    facet_wrap(
      vars(forcats::fct_relevel(as.factor(Cohort), var1, var2)),
      ncol = 2,
      scales = "free_x"
    )
  print(p_enrich)
  dev.off()
  
  # ----------------------------------------------------------------------------
  # alpha diversity
  # ----------------------------------------------------------------------------
  carry_cols  <- c("Sex", "Age")
  if (!is.null(paired_col)) carry_cols  <-  c(carry_cols, paired_col)
  alpha_div <- phiper::compute_alpha(ps_cmp,
                             group_cols = "group_char",
                             carry_cols = carry_cols
                             )
  
  dir.create(file.path(out_dir, "alpha_diversity"), 
             recursive = T, showWarnings = F)
  write.xlsx(alpha_div$group_char, file.path(out_dir, 
                                             "alpha_diversity", "table.xlsx"))
  
  pv <- compute_alpha_pval(
    data = alpha_div$group_char,
    metric = "richness",
    group_col = "group_char",
    comparisons = list(cmp),
    paired_id_col = paired_col
  )
  
  svglite::svglite(
    file.path(out_dir, "alpha_diversity", "plot.svg"),
    width = 20/2.54,
    height = 13/2.54,
    bg = "white"
  )
  p_alpha <- plot_alpha(alpha_div,
                        metric = "richness",
                        group_col = "group_char",
                        x_order = cmp,
                        custom_colors=group_palette,
                        text_size = 14) +
    stat_pvalue_manual(
      pv,
      label = "p",
      y.position = "y.position",
      tip.length = 0.01,
      label.size = 4.5,
      family = "Montserrat"
    ) +
    coord_cartesian(clip = "off") +
    theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )
    )
  p_alpha$layers$geom_boxplot$stat_params$width <- 0.4
  p_alpha$layers$geom_point$position$width <- 0.1
  print(p_alpha)
  dev.off()
  
  # species alpha diversity
  alpha_div <- compute_alpha(ps_cmp,
                             group_cols = "group_char",
                             carry_cols = carry_cols,
                             ranks = "species")
  write.xlsx(alpha_div$group_char, file.path(out_dir, 
                                             "alpha_diversity", "table_species.xlsx"))
  # Calculate and format p-values first
  pv <- compute_alpha_pval(
    data = alpha_div$group_char,
    metric = "richness",
    group_col = "group_char",
    comparisons = list(cmp),
    paired_id_col = paired_col
  )
  
  svglite::svglite(
    file.path(out_dir, "alpha_diversity", "plot_species.svg"),
    width = 20/2.54,
    height = 13/2.54,
    bg = "white"
  )
  p_alpha_sh <- plot_alpha(
    alpha_div,
    #metric = "shannon_diversity",
    group_col = "group_char",
    custom_colors=group_palette,
    x_order = cmp,
    text_size = 14
  )
  # if (!is.null(paired_col)) {
  #   p_alpha_sh <- p_alpha_sh +
  #     geom_line(aes(group = .data[[paired_col]]), colour = "grey60", linewidth = 0.35, alpha = 0.6)
  # }
  p_alpha_sh <- p_alpha_sh +
    stat_pvalue_manual(
      pv,
      label = "p",
      y.position = "y.position",
      tip.length = 0.01,
      label.size = 4.5,
      family = "Montserrat"
    ) +
    coord_cartesian(clip = "off") +
    theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )
    ) 
  p_alpha_sh$layers$geom_boxplot$stat_params$width <- 0.4
  p_alpha_sh$layers$geom_point$position$width <- 0.1
  print(p_alpha_sh)
  dev.off()
  
  # ----------------------------------------------------------------------------
  # beta diversity
  # ----------------------------------------------------------------------------
  group_n <- ps_cmp$data_long%>%
    dplyr::filter(group_char %in% c(var1, var2)) %>%
    dplyr::distinct(sample_id, group_char) %>%
    dplyr::count(group_char, name = "n") %>% 
    collect() %>% 
    dplyr::mutate(label = paste0(group_char, "\n(n = ", n, ")")) %>%
    dplyr::select(group_char, label) %>%
    tibble::deframe()
  
  dir.create(file.path(out_dir, "beta_diversity"), 
             recursive = T, showWarnings = F)
  
  dist_bc <- phiper:::compute_distance(ps_cmp, value_col = "exist",
                                       method_normalization = "auto",
                                       distance = "kulczynski", n_threads = DIST_CORES)     
  dist_mat <- as.matrix(dist_bc)
  openxlsx::write.xlsx(
    dist_mat,
    file = file.path(out_dir, "beta_diversity", "distance_matrix.xlsx"),
    rowNames = TRUE
  )

  
  # PCOA ------------
  if (!is.null(paired_col)) {
    
    beta_cols <- parse_longitudinal_group_col(
      group_col = group_col,
      paired_col = paired_col
    )
    
    permanova_res <- phiper:::compute_permanova(
      dist_bc,
      ps = ps_cmp,
      group_col = beta_cols$group_col,
      time_col = beta_cols$time_col,
      subject_col = paired_col
    )
    
    disp_res <- phiper:::compute_dispersion(
      dist_bc,
      ps = ps_cmp,
      group_col = beta_cols$group_col,
      time_col = beta_cols$time_col,
      subject_col = paired_col,
      permutations = 9999
    )
    
    disp_res$distances <- disp_res$distances %>%
      dplyr::left_join(
        ps_cmp$data_long %>%
          dplyr::select(sample_id, group_char) %>%
          dplyr::distinct() %>%
          dplyr::collect(),
        by = "sample_id"
      ) %>% 
      dplyr::mutate(
        level_original = level,
        level = group_char
      )
    
    perplex <- length(disp_res$distances$sample_id) / 2 - 1
    
    beta_meta_cols <- c("group_char", paired_col)
    
    
  } else {
    permanova_res <- phiper:::compute_permanova(
      dist_bc,
      ps = ps_cmp,
      group_col = "group_char"
    )
    disp_res <- phiper:::compute_dispersion(
      dist_bc,
      ps = ps_cmp,
      group_col = "group_char",
      permutations = 9999
    )
    perplex <- length(disp_res$distances$sample_id) - 1
    
    beta_meta_cols <- c("group_char")
  }
  saveRDS(permanova_res, file.path(out_dir, 
                                   "beta_diversity", "permanova_results.rds"))
  saveRDS(disp_res, file.path(out_dir, 
                              "beta_diversity", "dispersion_results.rds"))
  
  n_axes_pcoa <- max(2L, min(109L, length(disp_res$distances$sample_id) - 1L))
  pcoa_res <- phiper:::compute_pcoa(dist_bc,
                                    neg_correction = "none",
                                    n_axes = n_axes_pcoa)
  # add group information to PCoA sample coordinates
  pcoa_res$sample_coords <- pcoa_res$sample_coords %>%
    dplyr::left_join(
      ps_cmp$data_long %>%
        dplyr::select(dplyr::all_of(c("sample_id", beta_meta_cols))) %>%
        dplyr::distinct(),
      by = "sample_id",
      copy = TRUE
    )
  saveRDS(pcoa_res, file.path(out_dir, "beta_diversity", "pcoa_results.rds"))
  
  lab_perm <- paste0("PERMANOVA p = ", format_pval(permanova_res$p_adjust),
                     ", R² = ", round(permanova_res$R2, 3)
                     )
  #lab_perm <- paste0("PERMANOVA p = ", format_pval(permanova_res$p_adjust))
  lab_disp <- paste0("Dispersion p = ", format_pval(disp_res$tests$p_adjust))
  
  # PCoA plot with group centroids and ellipses
  svglite::svglite(
    file.path(out_dir, "beta_diversity", "pcoa_plot.svg"),
    width = 20/2.54,
    height = 20/2.54,
    bg = "white"
  )
  p_pcoa <- phiper:::plot_pcoa(pcoa_res,
                               axes = c(1, 2),
                               group_col = "group_char",
                               ellipse_by = "group",
                               show_centroids = TRUE,
                               point_size = 2,
                               centroid_size=3.5) +
    theme(
      axis.title = element_text(size =14),
      axis.text = element_text(size =14),
      legend.title = element_blank(),
      legend.text = element_text(size = 13)
    ) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = paste(lab_perm, lab_disp, sep = "\n"),
      hjust = 1.1,  # nudge left a bit
      vjust = -0.1,  # nudge up a bit
      size = 4
    ) +
    scale_color_manual(
      values = group_palette,
      breaks = c(var1, var2),
      labels = group_n,
      name = NULL
    )
  
  if (!is.null(paired_col)) {
    
    line_data <- pcoa_res$sample_coords %>%
      dplyr::filter(!is.na(.data[[paired_col]])) %>%
      dplyr::group_by(.data[[paired_col]]) %>%
      dplyr::filter(dplyr::n() >= 2) %>%
      dplyr::ungroup()
    
    p_pcoa <- p_pcoa +
      geom_line(
        data = line_data,
        aes(
          x = PCoA1,
          y = PCoA2,
          group = .data[[paired_col]]
        ),
        inherit.aes = FALSE,
        linewidth = 0.3,
        alpha = 0.35
      )
  }
  # if (!is.null(paired_col)) {
  #   
  #   label_data <- pcoa_res$sample_coords %>%
  #     dplyr::filter(!is.na(.data[[paired_col]])) %>%
  #     dplyr::group_by(.data[[paired_col]]) %>%
  #     dplyr::summarise(
  #       PCoA1 = mean(PCoA1, na.rm = TRUE),
  #       PCoA2 = mean(PCoA2, na.rm = TRUE),
  #       label = dplyr::first(.data[[paired_col]]),
  #       .groups = "drop"
  #     )
  #   
  #   p_pcoa <- p_pcoa +
  #     ggrepel::geom_text_repel(
  #       data = label_data,
  #       aes(
  #         x = PCoA1,
  #         y = PCoA2,
  #         label = label
  #       ),
  #       inherit.aes = FALSE,
  #       size = 3,
  #       max.overlaps = 50
  #     )
  # }
  tmp <- c(var1, var2)
  ellipse_idx <- which(vapply(p_pcoa$layers, function(x) inherits(x$stat, "StatEllipse"), logical(1)))
  ellipse_groups <- sapply(ellipse_idx, function(i) {
    dat_i <- p_pcoa$layers[[i]]$data
    if (is.null(dat_i) || !"group_char" %in% names(dat_i)) return(NA_character_)
    unique(as.character(dat_i[["group_char"]]))[1]
  })
  for (j in seq_along(ellipse_idx)) {
    g <- ellipse_groups[j]
    if (!is.na(g) && g %in% tmp) {
      p_pcoa$layers[[ellipse_idx[j]]]$aes_params$colour <- unname(group_palette[g])
    }
  }
  
  print(p_pcoa)
  dev.off()


  # scree plot for first 15 axes of PCoA
  svglite::svglite(
    file.path(out_dir, "beta_diversity", "scree_plot.svg"),
    width = 20/2.54,
    height = 20/2.54,
    bg = "white"
  )
  
  n_axes_scree <- max(2L, min(15L, length(disp_res$distances$sample_id) - 1L))
  p_scree <- phiper:::plot_scree(pcoa_res,
                                 n_axes = n_axes_scree,
                                 type = "line")  +
    theme(
      axis.text  = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13)
    )
  print(p_scree)
  dev.off()
  
  # dispersion plot ----
  # determine which contrast label actually exists in the dispersion object
  available_contrasts <- unique(disp_res$distances$contrast)
  available_scope <- unique(disp_res$distances$scope)
  # preferred pairwise label: var1 vs var2 (e.g. "control vs MCI")
  pair_contrast <- paste(var1, "vs", var2)
  
  if (pair_contrast %in% available_contrasts) {
    contrast_to_use <- pair_contrast
  } else if ("<global>" %in% available_contrasts) {
    # fallback: use global dispersion if pairwise distances are not stored
    contrast_to_use <- "<global>"
  } else {
    # last resort: just take the first available contrast and warn
    contrast_to_use <- available_contrasts[1]
    message(
      "Warning: requested contrast '", pair_contrast,
      "' not found in disp_res$distances$contrast. Using '",
      contrast_to_use, "' instead."
    )
  }
  
  svglite::svglite(
    file.path(out_dir, "beta_diversity", "dispersion_plot.svg"),
    width = 20/2.54,
    height = 20/2.54,
    bg = "white"
  )
  p_disp <- phiper:::plot_dispersion(
    disp_res,
    scope        = available_scope,
    contrast     = contrast_to_use,
    show_violin  = TRUE,
    show_box     = TRUE,
    show_points  = TRUE
  )
  p_disp$data$level <- factor(p_disp$data$level, levels = c(var1, var2))
  p_disp <- p_disp +
    labs(x = NULL) +
    scale_colour_manual(values = group_palette,
                        breaks = c(var1, var2),
                        name = "Groups") +
    scale_fill_manual(values = group_palette,
                      breaks = c(var1, var2),
                      name = "Groups") + 
    scale_x_discrete(limits = c(var1, var2), labels = group_n) +
    theme(
      axis.text  = element_text(size = 14),
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      ),
      axis.title = element_text(size = 14),
      legend.position = "none",
      legend.text = element_text(size = 13)
    ) 
  p_disp$layers$geom_boxplot$stat_params$width <- 0.4
  print(p_disp)
  dev.off()
  
  # tSNE -----
  tsne_perplexity <- max(2L, min(10L, floor(perplex)))
  tsne_res <- phiper:::compute_tsne(ps = ps_cmp,
                                    dist_obj = dist_bc,
                                    dims = 2L,
                                    max_iter = 4000L,
                                    perplexity = tsne_perplexity,
                                    meta_cols = beta_meta_cols)
  openxlsx::write.xlsx(tsne_res, 
                       file = file.path(out_dir, "beta_diversity", "tsne2d_results.xlsx"),
                       rowNames = TRUE)
  
  
  svglite::svglite(
    file.path(out_dir, "beta_diversity", "tsne2d_plot.svg"),
    width = 20/2.54,
    height = 20/2.54,
    bg = "white"
  )
  
  p_tsne2d <- phiper:::plot_tsne(tsne_res,
                                 view = "2d",
                                 colour = "group_char",
                                 palette = group_palette) +
    scale_color_manual(
      values = group_palette,
      breaks = c(var1, var2),
      labels = group_n[c(var1, var2)],
      name = NULL
    ) +
    theme(
      axis.text  = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13)
    )
  
  if (!is.null(paired_col)) {
    line_data <- tsne_res %>%
      dplyr::filter(!is.na(.data[[paired_col]])) %>%
      dplyr::group_by(.data[[paired_col]]) %>%
      dplyr::filter(dplyr::n() >= 2) %>%
      dplyr::arrange(.data[[paired_col]], group_char) %>%
      dplyr::ungroup()
    
    p_tsne2d <- p_tsne2d +
      geom_line(
        data = line_data,
        aes(
          x = tSNE1,
          y = tSNE2,
          group = .data[[paired_col]]
        ),
        inherit.aes = FALSE,
        linewidth = 0.3,
        alpha = 0.35
      )
  }
  print(p_tsne2d)
  dev.off()
  
  tsne_res <- phiper:::compute_tsne(ps = ps_cmp,
                                    dist_obj = dist_bc,
                                    dims = 3L,
                                    max_iter = 4000L,
                                    perplexity =  tsne_perplexity,
                                    meta_cols = beta_meta_cols)
  openxlsx::write.xlsx(tsne_res, file = file.path(out_dir,
                                                  "beta_diversity",
                                                  "tsne3d_results.xlsx"),
                       rowNames = TRUE)
  
  p3d <- phiper:::plot_tsne(tsne_res,
                            view = "3d",
                            colour = "group_char",
                            palette = group_palette)
  htmlwidgets::saveWidget(p3d, file = file.path(out_dir, "beta_diversity",
                                                "tsne3d_plot.html"),
                          selfcontained = TRUE)
  
  # CAPSCALE -----
  if (!is.null(paired_col)) {
    #library(permute)
    cap_formula <- stats::as.formula(paste0("~ ", beta_cols$time_col, " + Condition(", paired_col, ")"))
    cap_meta_cols <- c("sample_id", "group_char", beta_cols$time_col, paired_col)
    cap_term <- beta_cols$time_col
  } else {
    cap_formula <- stats::as.formula(paste0("~ ", "group_char"))
    cap_meta_cols <- c("sample_id", "group_char")
    cap_term <- "group_char"
    perm_ctrl <- 9999
  }
  
  message("CAP formula used: ", deparse(cap_formula))
  cap_res <- phiper:::compute_capscale(dist_bc,
                                       ps = ps_cmp,
                                       formula = cap_formula,
                                       permutations = 9999)
  saveRDS(cap_res, file.path(out_dir, "beta_diversity", "capscale_results.rds"))
  
  cap_meta <- ps_cmp$data_long %>%
    dplyr::select(dplyr::all_of(cap_meta_cols)) %>%
    dplyr::distinct() %>%
    dplyr::collect()
  cap_plot_df <- cap_res$sample_coords %>%
    dplyr::left_join(
      cap_meta,
      by = "sample_id"
    )
  if (!is.null(paired_col)) perm_ctrl <- permute::how(nperm = 9999, blocks = cap_plot_df[[paired_col]])
  anova_cap <- vegan::anova.cca(
    cap_res$cap_model,
    permutations = perm_ctrl,
    by = "terms"
  )
  cap_p <- anova_cap %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::filter(term == cap_term) %>%
    dplyr::pull(`Pr(>F)`)
  if (length(cap_p) == 0) cap_p <- NA_real_
  # cap_p <- cap_res$perm_terms %>%
  #   dplyr::filter(term == cap_term) %>%
  #   dplyr::pull(`Pr(>F)`)
  # if (length(cap_p) == 0) cap_p <- NA_real_
  if(is.null(paired_col)) cap_term <- group_col
  cap_label <- paste0(
    "CAP ", cap_term, " p = ", format_pval(cap_p),
    "\nR² = ", round(cap_res$r2, 3),
    "; adj. R² = ", round(cap_res$r2_adj, 3)
  )
  
  p_cap <- ggplot(cap_plot_df, aes(x = group_char, y = CAP1, colour = group_char, fill = group_char)) +
    geom_boxplot(width = 0.4, alpha = 1, outlier.shape = NA, linewidth = 0.5, colour = "black")
  if (!is.null(paired_col)) {
    p_cap <- p_cap +
      geom_line(aes(group = .data[[paired_col]]), colour = "grey60", linewidth = 0.35, alpha = 0.6)
  }
  
  p_cap <- p_cap +  
    geom_point(size = 1.8, alpha = 1, stroke = 0.4, shape = 21, colour = "black") +
               #position = position_jitter(width = 0.08, height = 0)) +
    annotate("text", x = Inf, y = Inf, label = cap_label, hjust = 1.05, vjust = 1.2, size = 4) +
    scale_fill_manual(
      values = group_palette,
      breaks = c(var1, var2),
      labels = group_n,
      name = NULL
    ) +
    labs(x = NULL, y = "CAP1",
         title = if (!is.null(paired_col)) {
           paste0("CAP1 constrained by ", beta_cols$time_col, ", conditioned on subject")
         } else {
           "Constrained ordination by group"
         }
    ) +
    scale_x_discrete(limits = c(var1, var2), labels = group_n) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      plot.margin = margin(10, 20, 10, 10)
    )
  ggsave(
    file.path(out_dir, "beta_diversity", "capscale_plot.svg"),
    plot = p_cap,
    bg = "white",
    device = svglite::svglite,
    width = 20/2.54,
    height = 20/2.54,
    units = "in"
  )
  
  
  if (is_longitudinal){
    # ----------------------------------------------------------------------------
    # Stability analysis
    # ----------------------------------------------------------------------------
    dir.create(file.path(out_dir, "similarity"), 
               recursive = T, showWarnings = F)
    
    if (identical(beta_cols$time_col, beta_cols$group_col)) {
      stop("time_col and group_col must be different columns.")
    }
    df_long <- ps_cmp$data_long %>%
      select(peptide_id, subject_id, sample_id, exist, 
             Timepoint_tmp = all_of(beta_cols$time_col),
             Group_tmp     = all_of(beta_cols$group_col)
      ) %>%
      filter(!is.na(subject_id), !is.na(sample_id), !is.na(peptide_id)) %>%
      mutate(
        Timepoint_tmp = as.character(Timepoint_tmp),
        Group_tmp     = as.character(Group_tmp),
        exist         = as.integer(exist)
      ) %>%
      collect()
    
    # Detect the two compared groups/timepoints
    comparison_levels <- df_long %>%
      distinct(Group_tmp, Timepoint_tmp) %>%
      arrange(Timepoint_tmp, Group_tmp) %>%
      mutate(comp_label = paste0(Group_tmp, " at ", Timepoint_tmp))
    if (nrow(comparison_levels) != 2) {
      stop(
        "Expected exactly two group/timepoint combinations, but found: ",
        paste(comparison_levels$comp_label, collapse = ", ")
      )
    }
    x_group <- comparison_levels$Group_tmp[1]
    x_time  <- comparison_levels$Timepoint_tmp[1]
    x_title_raw <- comparison_levels$comp_label[1]
    
    y_group <- comparison_levels$Group_tmp[2]
    y_time  <- comparison_levels$Timepoint_tmp[2]
    y_title_raw <- comparison_levels$comp_label[2]
    # Sample annotation and ordering
    ann_df <- df_long %>%
      distinct(subject_id, sample_id, Group_tmp, Timepoint_tmp) %>%
      mutate(
        comp_label = paste0(Group_tmp, " at ", Timepoint_tmp)
      ) %>% 
      as.data.frame()
    x_samples_df <- ann_df %>%
      filter(Group_tmp == x_group, Timepoint_tmp == x_time) %>%
      arrange(subject_id)
    y_samples_df <- ann_df %>%
      filter(Group_tmp == y_group, Timepoint_tmp == y_time) %>%
      arrange(subject_id)
    
    if (!identical(x_samples_df$subject_id, y_samples_df$subject_id)) {
      stop("Subject IDs are not in the same order between the two paired groups.")
    }
    all_sample_order <- c(x_samples_df$sample_id, y_samples_df$sample_id)
    x_title <- paste0(x_title_raw, " (n = ", nrow(x_samples_df), ")")
    y_title <- paste0(y_title_raw, " (n = ", nrow(y_samples_df), ")")
    
    # Build peptide x sample binary matrix
    exist_mat <- df_long %>%
      group_by(peptide_id, sample_id) %>%
      summarise(exist = max(exist, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(
        names_from = sample_id,
        values_from = exist,
        values_fill = 0
      )
    peptide_ids <- exist_mat$peptide_id
    exist_mat <- exist_mat %>%
      select(-peptide_id) %>%
      as.matrix()
    rownames(exist_mat) <- peptide_ids
    storage.mode(exist_mat) <- "numeric"
    exist_mat <- exist_mat[, all_sample_order, drop = FALSE]
    
    # Kulczynski sample-sample similarity
    kulczynski_dist <- vegan::vegdist(
      t(exist_mat),
      method = "kulczynski",
      binary = TRUE
    )
    kulczynski_sim <- 1 - as.matrix(kulczynski_dist)
    kulczynski_sim <- kulczynski_sim[y_samples_df$sample_id, x_samples_df$sample_id, drop = FALSE]
    
    
    
    plot_df <- as.data.frame(as.table(kulczynski_sim)) %>%
      rename(
        y_sample   = Var1,
        x_sample   = Var2,
        similarity = Freq
      ) %>%
      mutate(
        x_subject = x_samples_df$subject_id[match(x_sample, x_samples_df$sample_id)],
        y_subject = y_samples_df$subject_id[match(y_sample, y_samples_df$sample_id)],
        
        # same order for x and y gives diagonal bottom-left -> top-right
        x_subject = factor(x_subject, levels = x_samples_df$subject_id),
        y_subject = factor(y_subject, levels = y_samples_df$subject_id),
        comparison = if_else(
          as.character(x_subject) == as.character(y_subject),
          "Within individual",
          "Between individuals"
        ),
        comparison = factor(
          comparison,
          levels = c("Within individual", "Between individuals")
        )
      )
    saveRDS(plot_df, file.path(out_dir, "similarity", "similarity_results.rds"))
    # Full similarity heatmap
    # p <- ggplot(plot_df, aes(x = x_subject, y = y_subject, fill = similarity)) +
    #   geom_tile(color = "grey85", linewidth = 0.15) +
    #   coord_fixed() +
    #   scale_fill_viridis_c(
    #     name = "Kulczynski similarity coefficient",
    #     limits = c(0, 1),
    #     guide = guide_colorbar(
    #       title.position = "right",
    #       title.vjust = 0.5,
    #       barheight = grid::unit(35, "mm"),
    #       barwidth  = grid::unit(5, "mm")
    #     ),
    #     option = "B",
    #     direction = -1
    #   ) +
    #   labs(
    #     x = x_title,
    #     y = y_title
    #   ) +
    #   theme_minimal(base_size = 13) +
    #   theme(
    #     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    #     axis.text.y = element_text(size = 9),
    #     axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    #     axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    #     panel.grid = element_blank(),
    #     
    #     legend.position = "right",
    #     legend.title = element_text(
    #       face = "bold",
    #       angle = 90,
    #       hjust = 0.5,
    #       vjust = 0.5
    #     ),
    #     legend.text = element_text(size = 10),
    #     plot.margin = margin(2, 2, 2, 2)
    #   )
    # p
    comp_palette_df <- ann_df %>%
      distinct(comp_label, Group_tmp, Timepoint_tmp) %>%
      mutate(group_key = paste(Group_tmp, Timepoint_tmp, sep = "_"),
             color = group_palette[group_key])
    comp_palette <- setNames(comp_palette_df$color, comp_palette_df$comp_label)
    
    x_col <- comp_palette[[x_title_raw]]
    y_col <- comp_palette[[y_title_raw]]
    x_title_md <- make_colored_axis_title(label = x_title, color = x_col, square_size = 14,text_size = 14, baseline_shift = "-1pt")
    y_title_md <- make_colored_axis_title(label = y_title, color = y_col, square_size = 14, text_size = 14, baseline_shift = "-1pt")
    
    p_heatmap <- ggplot(plot_df, aes(x = x_subject, y = y_subject, fill = similarity)) +
      geom_tile(color = "grey85", linewidth = 0.15) +
      coord_fixed() +
      scale_fill_viridis_c(
        option = "inferno",
        direction = -1,
        name = "Kulczynski similarity coefficient",
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = seq(0, 1, by = 0.2),
        guide = guide_colorbar(
          title.position = "right",
          title.vjust = 0.5,
          barheight = grid::unit(75, "mm"),   # make legend longer
          barwidth  = grid::unit(6, "mm")
        )
      ) +
      labs(
        x = x_title_md,
        y = y_title_md
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.title.x = ggtext::element_markdown(face = "bold", margin = margin(t = 10), lineheight = 1, hjust = 0.5),
        axis.title.y = ggtext::element_markdown(face = "bold", margin = margin(r = 10), lineheight = 1, angle = 90),
        
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 9),
        
        panel.grid = element_blank(),
        
        legend.position = "right",
        legend.title = element_text(
          face = "bold",
          angle = 90,
          hjust = 0.5,
          vjust = 0.5
        ),
        legend.text = element_text(size = 10),
        
        plot.margin = margin(1, 1, 1, 1)
      )
    ggsave(
      file.path(out_dir, "similarity", "kulczynski_heatmap.svg"),
      plot = p_heatmap,
      bg = "white",
      device = svglite::svglite,
      width = 9,
      height = 7.5,
      units = "in"
    )
    
    p_box <- ggplot(plot_df, aes(x = comparison, y = similarity, fill = comparison)) +
      geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.95, colour = "black") +
      geom_jitter(width = 0.1, size = 2, alpha = 1, shape = 21, color = "black", stroke = 0.35) +
      scale_fill_manual(
        values = c(
          "Within individual" = "#420A68FF",
          "Between individuals" = "#FCA50AFF"
        )
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = seq(0, 1, by = 0.2),
      ) +
      labs(
        x = NULL,
        y = "Kulczynski similarity coefficient"
      ) +
      coord_cartesian(clip = "off") +
      theme_classic(base_size = 14) +
      theme(
        legend.position = "none",
        axis.text.x = element_text( size = 14),
        axis.text.y = element_text( size = 14),
        axis.title.y = element_text(face = "bold"),
        panel.grid.major.y = element_line(
          color = "grey85",
          linewidth = 0.35
        )
        #panel.grid.minor = element_blank()
      )
    ggsave(
      file.path(out_dir, "similarity", "kulczynski_boxplot.svg"),
      plot = p_box,
      bg = "white",
      device = svglite::svglite,
      width = 7.5,
      height = 7.5,
      units = "in"
    )
    
    p_box_clean <- p_box +
      labs(y = NULL) +
      scale_x_discrete(
        labels = c(
          "Within individual" = "Within\nindividual",
          "Between individuals" = "Between\nindividuals"
        )
      ) +
      theme(
        axis.title.y = element_blank(),
        aspect.ratio = 1.45,   # important: prevents boxplot from becoming too tall
        plot.margin = margin(1, 1, 1, 1)
      )
    p_combined <- p_heatmap + p_box_clean +
      patchwork::plot_layout(widths = c(1.5, 1))
    ggsave(
      file.path(out_dir, "similarity", "kulczynski_combined.svg"),
      plot = p_combined,
      bg = "white",
      device = svglite::svglite,
      width = 11,
      height = 7.5,
      units = "in"
    )
  }
  
  # ----------------------------------------------------------------------------
  # POP framework
  # ----------------------------------------------------------------------------
  data_frameworks <- readRDS(file.path(out_dir, paste0(label_dir, "_data.rds")))
  if(is.null(paired_col)){
    data_frameworks$group_char <- factor(data_frameworks$group_char, levels = rev(cmp))
  }  
  
  dir.create(file.path(out_dir, "POP_framework"), recursive = TRUE,
             showWarnings = FALSE)
  paired_prev <- if (is_longitudinal && !is.null(paired_col)) paired_col else FALSE
  
  prev_res_pep <- phiper::compute_pop(
    x                 = data_frameworks,
    exist_col         = "exist",
    group_cols        = "group_char",
    rank_cols         = "peptide_id",
    paired            = paired_prev,
    peptide_library   = peplib  
  )
  pep_tbl <- extract_tbl(prev_res_pep)
  write.csv(pep_tbl, file.path(out_dir, "POP_framework", "single_peptide.csv"))
  
  ranks_tax <- c("phylum", "class", "order", "family", "genus", "species")
  prev_res_rank <- phiper::compute_pop(
    x                 = data_frameworks,
    exist_col         = "exist",
    group_cols        = "group_char",
    rank_cols         = ranks_tax,
    paired            = paired_prev,
    peptide_library   = peplib
  )
  rank_tbl <- extract_tbl(prev_res_rank)
  write.csv(rank_tbl, file.path(out_dir, "POP_framework", "taxa_ranks.csv"))
  
  ranks_combined <- c(ranks_tax, "peptide_id")
  plots_dir <- file.path(out_dir, "POP_framework", "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (rank_name in ranks_combined) {
    #rank_name <- "species"
    rank_chr <- as.character(rank_name)
    out_name <- file.path(plots_dir, rank_chr)
    df_rank <- if (rank_chr == "peptide_id") {
      pep_tbl
    } else {
      rank_tbl %>% filter(rank == rank_chr)
    }
    
    df_rank_static <- downsample_for_static(df_rank, prop = 1, seed = 1L)
    p_static <- scatter_static(
      df   = df_rank_static,
      rank = rank_chr,
      xlab = group_n[df_rank$group1[1]],
      ylab = group_n[df_rank$group2[1]],
      point_size       = 2,
      jitter_width_pp  = 0.15,
      jitter_height_pp = 0.15,
      point_alpha      = 0.85,
      font_size        = 14
    ) +
      ggplot2::coord_cartesian(xlim = c(-2, 102), ylim = c(-2, 102),
                               expand = TRUE) +
      ggplot2::theme(
        plot.margin = grid::unit(c(12, 12, 12, 12), "pt"),
        text        = ggplot2::element_text(family = "Montserrat"),
        axis.text  = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 13)
      )
    
    svglite::svglite(
      file.path(paste0(out_name, "_static.svg")),
      width = 20/2.54,
      height = 20/2.54,
      bg = "white"
    )
    print(p_static)
    dev.off()
    
    
    p_inter <- scatter_interactive(
      df   = df_rank,
      rank = rank_chr,
      xlab = group_n[df_rank$group1[1]],
      ylab = group_n[df_rank$group2[1]],
      peplib = peplib,
      point_size = 10,
      jitter_width_pp  = 0.25,
      jitter_height_pp = 0.25,
      point_alpha = 0.85,
      font_size = 14
    )
    
    p_inter <- plotly::layout(
      p_inter,
      autosize = TRUE,
      margin   = list(l = 70, r = 30, t = 10, b = 70),
      xaxis    = list(range = c(-2, 102), automargin = TRUE),
      yaxis    = list(range = c(-2, 102), automargin = TRUE)
    )
    htmlwidgets::saveWidget(
      p_inter,
      file = paste0(out_name, "_interactive.html"),
      selfcontained = TRUE
    )
  }
  
  # ----------------------------------------------------------------------------
  # DELTA framework
  # ----------------------------------------------------------------------------
  dir.create(file.path(out_dir, "DELTA_framework"), recursive = TRUE,
             showWarnings = FALSE)
  
  peplib[] <- lapply(peplib, as.character)
  log_file_current <- if (is.null(LOG_FILE)) {
    file.path(out_dir, "DELTA_framework", "log.txt")
  } else {
    LOG_FILE
  }
  
  if(!is.null(paired_col)){
    stat_mode <- "srlr_paired"
    paired_by <- paired_col
  } else{
    stat_mode <- "srlr"
    paired_by <- NULL
    data_frameworks$subject_id <- data_frameworks$sample_id
  }
  
  #data_frameworks <- readRDS("/home/creyna/Vogl-lab_Projects_git/BC-Engl/Data/allBC/noSARs/results/VitalStatus_stratified_Timepoint/Alive_BL_vs._Alive_FU/Dead_BL_vs._Alive_BL_data.rds")
  
  res <- phiper::compute_delta(
    x                  = data_frameworks,
    exist_col          = "exist",
    rank_cols          = c("class", "order", "family", "genus", "species"),
    group_cols         = "group_char",
    peptide_library    = peplib,
    B_permutations     = 150000L,
    weight_mode        = "equal", #"n_eff_sqrt",
    stat_mode          = stat_mode,
    aggregate_stat     = "af",
    winsor_z           = Inf,
    # rank_feature_keep   = list(
    #   phylum  = NULL, class = NULL, order = NULL, family = NULL, genus = NULL,
    #   species = NULL
    # ),
    log                = LOG,
    log_file           = log_file_current,
    paired_by          = paired_by
  )
  res <- as.data.frame(res)
  write.csv(res, file = file.path(out_dir, "DELTA_framework",
                                  "delta_table.csv"))
  
  # optional: save filtered table too
  res_plot <- res[!is.na(res$m_eff) & res$m_eff > 5, ]
  write.csv(res_plot, 
            file = file.path(out_dir, "DELTA_framework", "delta_table_m_eff_gt5.csv"))
  
  
  tax_ranks <- c("domain", "kingdom", "phylum", "class", "order", "family",
                 "genus", "species")
  # for (taxon in c("order", "family", "genus", "species", "all")){
  #   if(taxon %in% "all"){
  #       std_idx <- res$rank %in% tax_ranks
  #       res$feature[!std_idx] <- as.character(res$rank[!std_idx])
  #       res$rank <- taxon
  #   }
  
  #   svglite::svglite(
  #     file.path(out_dir, "DELTA_framework",
  #               paste0("significant_static_", taxon, "_forestplot.svg")),
  #     width = 20/2.54,
  #     height = 20/2.54,
  #     bg = "white"
  #   )
  
  #   message("Creating forest plot for taxon: ", taxon, ", group1 = ", df_rank$group1[1], ", group2 = ", df_rank$group2[1])
  
  #   p_forest_unc <- phiper::forestplot(
  #     results_tbl          = res,
  #     rank_of_interest     = taxon,
  #     use_diverging_colors = TRUE,
  #     filter_significant   = "p_perm",
  #     left_label           = paste0("More in ", as.character(df_rank$group1[1])),#df_rank$group1[1]),
  #     right_label          = paste0("More in ", as.character(df_rank$group2[1])),#df_rank$group2[1]),
  #     label_vjust           = -0.9,
  #     y_pad                 = 0.3,
  #     label_x_gap_frac      = -0.3,
  #     statistic_to_plot     = "T_stand"
  #   )
  #   print(p_forest_unc)
  #   dev.off()
  
  #   p_inter <- phiper::forestplot_interactive(
  #     results_tbl            = res,
  #     rank_of_interest       = taxon,
  #     statistic_to_plot      = "T_stand",
  #     use_diverging_colors   = TRUE,
  #     filter_significant     = "p_perm",
  #     left_label           = paste0("More in ", as.character(df_rank$group1[1])),#df_rank$group1[1]),
  #     right_label          = paste0("More in ", as.character(df_rank$group2[1])),#df_rank$group2[1]),
  #     arrow_length_frac      = 0.35,
  #     label_x_gap_frac       = -0.3,
  #     label_y_offset        = -0.9
  #   )$plot
  #   htmlwidgets::saveWidget(
  #     p_inter,
  #     file = file.path(out_dir, "DELTA_framework",
  #                     paste0("significant_interactive_", taxon, "_forestplot.html")),
  #     selfcontained = TRUE
  #   )
  # }
  
  for (taxon in c("order", "family", "genus", "species", "all")){
    if(taxon %in% "all"){
      std_idx <- res_plot$rank %in% tax_ranks
      res_plot$feature[!std_idx] <- as.character(res_plot$rank[!std_idx])
      res_plot$rank <- taxon
    }
    
    svglite::svglite(
      file.path(out_dir, "DELTA_framework",
                paste0("significant_static_", taxon, "_forestplot.svg")),
      width = 20/2.54,
      height = 20/2.54,
      bg = "white"
    )
    
    message("Creating forest plot for taxon: ", taxon, ", group1 = ", df_rank$group1[1], ", group2 = ", df_rank$group2[1])
    
    p_forest_unc <- phiper::forestplot(
      results_tbl          = res_plot,
      rank_of_interest     = taxon,
      use_diverging_colors = TRUE,
      filter_significant   = "p_perm",
      left_label           = paste0("More in ", as.character(df_rank$group1[1])),#df_rank$group1[1]),
      right_label          = paste0("More in ", as.character(df_rank$group2[1])),#df_rank$group2[1]),
      label_vjust           = -0.9,
      y_pad                 = 0.3,
      label_x_gap_frac      = -0.3,
      statistic_to_plot     = "T_stand"
    )
    print(p_forest_unc)
    dev.off()
    
    p_inter <- phiper::forestplot_interactive(
      results_tbl            = res_plot,
      rank_of_interest       = taxon,
      statistic_to_plot      = "T_stand",
      use_diverging_colors   = TRUE,
      filter_significant     = "p_perm",
      left_label           = paste0("More in ", as.character(df_rank$group1[1])),#df_rank$group1[1]),
      right_label          = paste0("More in ", as.character(df_rank$group2[1])),#df_rank$group2[1]),
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
  
  # ----------------------------------------------------------------------------
  # interesting features: export + per-feature plots
  # ----------------------------------------------------------------------------
  #always_keep <- c("Staphylococcus aureus", "Norwalk virus")
  always_keep <- c()
  
  res_filtered <- res_plot %>%
    dplyr::mutate(
      .force_keep = (.data$feature %in% always_keep) | (.data$rank %in% always_keep)
    ) %>%
    dplyr::filter(
      .force_keep | (.data$p_perm < 0.05)
    ) %>%
    dplyr::arrange(
      dplyr::desc(.force_keep),
      dplyr::desc(.data$T_obs_stand)
      #dplyr::desc(.data$cross_prev_mean)
    ) %>%
    dplyr::select(-.force_keep)
  
  write.csv(
    res_filtered,
    file      = file.path(out_dir, "DELTA_framework",
                          "delta_table_interesting.csv"),
    row.names = FALSE
  )
  
  tax_cols <- intersect(c(tax_ranks), names(peplib))
  
  res_with_pep <- res_filtered %>%
    mutate(
      match_info      = map(feature, ~ get_binary_and_ids(.x, peplib,
                                                          tax_cols)),
      binary_in_peplib= map(match_info, "present"),
      peptide_ids     = map(match_info, "peptide_ids"),
      pep_tbl_subset  = map(peptide_ids, ~ pep_tbl %>% filter(feature %in% .x))
    ) %>%
    select(-match_info)
  
  dir.create(file.path(out_dir, "DELTA_framework", "interesting_features"),
             recursive = TRUE, showWarnings = FALSE)
  
  
  plot_feature_all <- function(feature_name,
                               group1,
                               group2,
                               feature_data,
                               out_dir) {
    if (is.null(feature_data) || nrow(feature_data) == 0L) {
      message("Skipping ", feature_name, " (no peptide data)")
      return(invisible(NULL))
    }
    safe_name <- gsub("[^A-Za-z0-9_-]+", "_", as.character(feature_name))
    base_dir       <- file.path(out_dir, "DELTA_framework",
                                "interesting_features")
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
    #scatter_dir    <- file.path(base_dir, "scatter")
    #dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)
    file_prefix    <- file.path(base_dir, safe_name)
    #scatter_prefix <- file.path(scatter_dir, safe_name)
    
    SHOW_BG <- TRUE
    BG_MAX_N_INTERACTIVE <- Inf
    BG_SEED <- 1L
    
    bg_df <- NULL
    if (isTRUE(SHOW_BG)) {
      bg_df <- pep_tbl
      
      if ("feature" %in% names(bg_df) && "feature" %in% names(feature_data)) {
        bg_df <- bg_df %>% dplyr::filter(!(feature %in% feature_data$feature))
      }
      
      # keep only unique (percent1, percent2) combos; drop the rest at random
      if (all(c("percent1", "percent2") %in% names(bg_df))) {
        set.seed(BG_SEED)  # makes the random pick reproducible
        bg_df <- bg_df %>%
          dplyr::slice_sample(prop = 1) %>%  # shuffle rows
          dplyr::distinct(percent1, percent2, .keep_all = TRUE)
      }
      
      if (nrow(bg_df) > BG_MAX_N_INTERACTIVE) {
        set.seed(BG_SEED)
        bg_df <- bg_df %>% dplyr::slice_sample(n = BG_MAX_N_INTERACTIVE)
      }
    }
    
    feature_data_static <- downsample_for_static(feature_data, prop = 1, seed = 1L)
    bg_df_static <- downsample_for_static(bg_df, prop = 1, seed = BG_SEED)
    
    ## ---------------- SCATTER STATIC ----------------
    svglite::svglite(
      paste0(file_prefix, "_scatter_static.svg"),
      width = 20/2.54,
      height = 20/2.54,
      bg = "white"
    )
    
    p_scatter <- scatter_static(
      df   = feature_data_static,
      xlab = group1,
      ylab = group2,
      point_size       = 2,
      point_alpha      = 0.85,
      jitter_width_pp  = 0.15,
      jitter_height_pp = 0.15,
      font_family      = "Montserrat",
      font_size        = 14
    ) +
      ggplot2::coord_cartesian(xlim = c(-2, 102),
                               ylim = c(-2, 102),
                               expand = TRUE) +
      ggplot2::theme(
        plot.margin = grid::unit(c(12, 12, 12, 12), "pt"),
        text        = ggplot2::element_text(family = "Montserrat"),
        axis.text  = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 13)
      )  +
      ggplot2::ggtitle(feature_name)
    
    p_scatter <- add_background_static(
      p_scatter,
      bg = bg_df_static,
      size = 1,
      alpha = 0.35
    )
    
    print(p_scatter)
    dev.off()
    
    ## ---------------- SCATTER INTERACTIVE ----------------
    # make sure background cannot “inherit” any categories
    if (!is.null(bg_df) && nrow(bg_df)) {
      bg_df <- bg_df %>%
        dplyr::select(-dplyr::any_of(c("category"))) %>%
        dplyr::mutate(category = "all peptides")
    }
    
    cat_cols <- c(
      "significant (wBH, per rank)" = "#FF1744",
      "nominal only"                = "#00E676",
      "not significant"             = "#2979FF",
      "all peptides"                = "#7A7A7A"
    )
    
    p_inter <- scatter_interactive(
      df = feature_data,
      xlab = group1,
      ylab = group2,
      peplib = peplib,
      show_background   = TRUE,
      background_df     = bg_df,
      background_name   = "all peptides",
      background_color  = "#808080",
      background_alpha  = 0.40,
      background_size   = 7,
      background_max_n  = Inf,
      category_colors = cat_cols,
      point_size  = 11,
      point_alpha = 0.95,
      jitter_width_pp  = 0.05,
      jitter_height_pp = 0.05,
      font_size = 14
    )
    
    p_inter <- plotly::layout(
      p_inter,
      autosize = TRUE,
      margin   = list(l = 70, r = 30, t = 10, b = 70),
      xaxis    = list(range = c(-2, 102), automargin = TRUE),
      yaxis    = list(range = c(-2, 102), automargin = TRUE)
    )
    
    htmlwidgets::saveWidget(
      p_inter,
      file = paste0(file_prefix, "_scatter_interactive.html"),
      selfcontained = TRUE
    )
    
    ## ---------------- DELTA PLOT: conditional smooth ----------------
    use_smooth <- nrow(feature_data) >= 7
    smooth_k   <- if (use_smooth) 3L else 1L
    
    ## ---------------- DELTA PLOT STATIC ----------------
    svglite::svglite(
      paste0(file_prefix, "_deltaplot_static.svg"),
      width = 20/2.54,
      height = 20/2.54,
      bg = "white"
    )

    p_delta_static <- tryCatch(
      {
        deltaplot(
          prev_tbl              = feature_data,
          group_pair_values     = as.character(c(group1, group2)),
          group_labels          = as.character(c(group1, group2)),
          point_jitter_width    = 0.01,
          point_jitter_height   = 0.01,
          point_alpha           = 0.6,
          point_size            = 6,
          add_smooth            = use_smooth,
          smooth_k              = smooth_k,
          arrow_head_length_mm  = 4
        ) +
          ggplot2::ggtitle(
            paste0("Per-peptide shift vs. pooled prevalence in ",feature_name)
          ) +
          ggplot2::theme(
            text        = ggplot2::element_text(family = "Montserrat"),
            title = element_text(size =12),
            axis.text  = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 12)
          ) 
      },
      error = function(e) {
        message(
          "Delta static plot failed for ", feature_name,
          " (", group1, " vs ", group2, "): ", conditionMessage(e)
        )
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::ggtitle(
            paste0("No delta plot for ", feature_name, "\n(", group1, " vs ",
                   group2, ")")
          ) +
          ggplot2::theme(text = ggplot2::element_text(family = "Montserrat"))
      }
    )
    print(p_delta_static)
    dev.off()
    
    ## ---------------- DELTA PLOT INTERACTIVE ----------------
    # p_delta <- tryCatch(
    #   deltaplot_interactive(
    #     prev_tbl            = feature_data,
    #     group_pair_values     = as.character(c(group1, group2)),
    #     group_labels          = as.character(c(group1, group2)),
    #     point_alpha         = 0.6,
    #     point_size          = 6,
    #     add_smooth          = use_smooth,
    #     smooth_k            = smooth_k,
    #     arrow_length_frac   = 0.35,   # old arrow_frac_h
    #     point_jitter_width  = 0.01,
    #     point_jitter_height = 0.01
    #   ),
    #   error = function(e) {
    #     message(
    #       "Delta interactive plot failed for ", feature_name,
    #       " (", group1, " vs ", group2, "): ", conditionMessage(e)
    #     )
    #     NULL
    #   }
    # )
    # if (!is.null(p_delta)) {
    #   htmlwidgets::saveWidget(
    #     p_delta,
    #     file = paste0(file_prefix, "_deltaplot_interactive.html"),
    #     selfcontained = TRUE
    #   )
    # }
    
    ## ---------------- ECDF STATIC ----------------
    svglite::svglite(
      paste0(file_prefix, "_ecdfplot_static.svg"),
      width = 20/2.54,
      height = 20/2.54,
      bg = "white"
    )

    p_ecdf_static <- tryCatch(
      {
        ecdf_plot(
          prev_tbl            = feature_data,
          group_pair_values     = as.character(c(group1, group2)),
          group_labels          = as.character(c(group1, group2)),
          group1_line_color = group_palette[[as.character(group1)]],
          group2_line_color = group_palette[[as.character(group2)]],
          line_width_pt       = 1,
          line_alpha          = 1,
          show_median_lines   = TRUE,
          show_ks_test        = TRUE
        ) +
          ggplot2::ggtitle(
            paste0("ECDF of per-peptide prevalence in ",feature_name)
          ) +
          ggplot2::theme(text = ggplot2::element_text(family = "Montserrat",
                                                      size = 12))
      },
      error = function(e) {
        message(
          "ECDF static plot failed for ", feature_name,
          " (", group1, " vs ", group2, "): ", conditionMessage(e)
        )
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::ggtitle(
            paste0("No ECDF plot for ", feature_name, "\n(", group1, " vs ",
                   group2, ")")
          ) +
          ggplot2::theme(
            text        = ggplot2::element_text(family = "Montserrat"),
            title = element_text(size =12),
            axis.text  = element_text(size = 14),
            axis.title = element_text(size = 14),
            legend.text = element_text(size = 13)
          )
      }
    )
    print(p_ecdf_static)
    dev.off()
    
    ## ---------------- ECDF INTERACTIVE ----------------
    #   p_ecdf <- tryCatch(
    #     ecdf_plot_interactive(
    #       prev_tbl            = feature_data,
    #       group_pair_values     = as.character(c(group1, group2)),
    #       group_labels          = as.character(c(group1, group2)),
    #       group1_line_color = group_palette[[as.character(group1)]],
    #       group2_line_color = group_palette[[as.character(group2)]],
    #       line_width_px       = 2,
    #       line_alpha          = 1,
    #       show_median_lines   = TRUE,
    #       show_ks_test        = TRUE
    #     )
    #     ,
    #     error = function(e) {
    #       message(
    #         "ECDF interactive plot failed for ", feature_name,
    #         " (", group1, " vs ", group2, "): ", conditionMessage(e)
    #       )
    #       NULL
    #     }
    #   )
    #   if (!is.null(p_ecdf)) {
    #     htmlwidgets::saveWidget(
    #       p_ecdf,
    #       file = paste0(file_prefix, "_ecdfplot_interactive.html"),
    #       selfcontained = TRUE
    #     )
    #   }
    
    #   invisible(NULL)
  }
  
  # ------------------------------------------------------------------------------
  # create the plots for interesting_features
  # ------------------------------------------------------------------------------
  n_features <- nrow(res_with_pep)
  for (i in seq_len(n_features)) {
    feature_name <- res_with_pep$feature[i]
    group1       <- res_with_pep$pep_tbl_subset[[i]]$group1[1]
    group2       <- res_with_pep$pep_tbl_subset[[i]]$group2[1]
    feature_data <- res_with_pep$pep_tbl_subset[[i]]
    message("Plotting [", i, "/", n_features, "]: ", feature_name)
    plot_feature_all(feature_name, group1, group2, feature_data, out_dir)
  }
}



# ------------------------------------------------------------------------------
# restore original future plan
# ------------------------------------------------------------------------------
future::plan(original_plan)