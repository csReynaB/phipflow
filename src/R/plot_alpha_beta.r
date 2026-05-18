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
N_CORES <- 30
LOG <- TRUE
LOG_FILE <- NULL
MAX_GB <- 40
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
                                        distance = "kulczynski", n_threads = 8)                                      

    dir.create(file.path(out_dir, "beta_diversity"), 
            recursive = T, showWarnings = F)


    dist_mat <- as.matrix(dist_bc)
    openxlsx::write.xlsx(
    dist_mat,
    file = file.path(out_dir, "beta_diversity", "distance_matrix.xlsx"),
    rowNames = TRUE
    )


    pcoa_res <- phiper:::compute_pcoa(dist_bc,
                                    neg_correction = "none",
                                    n_axes = 109)
    saveRDS(pcoa_res, file.path(out_dir, "beta_diversity", "pcoa_results.rds"))

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

    tsne_res <- phiper:::compute_tsne(ps = ps,
                                    dist_obj = dist_bc,
                                    dims = 2L,
                                    perplexity = min(15, length(disp_res$distances$sample_id) - 1),
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

    tsne_res <- phiper:::compute_tsne(ps = ps,
                                    dist_obj = dist_bc,
                                    dims = 3L,
                                    perplexity = 20,
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
                                n_axes = 15,
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
  alpha_div <- compute_alpha(ps_cmp,
                                       group_cols = "group_char",
                                       carry_cols = c("Sex", "Age"))
  
  dir.create(file.path(out_dir, "alpha_diversity"), 
             recursive = T, showWarnings = F)
  write.xlsx(alpha_div$group_char, file.path(out_dir, 
                                             "alpha_diversity", "table.xlsx"))
  
  pv <- compare_means(
    formula = reformulate("group_char", response = "richness"),
    data = alpha_div$group_char,
    method = "wilcox.test",
    comparisons = list(cmp)
  ) %>%
    mutate(
      p = format_pval(p),
      y.position = max(alpha_div$group_char$richness, na.rm = TRUE) * 1.05
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
                                       carry_cols = c("Sex", "Age"),
                                       ranks = "species")
  write.xlsx(alpha_div$group_char, file.path(out_dir, 
                                             "alpha_diversity", "table_species.xlsx"))
  # Calculate and format p-values first
  pv <- compare_means(
    formula = reformulate("group_char", response = "richness"),
    data = alpha_div$group_char,
    method = "wilcox.test",
    comparisons = list(cmp)
  ) %>%
    mutate(
      p = format_pval(p),
      y.position = max(alpha_div$group_char$richness, na.rm = T) * 1.05
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
  ) +
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
  dist_bc <- phiper:::compute_distance(ps_cmp, value_col = "exist",
                                       method_normalization = "auto",
                                       distance = "kulczynski", n_threads = 8)                                      
  
  dir.create(file.path(out_dir, "beta_diversity"), 
             recursive = T, showWarnings = F)
  
  
  dist_mat <- as.matrix(dist_bc)
  openxlsx::write.xlsx(
    dist_mat,
    file = file.path(out_dir, "beta_diversity", "distance_matrix.xlsx"),
    rowNames = TRUE
  )
  
  
  pcoa_res <- phiper:::compute_pcoa(dist_bc,
                                    neg_correction = "none",
                                    n_axes = 109)
  saveRDS(pcoa_res, file.path(out_dir, "beta_diversity", "pcoa_results.rds"))
  
  cap_res <- phiper:::compute_capscale(dist_bc,
                                       ps = ps_cmp,
                                       formula = ~ group_char)
  saveRDS(cap_res, file.path(out_dir, "beta_diversity", "capscale_results.rds"))
  
  permanova_res <- phiper:::compute_permanova(dist_bc,
                                              ps = ps_cmp,
                                              group_col = "group_char")
  saveRDS(permanova_res, file.path(out_dir, 
                                   "beta_diversity", "permanova_results.rds"))
  
  
  disp_res <- phiper:::compute_dispersion(dist_bc,
                                          ps = ps_cmp,
                                          group_col = "group_char")
  saveRDS(disp_res, file.path(out_dir, 
                              "beta_diversity", "dispersion_results.rds"))
  
  #print(disp_res)
  
  tsne_res <- phiper:::compute_tsne(ps = ps_cmp,
                                    dist_obj = dist_bc,
                                    dims = 2L,
                                    perplexity = min(15, length(disp_res$distances$sample_id) - 1),
                                    meta_cols = c("group_char"))
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
      name = "Groups"
    ) +
    theme(
      axis.text  = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13)
    )
  print(p_tsne2d)
  dev.off()
  
  tsne_res <- phiper:::compute_tsne(ps = ps_cmp,
                                    dist_obj = dist_bc,
                                    dims = 3L,
                                    perplexity =  min(15, length(disp_res$distances$sample_id) - 1),
                                    meta_cols = c("group_char"))
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
  
  # add group information to PCoA sample coordinates
  pcoa_res$sample_coords <- pcoa_res$sample_coords %>%
    dplyr::left_join(
      ps_cmp$data_long %>%
        dplyr::select(sample_id, group_char) %>%
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
      name = "Groups"
    )
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
  
  p_scree <- phiper:::plot_scree(pcoa_res,
                                 n_axes =  min(15, length(disp_res$distances$sample_id) - 1),
                                 type = "line")  +
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
    scope        = "group",
    contrast     = contrast_to_use,
    show_violin  = TRUE,
    show_box     = TRUE,
    show_points  = TRUE
  )
  p_disp$data$level <- factor(p_disp$data$level, levels = c(var1, var2))
  p_disp <- p_disp +
    labs(x = "Groups") +
    scale_colour_manual(values = group_palette,
                        breaks = c(var1, var2),
                        name = "Groups") +
    scale_fill_manual(values = group_palette,
                      breaks = c(var1, var2),
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
    )
  p_disp$layers$geom_boxplot$stat_params$width <- 0.4
  print(p_disp)
  dev.off()
  
}
# ------------------------------------------------------------------------------
# restore original future plan
# ------------------------------------------------------------------------------
future::plan(original_plan)
