#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

extend_palette_distinct <- function(palette, n = 50, n_candidates = 1000) {
  palette <- unique(palette)

  if (length(palette) >= n) {
    return(palette[seq_len(n)])
  }

  candidates <- grDevices::hcl(
    h = seq(0, 360, length.out = n_candidates),
    c = 80,
    l = rep(c(35, 50, 65), length.out = n_candidates)
  )

  candidates <- unique(candidates)
  candidates <- setdiff(candidates, palette)

  if (length(candidates) == 0L) {
    stop("No candidate colors available to extend palette.", call. = FALSE)
  }

  selected <- palette

  if (length(selected) == 0L) {
    selected <- candidates[1]
    candidates <- candidates[-1]
  }

  while (length(selected) < n) {
    if (length(candidates) == 0L) {
      stop("Ran out of candidate colors before reaching n = ", n, call. = FALSE)
    }

    rgb_selected <- grDevices::col2rgb(selected) / 255
    rgb_candidates <- grDevices::col2rgb(candidates) / 255

    min_dist <- apply(rgb_candidates, 2, function(candidate_rgb) {
      min(colSums((rgb_selected - candidate_rgb)^2))
    })

    best_idx <- which.max(min_dist)
    selected <- c(selected, candidates[best_idx])
    candidates <- candidates[-best_idx]
  }

  selected
}



# ------------------------------------------------------------------------------
# Build active configuration
# ------------------------------------------------------------------------------

build_group_config <- function(active_group_name,
                               default_longitudinal = FALSE,
                               group_definitions,
                               fallback_palette = phip_palette,
                               manual_comparisons = NULL,
                               manual_longitudinal = NULL) {
  
  if (!active_group_name %in% names(group_definitions)) {
    stop(
      "Unknown active_group_name: ", active_group_name,
      ". Available: ", paste(names(group_definitions), collapse = ", ")
    )
  }
  
  group_palette_full <- make_group_palette_full(
    group_definitions = group_definitions,
    fallback_palette = fallback_palette
  )
  
  cfg <- group_definitions[[active_group_name]]
  group_col <- cfg$group_col
  groups <- cfg$groups
  
  if (length(groups) == 0) {
    stop("No groups defined for active_group_name = ", active_group_name)
  }
  
  missing_labels <- setdiff(groups, names(group_palette_full))
  if (length(missing_labels) > 0) {
    stop(
      "These labels are missing from group_palette_full: ",
      paste(missing_labels, collapse = ", ")
    )
  }
  
  # --------------------------------------------------------------------------
  # Decide comparison source
  # Priority:
  # 1. manual_comparisons passed externally
  # 2. comparisons defined inside group_definitions[[active_group_name]]
  # 3. automatic pairwise comparisons
  # --------------------------------------------------------------------------
  
  if (!is.null(manual_comparisons)) {
    
    comparisons <- manual_comparisons
    
    if (!is.null(manual_longitudinal)) {
      longitudinal <- manual_longitudinal
    } else {
      longitudinal <- rep(default_longitudinal, length(comparisons))
    }
    
  } else if (!is.null(cfg$comparisons)) {
    
    comparisons <- cfg$comparisons
    
    if (!is.null(cfg$longitudinal)) {
      longitudinal <- cfg$longitudinal
    } else {
      longitudinal <- rep(default_longitudinal, length(comparisons))
    }
    
  } else {
    
    comparisons <- combn(groups, 2, simplify = FALSE)
    longitudinal <- rep(default_longitudinal, length(comparisons))
  }

  # --------------------------------------------------------------------------
  # Validate comparisons
  # --------------------------------------------------------------------------

  bad_cmp_length <- which(lengths(comparisons) != 2)

  if (length(bad_cmp_length) > 0) {
    stop(
      "Each comparison must contain exactly two group labels. Problem in comparison(s): ",
      paste(bad_cmp_length, collapse = ", ")
    )
  }

  cmp_labels <- unique(unlist(comparisons))
  bad_labels <- setdiff(cmp_labels, groups)

  if (length(bad_labels) > 0) {
    stop(
      "These labels in comparisons are not present in active group '",
      active_group_name, "': ",
      paste(bad_labels, collapse = ", ")
    )
  }

  if (!is.logical(longitudinal)) {
    stop("longitudinal must be a logical vector: TRUE/FALSE.")
  }

  if (length(longitudinal) != length(comparisons)) {
    stop(
      "longitudinal must have the same length as comparisons for active group '",
      active_group_name, "'."
    )
  }

  list(
    active_group_name = active_group_name,
    default_longitudinal = default_longitudinal,
    group_col = group_col,
    groups = groups,
    comparisons = comparisons,
    longitudinal = longitudinal,
    group_palette = group_palette_full[groups],
    group_palette_full = group_palette_full
  )
}

# ------------------------------------------------------------------------------
# Build global palette from all labels across all group definitions
# ------------------------------------------------------------------------------
make_group_palette_full <- function(group_definitions,
                                    fallback_palette = phip_palette) {
  all_labels <- unique(unlist(lapply(group_definitions, `[[`, "groups")))
  n_labels <- length(all_labels)

  if (n_labels == 0L) {
    stop("No group labels found in group_definitions.", call. = FALSE)
  }

  if (n_labels >= 12) {
    palette_idx <- c(seq_len(n_labels)[-12], n_labels + 1L)
  } else {
    palette_idx <- seq_len(n_labels)
  }

  if (max(palette_idx) > length(fallback_palette)) {
    stop(
      "Not enough colors in palette: requested color index ",
      max(palette_idx), " but palette has only ",
      length(fallback_palette), " colors. Extend phip_palette.",
      call. = FALSE
    )
  }

  stats::setNames(
    fallback_palette[palette_idx],
    all_labels
  )
}

# ------------------------------------------------------------------------------
# Other helpers for phiper analysis
# ------------------------------------------------------------------------------


# columns that are always retained when exporting data
base_cols <- c("sample_id", "peptide_id", "group_char", "exist")
special_features <- c(
    "is_IEDB_or_cntrl", "is_auto", "is_infect", "is_EBV",
    "is_toxin", "is_PNP", "is_EM", "is_MPA", "is_patho",
    "is_probio", "is_IgA", "is_flagellum", "signalp6_slow",
    "is_topgraph_new", "is_allergens", "anno_is_fungi", "anno_is_food",
    "anno_is_homo_sapiens", "anno_is_lacto_phage"
)

# Function to format p-values nicely -------------------------------------------
format_pval <- function(p, alpha = 0.05) {
  vapply(p, function(p_one) {
    drop_zeros <- function(x) sub("\\.?0+$", "", x)

    if (is.na(p_one)) return("NA")

    if (p_one > alpha) {
      raw <- formatC(p_one, digits = 2, format = "f")
      raw <- drop_zeros(raw)
      return(paste0("ns [", raw, "]"))
    }

    if (p_one >= 0.001) {
      raw <- formatC(p_one, digits = 3, format = "f")
      raw <- drop_zeros(raw)
      return(raw)
    }

    raw <- formatC(p_one, digits = 2, format = "e")
    raw <- sub("([0-9]+)\\.0+e", "\\1e", raw)
    raw <- sub("([0-9]+\\.[0-9]*[1-9])0+e", "\\1e", raw)

    raw
  }, character(1))
}

# ------------------------------------------------------------------------------
# I/O helpers
# ------------------------------------------------------------------------------
# save an R object to an RDS file, creating parent directories if needed
save_rds_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = path)
}

# select requested columns (base + extras), optionally drop rows with NA in
# selected columns, collect to R, and save to RDS
make_and_save <- function(data, out_path, extra_vars = NULL, drop_na = TRUE) {
  cols_to_select <- unique(c(base_cols, extra_vars %||% character(0)))
  
  avail_cols <- colnames(data$data_long)
  missing_cols <- setdiff(cols_to_select, avail_cols)
  if (length(missing_cols) > 0) {
    message("Skipping missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  df <- data %>%
    dplyr::select(dplyr::any_of(cols_to_select))
  
  if (isTRUE(drop_na)) {
    cols_present <- intersect(cols_to_select, colnames(df$data_long))
    df <- df %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(cols_present), ~ !is.na(.x)))
  }
  
  df <- df %>% dplyr::collect()
  save_rds_safe(df, out_path)
  invisible(df)
}


extract_tbl <- function(obj) {
  if (is.data.frame(obj)) {
    return(tibble::as_tibble(obj))
  }
  for (nm in c("data", "table", "tbl", "df", "result", "results")) {
    if (!is.null(obj[[nm]])) {
      return(tibble::as_tibble(obj[[nm]]))
    }
  }
  out <- try(tibble::as_tibble(obj), silent = TRUE)
  if (!inherits(out, "try-error")) {
    return(out)
  }
  out <- try(as.data.frame(obj), silent = TRUE)
  if (!inherits(out, "try-error")) {
    return(tibble::as_tibble(out))
  }
  stop("Cannot extract a data table from the POP result object.")
}


    
downsample_for_static <- function(df, prop = 0.1, seed = 1L) {
  if (is.null(df) || !nrow(df)) {
    return(df)
  }
  n <- nrow(df)
  size <- max(1L, floor(n * prop))
  if (size >= n) {
    return(df)
  }
  set.seed(seed)
  dplyr::slice_sample(df, n = size)
}


get_binary_and_ids <- function(feature, peplib, tax_cols,
                               peptide_col = "peptide_id") {

  if (!peptide_col %in% names(peplib)) {
    stop("peptide_col not found in peptide library: ", peptide_col, call. = FALSE)
  }

  if (feature %in% special_features && feature %in% names(peplib)) {
    vals <- peplib[[feature]]
    present <- !is.na(vals) & as.logical(vals)

  } else if (feature %in% names(peplib)) {
    vals <- peplib[[feature]]
    present <- as.logical(vals)
    present[is.na(present)] <- FALSE

  } else {
    if (length(tax_cols) == 0L) {
      stop("No taxonomic columns found in peptide library.", call. = FALSE)
    }

    matches <- lapply(tax_cols, function(col) {
      vals <- peplib[[col]]
      !is.na(vals) & vals == feature
    })

    present <- Reduce(`|`, matches)
  }

  peptide_ids <- as.character(peplib[[peptide_col]][present])
  peptide_ids <- unique(peptide_ids)

  list(present = present, peptide_ids = peptide_ids)
}


add_background_static <- function(p, bg,
                                    size = 0.8,
                                    alpha = 0.12,
                                    color = "#808080") {
  if (is.null(bg) || !nrow(bg)) return(p)
    
  bg_layer <- ggplot2::geom_point(
    data = bg,
    mapping = ggplot2::aes(x = percent1, y = percent2),
    inherit.aes = FALSE,
    color = color,
    size = size,
    alpha = alpha,
    show.legend = FALSE
  )
    
  p$layers <- c(list(bg_layer), p$layers)
  p
}




# ------------------------------------------------------------------------------
# Load optional manual comparison configuration
# ------------------------------------------------------------------------------
load_manual_comparison_config <- function(manual_comparison_file = NULL) {
  if (is.null(manual_comparison_file) || !nzchar(manual_comparison_file)) {
    return(list(
      manual_comparisons = NULL,
      manual_longitudinal = NULL
    ))
  }

  if (!file.exists(manual_comparison_file)) {
    stop("Manual comparison file does not exist: ", manual_comparison_file)
  }

  env <- new.env(parent = parent.frame())
  source(manual_comparison_file, local = env)

  if (!exists("manual_comparisons", envir = env)) {
    stop("Manual comparison file must define `manual_comparisons`.")
  }

  manual_comparisons <- get("manual_comparisons", envir = env)

  manual_longitudinal <- if (exists("manual_longitudinal", envir = env)) {
    get("manual_longitudinal", envir = env)
  } else {
    NULL
  }

  if (!is.list(manual_comparisons)) {
    stop("`manual_comparisons` must be a list.")
  }

  invalid_comparisons <- purrr::keep(
    manual_comparisons,
    ~ !is.character(.x) || length(.x) != 2
  )

  if (length(invalid_comparisons) > 0) {
    stop("Each manual comparison must be a character vector of length 2.")
  }

  if (!is.null(manual_longitudinal)) {
    if (!is.logical(manual_longitudinal)) {
      stop("`manual_longitudinal` must be logical: TRUE/FALSE.")
    }

    if (length(manual_longitudinal) != length(manual_comparisons)) {
      stop("`manual_longitudinal` must have the same length as `manual_comparisons`.")
    }
  }

  list(
    manual_comparisons = manual_comparisons,
    manual_longitudinal = manual_longitudinal
  )
}

  

compute_alpha_pval <- function(data,
                                 metric = "richness",
                                 group_col = "group_char",
                                 comparisons,
                                 paired_id_col = NULL,
                                 method = "wilcox.test") {
  purrr::map_dfr(comparisons, function(cmp) {
      
    g1 <- cmp[1]
    g2 <- cmp[2]
      
    df <- data %>%
      dplyr::filter(.data[[group_col]] %in% c(g1, g2)) %>%
      dplyr::mutate(
        !!group_col := factor(.data[[group_col]], levels = c(g1, g2))
      )
      
    if (is.null(paired_id_col)) {
        
      pv <- ggpubr::compare_means(
        formula = stats::reformulate(group_col, response = metric),
        data = df,
        method = method,
        comparisons = list(c(g1, g2))
      )
        
    } else {
        
      df_wide <- df %>%
        dplyr::select(
          dplyr::all_of(c(paired_id_col, group_col, metric))
        ) %>%
        tidyr::drop_na(
          dplyr::all_of(c(paired_id_col, group_col, metric))
        ) %>%
        tidyr::pivot_wider(
          names_from = dplyr::all_of(group_col),
          values_from = dplyr::all_of(metric)
        ) %>%
        tidyr::drop_na(
          dplyr::all_of(c(g1, g2))
        )
        
      if (nrow(df_wide) < 2) {
        pval <- NA_real_
      } else {
        pval <- stats::wilcox.test(
          df_wide[[g1]],
          df_wide[[g2]],
          paired = TRUE,
          exact = FALSE
        )$p.value
      }
        
      pv <- tibble::tibble(
        .y. = metric,
        group1 = g1,
        group2 = g2,
        p = pval,
        method = "Wilcoxon paired test"
      )
    }
      
    y_pos <- max(df[[metric]], na.rm = TRUE)
    if (!is.finite(y_pos)) y_pos <- 1

    pv %>%
      dplyr::mutate(
        p = format_pval(p),
        y.position = y_pos * 1.05
      )
  })
}
  
  
  
parse_longitudinal_group_col <- function(group_col, paired_col = NULL) {
      
  if (is.null(paired_col)) {
    return(list(
      group_col = "group_char",
      time_col = NULL
    ))
  }
    
  # Split only at the LAST underscore
  base_group_col <- sub("_[^_]+$", "", group_col)
  time_col <- sub("^.*_", "", group_col)
    
  list(
    group_col = base_group_col,
    time_col = time_col
  )
}
  
  
make_colored_axis_title <- function(label, color,
                                      square_size = 14,
                                      text_size = 13,
                                      baseline_shift = "-1pt") {
  paste0(
    "<span style='color:", color,
    "; font-size:", square_size, "pt;",
    " baseline-shift:", baseline_shift, ";'>&#9632;</span>",
    " ",
    "<span style='font-size:", text_size, "pt;'><b>", label, "</b></span>"
  )
}
  

safe_rescale01 <- function(x) {
    if (all(is.na(x))) return(rep(NA_real_, length(x)))
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(1, length(x)))
    scales::rescale(x, to = c(0, 1), from = rng)
}