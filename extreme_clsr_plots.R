# =============================================================================
# BLOW-UP ANALYSIS (BATCH MODE): per-scenario folders
# Reads from:   <scenario>/data_uniform/
# Writes to:    <scenario>/plots/
# =============================================================================

# -----------------------------
# 00) Packages
# -----------------------------
req_pkgs <- c("dplyr", "purrr", "tibble", "ggplot2", "readr", "stringr")
invisible(lapply(req_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install package: ", p)
  library(p, character.only = TRUE)
}))

# Optional packages used later (combined plot)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
has_cowplot   <- requireNamespace("cowplot", quietly = TRUE)
if (has_patchwork) library(patchwork)

# -----------------------------
# 01) ROOT PATH (EDIT THIS)
# -----------------------------

ROOT_DIR <- "/Users/nilson/Desktop/internship st jude/final_code/ALL_RESULTS/alpha_0_DGP/Uniform"


# Decide data subfolder name from ROOT_DIR
DATA_SUBDIR <- if (grepl("weib", tolower(ROOT_DIR))) "data_weibull" else "data_uniform"

# -----------------------------
# 02) User settings (global)
# -----------------------------


# Update this depending on result i am running 
FP_MODEL_NAME   <- "Z + D + X1 + X2"
# FP_MODEL_NAME   <- "Z + D"

# threshold_clsr  <- 100
threshold_high <- 100
threshold_low  <- 0.001

sample_ns       <- c(50, 100, 250, 500, 1500, 3000)

# Bin definitions (manual)
bin_breaks <- c(-Inf, 5, 10, 15, 20, 30, 40, 60, 100, 200, 400, 800, Inf)

# bin_labels <- c("0–5", "6–10", "11–15", "16–20", "21–30", "31–40",
#                 "41–60", "61–100", "101–200", "201–400", "401–800", ">800")



# Which event count for boxplots
EVENT_TYPE <- "post_events_upto_q99"  # or "events_post"

# Quantile index: q1/q2/q3 => q99 is typically row 3
q99_index <- 3L

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# 03) Helpers: robust structure access
# -----------------------------
detect_dgp_tag <- function(data_dir, out_dir) {
  x <- tolower(paste(data_dir, out_dir))
  if (grepl("unif", x) || grepl("uniform", x)) return("Uniform")
  if (grepl("weib", x) || grepl("weibull", x)) return("Weibull")
  return("DGP")
}

with_dgp_suffix <- function(filename, tag) {
  tools::file_path_sans_ext(filename) %>%
    paste0("_", tag, ".", tools::file_ext(filename))
}

get_post_upto_q <- function(rs) {
  if (!is.null(rs$post_events_upto_q)) return(rs$post_events_upto_q)
  if (!is.null(rs$post_upto_q))        return(rs$post_upto_q)
  stop("Neither post_events_upto_q nor post_upto_q found. Available names: ",
       paste(names(rs), collapse = ", "))
}

read_rep_summaries <- function(file) {
  x <- readRDS(file)
  
  # New structure: list(..., rep_results = list(...))
  if (is.list(x) && !is.null(x$rep_results)) return(x$rep_results)
  
  # Old structure: list of reps directly (heuristic)
  if (is.list(x) && length(x) > 0 && is.list(x[[1]]) && (!is.null(x[[1]]$b) || !is.null(x[[1]]$rep_id))) return(x)
  
  stop("Unrecognized rep summary structure in: ", file,
       "\nTop-level names: ", paste(names(x), collapse = ", "))
}

get_q99_count <- function(post_vec, q99_index = 3L) {
  if (is.null(post_vec)) return(NA_integer_)
  
  nm <- names(post_vec)
  if (!is.null(nm) && length(nm) > 0) {
    cand <- nm[grepl("99", nm)]
    if (length(cand) == 1L) return(as.integer(post_vec[[cand]]))
    return(as.integer(post_vec[[q99_index]]))
  } else {
    if (length(post_vec) < q99_index) return(NA_integer_)
    return(as.integer(post_vec[[q99_index]]))
  }
}

get_clsr_hat_q99 <- function(rs, fp_model = FP_MODEL_NAME, q_row = "q3", q_index = 3L) {
  m <- rs$est_CLSR_q
  if (is.null(m)) return(NA_real_)
  
  rn <- rownames(m); cn <- colnames(m)
  if (!is.null(rn) && !is.null(cn) && (q_row %in% rn) && (fp_model %in% cn)) {
    return(as.numeric(m[q_row, fp_model]))
  }
  
  j <- match(fp_model, cn)
  # if (is.na(j)) j <- 1L
  return(as.numeric(m[q_index, j]))
}

# Event-count getters
get_total_events <- function(rs) {
  if (!is.null(rs$ev_count)) return(as.integer(rs$ev_count))
  stop("rs$ev_count not found. Available names: ", paste(names(rs), collapse = ", "))
}
get_post_events <- function(rs) {
  if (!is.null(rs$ev_post)) return(as.integer(rs$ev_post))
  stop("rs$ev_post not found. Available names: ", paste(names(rs), collapse = ", "))
}

# -----------------------------
# 04) Plot helpers
# -----------------------------
make_events_boxplot <- function(df, title_suffix = NULL) {
  ggplot2::ggplot(df, ggplot2::aes(x = factor(n), y = events)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.5) +
    ggplot2::facet_wrap(~ dgp_label, scales = "free_y") +
    ggplot2::labs(
      x = "Sample size n",
      y = "Number of events",
      title = title_suffix
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text  = ggplot2::element_text(size = 10)
    )
}

create_blowup_combined_plot <- function(bin_summary, out_dir, tag,
                                        filename = "blowup_combined.png",
                                        width = 10, height = 8, dpi = 300) {
  if (!has_patchwork || !has_cowplot) {
    message("Skipping combined plot (need patchwork + cowplot).")
    return(invisible(NULL))
  }
  
  # Plot 1: COUNT
  p_count <- ggplot2::ggplot(
    bin_summary,
    ggplot2::aes(x = post_bin, y = n_blow, group = n_f, color = n_f)
  ) +
    ggplot2::geom_line(alpha = 0.9) +
    ggplot2::geom_point(size = 2, alpha = 0.9) +
    ggplot2::labs(y = "Number of extreme CLSR", color = "Sample size (n)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(5, 5, 0, 5)
    )
  
  # Plot 2: RATE
  p_rate <- ggplot2::ggplot(
    bin_summary,
    ggplot2::aes(x = post_bin, y = p_blow, group = n_f, color = n_f)
  ) +
    ggplot2::geom_line(alpha = 0.9) +
    ggplot2::geom_point(size = 2, alpha = 0.9) +
    ggplot2::labs(
      x = "Post-treatment event count",
      y = "Rate of extreme CLSR",
      color = "Sample size (n)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.margin = ggplot2::margin(0, 5, 5, 5)
    )
  
  # Legend extractor
  p_legend <- ggplot2::ggplot(
    bin_summary,
    ggplot2::aes(x = post_bin, y = n_blow, group = n_f, color = n_f)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(color = "Sample size (n)") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right")
  
  legend <- cowplot::get_legend(p_legend)
  
  p_combined <- (p_count / p_rate) | legend
  p_combined <- p_combined + patchwork::plot_layout(widths = c(4, 1))
  
  ggplot2::ggsave(
    filename = file.path(out_dir, with_dgp_suffix(filename, tag)),
    plot = p_combined,
    width = width, height = height, dpi = dpi
  )
  
  print(p_combined)
  message("Saved combined plot to: ", file.path(out_dir, with_dgp_suffix(filename, tag)))
  invisible(p_combined)
}

# =============================================================================
# 05) RUN ONE SCENARIO FOLDER
# =============================================================================
run_one_scenario <- function(data_dir, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  rds_summary <- file.path(data_dir, "fit_results_summary.rds")
  if (!file.exists(rds_summary)) {
    message("SKIP (missing fit_results_summary.rds): ", rds_summary)
    return(invisible(NULL))
  }
  
  tag <- detect_dgp_tag(data_dir, out_dir)
  treatment_time_type <- paste(tag, "treatment times")
  
  message("\n====================================================")
  message("Scenario DATA_DIR: ", data_dir)
  message("Scenario OUT_DIR : ", out_dir)
  message("DGP_TAG          : ", tag)
  message("====================================================")
  
  # Helper: per-row rep summaries
  .rep_summaries_file <- function(row_idx) {
    file.path(data_dir, sprintf("rep_summaries_scen_%04d.rds", row_idx))
  }
  
  # Load grid
  obj <- readRDS(rds_summary)
  if (is.null(obj$grid)) stop("fit_results_summary.rds missing $grid: ", rds_summary)
  param_grid <- obj$grid
  
  # Which grid rows to use (based on n)
  row_idxs <- which(param_grid$n %in% sample_ns)
  if (length(row_idxs) == 0) {
    message("SKIP (no grid rows with n in sample_ns for this scenario folder).")
    return(invisible(NULL))
  }
  
  # -----------------------------
  # Extract per-rep (n, post_count_q99, CLSR_hat_q99)
  # -----------------------------
  df_q99_raw <- purrr::map_dfr(row_idxs, function(row_idx) {
    f <- .rep_summaries_file(row_idx)
    if (!file.exists(f)) {
      message("  SKIP row_idx ", row_idx, " (missing rep_summaries file): ", f)
      return(tibble::tibble())
    }
    
    reps <- read_rep_summaries(f)
    
    purrr::map_dfr(reps, function(rs) {
      post_vec       <- get_post_upto_q(rs)
      post_count_q99 <- get_q99_count(post_vec, q99_index = q99_index)
      clsr_hat_q99   <- get_clsr_hat_q99(rs, fp_model = FP_MODEL_NAME, q_row = "q3", q_index = q99_index)
      
      tibble::tibble(
        scenario_id = row_idx,
        rep_id      = rs$b %||% rs$rep_id %||% NA_integer_,
        n           = as.integer(param_grid$n[row_idx]),
        post_count  = post_count_q99,
        CLSR_hat    = clsr_hat_q99
      )
    })
  })
  
  if (nrow(df_q99_raw) == 0) {
    message("SKIP (no rep rows loaded).")
    return(invisible(NULL))
  }
  
  # -----------------------------
  # Prepare + bin
  # -----------------------------
  df_q99 <- df_q99_raw %>%
    dplyr::filter(n %in% sample_ns) %>%
    dplyr::mutate(
      # blow = CLSR_hat > threshold_clsr,
      blow = (CLSR_hat > threshold_high) | (CLSR_hat < threshold_low),
      n_f  = factor(n, levels = sort(sample_ns))
    ) 
  
  
  message("DEBUG: min(post_count used for bins) = ", min(df_q99$post_count, na.rm = TRUE))
  

  
  #------------------------------------
  #------------------------------------
  # first bin label to start at the minimum x value actually present in the data (per scenario),
  #------------------------------------
  #------------------------------------

# 
#   min_x <- min(df_q99$post_count, na.rm = TRUE)
#   
#   # Replace -Inf with the observed minimum, so the first bin starts at min_x
#   bin_breaks_s <- c(min_x, 10, 15, 20, 30, 40, 60, 100, 200, 400, 800, Inf)
# 
#   
#   make_bin_labels_dyn <- function(breaks, min_x) {
#     labs <- character(length(breaks) - 1)
#     
#     for (i in seq_len(length(breaks) - 1)) {
#       lo <- breaks[i]
#       hi <- breaks[i + 1]
#       
#       if (is.infinite(hi)) {
#         labs[i] <- sprintf(">%s", lo)
#         next
#       }
#       
#       # If this is the first bin starting at min_x, label it as just "min_x"
#       if (i == 1L) {
#         labs[i] <- as.character(min_x)
#         next
#       }
#       
#       # Otherwise standard (lo,hi] shown as lo+1–hi
#       labs[i] <- sprintf("%s–%s", lo + 1, hi)
#     }
#     
#     labs
#   }
#   
#   bin_labels_dyn <- make_bin_labels_dyn(bin_breaks_s, min_x)
#   
#   df_binned <- df_q99 %>%
#     dplyr::mutate(
#       post_bin = cut(
#         post_count,
#         breaks = bin_breaks_s,
#         labels = bin_labels_dyn,
#         right = TRUE,
#         include.lowest = TRUE
#       )
#     ) %>%
#     dplyr::filter(!is.na(post_bin)) %>%
#     dplyr::mutate(post_bin = factor(post_bin, levels = bin_labels_dyn, ordered = TRUE))

  
  
  
  
  # --- Dynamic breaks + labels so first x tick is the observed minimum (e.g., 9 or 17) ---
  min_x <- min(df_q99$post_count, na.rm = TRUE)
  if (!is.finite(min_x)) stop("min_x not finite; check post_count.")
  
  # Take your base breaks but replace -Inf with min_x AND also include min_x as a breakpoint
  base_breaks <- c(5, 10, 15, 20, 30, 40, 60, 100, 200, 400, 800, Inf)
  
  # Find the "upper endpoint" of the bin that contains min_x under right=TRUE
  # Example: min_x=17 -> upper endpoint is 20; min_x=9 -> upper endpoint is 10
  upper_end <- base_breaks[which(min_x <= base_breaks)[1]]
  
  # Construct breaks so the first bin is [min_x, upper_end]
  # and then continue with the remaining base breaks after upper_end
  bin_breaks_s <- c(min_x, upper_end, base_breaks[base_breaks > upper_end])
  
  # Remove any accidental duplicates, keep sorted
  bin_breaks_s <- sort(unique(bin_breaks_s))
  
  # Label builder:
  #  - First bin label: "min_x" (just the number)
  #  - Subsequent bins: (lo, hi] displayed as (lo+1)–hi (for integer counts)
  make_bin_labels_from_breaks <- function(breaks, min_x) {
    labs <- character(length(breaks) - 1)
    for (i in seq_len(length(breaks) - 1)) {
      lo <- breaks[i]
      hi <- breaks[i + 1]
      
      if (i == 1L) {
        labs[i] <- as.character(min_x)
      } else if (is.infinite(hi)) {
        labs[i] <- sprintf(">%s", lo)
      } else {
        labs[i] <- sprintf("%s–%s", lo + 1, hi)
      }
    }
    labs
  }
  
  bin_labels_dyn <- make_bin_labels_from_breaks(bin_breaks_s, min_x)
  
  df_binned <- df_q99 %>%
    dplyr::mutate(
      post_bin = cut(
        post_count,
        breaks = bin_breaks_s,
        labels = bin_labels_dyn,
        right = TRUE,
        include.lowest = TRUE
      )
    ) %>%
    dplyr::filter(!is.na(post_bin)) %>%
    dplyr::mutate(post_bin = factor(post_bin, levels = bin_labels_dyn, ordered = TRUE))
  
  
  # df_binned <- df_q99 %>%
  #   dplyr::mutate(
  #     post_bin = cut(
  #       post_count,
  #       breaks = bin_breaks,
  #       labels = bin_labels,
  #       right = TRUE,
  #       include.lowest = TRUE
  #     )
  #   ) %>%
  #   dplyr::filter(!is.na(post_bin)) %>%
  #   dplyr::mutate(post_bin = factor(post_bin, levels = bin_labels, ordered = TRUE))
  # 
  # 
  # 
  
  
  
  bin_summary <- df_binned %>%
    dplyr::group_by(n_f, post_bin) %>%
    dplyr::summarise(
      n_reps = dplyr::n(),
      n_blow = sum(blow, na.rm = TRUE),
      p_blow = n_blow / n_reps,
      .groups = "drop"
    ) %>%
    dplyr::arrange(n_f, post_bin)
  
  # =============================================================================
  # PLOT 1: BINNED BLOW-UP COUNT
  # =============================================================================
  p_count_bin <- ggplot2::ggplot(
    bin_summary,
    ggplot2::aes(x = post_bin, y = n_blow, group = n_f, color = n_f)
  ) +
    ggplot2::geom_line(alpha = 0.9) +
    ggplot2::geom_point(size = 2, alpha = 0.9) +
    ggplot2::labs(
      x     = "Post-treatment event count",
      y     = "Number of extreme CLSR",
      color = "Sample size (n)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  ggplot2::ggsave(
    filename = file.path(out_dir, with_dgp_suffix("blowup_count_binned.png", tag)),
    plot = p_count_bin,
    width = 10, height = 6, dpi = 300
  )
  
  # =============================================================================
  # PLOT 2: BINNED WITHIN-BIN BLOW-UP RATE
  # =============================================================================
  p_rate_bin <- ggplot2::ggplot(
    bin_summary,
    ggplot2::aes(x = post_bin, y = p_blow, group = n_f, color = n_f)
  ) +
    ggplot2::geom_line(alpha = 0.9) +
    ggplot2::geom_point(size = 2, alpha = 0.9) +
    ggplot2::labs(
      x     = "Post-treatment event count",
      y     = "Rate of extreme CLSR",
      color = "Sample size (n)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  ggplot2::ggsave(
    filename = file.path(out_dir, with_dgp_suffix("blowup_rate_binned.png", tag)),
    plot = p_rate_bin,
    width = 10, height = 6, dpi = 300
  )
  
  # Combined plot 
  # create_blowup_combined_plot(bin_summary, out_dir = out_dir, tag = tag)
  
  # =============================================================================
  # EVENT BOXPLOTS (small + large)
  # =============================================================================
  make_dgp_label <- function(pg_row) paste0(treatment_time_type)
  
  ns_to_use <- sort(unique(param_grid$n))
  row_idxs_events <- which(param_grid$n %in% ns_to_use)
  
  df_events <- purrr::map_dfr(row_idxs_events, function(row_idx) {
    f <- .rep_summaries_file(row_idx)
    if (!file.exists(f)) return(tibble::tibble())
    
    reps <- read_rep_summaries(f)
    pg_row <- tibble::as_tibble(param_grid[row_idx, , drop = FALSE])
    dgp_lab <- make_dgp_label(pg_row)
    
    purrr::map_dfr(reps, function(rs) {
      ev <- NA_integer_
      if (EVENT_TYPE == "post_events_upto_q99") {
        post_vec <- get_post_upto_q(rs)
        ev <- get_q99_count(post_vec, q99_index = q99_index)
      } else if (EVENT_TYPE == "events_post") {
        ev <- get_post_events(rs)
      }
      
      tibble::tibble(
        row_idx     = row_idx,
        scenario_id = row_idx,
        rep_id      = rs$b %||% rs$rep_id %||% NA_integer_,
        n           = as.integer(param_grid$n[row_idx]),
        dgp_label   = dgp_lab,
        events      = ev
      )
    })
  })
  
  n_small <- c(50, 100)
  n_large <- c(250, 500, 1500, 3000, 5000)
  
  df_events_small <- dplyr::filter(df_events, n %in% n_small)
  df_events_large <- dplyr::filter(df_events, n %in% n_large)
  
  p_events_small <- make_events_boxplot(df_events_small, title_suffix = "Small sample sizes (n = 50, 100)")
  p_events_large <- make_events_boxplot(df_events_large, title_suffix = "Large sample sizes (n ≥ 250)")
  
  out_name_small <- if (EVENT_TYPE == "post_events_upto_q99") {
    "boxplot_events_post_upto_q99_n_50_100.png"
  } else {
    "boxplot_events_post_n_50_100.png"
  }
  
  out_name_large <- if (EVENT_TYPE == "post_events_upto_q99") {
    "boxplot_events_post_upto_q99_n_250_plus.png"
  } else {
    "boxplot_events_post_n_250_plus.png"
  }
  
  ggplot2::ggsave(
    filename = file.path(out_dir, with_dgp_suffix(out_name_small, tag)),
    plot = p_events_small,
    width = 13, height = 7, dpi = 300
  )
  
  ggplot2::ggsave(
    filename = file.path(out_dir, with_dgp_suffix(out_name_large, tag)),
    plot = p_events_large,
    width = 13, height = 7, dpi = 300
  )
  
  # =============================================================================
  # CSV TABLE: Total events + Post events by n
  # =============================================================================
  row_idxs_tbl <- which(param_grid$n %in% ns_to_use)
  
  df_events_counts <- purrr::map_dfr(row_idxs_tbl, function(row_idx) {
    f <- .rep_summaries_file(row_idx)
    if (!file.exists(f)) return(tibble::tibble())
    
    reps <- read_rep_summaries(f)
    
    purrr::map_dfr(reps, function(rs) {
      tibble::tibble(
        n            = as.integer(param_grid$n[row_idx]),
        scenario_id  = row_idx,
        rep_id       = rs$b %||% rs$rep_id %||% NA_integer_,
        total_events = get_total_events(rs),
        post_events  = get_post_events(rs)
      )
    })
  })
  
  if (nrow(df_events_counts) > 0) {
    events_summary_by_n <- df_events_counts %>%
      dplyr::group_by(n) %>%
      dplyr::summarise(
        n_replications = dplyr::n(),
        
        total_min    = min(total_events, na.rm = TRUE),
        total_median = stats::median(total_events, na.rm = TRUE),
        total_max    = max(total_events, na.rm = TRUE),
        
        post_min     = min(post_events, na.rm = TRUE),
        post_median  = stats::median(post_events, na.rm = TRUE),
        post_max     = max(post_events, na.rm = TRUE),
        
        .groups = "drop"
      ) %>%
      dplyr::arrange(n)
    
    out_csv <- file.path(out_dir, with_dgp_suffix("table_total_and_post_events_by_n.csv", tag))
    readr::write_csv(events_summary_by_n, out_csv)
    message("Saved table to: ", out_csv)
  } else {
    message("Skipping events summary CSV (no df_events_counts rows).")
  }
  
  message("DONE. Saved outputs to: ", out_dir)
  invisible(TRUE)
}

# =============================================================================
# 06) Batch runner over scenario directories
# =============================================================================
scenario_dirs <- list.dirs(ROOT_DIR, full.names = TRUE, recursive = FALSE)
# ignore folders like 0.logs or other dot/utility dirs
scenario_dirs <- scenario_dirs[!grepl("/0\\.", scenario_dirs)]

message("Found ", length(scenario_dirs), " scenario folders under ROOT_DIR.")

for (scen in scenario_dirs) {
  data_dir <- file.path(scen, DATA_SUBDIR)
  out_dir  <- file.path(scen, "plots")
  
  if (!dir.exists(data_dir)) {
    message("SKIP (missing ", DATA_SUBDIR, "): ", data_dir)
    next
  }
  
  tryCatch(
    run_one_scenario(data_dir, out_dir),
    error = function(e) {
      message("ERROR in scenario: ", scen)
      message("  ", conditionMessage(e))
      # continue to next scenario
      invisible(NULL)
    }
  )
}

message("\n========== ALL DONE ==========\n")


