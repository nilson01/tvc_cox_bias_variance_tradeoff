# ============================================================
# 00) Packages
# ============================================================
req_pkgs <- c("dplyr","tidyr","tibble","purrr","ggplot2","stringr","readr")
invisible(lapply(req_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install package: ", p)
  library(p, character.only = TRUE)
}))

# ============================================================
# 01) Paths (edit these two only)
# ============================================================
DATA_DIR <- "/Users/nilson/Desktop/internship st jude/final_code/data_uniform"
OUT_DIR  <- file.path("/Users/nilson/Desktop/internship st jude/final_code/1. Plots/plot_uniform")

# DATA_DIR <- "/Users/nilson/Desktop/internship st jude/final_code/data_weibull"



# ============================================================
# 02) USER SETTINGS: n groups (auto-generate BOTH)
# ============================================================
n_small <- c(50, 100)
n_large <- c(250, 500, 1500, 3000)


# Which models to plot (use the labels AFTER recode_model())
PLOT_MODELS <- c("CT", "TH")   # or c("FR","FP") depending on what exists


N_GROUPS <- list(
  small = n_small,
  large = n_large
)

# Which n values to HIDE LABELS for (data still plotted, just no label)
HIDE_LABEL_N_GAMMA_D <- c(100)
HIDE_LABEL_N_GAMMA_0 <- c(100)
HIDE_LABEL_N_CLSR    <- c()
HIDE_LABEL_N_L0      <- c()

# ============================================================
# 01b) DGP tag for filenames (Uniform vs Weibull)
# ============================================================

detect_dgp_tag <- function(data_dir, out_dir) {
  x <- tolower(paste(data_dir, out_dir))
  if (grepl("unif", x) || grepl("uniform", x)) return("Uniform")
  if (grepl("weib", x) || grepl("weibull", x)) return("Weibull")
  return("DGP")  # fallback instead of stop
}


with_dgp_suffix <- function(filename, tag) {
  tools::file_path_sans_ext(filename) %>%
    paste0("_", tag, ".", tools::file_ext(filename))
}

# cat("Using DGP_TAG =", DGP_TAG, "\n")

# ============================================================
# 03) Plot styling knobs
# ============================================================
BASE_SIZE        <- 14
AXIS_TITLE_SIZE  <- 18
AXIS_TEXT_SIZE   <- 16
LEGEND_TEXT_SIZE <- 18
X_ANGLE          <- 45

theme_paper <- function() {
  theme_bw(base_size = BASE_SIZE) +
    theme(
      axis.title  = element_text(size = AXIS_TITLE_SIZE),
      axis.text   = element_text(size = AXIS_TEXT_SIZE),
      axis.text.x = element_text(angle = X_ANGLE, hjust = 1, vjust = 1),
      legend.title = element_text(size = LEGEND_TEXT_SIZE),
      legend.text  = element_text(size = LEGEND_TEXT_SIZE),
      strip.text   = element_text(size = AXIS_TEXT_SIZE)
    )
}

# ============================================================
# 04) Helpers: read rds + model labels
# ============================================================
read_fit_obj <- function(rds_path) {
  obj <- readRDS(rds_path)
  if (!is.list(obj) || !("results" %in% names(obj))) {
    stop("fit_results_summary.rds must be a list with element $results.")
  }
  if (!is.data.frame(obj$results)) {
    stop("$results must be a data.frame / tibble inside fit_results_summary.rds.")
  }
  obj
}

recode_model <- function(x) {
  dplyr::case_when(
    x == "Z" ~ "CT",
    x == "Z + D" ~ "TH",
    x == "Z + X1 + X2" ~ "FR",
    x == "Z + D + X1 + X2" ~ "FP",
    TRUE ~ x
  )
}

add_model_label <- function(df) {
  df %>%
    mutate(
      model_label = recode_model(model),
      model_label = factor(model_label, levels = c("CT","TH","FR","FP"))
    )
}

# ============================================================
# 05) Build long data from obj$results (Bias + MSE + 95% interval)
# ============================================================

# ---- Bias/MSE for logHR_Z or logHR_D ----
build_bias_mse_long_from_results <- function(results_tbl,
                                             target = c("logHR_D","logHR_Z"),
                                             include_n = NULL,
                                             model_labels = PLOT_MODELS) {
  target <- match.arg(target)
  
  df <- results_tbl %>% as_tibble() %>% add_model_label()
  
  if (!is.null(include_n)) df <- df %>% filter(n %in% include_n)
  df <- df %>% filter(model_label %in% model_labels)
  
  if (nrow(df) == 0) return(tibble())
  
  # required columns check
  need_cols <- c("n","model_label","HR","gam_d",
                 if (target == "logHR_Z") "sim_Z" else "sim_D")
  miss <- setdiff(need_cols, names(df))
  if (length(miss) > 0) stop("Missing columns in results_tbl: ", paste(miss, collapse = ", "))
  
  if (target == "logHR_Z") {
    df <- df %>%
      transmute(
        n, model_label,
        sim   = sim_Z,
        truth = log(HR)
      )
  } else {
    df <- df %>%
      transmute(
        n, model_label,
        sim   = sim_D,
        truth = gam_d
      )
  }
  
  df_long <- df %>%
    mutate(diff = purrr::map2(sim, truth, ~ .x - .y)) %>%
    select(-sim, -truth) %>%
    tidyr::unnest_longer(diff, values_to = "diff") %>%
    mutate(
      Bias_val = diff,
      MSE_val  = diff^2
    ) %>%
    tidyr::pivot_longer(c(Bias_val, MSE_val), names_to = "Metric", values_to = "value") %>%
    mutate(Metric = dplyr::recode(Metric, Bias_val = "Bias", MSE_val = "MSE"))
  
  df_long %>%
    group_by(model_label, n, Metric) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      q025 = as.numeric(stats::quantile(value, 0.025, na.rm = TRUE)),
      q975 = as.numeric(stats::quantile(value, 0.975, na.rm = TRUE)),
      .groups = "drop"
    )
}

# ---- CLSR: FP vs FR with t-panels ----
build_clsr_long_from_results <- function(results_tbl,
                                         include_n = NULL,
                                         model_labels = PLOT_MODELS) {
  df <- results_tbl %>% as_tibble() %>% add_model_label()
  if (!is.null(include_n)) df <- df %>% filter(n %in% include_n)
  df <- df %>% filter(model_label %in% model_labels)
  
  if (nrow(df) == 0) return(tibble())
  
  # Expect list-cols for sim + truth
  t_info <- tibble::tibble(
    sim_col = c("sim_CLSR_q_c1","sim_CLSR_q_c2","sim_CLSR_q_c3"),
    tru_col = c("truth_CLSR_q_c1","truth_CLSR_q_c2","truth_CLSR_q_c3"),
    t_panel = c("75th Quantile","90th Quantile","99th Quantile")
  )
  
  purrr::pmap_dfr(t_info, function(sim_col, tru_col, t_panel) {
    miss <- setdiff(c("n","model_label", sim_col, tru_col), names(df))
    if (length(miss) > 0) stop("Missing columns in results_tbl for CLSR: ", paste(miss, collapse = ", "))
    
    tmp <- df %>%
      transmute(
        n, model_label,
        sim = .data[[sim_col]],
        tru = .data[[tru_col]]
      ) %>%
      mutate(diff = purrr::map2(sim, tru, ~ .x - .y)) %>%
      select(-sim, -tru) %>%
      tidyr::unnest_longer(diff, values_to = "diff") %>%
      mutate(
        t_panel = t_panel,
        Bias_val = diff,
        MSE_val  = diff^2
      ) %>%
      tidyr::pivot_longer(c(Bias_val, MSE_val), names_to = "Metric", values_to = "value") %>%
      mutate(Metric = dplyr::recode(Metric, Bias_val = "Bias", MSE_val = "MSE"))
    
    tmp %>%
      group_by(model_label, n, Metric, t_panel) %>%
      summarise(
        mean = mean(value, na.rm = TRUE),
        q025 = as.numeric(stats::quantile(value, 0.025, na.rm = TRUE)),
        q975 = as.numeric(stats::quantile(value, 0.975, na.rm = TRUE)),
        .groups = "drop"
      )
  })
}

# ---- L0: FR/FP with t-panels ----
build_L0_long_from_results <- function(results_tbl,
                                       include_n = NULL,
                                       model_labels = PLOT_MODELS) {
  df <- results_tbl %>% as_tibble() %>% add_model_label()
  if (!is.null(include_n)) df <- df %>% filter(n %in% include_n)
  df <- df %>% filter(model_label %in% model_labels)
  
  if (nrow(df) == 0) return(tibble())
  
  t_info <- tibble::tibble(
    sim_col = c("sim_L0_q_c1","sim_L0_q_c2","sim_L0_q_c3"),
    tru_col = c("truth_L0_q_c1","truth_L0_q_c2","truth_L0_q_c3"),
    t_panel = c("75th Quantile","90th Quantile","99th Quantile")
  )
  
  purrr::pmap_dfr(t_info, function(sim_col, tru_col, t_panel) {
    miss <- setdiff(c("n","model_label", sim_col, tru_col), names(df))
    if (length(miss) > 0) stop("Missing columns in results_tbl for L0: ", paste(miss, collapse = ", "))
    
    tmp <- df %>%
      transmute(
        n, model_label,
        sim = .data[[sim_col]],
        tru = .data[[tru_col]]
      ) %>%
      mutate(diff = purrr::map2(sim, tru, ~ .x - .y)) %>%
      select(-sim, -tru) %>%
      tidyr::unnest_longer(diff, values_to = "diff") %>%
      mutate(
        t_panel = t_panel,
        Bias_val = diff,
        MSE_val  = diff^2
      ) %>%
      tidyr::pivot_longer(c(Bias_val, MSE_val), names_to = "Metric", values_to = "value") %>%
      mutate(Metric = dplyr::recode(Metric, Bias_val = "Bias", MSE_val = "MSE"))
    
    tmp %>%
      group_by(model_label, n, Metric, t_panel) %>%
      summarise(
        mean = mean(value, na.rm = TRUE),
        q025 = as.numeric(stats::quantile(value, 0.025, na.rm = TRUE)),
        q975 = as.numeric(stats::quantile(value, 0.975, na.rm = TRUE)),
        .groups = "drop"
      )
  })
}

# ============================================================
# 06) Plot Bias + MSE stacked (Metric rows)
# ============================================================
plot_bias_mse_stacked <- function(df_long,
                                  y_lab,
                                  outfile, out_dir, tag,
                                  hide_label_n = NULL,
                                  width = 12.8, height = 6.2, dpi = 300) {
  
  if (is.null(df_long) || nrow(df_long) == 0) {
    warning("No data to plot for: ", outfile)
    return(invisible(NULL))
  }
  
  df_long <- df_long %>% mutate(model_label = droplevels(model_label))
  
  x_breaks <- sort(unique(df_long$n))
  x_labels <- ifelse(x_breaks %in% hide_label_n, "", as.character(x_breaks))
  
  p <- ggplot(df_long, aes(x = n, group = model_label, color = model_label)) +
    geom_ribbon(
      aes(ymin = q025, ymax = q975, fill = model_label),
      alpha = 0.15, color = NA
    ) +
    geom_line(aes(y = mean), linewidth = 1.0) +
    geom_point(aes(y = mean, shape = model_label), size = 2.7) +
    facet_grid(Metric ~ ., scales = "free_y") +
    guides(fill = "none") +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(mult = c(0.06, 0.12))
    ) +
    labs(x = "Sample size", y = y_lab, color = "Model", shape = "Model") +
    scale_color_brewer(palette = "Dark2") +
    theme_paper()
  
  print(p)
  # ggsave(file.path(OUT_DIR, with_dgp_suffix(outfile)), p,
  #        width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, with_dgp_suffix(outfile, tag)), p,
         width = width, height = height, dpi = dpi)
  invisible(p)
}

# ============================================================
# 07) CLSR plot
# ============================================================
plot_clsr_bias_mse <- function(df_clsr_long,
                               outfile, out_dir, tag,
                               hide_label_n = NULL, errbar_width = 150, 
                               width = 13.8, height = 7.2, dpi = 300) {
  if (is.null(df_clsr_long) || nrow(df_clsr_long) == 0) {
    warning("No CLSR long data to plot; skipping: ", outfile)
    return(invisible(NULL))
  }
  
  x_breaks <- sort(unique(df_clsr_long$n))
  x_labels <- ifelse(x_breaks %in% hide_label_n, "", as.character(x_breaks))
  
  p <- ggplot(df_clsr_long, aes(x = n, color = model_label, group = model_label)) +
    geom_hline(
      data = subset(df_clsr_long, Metric == "Bias"),
      aes(yintercept = 0),
      linewidth = 0.4, linetype = "dotted", inherit.aes = FALSE
    ) +
    # geom_errorbar(aes(ymin = q025, ymax = q975), width = 0, linewidth = 0.6) +
    # geom_point(aes(y = mean, shape = model_label), size = 2.7) +
    # geom_line(aes(y = mean), linewidth = 0.95) +
    { pd <- position_dodge(width = 0.25); NULL } +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = errbar_width, linewidth = 0.6, position = pd) +
    geom_point(aes(y = mean, shape = model_label), size = 2.7, position = pd) +
    geom_line(aes(y = mean), linewidth = 0.95) +
    facet_grid(Metric ~ t_panel, scales = "free_y") +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(mult = c(0.06, 0.14))
    ) +
    labs(x = "Sample size", y = "Estimation error (CLSR)", color = "Model", shape = "Model") +
    # labs(x = "Sample size", y = "Log Estimation error (CLSR)", color = "Model", shape = "Model") +
    scale_color_brewer(palette = "Dark2") +
    theme_paper()
  
  print(p)
  # ggsave(file.path(OUT_DIR, with_dgp_suffix(outfile)), p,
  #        width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, with_dgp_suffix(outfile, tag)), p,
         width = width, height = height, dpi = dpi)
  invisible(p)
}

# ============================================================
# 07b) L0 plot
# ============================================================
plot_L0_bias_mse <- function(df_L0_long,
                             outfile, out_dir, tag,
                             hide_label_n = NULL,errbar_width = 150,
                             width = 13.8, height = 7.2, dpi = 300) {
  if (is.null(df_L0_long) || nrow(df_L0_long) == 0) {
    warning("No L0 long data to plot; skipping: ", outfile)
    return(invisible(NULL))
  }
  
  x_breaks <- sort(unique(df_L0_long$n))
  x_labels <- ifelse(x_breaks %in% hide_label_n, "", as.character(x_breaks))
  
  p <- ggplot(df_L0_long, aes(x = n, color = model_label, group = model_label)) +
    geom_hline(
      data = subset(df_L0_long, Metric == "Bias"),
      aes(yintercept = 0),
      linewidth = 0.4, linetype = "dotted", inherit.aes = FALSE
    ) +
    # geom_errorbar(aes(ymin = q025, ymax = q975), width = 0, linewidth = 0.6) +
    # geom_point(aes(y = mean, shape = model_label), size = 2.7) +
    # geom_line(aes(y = mean), linewidth = 0.95) +
    { pd <- position_dodge(width = 0.25); NULL } +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = errbar_width, linewidth = 0.6, position = pd) +
    geom_point(aes(y = mean, shape = model_label), size = 2.7, position = pd) +
    geom_line(aes(y = mean), linewidth = 0.95) +
    facet_grid(Metric ~ t_panel, scales = "free_y") +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(mult = c(0.06, 0.14))
    ) +
    labs(
      x = "Sample size",
      y = expression(paste("Estimation error (", Lambda[0], ")")),
      color = "Model",
      shape = "Model"
    ) +
    scale_color_brewer(palette = "Dark2") +
    theme_paper()
  
  print(p)
  # ggsave(file.path(OUT_DIR, with_dgp_suffix(outfile)), p,
  #        width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, with_dgp_suffix(outfile, tag)), p,
         width = width, height = height, dpi = dpi)
  invisible(p)
}

# ============================================================
# 08) Run (auto-generate small + large plots)
# ============================================================

run_one_scenario <- function(data_dir, out_dir) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  rds_path <- file.path(data_dir, "fit_results_summary.rds")
  if (!file.exists(rds_path)) {
    message("SKIP (missing RDS): ", rds_path)
    return(invisible(NULL))
  }
  
  tag <- detect_dgp_tag(data_dir, out_dir)
  message("---- Running scenario ----")
  message("DATA_DIR: ", data_dir)
  message("OUT_DIR : ", out_dir)
  message("DGP_TAG : ", tag)
  
  fit_obj <- read_fit_obj(rds_path)
  res_tbl <- as_tibble(fit_obj$results)
  
  all_n <- sort(unique(res_tbl$n))
  message("Available n in results: ", paste(all_n, collapse = ", "))
  
  for (grp_name in names(N_GROUPS)) {
    
    include_n <- N_GROUPS[[grp_name]]
    err_w <- if (grp_name == "small") 20 else 150
    
    hide_d <- intersect(HIDE_LABEL_N_GAMMA_D, include_n)
    hide_0 <- intersect(HIDE_LABEL_N_GAMMA_0, include_n)
    hide_c <- intersect(HIDE_LABEL_N_CLSR,    include_n)
    hide_l <- intersect(HIDE_LABEL_N_L0,      include_n)
    
    message("\n=============================")
    message("Plot group: ", grp_name)
    message("Including n: ", paste(include_n, collapse = ", "))
    message("=============================")
    
    df_long_d <- build_bias_mse_long_from_results(
      res_tbl, target = "logHR_D", include_n = include_n, model_labels = PLOT_MODELS
    )
    
    plot_bias_mse_stacked(
      df_long_d,
      y_lab        = bquote("Estimation error (" * gamma[d] * ")"),
      outfile      = paste0("fig_bias_mse_logHRD_vs_n_", grp_name, ".png"),
      hide_label_n = hide_d,
      out_dir      = out_dir,
      tag          = tag,
      width        = 13.2, height = 6.4
    )
    
    df_long_z <- build_bias_mse_long_from_results(
      res_tbl, target = "logHR_Z", include_n = include_n, model_labels = PLOT_MODELS
    )
    
    plot_bias_mse_stacked(
      df_long_z,
      y_lab        = bquote("Estimation error (" * gamma[0] * ")"),
      outfile      = paste0("fig_bias_mse_logHRZ_vs_n_", grp_name, ".png"),
      hide_label_n = hide_0,
      out_dir      = out_dir,
      tag          = tag,
      width        = 13.2, height = 6.4
    )
    
    
    
    df_clsr_long_ci <- build_clsr_long_from_results(
      res_tbl, include_n = include_n, model_labels = PLOT_MODELS
    )
    
    
    # =============================================================================
    # =============================================================================
    # plot log-transformed magnitudes for the small-n group (50, 100) 
    # =============================================================================
    # =============================================================================
    
    # 
    # if (grp_name == "small") df_clsr_long_ci <- dplyr::mutate(df_clsr_long_ci,
    #                                                           mean = ifelse(Metric=="MSE", log1p(mean), sign(mean)*log1p(abs(mean))),
    #                                                           q025 = ifelse(Metric=="MSE", log1p(q025), sign(q025)*log1p(abs(q025))),
    #                                                           q975 = ifelse(Metric=="MSE", log1p(q975), sign(q975)*log1p(abs(q975)))
    # )
    # 
    # 
    
    plot_clsr_bias_mse(
      df_clsr_long_ci,
      outfile      = paste0("fig_clsr_bias_mse_fp_vs_fr_", grp_name, ".png"),
      hide_label_n = hide_c,
      errbar_width = err_w,
      out_dir      = out_dir,
      tag          = tag,
      width        = 14.6, height = 7.4
    )
    
    df_L0_long <- build_L0_long_from_results(
      res_tbl, include_n = include_n, model_labels = PLOT_MODELS
    )
    
    plot_L0_bias_mse(
      df_L0_long,
      outfile      = paste0("fig_L0_bias_mse_vs_n_", grp_name, ".png"),
      hide_label_n = hide_l,
      errbar_width = err_w,
      out_dir      = out_dir,
      tag          = tag,
      width        = 14.6, height = 7.4
    )
  }
  
  message("Saved plots to: ", out_dir)
  invisible(TRUE)
}


# ============================================================
# 09) Loop over each scenario folder under ROOT_DIR
# ============================================================
scenario_dirs <- list.dirs(ROOT_DIR, full.names = TRUE, recursive = FALSE)

# ignore folders that start with "0." (like 0.logs)
scenario_dirs <- scenario_dirs[!grepl("/0\\.", scenario_dirs)]

for (scen in scenario_dirs) {
  data_dir <- file.path(scen, "data_uniform")
  out_dir  <- file.path(scen, "plots")
  
  if (!dir.exists(data_dir)) {
    message("SKIP (missing data_uniform): ", data_dir)
    next
  }
  
  tryCatch(
    run_one_scenario(data_dir, out_dir),
    error = function(e) {
      message("ERROR in scenario folder: ", scen)
      message("  ", conditionMessage(e))
      # keep going
      invisible(NULL)
    }
  )
}







