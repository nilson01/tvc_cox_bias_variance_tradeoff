library(survival)
library(dplyr)
library(tibble)

# ============================================================
# 1) Baseline hazard from tt() fit (Breslow-style), robust to NA
# ============================================================
bh_from_tt_safe <- function(fit, dat) {
  cf <- coef(fit)
  if (is.null(cf) || !length(cf)) {
    return(data.frame(time = numeric(0), hazard = numeric(0)))
  }
  cf[!is.finite(cf)] <- 0
  nm <- names(cf)
  
  t_ev <- sort(unique(dat$futime[dat$status == 1]))
  if (!length(t_ev)) return(data.frame(time = numeric(0), hazard = numeric(0)))
  
  lp_at_time <- function(t) {
    lp <- rep(0, nrow(dat))
    
    if ("W" %in% nm)  lp <- lp + cf["W"]  * dat$W
    if ("X1" %in% nm) lp <- lp + cf["X1"] * dat$X1
    if ("X2" %in% nm) lp <- lp + cf["X2"] * dat$X2
    
    if ("tt(T_treat)" %in% nm) {
      Zt <- as.integer(t > dat$T_treat)
      lp <- lp + cf["tt(T_treat)"] * Zt
    }
    if ("tt(D_treat)" %in% nm) {
      Dt <- pmax(0, t - dat$T_treat)
      lp <- lp + cf["tt(D_treat)"] * Dt
    }
    lp
  }
  
  Hcum <- numeric(length(t_ev))
  for (k in seq_along(t_ev)) {
    t  <- t_ev[k]
    dk <- sum(dat$status == 1 & dat$futime == t)
    
    at_risk <- dat$futime >= t
    lp_t <- lp_at_time(t)
    denom <- sum(exp(lp_t[at_risk]))
    
    jump <- if (is.finite(denom) && denom > 0) dk / denom else NA_real_
    Hcum[k] <- if (k == 1) jump else Hcum[k - 1] + jump
  }
  
  data.frame(time = t_ev, hazard = Hcum)
}

# ============================================================
# 2) CLSR, L0, h0 at horizon
# ============================================================
clsr_L0_h0_at_t <- function(bh, gamma0, gamma_d, t_lim) {
  if (is.null(bh) || !nrow(bh) || !is.finite(t_lim)) {
    return(c(CLSR = NA_real_, L0 = NA_real_, h0 = NA_real_))
  }
  
  eps <- 1e-12
  idx <- which(bh$time <= t_lim + eps)
  if (!length(idx)) return(c(CLSR = NA_real_, L0 = 0, h0 = NA_real_))
  
  times <- bh$time[idx]
  Hcum  <- bh$hazard[idx]
  dLam  <- diff(c(0, Hcum))
  
  ok <- is.finite(times) & is.finite(dLam) & (dLam > 0)
  times <- times[ok]
  dLam  <- dLam[ok]
  if (!length(dLam)) return(c(CLSR = NA_real_, L0 = 0, h0 = NA_real_))
  
  L0 <- sum(dLam)
  if (!is.finite(L0) || L0 <= 0) return(c(CLSR = NA_real_, L0 = L0, h0 = NA_real_))
  
  # crude last-interval hazard proxy near t_lim
  if (length(times) >= 2) {
    dt_last <- times[length(times)] - times[length(times) - 1]
    h0 <- if (is.finite(dt_last) && dt_last > 0) dLam[length(dLam)] / dt_last else NA_real_
  } else {
    h0 <- NA_real_
  }
  
  # if no D term, CLSR reduces to exp(gamma0)
  if (!is.finite(gamma_d) || abs(gamma_d) < 1e-12) {
    CLSR <- if (is.finite(gamma0)) exp(gamma0) else NA_real_
    return(c(CLSR = CLSR, L0 = L0, h0 = h0))
  }
  if (!is.finite(gamma0)) return(c(CLSR = NA_real_, L0 = L0, h0 = h0))
  
  # CLSR = exp(gamma0) * (sum dLam * exp(gamma_d * time)) / (sum dLam)
  v <- log(dLam) + gamma_d * times
  m <- max(v)
  log_num  <- m + log(sum(exp(v - m)))
  log_den  <- log(L0)
  log_clsr <- gamma0 + (log_num - log_den)
  
  CLSR <- if (log_clsr > 700) Inf else exp(log_clsr)
  c(CLSR = CLSR, L0 = L0, h0 = h0)
}

# ============================================================
# 3) Model specs (4 models)
#    X_i(t) = Z_i(t) = 1(t > T_i)
#    D_i(t) = max(0, t - T_i)
#    W_i    = baseline covariate
# ============================================================

model_specs <- list(
  Z = list(
    pretty = "Z(t) only: exp{ gamma0 * Z_i(t) }",
    add_D  = FALSE,
    add_W  = FALSE
  ),
  ZD = list(
    pretty = "Z(t) + D(t): exp{ gamma0 * Z_i(t) + gamma_d * D_i(t) }",
    add_D  = TRUE,
    add_W  = FALSE
  ),
  ZW = list(
    pretty = "Z(t) + W: exp{ gamma0 * Z_i(t) + delta_Y * W_i }",
    add_D  = FALSE,
    add_W  = TRUE
  ),
  ZDW = list(
    pretty = "Z(t) + D(t) + W: exp{ gamma0 * Z_i(t) + gamma_d * D_i(t) + delta_Y * W_i }",
    add_D  = TRUE,
    add_W  = TRUE
  )
)

# ============================================================
# 4) Fit a selected model
# ============================================================
fit_one_model <- function(dat, model_name) {
  if (!model_name %in% names(model_specs)) {
    stop("Unknown model_name. Use one of: ", paste(names(model_specs), collapse = ", "))
  }
  spec <- model_specs[[model_name]]
  
  # Build formula + tt()
  if (spec$add_D) {
    base_fml <- Surv(futime, status) ~ tt(T_treat) + tt(D_treat)
    tt_funs <- list(
      function(x, t, ...) as.integer(t > x),  # Z(t)
      function(x, t, ...) pmax(0, t - x)      # D(t)
    )
  } else {
    base_fml <- Surv(futime, status) ~ tt(T_treat)
    tt_funs <- function(x, t, ...) as.integer(t > x)  # Z(t)
  }
  
  if (spec$add_W) {
    fml <- update(base_fml, . ~ . + W)
  } else {
    fml <- base_fml
  }
  
  fit <- survival::coxph(
    fml,
    data    = dat,
    ties    = "breslow",
    control = survival::coxph.control(timefix = FALSE),
    tt      = tt_funs
  )
  
  cf <- coef(fit)
  gamma0_hat <- if (!is.null(cf) && "tt(T_treat)" %in% names(cf)) unname(cf["tt(T_treat)"]) else NA_real_
  gammaD_hat <- if (!is.null(cf) && "tt(D_treat)" %in% names(cf)) unname(cf["tt(D_treat)"]) else NA_real_
  deltaY_hat <- if (!is.null(cf) && "W" %in% names(cf))           unname(cf["W"])           else NA_real_
  
  bh <- bh_from_tt_safe(fit, dat)
  
  list(
    model_name = model_name,
    pretty     = spec$pretty,
    fit        = fit,
    gamma0_hat = gamma0_hat,
    gammaD_hat = gammaD_hat,
    deltaY_hat = deltaY_hat,
    bh         = bh
  )
}

# ============================================================
# 5) Fit many models + produce CLSR comparison table
# ============================================================
run_models_and_clsr <- function(dat, models = c("Z", "ZD"), horizons = NULL, print_fits = TRUE) {
  models <- unique(models)
  bad <- setdiff(models, names(model_specs))
  if (length(bad)) stop("Unknown model(s): ", paste(bad, collapse = ", "))
  
  # default horizons: event-time quantiles
  if (is.null(horizons)) {
    horizons <- as.numeric(stats::quantile(
      dat$futime[dat$status == 1],
      probs = c(0.75, 0.8, 0.85, 0.90, 0.95, 0.99),
      na.rm = TRUE,
      names = FALSE,
      type = 7
    ))
  } else {
    horizons <- as.numeric(horizons)
  }
  
  res_list <- lapply(models, function(m) fit_one_model(dat, m))
  names(res_list) <- models
  
  if (isTRUE(print_fits)) {
    for (m in models) {
      cat("\n============================================================\n")
      cat("Model:", m, "\n")
      cat(model_specs[[m]]$pretty, "\n")
      cat("------------------------------------------------------------\n")
      print(res_list[[m]]$fit)
      cat("gamma0_hat =", res_list[[m]]$gamma0_hat, "\n")
      cat("gammaD_hat =", res_list[[m]]$gammaD_hat, "\n")
      cat("deltaY_hat =", res_list[[m]]$deltaY_hat, "\n")
      cat("bh head:\n"); print(head(res_list[[m]]$bh, 6))
    }
  }
  
  # CLSR table for each model at each horizon
  clsr_tbl <- dplyr::bind_rows(lapply(models, function(m) {
    r <- res_list[[m]]
    dplyr::bind_rows(lapply(horizons, function(tlim) {
      mm <- clsr_L0_h0_at_t(r$bh, r$gamma0_hat, r$gammaD_hat, tlim)
      tibble(
        model = m,
        model_description = model_specs[[m]]$pretty,
        horizon_time = as.numeric(tlim),
        CLSR = unname(mm["CLSR"]),
        baseline_cumhaz_L0 = unname(mm["L0"]),
        baseline_hazard_h0 = unname(mm["h0"]),
        gamma0_hat = r$gamma0_hat,
        gammaD_hat = r$gammaD_hat,
        deltaY_hat = r$deltaY_hat
      )
    }))
  }))
  
  list(results = res_list, horizons = horizons, clsr_table = clsr_tbl)
}
