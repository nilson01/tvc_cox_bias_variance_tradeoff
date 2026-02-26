
# Sys.getenv("SLURM_SUBMIT_DIR", unset = getwd())

# On cluster: use the directory 
# setwd(Sys.getenv("SLURM_SUBMIT_DIR", unset = getwd()))

# setwd("/scratch/user/nchapagain/tvccsb")





#..............................................................................


#..............................................................................
#..............................................................................
## ---- Packages ----
#..............................................................................
#..............................................................................
# 


pkgs <- c(
  "data.table","survival","numDeriv","coxphf","knitr",
  "ggplot2", "pracma","nleqslv",
  "future","furrr","progressr",
  "dplyr","purrr","tibble","tidyr","parallel"
)

suppressPackageStartupMessages({
  ok <- vapply(pkgs, function(p)
    require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE),
    logical(1)
  )
})
if (any(!ok)) stop("Missing packages: ", paste(pkgs[!ok], collapse = ", "))

options(repr.plot.width = 28, repr.plot.height = 10)




# 00_setup.R

# ----TOP----



#..............................................................................
#..............................................................................
##----Global variables----
#..............................................................................
#..............................................................................


## Parameters


B  <- 1000




# ---- #

# shape_vals_beta_E <- c(1)
# shape_vals_beta_E <- c(1.5)
shape_vals_beta_E <- c(0.5) # 1, 1,5

# HR_vals    <- c(exp(1)) # gamma_0 = 1
# HR_vals    <- c(exp(-1)) # gamma_0 = -1
HR_vals    <- c(1) # gamma_0 = 0

# gam_ds_log <- c(1)
gam_ds_log <- c(0.5)
# gam_ds_log <- c(0)


a_l <- 0; a_h <- 1;

alpha0s  <- c(a_l)
a_x1s    <- c(a_l)
a_x2s    <- c(a_l)


# alpha0s  <- c(a_h)
# a_x1s    <- c(a_h)
# a_x2s    <- c(a_h)


# weibull
shape_vals_beta_T <- c(1)

n_vals     <- c(50, 100, 250, 500, 1500, 3000) # 100, 512
q_probs  <- c(0.75, 0.90, 0.99)
q_names  <- paste0("q_", as.integer(round(100 * q_probs)))


# lambda_E   <- 3
# lambda_T   <- 6


lambda_E   <- 3
lambda_T   <- 6

study_time <- 1
tol        <- 1e-5
max_it     <- 50
quantile   <- 0.75





# Model specs and grid


# model_specs <- list(
#   "Z + D + X1 + X2" = list(
#     fml     = Surv(futime, status) ~ X1 + X2 + tt(T_treat) + tt(D_treat),
#     tt_vars = c("Z0","D0")
#   ),
#   "Z + X1 + X2" = list(
#     fml     = Surv(futime, status) ~ X1 + X2 + tt(T_treat),
#     tt_vars = c("Z0")
#   )
# )





model_specs <- list(
  "Z + D" = list(
    fml     = Surv(futime, status) ~ tt(T_treat) + tt(D_treat),
    tt_vars = c("Z0","D0")
  ),
  "Z" = list(
    fml     = Surv(futime, status) ~ tt(T_treat),
    tt_vars = c("Z0")
  )
)





# -------------------------------
# Run ID builder (jobid + params + B)
# -------------------------------
make_run_id <- function(B,
                        beta_E, lambda_E, beta_T, lambda_T,
                        gamma_0, gamma_d,
                        alpha0, alpha_1, alpha_2,
                        n_vals) {
  
  # --- SLURM identity (ONLY job id) ---
  job_id <- Sys.getenv("SLURM_JOB_ID", unset = "")
  
  base_id <- if (nzchar(job_id)) {
    paste0("job_", job_id)
  } else {
    paste0("local_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  # --- compact formatting helpers ---
  fmt <- function(x) {
    x <- as.numeric(x)
    s <- format(x, scientific = FALSE, trim = TRUE)
    s <- sub("\\.0+$", "", s)
    s <- sub("(\\.[0-9]*?)0+$", "\\1", s)
    s
  }
  fmt_vec <- function(x) paste(fmt(x), collapse = "-")
  
  # --- alpha collapse rule: only write alpha_i_0 or alpha_i_1 ---
  alpha_i_tag <- if (all(c(alpha0, alpha_1, alpha_2) == 0)) "alpha_i_0" else "alpha_i_1"
  
  # --- build param tag string ---
  param_tag <- paste(
    paste0("B_",        as.integer(B)),
    paste0("beta_E_",   fmt(beta_E)),
    paste0("lambda_E_", fmt(lambda_E)),
    paste0("beta_T_",   fmt(beta_T)),
    paste0("lambda_T_", fmt(lambda_T)),
    paste0("gamma_0_",  fmt(gamma_0)),
    paste0("gamma_d_",  fmt(gamma_d)),
    alpha_i_tag,
    paste0("n_",        fmt_vec(n_vals)),
    sep = "_"
  )
  
  paste0(param_tag, "__", base_id)
  
}


PREPARE_DATA <- TRUE



RUN_ID <- make_run_id(
  B        = B,
  beta_E   = shape_vals_beta_E[1],
  lambda_E = lambda_E,
  beta_T   = shape_vals_beta_T[1],
  lambda_T = lambda_T,
  gamma_0  = log(HR_vals)[1], 
  gamma_d  = gam_ds_log[1], 
  alpha0   = alpha0s[1],
  alpha_1  = a_x1s[1],
  alpha_2  = a_x2s[1],
  n_vals   = n_vals
)








BASE_DIR <- "/scratch/user/nchapagain/tvccsb"

# One folder per run
RUN_DIR <- file.path(BASE_DIR, "tvccsb_runs", RUN_ID)

# Dataset folder inside the run folder ( switch depending on DGP type)
DATA_FOLDER_NAME <- "data_weibull"   # or "data_uniform" "data_weibull"


DATA_DIR <- file.path(RUN_DIR, DATA_FOLDER_NAME)

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

cat("RUN_ID           =", RUN_ID, "\n")
cat("RUN_DIR          =", RUN_DIR, "\n")
cat("DATA_FOLDER_NAME =", DATA_FOLDER_NAME, "\n")
cat("DATA_DIR         =", DATA_DIR, "\n")

# Fail fast if not writable
tf <- file.path(DATA_DIR, sprintf(".write_test_%d", Sys.getpid()))
tryCatch({
  writeLines("ok", tf)
  unlink(tf)
}, error = function(e) {
  stop("DATA_DIR is not writable: ", DATA_DIR, "\n", conditionMessage(e))
})




SEED_BASE    <- 1234
if (PREPARE_DATA) dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

## Threads
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}






param_grid <- expand.grid(
  beta_T = shape_vals_beta_T,
  beta_E = shape_vals_beta_E,
  HR     = HR_vals,
  gam_d  = gam_ds_log,
  n      = n_vals,
  alpha0 = alpha0s,
  a_x1   = a_x1s,
  a_x2   = a_x2s,
  KEEP.OUT.ATTRS = FALSE
)

cat("<----------------\n",
    "  Monte-Carlo replicates (B) = ", B, "\n",
    "  #Scenarios (grid rows)     = ", nrow(param_grid), "\n",
    "------------------>\n", sep = "")

## Progress handlers
options(progressr.enable = TRUE, progressr.clear = FALSE)
progressr::handlers(global = TRUE)
progressr::handlers(progressr::handler_txtprogressbar)
progress_every <- max(1L, floor(B / 20L))
progress_steps <- nrow(param_grid) * ceiling(B / progress_every)




# ----EOF----




# 01_helper functions

# ----TOP----


.rep_audit_counts <- function(dat) {
  n <- nrow(dat)
  n_events   <- sum(dat$status == 1)
  n_censored <- n - n_events
  event_rate <- if (n > 0) n_events / n else NA_real_
  list(
    n = n,
    n_events = n_events,
    n_censored = n_censored,
    event_rate = event_rate
  )
}



# Inversion helpers

#..............................................................................
#..............................................................................
# ----Helper for β_E = 1----
#..............................................................................
#..............................................................................

invert_branchB_beta1A <- function(y, T0, eta, lambda_E, gamma0, gam_d) {
  # constants
  c0    <- exp(eta)
  c1    <- exp(eta + gamma0 - gam_d * T0)
  gamma <- gam_d
  lam   <- lambda_E
  w0    <- lam * T0
  
  # analytic inversion: (4.1) with gamma != 0, n=1
  term <- exp(gamma * T0) +
    gamma * (y - c0 * w0) / (lam * c1)
  t    <- log(term) / gamma
  t
}







# Numeric inverse




# =========================
# Safety helpers
# =========================
.safe_exp <- function(x) {
  # exp(709) is near the largest finite double
  x <- pmin(pmax(x, -745), 700)
  exp(x)
}

.safe_finite <- function(x, fallback) {
  if (!is.finite(x)) fallback else x
}

.clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)


# =========================
# Branch B inversion (gamma_d != 0) via NR + fallbacks
# =========================
invert_branchB_NR <- function(y, T_treat, eta,
                              beta_E, lambda_E,
                              gamma0, gam_d,
                              study_time,
                              q_cap  = 0.999,
                              tol    = 1e-9,
                              max_it = 50) {
  
  if (gam_d == 0) stop("gam_d == 0: use closed-form branch.")
  n <- length(y)
  stopifnot(length(T_treat) == n, length(eta) == n)
  
  # Weibull-based time scale cap (baseline event-time quantile)
  cap_weib <- stats::qweibull(q_cap, shape = beta_E, scale = 1 / lambda_E)
  

  H_int <- function(t, T0) {
    if (t <= T0) return(0)
    
    lower <- max(T0, 1e-12)  # avoid log(0) & u^(beta-1) singular at 0
    upper <- t
    
    res <- tryCatch(
      stats::integrate(
        f = function(u) {
          u2 <- pmax(u, 1e-12)
          # log integrand = (beta_E-1)log(u) + gam_d*u
          log_val <- (beta_E - 1) * log(u2) + gam_d * u2
          .safe_exp(log_val)
        },
        lower        = lower,
        upper        = upper,
        rel.tol      = 1e-8,
        subdivisions = 100
      ),
      error = function(e) list(value = NA_real_)
    )
    res$value
  }
  
  out <- numeric(n)
  res_vals <- numeric(n)
  
  for (i in seq_len(n)) {
    
    T0   <- T_treat[i]
    yi   <- y[i]
    etai <- eta[i]
    
    # constants A and B 
    A  <- exp(etai) * (lambda_E * T0)^beta_E
    Bc <- exp(etai + gamma0 - gam_d * T0) * lambda_E^beta_E * beta_E
    
    # Region A (event before treatment)
    if (yi <= A) {
      out[i] <- ((yi / exp(etai))^(1 / beta_E)) / lambda_E
      res_vals[i] <- 0
      next
    }
    
    # -------------------------
    # Safe search interval (bounded)
    # -------------------------
    eps  <- 1e-10
    t_lo <- max(T0 + eps, 1e-12)
    

    t_hi <- T0 + max(
      study_time,
      10 * cap_weib,
      20 / max(abs(gam_d), 1e-3)
    )
    
    if (!is.finite(t_hi) || t_hi <= t_lo) {
      t_hi <- T0 + max(10 * cap_weib, 1)
    }
    
    # Hard cap for any bracket expansion (no 1e6 wandering)
    t_hard <- t_hi
    
    # -------------------------
    # f(t) and f'(t)
    # -------------------------
    f_fun <- function(t) {
      t <- .clamp(t, t_lo, t_hard)
      H <- H_int(t, T0)
      if (!is.finite(H)) return(1e100)   # keep solver moving
      val <- A + Bc * H - yi
      .safe_finite(val, 1e100)
    }
    
    # NOTE: derivative of A + Bc*∫(...) is Bc * t^(beta-1) * exp(gam_d t)
    # Compute in log-space to avoid overflow.
    df_fun <- function(t) {
      t <- .clamp(t, t_lo, t_hard)
      # log(df) = log(Bc) + (beta_E-1)log(t) + gam_d*t
      log_df <- log(Bc) + (beta_E - 1) * log(pmax(t, 1e-12)) + gam_d * t
      val <- .safe_exp(log_df)
      .safe_finite(val, 1e-12)
    }
    
    # -------------------------
    # Initial guess: Newton-style increment but bounded
    # -------------------------
    denom0  <- df_fun(t_lo)
    raw_inc <- (yi - A) / .safe_finite(denom0, 1e-8)
    
    max_inc <- 10 / max(abs(gam_d), 1e-3)
    inc <- pmin(pmax(raw_inc, -max_inc), max_inc)
    
    if (!is.finite(inc) || abs(inc) < 1e-6) inc <- sign(raw_inc) * (1 / lambda_E)
    
    t0 <- .clamp(T0 + inc, t_lo, t_hard)
    
    # =========================================================
    # Attempt 1: nleqslv with jac
    # Attempt 2: nleqslv Broyden
    # Attempt 3: uniroot bracket (monotone increasing f)
    # Final: return t_hard
    # =========================================================
    
    sol <- tryCatch(
      nleqslv::nleqslv(
        x       = t0,
        fn      = f_fun,
        jac     = df_fun,
        control = list(ftol = tol, maxit = max_it, allowSingular = TRUE)
      ),
      error = function(e) NULL
    )
    
    ok1 <- !is.null(sol) &&
      is.finite(sol$x) &&
      sol$x > T0 &&
      (sol$termcd %in% c(1, 2)) &&
      is.finite(f_fun(sol$x)) &&
      abs(f_fun(sol$x)) <= 1e-4
    
    if (!ok1) {
      sol2 <- tryCatch(
        nleqslv::nleqslv(
          x       = t0,
          fn      = f_fun,
          method  = "Broyden",
          control = list(ftol = tol, maxit = max_it)
        ),
        error = function(e) NULL
      )
      
      ok2 <- !is.null(sol2) &&
        is.finite(sol2$x) &&
        sol2$x > T0 &&
        (sol2$termcd %in% c(1, 2)) &&
        is.finite(f_fun(sol2$x)) &&
        abs(f_fun(sol2$x)) <= 1e-4
      
      if (ok2) sol <- sol2
    }
    
    ok_final <- !is.null(sol) &&
      is.finite(sol$x) &&
      sol$x > T0 &&
      is.finite(f_fun(sol$x)) &&
      abs(f_fun(sol$x)) <= 1e-4
    
    if (ok_final) {
      out[i] <- sol$x
      res_vals[i] <- f_fun(sol$x)
      next
    }
    
    # -------------------------
    # Attempt 3: uniroot bracket
    # -------------------------
    f_lo <- f_fun(t_lo)        # should be negative when yi > A
    f_hi <- f_fun(t_hard)
    
    # Expand upper bound ONLY up to t_hard (which is already bounded)
    if (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0) {
      t_try <- t_hard
      # Try a few increasing points between t_lo and t_hard
      # (monotone f => if no sign change by t_hard, no root in [t_lo,t_hard])
      for (k in 1:10) {
        t_mid <- t_lo + (k / 10) * (t_hard - t_lo)
        f_mid <- f_fun(t_mid)
        if (is.finite(f_mid) && f_lo * f_mid <= 0) {
          t_hi <- t_mid
          f_hi <- f_mid
          break
        }
      }
    } else {
      t_hi <- t_hard
    }
    
    root <- NULL
    if (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi <= 0) {
      root <- tryCatch(
        stats::uniroot(f_fun, lower = t_lo, upper = t_hi, tol = tol)$root,
        error = function(e) NULL
      )
    }
    
    if (!is.null(root) && is.finite(root) && root > T0) {
      out[i] <- root
      res_vals[i] <- f_fun(root)
    } else {
      # Final placeholder: bounded (not astronomically huge)
      out[i] <- t_hard
      res_vals[i] <- f_fun(t_hard)
    }
  }
  
  list(times = out, residuals = res_vals)
}




#..............................................................................
#..............................................................................
# ---- Breslow baseline from coxph with tt() ----
#..............................................................................
#..............................................................................


# breslow estimator with different model case handling

bh_from_tt <- function(fit, dat) {
  # Coefficients actually estimated in this fit
  cf <- coef(fit)
  if (is.null(cf) || !length(cf)) {
    return(data.frame(time = numeric(0), hazard = numeric(0)))
  }
  nm <- names(cf)
  
  # Event times
  t_ev <- sort(unique(dat$futime[dat$status == 1]))
  if (!length(t_ev)) return(data.frame(time = numeric(0), hazard = numeric(0)))
  
  # Helper to build the time-specific linear predictor for all subjects at time t
  lp_at_time <- function(t) {
    lp <- rep(0, nrow(dat))
    
    # Static covariates if present in the model
    if ("X1" %in% nm) lp <- lp + cf["X1"] * dat$X1
    if ("X2" %in% nm) lp <- lp + cf["X2"] * dat$X2
    
    # Time-dependent terms created via tt()
    # Name matching follows survival::coxph default: "tt(varname)"
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
    
    # Risk set just before t
    at_risk <- dat$futime >= t
    
    # Time-specific linear predictor and denominator
    lp_t  <- lp_at_time(t)
    denom <- sum(exp(lp_t[at_risk]))
    
    jump <- if (denom > 0) dk / denom else 0
    Hcum[k] <- if (k == 1) jump else Hcum[k - 1] + jump
  }
  
  data.frame(time = t_ev, hazard = Hcum)
}







#..............................................................................
#..............................................................................
# ----True baseline hazard pieces----
#..............................................................................
#..............................................................................

# TRUE baseline (Weibull) 
true_Lambda0 <- function(t, beta_E, lambda_E) (lambda_E * t)^beta_E # cumulative baseline hazard
true_h0      <- function(t, beta_E, lambda_E) (lambda_E^beta_E) * beta_E * t^(beta_E - 1) # Baseline hazard




#..............................................................................
#..............................................................................
# ---- CLSR/L0/h0 at a given horizon t_lim from basehaz jumps ----
#..............................................................................
#..............................................................................


.clsr_L0_h0_at_t <- function(bh, gamma0, gamma_d, t_lim) {
  if (is.null(bh) || !is.finite(t_lim)) return(c(CLSR = NA_real_, L0 = NA_real_, h0 = NA_real_))
  
  eps <- 1e-12
  idx <- which(bh$time <= t_lim + eps)
  if (length(idx) == 0) return(c(CLSR = NA_real_, L0 = 0, h0 = NA_real_))
  
  times <- bh$time[idx]
  Hcum  <- bh$hazard[idx]
  dLam  <- diff(c(0, Hcum))
  
  # guard against nonpositive / nonfinite jumps
  ok <- is.finite(times) & is.finite(dLam) & (dLam > 0)
  times <- times[ok]
  dLam  <- dLam[ok]
  
  L0 <- sum(dLam)
  if (!is.finite(L0) || L0 <= 0) {
    out <- c(CLSR = NA_real_, L0 = L0, h0 = NA_real_)
    names(out) <- c("CLSR","L0","h0")
    return(out)
  }
  
  # last-interval proxy for baseline hazard at t_lim
  dt_last <- if (length(times) >= 2) times[length(times)] - times[length(times)-1] else times[1]
  dlast   <- dLam[length(dLam)]
  h0      <- if (is.finite(dt_last) && dt_last > 0) dlast / dt_last else NA_real_
  
  # === Fast path when D is absent or gamma_d ~ 0 ===
  if (!is.finite(gamma_d) || abs(gamma_d) < 1e-12) {
    CLSR <- if (is.finite(gamma0)) exp(gamma0) else NA_real_
    out <- c(CLSR = CLSR, L0 = L0, h0 = h0)
    names(out) <- c("CLSR","L0","h0")
    return(out)
  }
  
  # === Numerically stable CLSR via log-sum-exp ===
  # CLSR = exp(gamma0) * (sum dLam * exp(gamma_d * times)) / L0
  if (!is.finite(gamma0)) {
    out <- c(CLSR = NA_real_, L0 = L0, h0 = h0)
    names(out) <- c("CLSR","L0","h0")
    return(out)
  }
  
  v <- log(dLam) + gamma_d * times   # log-weights
  m <- max(v)
  # numerator = exp(m) * sum(exp(v-m))
  log_num <- m + log(sum(exp(v - m)))
  log_den <- log(L0)
  log_clsr <- gamma0 + (log_num - log_den)
  
  # avoid exp overflow; R overflows around ~709
  CLSR <- if (log_clsr > 700) Inf else exp(log_clsr)
  
  out <- c(CLSR = CLSR, L0 = L0, h0 = h0)
  names(out) <- c("CLSR","L0","h0")
  out
}



# ----EOF----



# 03_cache_and_run.R


# ----TOP----


#..............................................................................
#..............................................................................
# File helpers
#..............................................................................
#..............................................................................

.scenario_file <- function(row_idx, B, data_dir = DATA_DIR) {
  file.path(data_dir, sprintf("scen_%04d__B_%d.rds", row_idx, B))
}



.rep_summaries_file <- function(row_idx, data_dir = DATA_DIR) {
  file.path(data_dir, sprintf("rep_summaries_scen_%04d.rds", row_idx))
}



.rep_done_dir <- function(row_idx, data_dir = DATA_DIR) {
  d <- file.path(data_dir, sprintf("rep_done_scen_%04d", row_idx))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

.rep_done_file <- function(row_idx, b, data_dir = DATA_DIR) {
  file.path(.rep_done_dir(row_idx, data_dir), sprintf("rep_%05d.done", b))
}

.count_done_reps <- function(row_idx, data_dir = DATA_DIR) {
  d <- .rep_done_dir(row_idx, data_dir)
  if (!dir.exists(d)) return(0L)
  length(list.files(d, pattern = "\\.done$", full.names = FALSE))
}



#..............................................................................
#..............................................................................
# ----Precompute all datasets for each scenario----
#..............................................................................
#..............................................................................

precompute_all_datasets <- function() {
  cat("Preparing datasets into: ", normalizePath(DATA_DIR, winslash = "/"), "\n", sep = "")
  pb_total <- nrow(param_grid) * B
  progressr::with_progress({
    p <- progressr::progressor(steps = pb_total)
    for (row_idx in seq_len(nrow(param_grid))) {
      row <- param_grid[row_idx, ]
      be  <- row$beta_E
      n   <- row$n
      alpha0 <- row$alpha0
      hrZ <- row$HR
      lgD <- row$gam_d
      gamma0 <- log(hrZ)

      # bound  <- stats::qweibull(p = quantile, shape = be, scale = 1 / lambda_E) # UNCOMMENT THIS FOR Uniform DGP
      bound  <- Inf 

      
      # Build list of B datasets for this scenario
      datasets <- vector("list", length = B)
      for (b in seq_len(B)) {
        seed <- SEED_BASE * row_idx + b
        datasets[[b]] <- generate_cohort(
          n = n,
          beta_T = row$beta_T, lambda_T = lambda_T,
          beta_E = be,         lambda_E = lambda_E,
          gamma0 = gamma0,     gam_d    = lgD,
          study_time = study_time,
          alpha0 = alpha0, a_x1 = row$a_x1, a_x2 = row$a_x2,
          tol = tol, treat_bound = bound, max_it = max_it,
          seed = seed
        )
        p(message = sprintf("prep scen %d/%d rep %d/%d", row_idx, nrow(param_grid), b, B), amount = 1)
      }
      
      # Save exactly the list plus a small header for audit
      meta <- list(
        param_row   = row,
        lambda_E    = lambda_E,
        lambda_T    = lambda_T,
        study_time  = study_time,
        alpha0      = alpha0,
        tol         = tol,
        max_it      = max_it,
        quantile    = quantile,
        seed_base   = SEED_BASE,
        created_at  = Sys.time()
      )
      saveRDS(list(meta = meta, datasets = datasets), .scenario_file(row_idx, B))
    }
  })
  cat("Precompute complete.\n")
}




#..............................................................................
#..............................................................................
# ----Scenario runner----
#..............................................................................
#..............................................................................


run_scenario <- function(row_idx) {
  row    <- param_grid[row_idx, ]
  bt     <- row$beta_T
  be     <- row$beta_E
  hr_Z   <- row$HR
  loghr_D<- row$gam_d
  n      <- row$n
  a_x1   <- row$a_x1
  a_x2   <- row$a_x2
  gamma0 <- log(hr_Z)
  
  # Load datasets file
  fpath <- .scenario_file(row_idx, B)
  
  if (!file.exists(fpath)) {
    stop("Cached data not found for scenario ", row_idx,
         ". Expected file: ", fpath,
         ". Set PREPARE_DATA <- TRUE to generate.")
  }
  payload <- readRDS(fpath)
  
  # CHECK: structure and length
  if (!is.list(payload) || is.null(payload$datasets))
    stop("Malformed cache file: missing $datasets in ", fpath)
  if (length(payload$datasets) != B)
    stop("Cache file has ", length(payload$datasets), " reps but B = ", B, " in memory.")
  
  # pooled quantiles per scenario (outside replication loop)
  pooled_ev_times <- sort(unlist(lapply(
    payload$datasets,
    function(d) d$futime[d$status == 1 & d$futime <= study_time]
  ), use.names = FALSE))
  
  q_times_all <- if (length(pooled_ev_times) > 0L) {
    as.numeric(stats::quantile(pooled_ev_times,
                               probs = q_probs,
                               type = 7, names = FALSE))
  } else rep(NA_real_, 3)
  
  
  # ---- Scenario-level audit metadata ----
  scenario_id <- row_idx
  scenario_created_at <- Sys.time()
  
  scenario_q_times_all <- q_times_all
  names(scenario_q_times_all) <- q_names
  
  
  
  M <- length(model_specs)
  
  scen_lab <- sprintf("Scenario %d/%d | n=%d HR(Z)=%.3f logHR(D)=%.3f",
                      row_idx, nrow(param_grid), n, hr_Z, loghr_D)
  cat("\n", scen_lab, " - Starting parallel replications...\n", sep = "")
  
  
  already <- .count_done_reps(row_idx)
  cat("Already completed reps on disk:", already, "out of", B, "\n")
  
  
  
  # =========================================================================
  # Define a function to process ONE replication
  # =========================================================================
  process_one_rep <- function(b) {
    dat <- payload$datasets[[b]]
    if (!is.data.frame(dat)) stop("Cached object for scen ", row_idx, " rep ", b, " is not a data.frame")
    
    
    # ---- Rep-level audit metadata ----
    rep_id <- b
    seed   <- SEED_BASE * row_idx + b   # matches precompute seed logic
    created_at <- Sys.time()
    
    audit <- .rep_audit_counts(dat)
    
    
    # Basic counts
    ev_count <- sum(dat$status)
    ev_pre_val   <- sum(dat$status == 1 & dat$futime <  dat$T_treat)
    ev_post_val  <- sum(dat$status == 1 & dat$futime >  dat$T_treat)
    ev_tie_val   <- sum(dat$status == 1 & dat$futime == dat$T_treat)
    prop_post_val <- if (ev_count > 0) ev_post_val / ev_count else NA_real_
    
    # D(t) among post-treatment failures
    D_fail <- dat$futime - dat$T_treat
    D_post_fail <- D_fail[dat$status == 1 & dat$futime > dat$T_treat]
    Dpost_n_val   <- length(D_post_fail)
    Dpost_sd_val  <- if (Dpost_n_val > 1) stats::sd(D_post_fail) else NA_real_
    Dpost_med_val <- if (Dpost_n_val > 0) stats::median(D_post_fail) else NA_real_
    Dpost_max_val <- if (Dpost_n_val > 0) max(D_post_fail) else NA_real_
    
    # horizon-specific post-treatment failures
    post_upto_q_vec <- vapply(q_times_all, function(tlim) {
      sum(dat$status == 1 & dat$futime <= tlim & dat$futime > dat$T_treat)
    }, integer(1))
    names(post_upto_q_vec) <- q_names
    
    # True CLSR computation
    q_times <- q_times_all
    disc_true_clsr <- function(tlim) {
      if (!is.finite(tlim)) return(NA_real_)
      ev_times_rep <- sort(dat$futime[dat$status == 1 & dat$futime <= study_time])
      tt <- ev_times_rep[ev_times_rep <= tlim]
      if (length(tt) == 0) return(NA_real_)
      L0_true_at <- true_Lambda0(tt, be, lambda_E)
      dLam_true  <- diff(c(0, L0_true_at))
      exp(gamma0) * (sum(dLam_true * exp(loghr_D * tt)) / sum(dLam_true))
    }
    clsr_true_vec <- vapply(q_times, disc_true_clsr, numeric(1))
    
    L0_true_vec <- ifelse(is.finite(q_times), true_Lambda0(q_times, be, lambda_E), NA_real_)
    h0_true_vec <- ifelse(is.finite(q_times), true_h0(q_times, be, lambda_E), NA_real_)
    
    # Prepare D_treat column
    dat$D_treat <- dat$T_treat
    
    # Storage for this replication's model results
    est_Z_vec <- rep(NA_real_, M)
    est_D_vec <- rep(NA_real_, M)
    bh_rmse_vec <- rep(NA_real_, M)
    bh_maxae_vec <- rep(NA_real_, M)
    bh_badjump_vec <- rep(NA_real_, M)
    
    est_CLSR_q <- matrix(NA_real_, nrow = 3, ncol = M)
    est_L0_q   <- matrix(NA_real_, nrow = 3, ncol = M)
    est_h0_q   <- matrix(NA_real_, nrow = 3, ncol = M)
    
    
    # ---- Add names so we can extract by model name later ----
    names(est_Z_vec) <- names(model_specs)
    names(est_D_vec) <- names(model_specs)
    names(bh_rmse_vec) <- names(model_specs)
    names(bh_maxae_vec) <- names(model_specs)
    names(bh_badjump_vec) <- names(model_specs)
    
    colnames(est_CLSR_q) <- names(model_specs)
    colnames(est_L0_q)   <- names(model_specs)
    colnames(est_h0_q)   <- names(model_specs)
    
    rownames(est_CLSR_q) <- c("q1","q2","q3")
    rownames(est_L0_q)   <- c("q1","q2","q3")
    rownames(est_h0_q)   <- c("q1","q2","q3")
    
    
    
    j <- 0L
    for (mdl_name in names(model_specs)) {
      j <- j + 1L
      spec <- model_specs[[mdl_name]]
      
      tt_funs <- switch(paste(spec$tt_vars, collapse = ","),
                        "Z0,D0" = list(
                          function(x, t, ...) as.integer(t > x),
                          function(x, t, ...) pmax(0, t - x)
                        ),
                        "Z0" = function(x, t, ...) as.integer(t > x),
                        "D0" = function(x, t, ...) pmax(0, t - x),
                        NULL
      )
      
      fit <- NULL
      vals <- tryCatch({
        fit <- survival::coxph(
          spec$fml,
          data    = dat,
          ties    = "breslow",
          control = survival::coxph.control(timefix = FALSE),
          tt      = tt_funs
        )
        cf <- coef(fit)
        vZ <- if (!is.null(cf) && "tt(T_treat)" %in% names(cf)) unname(cf["tt(T_treat)"]) else NA_real_
        vD <- if (!is.null(cf) && "tt(D_treat)" %in% names(cf)) unname(cf["tt(D_treat)"]) else NA_real_
        c(Z = vZ, D = vD)
      }, error = function(e) c(Z = NA_real_, D = NA_real_))
      
      bh <- if (is.null(fit)) NULL else tryCatch(bh_from_tt(fit, dat), error = function(e) NULL)
      
      # Baseline hazard diagnostics
      if (!is.null(bh) && nrow(bh) > 0) {
        t_eval <- bh$time
        H_est  <- bh$hazard
        H_true <- true_Lambda0(t_eval, be, lambda_E)
        diffs  <- H_est - H_true
        bh_rmse_vec[j]    <- sqrt(mean(diffs^2))
        bh_maxae_vec[j]   <- max(abs(diffs))
        dH <- diff(c(0, H_est))
        bh_badjump_vec[j] <- as.numeric(any(dH < -1e-10))
      }
      
      # CLSR estimates
      has_D <- !is.null(fit) && "tt(D_treat)" %in% names(coef(fit))
      gamma0_hat <- vals["Z"]
      gamma_d_hat <- if (has_D) vals["D"] else 0
      
      m1 <- .clsr_L0_h0_at_t(bh, gamma0_hat, gamma_d_hat, q_times[1])
      m2 <- .clsr_L0_h0_at_t(bh, gamma0_hat, gamma_d_hat, q_times[2])
      m3 <- .clsr_L0_h0_at_t(bh, gamma0_hat, gamma_d_hat, q_times[3])
      
      est_CLSR_q[1, j] <- m1["CLSR"]; est_L0_q[1, j] <- m1["L0"]; est_h0_q[1, j] <- m1["h0"]
      est_CLSR_q[2, j] <- m2["CLSR"]; est_L0_q[2, j] <- m2["L0"]; est_h0_q[2, j] <- m2["h0"]
      est_CLSR_q[3, j] <- m3["CLSR"]; est_L0_q[3, j] <- m3["L0"]; est_h0_q[3, j] <- m3["h0"]
      
      est_Z_vec[j] <- vals["Z"]
      est_D_vec[j] <- if (has_D) vals["D"] else NA_real_
    }
    
    
    # mark replication as completed (durable progress marker)
    done_path <- .rep_done_file(row_idx, b)
    writeLines(
      sprintf("done b=%d time=%s pid=%d", b, as.character(Sys.time()), Sys.getpid()),
      done_path
    )
    
    
    
    # Return all results for this replication as a list
    list(
      # ---- audit / metadata (NEW, additive) ----
      scenario_id = scenario_id,
      rep_id      = rep_id,
      seed        = seed,
      created_at  = created_at,
      q_times_all = scenario_q_times_all,
      
      n           = audit$n,
      n_events    = audit$n_events,
      n_censored  = audit$n_censored,
      event_rate  = audit$event_rate,
      
      b = b,
      ev_count = ev_count,
      ev_pre = ev_pre_val,
      ev_post = ev_post_val,
      ev_tie = ev_tie_val,
      prop_post = prop_post_val,
      Dpost_n = Dpost_n_val,
      Dpost_sd = Dpost_sd_val,
      Dpost_med = Dpost_med_val,
      Dpost_max = Dpost_max_val,
      post_upto_q = post_upto_q_vec,
      truth_CLSR = clsr_true_vec,
      truth_L0 = L0_true_vec,
      truth_h0 = h0_true_vec,
      est_Z = est_Z_vec,
      est_D = est_D_vec,
      bh_rmse = bh_rmse_vec,
      bh_maxae = bh_maxae_vec,
      bh_badjump = bh_badjump_vec,
      est_CLSR_q = est_CLSR_q,
      est_L0_q = est_L0_q,
      est_h0_q = est_h0_q
    )
  }
  
  # =========================================================================
  # RUN REPLICATIONS IN PARALLEL
  # =========================================================================
  # rep_results <- furrr::future_map(
  #   seq_len(B),
  #   process_one_rep,
  #   .options = furrr::furrr_options(
  #     seed = TRUE,
  #     packages = c("survival", "stats")
  #   ),
  #   .progress = TRUE
  # )
  
  rep_results <- furrr::future_map(
    seq_len(B),
    process_one_rep,
    .options = furrr::furrr_options(
      seed = TRUE,
      packages = c("survival", "stats", "dplyr", "tibble", "purrr"),   
      stdout = TRUE
    ),
    .progress = TRUE
  )
  
  
  saveRDS(
    list(
      scenario_id = scenario_id,
      param_row   = row,
      meta        = payload$meta,
      q_times_all = scenario_q_times_all,
      created_at  = scenario_created_at,
      rep_results = rep_results
    ),
    .rep_summaries_file(row_idx)
  )
  
  
  
  cat(scen_lab, " - Done.\n", sep = "")
  
  # =========================================================================
  # AGGREGATE RESULTS FROM ALL REPLICATIONS
  # =========================================================================
  ev_counts <- vapply(rep_results, `[[`, integer(1), "ev_count")
  ev_pre    <- vapply(rep_results, `[[`, integer(1), "ev_pre")
  ev_post   <- vapply(rep_results, `[[`, integer(1), "ev_post")
  ev_tie    <- vapply(rep_results, `[[`, integer(1), "ev_tie")
  prop_post <- vapply(rep_results, `[[`, numeric(1), "prop_post")
  Dpost_n   <- vapply(rep_results, `[[`, integer(1), "Dpost_n")
  Dpost_sd  <- vapply(rep_results, `[[`, numeric(1), "Dpost_sd")
  Dpost_med <- vapply(rep_results, `[[`, numeric(1), "Dpost_med")
  Dpost_max <- vapply(rep_results, `[[`, numeric(1), "Dpost_max")
  
  # Post-treatment events up to q horizons: B x K matrix
  post_upto_q_mat <- do.call(rbind, lapply(rep_results, `[[`, "post_upto_q"))
  
  # Truth vectors (length B each)
  truth_CLSR_q_c1 <- vapply(rep_results, function(r) r$truth_CLSR[1], numeric(1))
  truth_CLSR_q_c2 <- vapply(rep_results, function(r) r$truth_CLSR[2], numeric(1))
  truth_CLSR_q_c3 <- vapply(rep_results, function(r) r$truth_CLSR[3], numeric(1))
  truth_L0_q_c1   <- vapply(rep_results, function(r) r$truth_L0[1], numeric(1))
  truth_L0_q_c2   <- vapply(rep_results, function(r) r$truth_L0[2], numeric(1))
  truth_L0_q_c3   <- vapply(rep_results, function(r) r$truth_L0[3], numeric(1))
  truth_h0_q_c1   <- vapply(rep_results, function(r) r$truth_h0[1], numeric(1))
  truth_h0_q_c2   <- vapply(rep_results, function(r) r$truth_h0[2], numeric(1))
  truth_h0_q_c3   <- vapply(rep_results, function(r) r$truth_h0[3], numeric(1))
  
  # Model estimates: B x M matrices
  est_Z <- do.call(rbind, lapply(rep_results, `[[`, "est_Z"))
  est_D <- do.call(rbind, lapply(rep_results, `[[`, "est_D"))
  bh_rmse_mat    <- do.call(rbind, lapply(rep_results, `[[`, "bh_rmse"))
  bh_maxae_mat   <- do.call(rbind, lapply(rep_results, `[[`, "bh_maxae"))
  bh_badjump_mat <- do.call(rbind, lapply(rep_results, `[[`, "bh_badjump"))
  
  colnames(est_Z) <- names(model_specs)
  colnames(est_D) <- names(model_specs)
  colnames(bh_rmse_mat)    <- names(model_specs)
  colnames(bh_maxae_mat)   <- names(model_specs)
  colnames(bh_badjump_mat) <- names(model_specs)
  
  # CLSR/L0/h0 estimates: need to extract per model
  # est_CLSR_q is 3 x M for each rep
  # CLSR/L0/h0 estimates: need to extract per model
  # est_CLSR_q is 3 x M for each rep - NOTE THE COMMA AFTER THE INDEX!
  est_CLSR_q_c1 <- do.call(rbind, lapply(rep_results, function(r) r$est_CLSR_q[1, , drop = FALSE]))
  est_CLSR_q_c2 <- do.call(rbind, lapply(rep_results, function(r) r$est_CLSR_q[2, , drop = FALSE]))
  est_CLSR_q_c3 <- do.call(rbind, lapply(rep_results, function(r) r$est_CLSR_q[3, , drop = FALSE]))
  est_L0_q_c1   <- do.call(rbind, lapply(rep_results, function(r) r$est_L0_q[1, , drop = FALSE]))
  est_L0_q_c2   <- do.call(rbind, lapply(rep_results, function(r) r$est_L0_q[2, , drop = FALSE]))
  est_L0_q_c3   <- do.call(rbind, lapply(rep_results, function(r) r$est_L0_q[3, , drop = FALSE]))
  est_h0_q_c1   <- do.call(rbind, lapply(rep_results, function(r) r$est_h0_q[1, , drop = FALSE]))
  est_h0_q_c2   <- do.call(rbind, lapply(rep_results, function(r) r$est_h0_q[2, , drop = FALSE]))
  est_h0_q_c3   <- do.call(rbind, lapply(rep_results, function(r) r$est_h0_q[3, , drop = FALSE]))
  
  colnames(est_CLSR_q_c1) <- names(model_specs)
  colnames(est_CLSR_q_c2) <- names(model_specs)
  colnames(est_CLSR_q_c3) <- names(model_specs)
  colnames(est_L0_q_c1)   <- names(model_specs)
  colnames(est_L0_q_c2)   <- names(model_specs)
  colnames(est_L0_q_c3)   <- names(model_specs)
  colnames(est_h0_q_c1)   <- names(model_specs)
  colnames(est_h0_q_c2)   <- names(model_specs)
  colnames(est_h0_q_c3)   <- names(model_specs)
  
  # =========================================================================
  # BUILD OUTPUT TIBBLE (same structure as before)
  # =========================================================================
  out <- lapply(seq_len(M), function(j) {
    tibble::tibble(
      n = n, beta_T = bt, beta_E = be, a_x1 = a_x1, a_x2 = a_x2,
      gam_d = loghr_D, HR = hr_Z,
      model = names(model_specs)[j],
      sim_Z = list(est_Z[, j]),
      sim_D = list(est_D[, j]),
      events = list(ev_counts),
      events_pre  = list(ev_pre),
      events_post = list(ev_post),
      events_tie  = list(ev_tie),
      prop_events_post = list(prop_post),
      D_post_n   = list(Dpost_n),
      D_post_sd  = list(Dpost_sd),
      D_post_med = list(Dpost_med),
      D_post_max = list(Dpost_max),
      post_events_upto_q = list(post_upto_q_mat),
      sim_CLSR_q_c1 = list(est_CLSR_q_c1[, j]),
      sim_CLSR_q_c2 = list(est_CLSR_q_c2[, j]),
      sim_CLSR_q_c3 = list(est_CLSR_q_c3[, j]),
      sim_L0_q_c1   = list(est_L0_q_c1[, j]),
      sim_L0_q_c2   = list(est_L0_q_c2[, j]),
      sim_L0_q_c3   = list(est_L0_q_c3[, j]),
      sim_h0_q_c1   = list(est_h0_q_c1[, j]),
      sim_h0_q_c2   = list(est_h0_q_c2[, j]),
      sim_h0_q_c3   = list(est_h0_q_c3[, j]),
      truth_CLSR_q_c1 = list(truth_CLSR_q_c1),
      truth_CLSR_q_c2 = list(truth_CLSR_q_c2),
      truth_CLSR_q_c3 = list(truth_CLSR_q_c3),
      truth_L0_q_c1   = list(truth_L0_q_c1),
      truth_L0_q_c2   = list(truth_L0_q_c2),
      truth_L0_q_c3   = list(truth_L0_q_c3),
      truth_h0_q_c1   = list(truth_h0_q_c1),
      truth_h0_q_c2   = list(truth_h0_q_c2),
      truth_h0_q_c3   = list(truth_h0_q_c3),
      sim_bh_rmse    = list(bh_rmse_mat[, j]),
      sim_bh_maxae   = list(bh_maxae_mat[, j]),
      sim_bh_badjmp  = list(bh_badjump_mat[, j])
    )
  })
  
  dplyr::bind_rows(out)
}

# ----EOF----


# 02_data.R


# ----TOP----


#..............................................................................
#..............................................................................
# ----Diagnostics container----
#..............................................................................
#..............................................................................

diagnostics <- data.frame(
  y        = numeric(),
  T_treat  = numeric(),
  eta      = numeric(),
  beta_E   = numeric(),
  lambda_E = numeric(),
  gamma0   = numeric(),
  gam_d    = numeric(),
  stringsAsFactors = FALSE
)






#..............................................................................
#..............................................................................
# ----Data generating process----
#..............................................................................
#..............................................................................

generate_cohort <- function(n,
                            beta_T ,
                            lambda_T ,
                            beta_E   ,
                            lambda_E ,
                            gamma0 ,
                            gam_d,
                            study_time,
                            alpha0 ,
                            a_x1 ,
                            a_x2  ,
                            tol  ,
                            treat_bound,
                            max_it,
                            
                            seed = NULL) {
  
  
  
  
  
  if (!is.null(seed)) set.seed(seed)
  
  ## 1. baseline covariates -------------------------------------------
  
  
  
  p_bin <- c(0.30)
  X <- cbind(
    sapply(p_bin, \(p) rbinom(n, 1, p)),
    matrix(rnorm(n), ncol = 1)
  )
  
  colnames(X) <- paste0("X", 1:2)
  
  
  ## 2. linear predictors ---------------------------------------------
  
  W_event <- cbind(1, X[, c(1,2)])
  
  # delta   <- c(alpha0, a_l, a_m, a_h, a_l, a_m, a_h)
  delta   <- c(alpha0, a_x1, a_x2)
  
  # eta_treat <- as.vector(W_treat %*% delta)
  eta_event <- as.vector(W_event %*% delta)
  
  ## 3. treatment time -------------------------------------------------
  uT       <- runif(n)
  
  
  ## treatment time: Uniform(0, treat_bound) -----------
  # T_treat <- runif(n, min = 0, max = treat_bound) # uncomment bound in the function precompute_all_datasets(.)
  
  
  # weibull treatment times
  T_treat   <- ( ( -log(1 - uT) )^(1 / beta_T) ) / lambda_T
  
  
  # purpose to avoid numerical problems for too small T_treat
  T_treat = pmax(T_treat, 1e-14)
  
  
  
  ## 4. event time -----------------------------------------------------
  uE   <- runif(n)
  y    <- -log(1 - uE)
  yThr <- exp(eta_event) * (lambda_E * T_treat)^beta_E
  
  T_event <- numeric(n)
  residuals <- numeric(n)
  
  idx_A <- y < yThr                                   # fail before Treatment
  
  
  
  
  ## ----- Region A: event before treatment ------------------------
  
  if (any(idx_A)) {
    T_event[idx_A] <- ( y[idx_A] / exp(eta_event[idx_A]) )^(1 / beta_E) / lambda_E
  }
  
  
  
  ## ----- 4B. Region B - event after treatment, gamma_d == 0 ----------------------------------
  idx_B0 <- (!idx_A) & (gam_d == 0)
  if (any(idx_B0)) {
    T_event[idx_B0] <-
      ((lambda_E * T_treat[idx_B0])^beta_E +
         (y[idx_B0] - yThr[idx_B0]) / exp(eta_event[idx_B0] + gamma0))^(1 / beta_E) /
      lambda_E
  }
  
  
  
  ## ----- 4C. Region B, γ_d ≠ 0 ----------------------------------
  idx_Bg <- (!idx_A) & (gam_d != 0)
  if (any(idx_Bg)) {
    # analytic vs numeric comparison for β_E = 1 or 2
    # if (beta_E %in% c(1, 2)) {
    if (beta_E %in% c(1)) {
      
      # 1) compute analytic times
      analytic_times <- invert_branchB_beta1A(
        y[idx_Bg],
        T_treat[idx_Bg],
        eta_event[idx_Bg],
        lambda_E,
        gamma0,
        gam_d
      )
      
      # 4) store the numeric solution 
      T_event[idx_Bg]   <- analytic_times

      
    } else {
      # fallback: pure numeric inversion
      res <- invert_branchB_NR(
        y         = y[idx_Bg],
        T_treat   = T_treat[idx_Bg],
        eta       = eta_event[idx_Bg],
        beta_E    = beta_E,
        lambda_E  = lambda_E,
        gamma0    = gamma0,
        gam_d     = gam_d,
        study_time= study_time,
        q_cap     = 0.999,
        tol       = tol,
        max_it    = max_it
      )
      T_event[idx_Bg]   <- res$times
      residuals[idx_Bg] <- res$residuals
    }
  }
  
  
  
  
  ## 5) Administrative censoring  ----------------------------
  status  <- as.integer(T_event <= study_time)
  T_obs   <- pmin(T_event, study_time)
  
  ## Ensure X1, X2 exist and are named
  X2cols <- as.data.frame(X[, 1:2, drop = FALSE])
  colnames(X2cols) <- c("X1","X2")
  
  ## 6) Return ONE ROW per subject (no splitting)
  base <- cbind(
    data.frame(
      id      = 1:n,
      futime  = T_obs,
      status  = status,
      T_treat = T_treat
    ),
    X2cols
  )
  
  base   # <- return 
  
  
}







# ----EOF----






# 04_main_and_summary.R

# ----TOP----



set_parallel_plan <- function() {
  workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  if (!is.finite(workers) || workers < 1) workers <- 1L
  
  if (.Platform$OS.type == "unix" && future::supportsMulticore()) {
    future::plan(future::multicore, workers = workers)
    plan_name <- "multicore"
  } else {
    future::plan(future::multisession, workers = workers)
    plan_name <- "multisession"
  }
  
  cat("Future plan:", plan_name, "| workers:", workers, "\n")
}



#..............................................................................
#..............................................................................
# Driver: build cache, run scenarios, gather results
#..............................................................................
#..............................................................................

run_pipeline <- function() {
  if (PREPARE_DATA) {
    total_start <- Sys.time()
    precompute_all_datasets()
    total_end <- Sys.time()
    cat(sprintf("Cache build time: %0.2f minutes\n",
                as.numeric(difftime(total_end, total_start, units = "mins"))))
  }
  
  fit_start <- Sys.time()
  
  # Loop over scenarios sequentially; parallelism happens INSIDE run_scenario
  results_list <- vector("list", nrow(param_grid))
  
  for (row_idx in seq_len(nrow(param_grid))) {
    results_list[[row_idx]] <- run_scenario(row_idx)
  }
  
  results <- dplyr::bind_rows(results_list)
  
  fit_end <- Sys.time()
  cat(sprintf("Fit time: %0.2f minutes\n",
              as.numeric(difftime(fit_end, fit_start, units = "mins"))))
  
  return(results)
}



#..............................................................................
#..............................................................................
# Build summaries and save
#..............................................................................
#..............................................................................


summarize_and_save <- function(results) {
  if (!is.null(results)) {
    summ_par <- results %>%
      dplyr::transmute(
        n, beta_T, beta_E, a_x1, a_x2, gam_d, HR, model,
        truth_logHR_Z = log(HR),      truth_HR_Z = HR,
        truth_logHR_D = gam_d,        truth_HR_D = exp(gam_d),
        
        # Z (log-HR)
        mean_logHR_Z = purrr::map_dbl(sim_Z, ~ mean(.x, na.rm = TRUE)),
        sd_logHR_Z   = purrr::map_dbl(sim_Z, ~ stats::sd(.x,  na.rm = TRUE)),
        bias_logHR_Z = purrr::map2_dbl(sim_Z, truth_logHR_Z, ~ mean(.x - .y, na.rm = TRUE)),
        mse_logHR_Z  = purrr::map2_dbl(sim_Z, truth_logHR_Z, ~ mean((.x - .y)^2, na.rm = TRUE)),
        
        # Z (HR)
        mean_HR_Z = purrr::map_dbl(sim_Z, ~ mean(exp(.x), na.rm = TRUE)),
        sd_HR_Z   = purrr::map_dbl(sim_Z, ~ stats::sd(exp(.x),  na.rm = TRUE)),
        bias_HR_Z = purrr::map2_dbl(sim_Z, truth_HR_Z, ~ mean(exp(.x) - .y, na.rm = TRUE)),
        mse_HR_Z  = purrr::map2_dbl(sim_Z, truth_HR_Z, ~ mean((exp(.x) - .y)^2, na.rm = TRUE)),
        na_frac_Z = purrr::map_dbl(sim_Z, ~ mean(is.na(.x))),
        
        # D (log-HR)
        mean_logHR_D = purrr::map_dbl(sim_D, ~ mean(.x, na.rm = TRUE)),
        sd_logHR_D   = purrr::map_dbl(sim_D, ~ stats::sd(.x,  na.rm = TRUE)),
        bias_logHR_D = purrr::map2_dbl(sim_D, truth_logHR_D, ~ mean(.x - .y, na.rm = TRUE)),
        mse_logHR_D  = purrr::map2_dbl(sim_D, truth_logHR_D, ~ mean((.x - .y)^2, na.rm = TRUE)),
        
        # D (HR)
        mean_HR_D = purrr::map_dbl(sim_D, ~ mean(exp(.x), na.rm = TRUE)),
        sd_HR_D   = purrr::map_dbl(sim_D, ~ stats::sd(exp(.x),  na.rm = TRUE)),
        bias_HR_D = purrr::map2_dbl(sim_D, truth_HR_D, ~ mean(exp(.x) - .y, na.rm = TRUE)),
        mse_HR_D  = purrr::map2_dbl(sim_D, truth_HR_D, ~ mean((exp(.x) - .y)^2, na.rm = TRUE)),
        na_frac_D = purrr::map_dbl(sim_D, ~ mean(is.na(.x))),
        
        # Events
        mean_events     = purrr::map_dbl(events, ~ mean(.x, na.rm = TRUE)),
        sd_events       = purrr::map_dbl(events, ~ stats::sd(.x, na.rm = TRUE)),
        min_events      = purrr::map_dbl(events, ~ min(.x, na.rm = TRUE)),
        max_events      = purrr::map_dbl(events, ~ max(.x, na.rm = TRUE)),
        mean_event_rate = mean_events / n,
        sd_event_rate   = sd_events   / n,
        
        
        
        
        
        # ============================================================
        # Post-treatment event timing diagnostics
        # ============================================================
        
        # --- Counts of events BEFORE treatment switch (per replication, then summarized) ---
        mean_events_pre  = purrr::map_dbl(events_pre,  ~ mean(.x, na.rm = TRUE)),       # avg # pre-treatment events across reps
        sd_events_pre    = purrr::map_dbl(events_pre,  ~ stats::sd(.x, na.rm = TRUE)),  # SD of pre-treatment event counts across reps
        min_events_pre   = purrr::map_dbl(events_pre,  ~ min(.x, na.rm = TRUE)),        # min pre-treatment events over reps
        max_events_pre   = purrr::map_dbl(events_pre,  ~ max(.x, na.rm = TRUE)),        # max pre-treatment events over reps
        
        # --- Counts of events AFTER treatment switch (per replication, then summarized) ---
        mean_events_post = purrr::map_dbl(events_post, ~ mean(.x, na.rm = TRUE)),       # avg # post-treatment events across reps
        sd_events_post   = purrr::map_dbl(events_post, ~ stats::sd(.x, na.rm = TRUE)),  # SD of post-treatment event counts across reps
        min_events_post  = purrr::map_dbl(events_post, ~ min(.x, na.rm = TRUE)),        # min post-treatment events over reps
        max_events_post  = purrr::map_dbl(events_post, ~ max(.x, na.rm = TRUE)),        # max post-treatment events over reps
        
        # --- Events exactly at treatment time (ties) ---
        # Interpretation: futime == T_treat. Ideally near 0 unless our generator creates ties.
        mean_events_tie  = purrr::map_dbl(events_tie,  ~ mean(.x, na.rm = TRUE)),       # avg # tie events across reps
        sd_events_tie    = purrr::map_dbl(events_tie,  ~ stats::sd(.x, na.rm = TRUE)),  # SD of tie counts across reps
        min_events_tie   = purrr::map_dbl(events_tie,  ~ min(.x, na.rm = TRUE)),        # min ties over reps
        max_events_tie   = purrr::map_dbl(events_tie,  ~ max(.x, na.rm = TRUE)),        # max ties over reps
        
        # --- Convert event counts into per-subject rates (divide by n) ---
        # Note: these are NOT hazard rates. Just “events per subject” for sanity checks.
        mean_event_pre_rate  = mean_events_pre  / n,  # avg (pre events)/n
        mean_event_post_rate = mean_events_post / n,  # avg (post events)/n
        
        # --- Proportion of events that happen post-treatment ---
        # Per replication: prop_events_post[b] = n_events_post[b] / n_events_total[b]
        # Here we summarize that proportion across reps.
        mean_prop_events_post = purrr::map_dbl(prop_events_post, ~ mean(.x, na.rm = TRUE)),                         # avg post-event proportion
        sd_prop_events_post   = purrr::map_dbl(prop_events_post, ~ stats::sd(.x, na.rm = TRUE)),                    # SD of post-event proportion
        q50_prop_events_post  = purrr::map_dbl(prop_events_post, ~ stats::quantile(.x, 0.50, na.rm = TRUE, type=7)),# median post-event proportion
        q90_prop_events_post  = purrr::map_dbl(prop_events_post, ~ stats::quantile(.x, 0.90, na.rm = TRUE, type=7)),# 90th %ile post-event proportion
        
        
        # ============================================================
        # D(t)=futime - T_treat among POST-treatment failures only
        # ============================================================
        
        mean_D_post_n   = purrr::map_dbl(D_post_n,   ~ mean(.x, na.rm = TRUE)),  # avg # post-treatment failures across reps
        mean_D_post_sd  = purrr::map_dbl(D_post_sd,  ~ mean(.x, na.rm = TRUE)),  # avg SD of (futime - T_treat) across reps
        mean_D_post_med = purrr::map_dbl(D_post_med, ~ mean(.x, na.rm = TRUE)),  # avg median delay after treatment across reps
        mean_D_post_max = purrr::map_dbl(D_post_max, ~ mean(.x, na.rm = TRUE)),  # avg max delay after treatment across reps
        
        # “distribution over replications” summaries (useful to see instability)
        q50_D_post_sd   = purrr::map_dbl(D_post_sd,  ~ stats::quantile(.x, 0.50, na.rm = TRUE, type = 7)), # median of per-rep SDs
        q90_D_post_sd   = purrr::map_dbl(D_post_sd,  ~ stats::quantile(.x, 0.90, na.rm = TRUE, type = 7)), # 90th %ile of per-rep SDs
        q50_D_post_max  = purrr::map_dbl(D_post_max, ~ stats::quantile(.x, 0.50, na.rm = TRUE, type = 7)), # median of per-rep max delays
        q90_D_post_max  = purrr::map_dbl(D_post_max, ~ stats::quantile(.x, 0.90, na.rm = TRUE, type = 7)), # 90th %ile of per-rep max delays
        
        
        # ============================================================
        # Horizon-specific post-treatment event counts (B x K matrix)
        # ============================================================
        
        mean_post_events_upto_q = purrr::map(post_events_upto_q, ~ colMeans(.x, na.rm = TRUE)), # named vector: mean count per horizon
        sd_post_events_upto_q   = purrr::map(post_events_upto_q, ~ apply(.x, 2, stats::sd, na.rm = TRUE)), # named vector: SD per horizon
        q50_post_events_upto_q  = purrr::map(post_events_upto_q, ~ apply(.x, 2, stats::quantile, probs = 0.50, na.rm = TRUE, type = 7)), # median per horizon
        q90_post_events_upto_q  = purrr::map(post_events_upto_q, ~ apply(.x, 2, stats::quantile, probs = 0.90, na.rm = TRUE, type = 7)), # 90th %ile per horizon
        
        
        
        
        # Q1
        mean_CLSR_q_c1 = purrr::map_dbl(sim_CLSR_q_c1, ~ mean(.x, na.rm = TRUE)),
        sd_CLSR_q_c1   = purrr::map_dbl(sim_CLSR_q_c1, ~ stats::sd(.x, na.rm = TRUE)),
        bias_CLSR_q_c1 = purrr::map2_dbl(sim_CLSR_q_c1, truth_CLSR_q_c1, ~ mean(.x - .y, na.rm = TRUE)),
        mse_CLSR_q_c1  = purrr::map2_dbl(sim_CLSR_q_c1, truth_CLSR_q_c1, ~ mean((.x - .y)^2, na.rm = TRUE)),
        
        mean_L0_q_c1 = purrr::map_dbl(sim_L0_q_c1, ~ mean(.x, na.rm = TRUE)),
        sd_L0_q_c1   = purrr::map_dbl(sim_L0_q_c1, ~ stats::sd(.x, na.rm = TRUE)),
        bias_L0_q_c1 = purrr::map2_dbl(sim_L0_q_c1, truth_L0_q_c1, ~ mean(.x - .y, na.rm = TRUE)),
        mse_L0_q_c1  = purrr::map2_dbl(sim_L0_q_c1, truth_L0_q_c1, ~ mean((.x - .y)^2, na.rm = TRUE)),
        
        mean_h0_q_c1 = purrr::map_dbl(sim_h0_q_c1, ~ mean(.x, na.rm = TRUE)),
        sd_h0_q_c1   = purrr::map_dbl(sim_h0_q_c1, ~ stats::sd(.x, na.rm = TRUE)),
        bias_h0_q_c1 = purrr::map2_dbl(sim_h0_q_c1, truth_h0_q_c1, ~ mean(.x - .y, na.rm = TRUE)),
        mse_h0_q_c1  = purrr::map2_dbl(sim_h0_q_c1, truth_h0_q_c1, ~ mean((.x - .y)^2, na.rm = TRUE)),
        
        # Q2
        mean_CLSR_q_c2 = purrr::map_dbl(sim_CLSR_q_c2, ~ mean(.x, na.rm = TRUE)),
        sd_CLSR_q_c2   = purrr::map_dbl(sim_CLSR_q_c2, ~ stats::sd(.x, na.rm = TRUE)),
        bias_CLSR_q_c2 = purrr::map2_dbl(sim_CLSR_q_c2, truth_CLSR_q_c2, ~ mean(.x - .y, na.rm = TRUE)),
        mse_CLSR_q_c2  = purrr::map2_dbl(sim_CLSR_q_c2, truth_CLSR_q_c2, ~ mean((.x - .y)^2, na.rm = TRUE)),
        
        mean_L0_q_c2 = purrr::map_dbl(sim_L0_q_c2, ~ mean(.x, na.rm = TRUE)),
        sd_L0_q_c2   = purrr::map_dbl(sim_L0_q_c2, ~ stats::sd(.x, na.rm = TRUE)),
        bias_L0_q_c2 = purrr::map2_dbl(sim_L0_q_c2, truth_L0_q_c2, ~ mean(.x - .y, na.rm = TRUE)),
        mse_L0_q_c2  = purrr::map2_dbl(sim_L0_q_c2, truth_L0_q_c2, ~ mean((.x - .y)^2, na.rm = TRUE)),
        
        mean_h0_q_c2 = purrr::map_dbl(sim_h0_q_c2, ~ mean(.x, na.rm = TRUE)),
        sd_h0_q_c2   = purrr::map_dbl(sim_h0_q_c2, ~ stats::sd(.x, na.rm = TRUE)),
        bias_h0_q_c2 = purrr::map2_dbl(sim_h0_q_c2, truth_h0_q_c2, ~ mean(.x - .y, na.rm = TRUE)),
        mse_h0_q_c2  = purrr::map2_dbl(sim_h0_q_c2, truth_h0_q_c2, ~ mean((.x - .y)^2, na.rm = TRUE)),
        
        # Q3
        mean_CLSR_q_c3 = purrr::map_dbl(sim_CLSR_q_c3, ~ mean(.x, na.rm = TRUE)),
        sd_CLSR_q_c3   = purrr::map_dbl(sim_CLSR_q_c3, ~ stats::sd(.x, na.rm = TRUE)),
        bias_CLSR_q_c3 = purrr::map2_dbl(sim_CLSR_q_c3, truth_CLSR_q_c3, ~ mean(.x - .y, na.rm = TRUE)),
        mse_CLSR_q_c3  = purrr::map2_dbl(sim_CLSR_q_c3, truth_CLSR_q_c3, ~ mean((.x - .y)^2, na.rm = TRUE)),
        
        
        
        mean_L0_q_c3 = purrr::map_dbl(sim_L0_q_c3, ~ mean(.x, na.rm = TRUE)),
        sd_L0_q_c3   = purrr::map_dbl(sim_L0_q_c3, ~ stats::sd(.x, na.rm = TRUE)),
        bias_L0_q_c3 = purrr::map2_dbl(sim_L0_q_c3, truth_L0_q_c3, ~ mean(.x - .y, na.rm = TRUE)),
        mse_L0_q_c3  = purrr::map2_dbl(sim_L0_q_c3, truth_L0_q_c3, ~ mean((.x - .y)^2, na.rm = TRUE)),
        
        mean_h0_q_c3 = purrr::map_dbl(sim_h0_q_c3, ~ mean(.x, na.rm = TRUE)),
        sd_h0_q_c3   = purrr::map_dbl(sim_h0_q_c3, ~ stats::sd(.x, na.rm = TRUE)),
        bias_h0_q_c3 = purrr::map2_dbl(sim_h0_q_c3, truth_h0_q_c3, ~ mean(.x - .y, na.rm = TRUE)),
        mse_h0_q_c3  = purrr::map2_dbl(sim_h0_q_c3, truth_h0_q_c3, ~ mean((.x - .y)^2, na.rm = TRUE)),
        
        # basehaz diagnostics
        mean_bh_rmse   = purrr::map_dbl(sim_bh_rmse,  ~ mean(.x, na.rm = TRUE)),
        sd_bh_rmse     = purrr::map_dbl(sim_bh_rmse,  ~ stats::sd(.x, na.rm = TRUE)),
        q50_bh_rmse    = purrr::map_dbl(sim_bh_rmse,  ~ stats::quantile(.x, 0.50, na.rm = TRUE, type = 7)),
        q90_bh_rmse    = purrr::map_dbl(sim_bh_rmse,  ~ stats::quantile(.x, 0.90, na.rm = TRUE, type = 7)),
        
        
        
      ) %>%
      dplyr::arrange(model, n, beta_T, beta_E, HR, gam_d, a_x1, a_x2)
    
    print(summ_par)
    
    # save
    saveRDS(list(results = results, summary = summ_par, grid = param_grid),
            file = file.path(DATA_DIR, "fit_results_summary.rds"))
  }
  
  
  
  # printout of sorted q names
  sorted_q_names <- q_names[order(as.numeric(sub("q_", "", q_names)))]
  cat(paste("q_", seq_along(sorted_q_names), ": ", sorted_q_names, sep = ""), sep = "\n")
  invisible(summ_par)
}




#..............................................................................
#..............................................................................
#---- Main ----
#..............................................................................
#..............................................................................

# ---- GLOBAL TIMER START ----
PIPELINE_START_TIME <- Sys.time()

set_parallel_plan()
results <- run_pipeline()

summ_par <- summarize_and_save(results)


print(as.data.frame(summ_par), row.names = FALSE, max = 1e6)


w <- warnings()           # returns an object of class "warnings"
# just a few:
if (!is.null(w)) print(utils::head(w, 10))



summ_par # this contains all things i need


# ---- GLOBAL TIMER END ----
PIPELINE_END_TIME <- Sys.time()

elapsed_secs <- as.numeric(difftime(PIPELINE_END_TIME, PIPELINE_START_TIME, units = "secs"))
elapsed_mins <- elapsed_secs / 60
elapsed_hours <- elapsed_mins / 60

cat("\n================ PIPELINE RUNTIME =================\n")
cat(sprintf("Start time : %s\n", PIPELINE_START_TIME))
cat(sprintf("End time   : %s\n", PIPELINE_END_TIME))
cat(sprintf("Elapsed    : %.2f seconds (%.2f minutes | %.2f hours)\n",
            elapsed_secs, elapsed_mins, elapsed_hours))
cat("==================================================\n")


# ----EOF----







