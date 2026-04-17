# run-sim-sbvar.R
# Simulation harness: H-SBVAR (NLL, D-trace) vs SB-VAR (3 BEA IC variants).
# Each method: CV-selected lambda + Stage 2 BEA pruning.
#
# Usage:
#   Rscript run-sim-sbvar.R --scenario=N [--nrep=N] [--outdir=PATH] [--sigscale=X]
#
#   --scenario=N   Required. 1, 2, or 3 (see generate-var-data.R).
#   --nrep=N       Number of replications (default 100).
#   --outdir=PATH  Output directory for .rds files (default ./sim).
#   --sigscale=X   Scale factor applied to all innovation covariance matrices (default 1).
#
# Methods compared (all use Stage 1 CV + Stage 2 BEA):
#   nll_s2      H-SBVAR NLL
#   dt_s2       H-SBVAR D-trace
#   sb_s2_rss   SB-VAR, BEA ic_type="rss"
#   sb_s2_sig   SB-VAR, BEA ic_type="sigma_scaled"
#   sb_s2_pro   SB-VAR, BEA ic_type="profile_lik"
#
# Stage 1 results (nll_s1, dt_s1, sb_s1) are also recorded.
# Results saved incrementally to run-sim-sbvar-<label>.rds.

# ============================================================
# CLI arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, key, default = NULL) {
  pat <- paste0("^--", key, "=(.+)$")
  m <- regmatches(args, regexpr(pat, args, perl = TRUE))
  if (length(m) == 0L) default else sub(pat, "\\1", m)
}

scenario_str <- parse_arg(args, "scenario")
if (is.null(scenario_str)) stop("--scenario=N is required (1, 2, or 3)")
scenario_arg <- as.integer(scenario_str)
if (!scenario_arg %in% 1:3) stop("--scenario must be 1, 2, or 3")

nrep_arg <- as.integer(parse_arg(args, "nrep", "100"))
outdir_arg <- parse_arg(args, "outdir", "./sim")
sigscale_arg <- as.numeric(parse_arg(args, "sigscale", "1"))

# ============================================================
# Dependencies
# ============================================================

source("generate-var-data.R")
source("cv-hsbvar-nll.R")
source("cv-hsbvar-dtrace.R")
source("hsbvar-bea.R")
source(file.path("ss", "cv_sb_stage1_clean.R"))
source(file.path("ss", "sb_stage2_bea_clean.R"))

# ============================================================
# Scenario dispatch
# ============================================================

GENERATE_FN <- switch(as.character(scenario_arg),
  "1" = generate_scenario1,
  "2" = generate_scenario2,
  "3" = generate_scenario3
)

meta <- GENERATE_FN(sigscale = sigscale_arg, seed = 0L)
N <- meta$n
P <- meta$p
D <- meta$d
TRUE_BREAKS <- meta$break_points
M0 <- length(TRUE_BREAKS)

# ============================================================
# Helpers
# ============================================================

#' Symmetric Hausdorff distance between two integer sets.
#' Returns n when exactly one set is empty.
hausdorff_dist <- function(est, tru, n) {
  if (length(est) == 0L && length(tru) == 0L) {
    return(0)
  }
  if (length(est) == 0L || length(tru) == 0L) {
    return(n)
  }
  max(
    max(vapply(tru, function(b) min(abs(est - b)), numeric(1L))),
    max(vapply(est, function(a) min(abs(tru - a)), numeric(1L)))
  )
}

#' Build n x (d*dp) matrix of true VAR coefficients (vec column-major per regime).
make_true_B_mat <- function(n, break_points, a_list, d, p) {
  dp <- d * p
  regime <- rep(1L, n)
  for (j in seq_along(break_points)) regime[seq(break_points[j], n)] <- j + 1L
  B_mat <- matrix(0, n, d * dp)
  for (r in seq_along(a_list)) {
    b_r <- unlist(lapply(a_list[[r]], as.vector)) # vec([A_1 | ... | A_p])
    rows <- which(regime == r)
    B_mat[rows, ] <- matrix(b_r, nrow = length(rows), ncol = d * dp, byrow = TRUE)
  }
  B_mat
}

#' Build n x d^2 matrix of true innovation covariances (vec per regime).
make_true_Sigma_mat <- function(n, break_points, sigma_list) {
  d2 <- length(as.vector(sigma_list[[1L]]))
  regime <- rep(1L, n)
  for (j in seq_along(break_points)) regime[seq(break_points[j], n)] <- j + 1L
  S_mat <- matrix(0, n, d2)
  for (r in seq_along(sigma_list)) {
    s_r <- as.vector(sigma_list[[r]])
    rows <- which(regime == r)
    S_mat[rows, ] <- matrix(s_r, nrow = length(rows), ncol = d2, byrow = TRUE)
  }
  S_mat
}

#' n x d matrix of one-step-ahead predictions from an n x (d*dp) B matrix.
compute_yhat_mv <- function(Y, B, p) {
  n <- nrow(Y)
  d <- ncol(Y)
  dp <- d * p
  Y_ext <- rbind(matrix(0, p, d), Y)
  X <- matrix(0, n, dp)
  for (lag in seq_len(p)) {
    cols <- ((lag - 1L) * d + 1L):(lag * d)
    X[, cols] <- Y_ext[(p - lag + 1L):(n + p - lag), , drop = FALSE]
  }
  Y_hat <- matrix(0, n, d)
  for (t in seq_len(n)) Y_hat[t, ] <- matrix(B[t, ], d, dp) %*% X[t, ]
  Y_hat
}

#' Compute metrics for one (stage, method) result.
compute_metrics <- function(cp_est, B_est, Sigma_est,
                            Y, true_B, true_Sigma, n, p, d) {
  ncp <- length(cp_est)
  Y_hat <- compute_yhat_mv(Y, B_est, p)
  list(
    ncp         = ncp,
    correct_ncp = (ncp == M0),
    hd          = hausdorff_dist(cp_est, TRUE_BREAKS, n),
    mse         = mean((Y - Y_hat)^2),
    B_err       = mean((B_est - true_B)^2),
    Sigma_err   = mean((Sigma_est - true_Sigma)^2),
    rel_cp      = if (ncp > 0L) sort(as.integer(cp_est)) / n else numeric(0),
    cp          = if (ncp > 0L) sort(as.integer(cp_est)) else integer(0)
  )
}

# ============================================================
# Single replication
# ============================================================

one_rep <- function(seed) {
  dat <- GENERATE_FN(sigscale = sigscale_arg, seed = seed)
  Y <- dat$y
  n <- dat$n
  p <- dat$p
  d <- dat$d

  true_B <- make_true_B_mat(n, TRUE_BREAKS, dat$a_list, d, p)
  true_Sigma <- make_true_Sigma_mat(n, TRUE_BREAKS, dat$sigma_list)

  n_eff <- n - p
  omega_author <- (1 / 1.75) * log(n_eff) * log(d)
  lambda2 <- log(n_eff) * log(d) / n_eff

  out <- list()

  # ---- H-SBVAR NLL -----------------------------------------------------------
  tryCatch(
    {
      cv_nll <- cv_hsbvar_nll(Y, p = p, verbose = FALSE)
      fit_nll <- hsbvar(Y, p = p, lambda = cv_nll$best$lambda)
      out$nll_s1 <- compute_metrics(
        fit_nll$cp, fit_nll$B_arr, fit_nll$Sigma_arr,
        Y, true_B, true_Sigma, n, p, d
      )
      bea_nll <- hsbvar_bea(fit_nll, Y, p = p)
      out$nll_s2 <- compute_metrics(
        bea_nll$cp, bea_nll$B, bea_nll$Sigma,
        Y, true_B, true_Sigma, n, p, d
      )
      out$nll_lambda <- cv_nll$best$lambda
    },
    error = function(e) {
      message(sprintf("  [seed %d] NLL error: %s", seed, conditionMessage(e)))
    }
  )

  # ---- H-SBVAR D-trace -------------------------------------------------------
  tryCatch(
    {
      cv_dt <- cv_hsbvar_dtrace(Y, p = p, verbose = FALSE)
      fit_dt <- hsbvar_dtrace(Y, p = p, lambda = cv_dt$best$lambda)
      out$dt_s1 <- compute_metrics(
        fit_dt$cp, fit_dt$B_arr, fit_dt$Sigma_arr,
        Y, true_B, true_Sigma, n, p, d
      )
      bea_dt <- hsbvar_bea(fit_dt, Y, p = p)
      out$dt_s2 <- compute_metrics(
        bea_dt$cp, bea_dt$B, bea_dt$Sigma,
        Y, true_B, true_Sigma, n, p, d
      )
      out$dt_lambda <- cv_dt$best$lambda
    },
    error = function(e) {
      message(sprintf("  [seed %d] D-trace error: %s", seed, conditionMessage(e)))
    }
  )

  # ---- SB-VAR (shared Stage 1; three BEA variants) ---------------------------
  tryCatch(
    {
      cv_sb <- cv_sbvar(Y, p = p, verbose = FALSE)
      fit_sb <- sbvar(Y, p = p, lambda = cv_sb$best$lambda)
      out$sb_s1 <- compute_metrics(
        fit_sb$cp, fit_sb$B_arr, fit_sb$Sigma_arr,
        Y, true_B, true_Sigma, n, p, d
      )

      bea_rss <- sbvar_bea(fit_sb, Y,
        p = p, lambda = lambda2,
        omega_n = omega_author, ic_type = "rss"
      )
      bea_sig <- sbvar_bea(fit_sb, Y,
        p = p,
        omega_n = omega_author, ic_type = "sigma_scaled"
      )
      bea_pro <- sbvar_bea(fit_sb, Y,
        p = p,
        omega_n = omega_author, ic_type = "profile_lik"
      )

      out$sb_s2_rss <- compute_metrics(
        bea_rss$cp, bea_rss$B, bea_rss$Sigma,
        Y, true_B, true_Sigma, n, p, d
      )
      out$sb_s2_sig <- compute_metrics(
        bea_sig$cp, bea_sig$B, bea_sig$Sigma,
        Y, true_B, true_Sigma, n, p, d
      )
      out$sb_s2_pro <- compute_metrics(
        bea_pro$cp, bea_pro$B, bea_pro$Sigma,
        Y, true_B, true_Sigma, n, p, d
      )
      out$sb_lambda <- cv_sb$best$lambda
    },
    error = function(e) {
      message(sprintf("  [seed %d] SB-VAR error: %s", seed, conditionMessage(e)))
    }
  )

  out
}

# ============================================================
# Main simulation loop
# ============================================================

run_label <- sprintf("scenario%d_sigscale%g", scenario_arg, sigscale_arg)

cat(sprintf(
  "\n=== SB-VAR simulation: Scenario %d  n=%d  d=%d  p=%d  m0=%d  true breaks: %s ===\n\n",
  scenario_arg, N, D, P, M0, paste(TRUE_BREAKS, collapse = ", ")
))
cat(sprintf("nrep=%d  outdir=%s  sigscale=%g\n\n", nrep_arg, outdir_arg, sigscale_arg))

dir.create(outdir_arg, showWarnings = FALSE, recursive = TRUE)
rds_path <- file.path(outdir_arg, sprintf("run-sim-sbvar-%s.rds", run_label))
cat(sprintf("Results saved incrementally to %s\n\n", rds_path))

t_start <- proc.time()
results <- vector("list", nrep_arg)

for (i in seq_len(nrep_arg)) {
  cat(sprintf("[%3d/%d]", i, nrep_arg))
  results[[i]] <- one_rep(i)
  local({
    r <- results[[i]]
    fcp <- function(key) {
      v <- r[[key]]
      if (is.null(v) || length(v$cp) == 0L) "none" else paste(v$cp, collapse = ",")
    }
    cat(sprintf(
      "  true=[%s]  nll=[%s]  dt=[%s]  sb_rss=[%s]  sb_sig=[%s]  sb_pro=[%s]",
      paste(TRUE_BREAKS, collapse = ","),
      fcp("nll_s2"), fcp("dt_s2"),
      fcp("sb_s2_rss"), fcp("sb_s2_sig"), fcp("sb_s2_pro")
    ))
  })
  cat("\n")
  saveRDS(results[seq_len(i)], rds_path)
}

elapsed <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("\nTotal elapsed: %.1f s  (%.1f s/rep)\n\n", elapsed, elapsed / nrep_arg))
cat(sprintf("Results saved to %s\n\n", rds_path))

# ============================================================
# Summarisation helpers
# ============================================================

extract_scalar <- function(results, variant, field) {
  vapply(results, function(r) {
    v <- r[[variant]]
    if (is.null(v)) NA_real_ else as.numeric(v[[field]])
  }, numeric(1L))
}

extract_rel_cp <- function(results, variant, m0) {
  t(vapply(results, function(r) {
    v <- r[[variant]]
    if (is.null(v) || length(v$rel_cp) != m0) rep(NA_real_, m0) else v$rel_cp
  }, numeric(m0)))
}

summarise_variant <- function(results, variant, m0) {
  ncp <- extract_scalar(results, variant, "ncp")
  correct <- extract_scalar(results, variant, "correct_ncp")
  hd <- extract_scalar(results, variant, "hd")
  mse <- extract_scalar(results, variant, "mse")
  B_err <- extract_scalar(results, variant, "B_err")
  Sig_err <- extract_scalar(results, variant, "Sigma_err")
  rel_cp <- extract_rel_cp(results, variant, m0)

  n_ok <- sum(!is.na(ncp))
  n_correct <- sum(correct, na.rm = TRUE)

  rel_ok <- rel_cp[!is.na(rel_cp[, 1L, drop = FALSE]), , drop = FALSE]
  n_rel <- nrow(rel_ok)
  mean_rel <- if (n_rel > 0L) colMeans(rel_ok) else rep(NA_real_, m0)
  se_rel <- if (n_rel > 1L) apply(rel_ok, 2L, sd) / sqrt(n_rel) else rep(NA_real_, m0)

  list(
    n_ok          = n_ok,
    avg_ncp       = mean(ncp, na.rm = TRUE),
    pct_correct   = 100 * n_correct / n_ok,
    mean_rel      = mean_rel,
    se_rel        = se_rel,
    mean_hd       = mean(hd, na.rm = TRUE),
    mean_mse      = mean(mse, na.rm = TRUE),
    mean_B_err    = mean(B_err, na.rm = TRUE),
    mean_Sig_err  = mean(Sig_err, na.rm = TRUE)
  )
}

# ============================================================
# Print results
# ============================================================

VARIANTS <- c(
  "nll_s1", "nll_s2",
  "dt_s1", "dt_s2",
  "sb_s1", "sb_s2_rss", "sb_s2_sig", "sb_s2_pro"
)
LABELS <- c(
  "NLL S1", "NLL S2",
  "Dtrace S1", "Dtrace S2",
  "SB-VAR S1", "SB-VAR rss", "SB-VAR sig", "SB-VAR pro"
)

sums <- setNames(
  lapply(VARIANTS, summarise_variant, results = results, m0 = M0),
  VARIANTS
)

# Build relative-location header dynamically
rel_hdr <- paste(sprintf(
  "%9s %7s", sprintf("Mb%d", seq_len(M0)),
  sprintf("SEb%d", seq_len(M0))
), collapse = " ")

# --- Table 1: break-point estimation ----------------------------------------
cat("==========================================================================\n")
cat(sprintf("Table 1: Break-point estimation  [Scenario %d]\n", scenario_arg))
cat(sprintf(
  "         True breaks: %s  |  True rel: %s\n",
  paste(TRUE_BREAKS, collapse = ", "),
  paste(sprintf("%.4f", TRUE_BREAKS / N), collapse = ", ")
))
cat("==========================================================================\n")

hdr <- sprintf("%-12s %8s %8s %s", "Method", "Avg ncp", "% m=m0", rel_hdr)
cat(hdr, "\n")
cat(strrep("-", nchar(hdr)), "\n")

for (k in seq_along(VARIANTS)) {
  s <- sums[[k]]
  locs <- paste(vapply(seq_len(M0), function(j) {
    sprintf("%9.4f %7.4f", s$mean_rel[j], s$se_rel[j])
  }, character(1L)), collapse = " ")
  cat(sprintf("%-12s %8.2f %8.1f %s\n", LABELS[k], s$avg_ncp, s$pct_correct, locs))
}

# --- Table 2: estimation accuracy -------------------------------------------
cat("\n==========================================================================\n")
cat(sprintf("Table 2: Estimation accuracy  [Scenario %d]\n", scenario_arg))
cat("==========================================================================\n")

hdr2 <- sprintf(
  "%-12s %10s %10s %10s %12s", "Method", "Mean HD", "Mean MSE",
  "Mean B-err", "Mean Sig-err"
)
cat(hdr2, "\n")
cat(strrep("-", nchar(hdr2)), "\n")

for (k in seq_along(VARIANTS)) {
  s <- sums[[k]]
  cat(sprintf(
    "%-12s %10.2f %10.5f %10.5f %12.5f\n",
    LABELS[k], s$mean_hd, s$mean_mse, s$mean_B_err, s$mean_Sig_err
  ))
}

# --- Lambda summary ----------------------------------------------------------
fmt_lam <- function(key) {
  v <- vapply(results, function(r) {
    x <- r[[key]]
    if (is.null(x)) NA_real_ else as.numeric(x)
  }, numeric(1L))
  sprintf(
    "median=%.4g  [%.4g, %.4g]",
    median(v, na.rm = TRUE), min(v, na.rm = TRUE), max(v, na.rm = TRUE)
  )
}
cat(sprintf("\nNLL lambda:    %s\n", fmt_lam("nll_lambda")))
cat(sprintf("D-trace lambda: %s\n", fmt_lam("dt_lambda")))
cat(sprintf("SB-VAR lambda:  %s\n", fmt_lam("sb_lambda")))

cat("\nDone.\n")
