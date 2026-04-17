# hsbvar-test.R
# Fits H-SBVAR on three simulation settings using both the FISTA solver
# (hsbvar-nll.R) and the exact CVXR solver (hsbvar-cvxr.R), then checks that
# the two objectives agree for each.
#
# Settings:
#   1. d=2, p=1, lambda=0.2  (baseline)
#   2. d=2, p=2, lambda=0.3  (higher lag order)
#   3. d=3, p=1, lambda=0.3  (higher dimension)

source("hsbvar-nll.R")
source("hsbvar-cvxr.R")

# -----------------------------------------------------------------------
# Helper: generate piecewise VAR(p) with one changepoint
# -----------------------------------------------------------------------
make_var <- function(n, d, p, cp_true, B1, B2, seed = 7) {
  set.seed(seed)
  Y <- matrix(0, n, d)
  for (t in 2:n) {
    B      <- if (t <= cp_true) B1 else B2
    idx    <- pmax(t - seq_len(p), 1L)
    x_t    <- as.vector(t(Y[idx, , drop = FALSE]))  # dp-vector
    Y[t, ] <- B %*% x_t + rnorm(d)
  }
  Y
}

# -----------------------------------------------------------------------
# Helper: run both solvers and print comparison
# -----------------------------------------------------------------------
run_comparison <- function(label, Y, p, lambda) {
  cat(strrep("-", 60), "\n")
  cat(label, "\n")
  cat(strrep("-", 60), "\n")

  t_nll  <- system.time(
    fit_nll  <- hsbvar(Y, p = p, lambda = lambda,
                       max_iter = 2000, verbose = FALSE)
  )["elapsed"]
  t_cvxr <- system.time(
    fit_cvxr <- hsbvar_cvxr(Y, p = p, lambda = lambda)
  )["elapsed"]

  cp_nll  <- if (length(fit_nll$cp)  == 0) "none" else paste(fit_nll$cp,  collapse = " ")
  cp_cvxr <- if (length(fit_cvxr$cp) == 0) "none" else paste(fit_cvxr$cp, collapse = " ")

  abs_diff <- abs(fit_nll$obj_val - fit_cvxr$obj_val)
  rel_diff <- abs_diff / (1.0 + abs(fit_cvxr$obj_val))

  cat(sprintf("  NLL FISTA  obj=%.8f  cp=%-20s  %s (%d iter)  %.2fs\n",
              fit_nll$obj_val,  cp_nll,
              fit_nll$stop_crit, fit_nll$n_iter, t_nll))
  cat(sprintf("  NLL CVXR   obj=%.8f  cp=%-20s  [%s]  %.2fs\n",
              fit_cvxr$obj_val, cp_cvxr,
              fit_cvxr$status, t_cvxr))
  cat(sprintf("  rel diff=%.2e  %s\n",
              rel_diff, if (rel_diff < 1e-3) "PASS" else "WARN"))
  cat("\n")
}

# -----------------------------------------------------------------------
# Setting 1: d=2, p=1, lambda=0.2
# -----------------------------------------------------------------------
n       <- 60
d       <- 2
p       <- 1
lambda  <- 0.2
cp_true <- 30
B1 <- matrix(c(0.5, 0.1, 0.1, 0.5),   d, d)
B2 <- matrix(c(-0.4, 0.1, 0.1, -0.4), d, d)
y1 <- make_var(n, d, p, cp_true, B1, B2)
run_comparison(
  sprintf("Setting 1: n=%d  d=%d  p=%d  lambda=%.2f  true cp=%d",
          n, d, p, lambda, cp_true),
  y1, p, lambda
)

# -----------------------------------------------------------------------
# Setting 2: d=2, p=2, lambda=0.3
# -----------------------------------------------------------------------
n       <- 60
d       <- 2
p       <- 2
lambda  <- 0.3
cp_true <- 30
B1 <- cbind(
  matrix(c(0.4, 0.1, 0.1, 0.4), d, d),
  matrix(c(0.1, 0.0, 0.0, 0.1), d, d)
)
B2 <- cbind(
  matrix(c(-0.3, 0.1, 0.1, -0.3), d, d),
  matrix(c(0.0, 0.1, 0.1, 0.0),   d, d)
)
y2 <- make_var(n, d, p, cp_true, B1, B2)
run_comparison(
  sprintf("Setting 2: n=%d  d=%d  p=%d  lambda=%.2f  true cp=%d",
          n, d, p, lambda, cp_true),
  y2, p, lambda
)

# -----------------------------------------------------------------------
# Setting 3: d=3, p=1, lambda=0.3
# -----------------------------------------------------------------------
n       <- 60
d       <- 3
p       <- 1
lambda  <- 0.3
cp_true <- 30
B1 <- matrix(c(0.5, 0.1, 0.0, 0.1, 0.4, 0.1, 0.0, 0.1, 0.5), d, d)
B2 <- matrix(c(-0.4, 0.1, 0.0, 0.1, -0.3, 0.1, 0.0, 0.1, -0.4), d, d)
y3 <- make_var(n, d, p, cp_true, B1, B2)
run_comparison(
  sprintf("Setting 3: n=%d  d=%d  p=%d  lambda=%.2f  true cp=%d",
          n, d, p, lambda, cp_true),
  y3, p, lambda
)
