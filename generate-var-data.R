# generate-var-data.R
# Simulation data generators for piecewise-stationary VAR(p) processes.

#' Generate a piecewise-stationary VAR(p) process (core generator)
#'
#' @param n            Total number of observations.
#' @param break_points Integer vector of break-point indices (1-based).
#'   Each index is the FIRST observation of the new regime.
#'   Length m0 for m0 breaks (m0 + 1 regimes).
#' @param a_list       List of coefficient-matrix lists, one per regime.
#'   Each element is a list of p d×d matrices [A_1, ..., A_p] such that
#'   Y_t = A_1 Y_{t-1} + ... + A_p Y_{t-p} + eps_t.
#' @param sigma_list   List of d×d positive-definite innovation covariance
#'   matrices, one per regime.
#' @param seed         RNG seed.
#' @param n_burnin     Number of burn-in observations to discard (default 200).
#' @return List with components:
#'   \item{y}{n × d matrix of observations.}
#'   \item{n}{Total observations.}
#'   \item{d}{Dimension of the process.}
#'   \item{p}{VAR lag order.}
#'   \item{break_points}{As supplied.}
#'   \item{a_list}{As supplied.}
#'   \item{sigma_list}{As supplied.}
#'   \item{true_sigma}{Length-n list of per-observation true covariance matrices.}
generate_var_piecewise <- function(n, break_points, a_list, sigma_list,
                                   seed = 42L, n_burnin = 200L) {
  n_regimes <- length(a_list)
  stopifnot(length(sigma_list) == n_regimes)
  stopifnot(length(break_points) == n_regimes - 1L)

  p <- length(a_list[[1L]])
  d <- nrow(a_list[[1L]][[1L]])

  for (r in seq_len(n_regimes)) {
    stopifnot(length(a_list[[r]]) == p)
    for (lag in seq_len(p)) {
      stopifnot(identical(dim(a_list[[r]][[lag]]), c(d, d)))
    }
    stopifnot(identical(dim(sigma_list[[r]]), c(d, d)))
  }

  n_total <- n + p + n_burnin

  set.seed(seed)

  # Regime assignment per observation index (1-based, length n)
  regime <- rep(1L, n)
  for (j in seq_along(break_points)) {
    regime[seq(break_points[j], n)] <- j + 1L
  }

  # Pre-compute upper Cholesky factors: chol(Sigma) = R s.t. t(R) %*% R = Sigma
  chol_list <- lapply(sigma_list, chol)

  # Draw all standard-normal innovations up front.
  # Draw exactly (n + n_burnin) rows — matching var.sim.break's nT = nobs + skip —
  # so that the RNG stream is identical when the same seed is used.
  # Noise for simulation step t (loop index p+1 .. n_total) lives at row t-p.
  n_noise <- n + n_burnin
  raw_noise <- matrix(rnorm(n_noise * d), nrow = n_noise, ncol = d)

  # Simulate the piecewise VAR with burn-in
  # t_obs <= 0 during burn-in (uses regime 1); t_obs in 1..n during observations
  y_mat <- matrix(0.0, nrow = n_total, ncol = d)

  for (t in (p + 1L):n_total) {
    t_obs <- t - p - n_burnin
    r <- if (t_obs >= 1L) regime[t_obs] else 1L
    coef_mats <- a_list[[r]]

    mu <- numeric(d)
    for (lag in seq_len(p)) {
      mu <- mu + coef_mats[[lag]] %*% y_mat[t - lag, ]
    }

    # innov = t(R) %*% z, z ~ N(0, I)  =>  innov ~ N(0, Sigma_r)
    # Index t-p maps the simulation step to the noise draw row (1-based).
    innov <- crossprod(chol_list[[r]], raw_noise[t - p, ])

    y_mat[t, ] <- mu + innov
  }

  y_mat <- y_mat[(p + n_burnin + 1L):n_total, , drop = FALSE]

  list(
    y            = y_mat,
    n            = n,
    d            = d,
    p            = p,
    break_points = break_points,
    a_list       = a_list,
    sigma_list   = sigma_list,
    true_sigma   = sigma_list[regime]
  )
}


#' Scenario 1: Small VAR(1), d=2, one coefficient break.
#'
#' Two regimes (n=200, break at t=100):
#'   1. A = [[0.60, 0.10], [0.10, 0.55]], Sigma = [[1.0, 0.2], [0.2, 1.0]]
#'   2. A = [[-0.50, 0.10], [0.10, -0.45]], Sigma = [[1.5, -0.3], [-0.3, 1.5]]
#'
#' Clear sign-flip in the diagonal of A at t=100; moderate covariance shift.
#'
#' @param seed RNG seed (default 42).
#' @return List from generate_var_piecewise.
generate_scenario1 <- function(sigscale = 1, seed = 42L) {
  d <- 2L
  A1 <- matrix(c(0.60, 0.10, 0.10, 0.55), d, d)
  Sigma1 <- matrix(c(1.0, 0.2, 0.2, 1.0), d, d)
  A2 <- matrix(c(-0.50, 0.10, 0.10, -0.45), d, d)
  Sigma2 <- matrix(c(1.5, -0.3, -0.3, 1.5), d, d)
  generate_var_piecewise(
    n            = 200L,
    break_points = 100L,
    a_list       = list(list(A1), list(A2)),
    sigma_list   = list(sigscale * Sigma1, sigscale * Sigma2),
    seed         = seed,
    n_burnin     = 200L
  )
}


#' Scenario 2: Small VAR(1), d=2, one coefficient break.
#'
#' Two regimes (n=200, break at t=100):
#'   1. A = [[0.60, 0.10], [0.10, 0.55]], Sigma = [[.1, 0], [0, .1]]
#'   2. A = [[-0.50, 0.10], [0.10, -0.45]], Sigma = [[.1, 0], [0, .1]]
#'
#' Clear sign-flip in the diagonal of A at t=100; moderate covariance shift.
#'
#' @param seed RNG seed (default 42).
#' @return List from generate_var_piecewise.
generate_scenario2 <- function(sigscale = 1, seed = 42L) {
  d <- 2L
  A1 <- matrix(c(0.60, 0.10, 0.10, 0.55), d, d)
  Sigma1 <- matrix(c(.1, 0, 0, .1), d, d)
  A2 <- matrix(c(-0.50, 0.10, 0.10, -0.45), d, d)
  Sigma2 <- matrix(c(.1, 0, 0, .1), d, d)
  generate_var_piecewise(
    n            = 200L,
    break_points = 100L,
    a_list       = list(list(A1), list(A2)),
    sigma_list   = list(sigscale * Sigma1, sigscale * Sigma2),
    seed         = seed,
    n_burnin     = 200L
  )
}


#' Scenario 3: Large sparse VAR(1), d=20, two breaks.
#'
#' Mirrors the DGP from Safikhani-Shojaie SBDetection_simulation_1.R.
#' Three equal regimes of 100 observations each (n=300, breaks at t=101, 201).
#' Each regime's 20×20 coefficient matrix has only a sparse superdiagonal:
#'   1. A[j, j+1] = -0.6  for j = 1..19
#'   2. A[j, j+1] =  0.75 for j = 1..19
#'   3. A[j, j+1] = -0.8  for j = 1..19
#' Innovation covariance: 0.01 * I_20 (identical across regimes).
#'
#' @param seed RNG seed (default 123456).
#' @return List from generate_var_piecewise.
generate_scenario3 <- function(sigscale = 1, seed = 123456L) {
  k <- 20L
  make_superdiag <- function(val) {
    A <- matrix(0.0, k, k)
    for (j in seq_len(k - 1L)) A[j, j + 1L] <- val
    A
  }
  Sigma <- 0.01 * diag(k)
  generate_var_piecewise(
    n = 300L,
    break_points = c(100L, 200L),
    a_list = list(
      list(make_superdiag(-0.6)),
      list(make_superdiag(0.75)),
      list(make_superdiag(-0.8))
    ),
    sigma_list = list(sigscale * Sigma, sigscale * Sigma, sigscale * Sigma),
    seed = seed,
    n_burnin = 200L
  )
}
