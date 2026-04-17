# hsbvar-bea.R
# Stage 2: Backward Elimination Algorithm for H-SBVAR (h-sbvar-bea.md).
#
# The Stage 1 joint group LASSO (hsbvar) over-selects candidate changepoints.
# This module prunes them using the natural-parameter IC from h-sbvar-bea.md:
# the total profiled NLL across segments plus m times a per-break penalty
# omega_n.  Removing a candidate tau means setting Theta_tau = Psi_tau = 0,
# collapsing both the coefficient and precision jumps at that time.
#
# The per-segment NLL is the Gaussian MLE: (n_j / 2) * log det(Sigma_hat_j).
# This is the multivariate analogue of the univariate (n_j / 2) * log(RSS/n_j)
# used in hsbar-bea.R, and equals the minimised natural-parameter loss on the
# segment up to the additive constant n * d * (1 + log(2*pi)) / 2 which is
# common to all segmentations and therefore cancels in IC comparisons.

source("hsbvar-nll.R")

# ------------------------------------------------------------------
# Internal: profiled NLL contribution for one segment (MLE Gaussian VAR)
# Returns (n_j / 2) * log det Sigma_hat, Inf if underdetermined.
# ------------------------------------------------------------------

seg_nll_var <- function(Y_seg, X_seg) {
  n_j <- nrow(Y_seg)
  dp <- ncol(X_seg)
  if (n_j <= dp) {
    return(Inf)
  }

  # B_hat: d x dp  (each row = VAR coefficients for one output series)
  B_hat <- tryCatch(
    t(solve(crossprod(X_seg), crossprod(X_seg, Y_seg))),
    error = function(e) {
      cf <- qr.coef(qr(X_seg), Y_seg) # dp x d
      cf[is.na(cf)] <- 0
      t(cf) # d x dp
    }
  )

  E <- Y_seg - X_seg %*% t(B_hat) # n_j x d  residuals
  Sigma_j <- crossprod(E) / n_j # d x d  MLE covariance
  Sigma_j <- (Sigma_j + t(Sigma_j)) * 0.5 # enforce symmetry

  ch <- tryCatch(chol(Sigma_j), error = function(e) NULL)
  if (is.null(ch)) {
    return(Inf)
  }

  (n_j / 2) * 2 * sum(log(diag(ch))) # (n_j/2) * log det Sigma_j
}

# ------------------------------------------------------------------
# Internal: total NLL = sum of seg_nll_var over all segments defined
# by the changepoint vector cps.
# ------------------------------------------------------------------

compute_g_var <- function(Y, Y_lag, cps, n) {
  breaks <- c(1L, cps, n + 1L)
  total <- 0.0
  for (b in seq_len(length(breaks) - 1L)) {
    rows <- breaks[b]:(breaks[b + 1L] - 1L)
    total <- total + seg_nll_var(
      Y[rows, , drop = FALSE],
      Y_lag[rows, , drop = FALSE]
    )
  }
  total
}

# ------------------------------------------------------------------
# Internal: segment-wise MLE B and Sigma for the pruned model.
# Returns B (n x d*dp) and Sigma (n x d^2) in the same vec'd layout
# as B_arr / Sigma_arr in hsbvar-nll.R.
# ------------------------------------------------------------------

refit_segments_var <- function(Y, Y_lag, cps, n, d, dp) {
  breaks <- c(1L, cps, n + 1L)
  B_out <- matrix(NA_real_, n, d * dp)
  Sig_out <- matrix(NA_real_, n, d * d)

  for (b in seq_len(length(breaks) - 1L)) {
    rows <- breaks[b]:(breaks[b + 1L] - 1L)
    n_j <- length(rows)
    Y_j <- Y[rows, , drop = FALSE]
    X_j <- Y_lag[rows, , drop = FALSE]

    B_j <- tryCatch(
      t(solve(crossprod(X_j), crossprod(X_j, Y_j))), # d x dp
      error = function(e) {
        cf <- qr.coef(qr(X_j), Y_j)
        cf[is.na(cf)] <- 0
        t(cf)
      }
    )

    E_j <- Y_j - X_j %*% t(B_j)
    Sigma_j <- crossprod(E_j) / n_j
    Sigma_j <- (Sigma_j + t(Sigma_j)) * 0.5

    B_out[rows, ] <- matrix(as.vector(B_j), nrow = n_j, ncol = d * dp, byrow = TRUE)
    Sig_out[rows, ] <- matrix(as.vector(Sigma_j), nrow = n_j, ncol = d * d, byrow = TRUE)
  }

  list(B = B_out, Sigma = Sig_out)
}

# ------------------------------------------------------------------
# Internal: data-driven omega_n selection (h-sbvar-bea.md, Section 5).
# Runs BEA to empty, records the loss increment Delta_i at each removal,
# clusters {Delta_i} by largest-gap rule, returns the minimum jump in
# the "large" (true-break) cluster.
# ------------------------------------------------------------------

select_omega_hsbvar <- function(Y, Y_lag, cps_init, n) {
  cps <- sort(cps_init)
  deltas <- numeric(length(cps))
  idx <- 0L

  while (length(cps) > 0L) {
    m <- length(cps)
    g_cur <- compute_g_var(Y, Y_lag, cps, n)

    g_rem <- vapply(seq_len(m), function(i) {
      compute_g_var(Y, Y_lag, cps[-i], n)
    }, numeric(1L))

    best_i <- which.min(g_rem - g_cur)
    idx <- idx + 1L
    deltas[idx] <- (g_rem - g_cur)[best_i]
    cps <- cps[-best_i]
  }

  deltas <- deltas[seq_len(idx)]
  if (idx < 2L) {
    return(max(deltas) * 1.01)
  }

  sorted <- sort(deltas)
  split <- which.max(diff(sorted))
  large <- sorted[(split + 1L):length(sorted)]
  if (length(large) == 0L) {
    return(max(sorted) * 1.01)
  }
  min(large)
}

# ------------------------------------------------------------------
# Public: BEA pruning of Stage 1 changepoints for H-SBVAR
# ------------------------------------------------------------------

#' Backward Elimination Algorithm for H-SBVAR (h-sbvar-bea.md)
#'
#' Prunes the over-selected candidate changepoints produced by Stage 1
#' (\code{hsbvar}) using the natural-parameter IC.  Removing a candidate
#' \eqn{\tau} sets \eqn{\Theta_\tau = \Psi_\tau = 0}, collapsing both
#' the coefficient and precision jumps at that time.
#'
#' @param fit     Output of \code{hsbvar()}.
#' @param Y       n x d numeric matrix (same data passed to Stage 1).
#' @param p       VAR lag order.  Inferred from the dimensions of
#'                \code{fit$Theta} when not supplied.
#' @param omega_n Penalty per break.  Default: BIC penalty
#'                \eqn{(d \cdot dp + d(d+1)/2)\log(n)/2}.
#'                Pass \code{"data"} to use the data-driven selection
#'                described in h-sbvar-bea.md Section 5.
#'
#' @return A list with:
#' \describe{
#'   \item{cp}{Integer vector of refined changepoint locations.}
#'   \item{B}{n x (d*dp) matrix of segment-wise MLE VAR coefficients
#'     (\code{vec(B_j)} stored column-major, constant within each segment).}
#'   \item{Sigma}{n x d^2 matrix of segment-wise MLE covariances.}
#'   \item{ic}{IC value at the final changepoint set.}
#'   \item{omega_n}{Penalty per break used.}
#' }
hsbvar_bea <- function(fit, Y, p = NULL, omega_n = NULL) {
  n <- nrow(Y)
  d <- ncol(Y)

  if (is.null(p)) {
    # infer p from the number of columns of Theta: ncol = d * d * p
    p <- ncol(fit$Theta) / (d * d)
    if (p != round(p) || p < 1L) stop("Cannot infer p from fit; supply p explicitly.")
    p <- as.integer(round(p))
  }

  dp <- d * p

  # ---- lagged regressor matrix (same convention as hsbvar-nll.R) -----------
  Y_ext <- rbind(matrix(0, p, d), Y) # (n+p) x d
  Y_lag <- matrix(0.0, n, dp)
  for (k in seq_len(p)) {
    cols <- ((k - 1L) * d + 1L):(k * d)
    Y_lag[, cols] <- Y_ext[(p - k + 1L):(n + p - k), , drop = FALSE]
  }

  # ---- resolve omega_n -------------------------------------------------------
  if (identical(omega_n, "data")) {
    cps_init <- sort(fit$cp)
    omega_n <- select_omega_hsbvar(Y, Y_lag, cps_init, n)
  } else if (is.null(omega_n)) {
    omega_n <- (d * dp + d * (d + 1L) / 2) * log(n) / 2
  }

  # ---- initialise ------------------------------------------------------------
  cps <- sort(fit$cp)

  if (length(cps) == 0L) {
    seg <- refit_segments_var(Y, Y_lag, cps, n, d, dp)
    return(list(
      cp      = integer(0L),
      B       = seg$B,
      Sigma   = seg$Sigma,
      ic      = compute_g_var(Y, Y_lag, cps, n),
      omega_n = omega_n
    ))
  }

  current_ic <- compute_g_var(Y, Y_lag, cps, n) + length(cps) * omega_n

  # ---- iterative backward elimination ----------------------------------------
  repeat {
    m <- length(cps)
    if (m == 0L) break

    ic_without <- vapply(seq_len(m), function(i) {
      compute_g_var(Y, Y_lag, cps[-i], n) + (m - 1L) * omega_n
    }, numeric(1L))

    best_i <- which.min(ic_without)

    if (ic_without[best_i] <= current_ic) {
      cps <- cps[-best_i]
      current_ic <- ic_without[best_i]
    } else {
      break
    }
  }

  # ---- refit pruned segments -------------------------------------------------
  seg <- refit_segments_var(Y, Y_lag, cps, n, d, dp)

  list(
    cp      = cps,
    B       = seg$B,
    Sigma   = seg$Sigma,
    ic      = current_ic,
    omega_n = omega_n
  )
}
