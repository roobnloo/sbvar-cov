# sb_stage2_bea_clean.R
# Stage 2: Backward Elimination Algorithm for SB-VAR (safikhani-supplemental.md, Appendix C).
#
# Faithfully implements the BEA of Safikhani & Shojaie (2020): the segment loss
# L_n is the LASSO objective (RSS + lambda * ||Phi||_1) returned by var.lasso.brk.
# This captures only coefficient changes -- covariance shifts are invisible to L_n.
#
# Interface and structure mirror hsbvar-bea.R for direct method comparison.
# The segment loss is the key methodological difference: hsbvar uses the Gaussian
# MLE NLL (n_j/2)*log(det(Sigma_j)), which responds to both coefficient and
# covariance changes; sbvar uses the LASSO objective, which does not.
#
# Omega selection uses the same largest-gap rule as hsbvar-bea.R (instead of the
# k-means described in Appendix G) so that omega calibration is comparable.
#
# B and Sigma outputs are OLS refits (not LASSO), matching hsbvar-bea.R convention.

source(file.path("ss", "sb_estimation.R"))

# ------------------------------------------------------------------
# Internal: OLS Frobenius RSS for one segment  (sigma_scaled base loss).
# Returns tr(E_j'E_j) via QR.  Returns Inf if underdetermined.
# ------------------------------------------------------------------

seg_ols_loss_sb <- function(Y_seg, p) {
  n_j <- nrow(Y_seg)
  d   <- ncol(Y_seg)
  if (n_j <= p) return(Inf)

  Y_ext <- rbind(matrix(0, p, d), Y_seg)
  X_j   <- matrix(0.0, n_j, d * p)
  for (lag in seq_len(p)) {
    cols <- ((lag - 1L) * d + 1L):(lag * d)
    X_j[, cols] <- Y_ext[(p - lag + 1L):(n_j + p - lag), , drop = FALSE]
  }

  B_j <- tryCatch({
    cf <- qr.coef(qr(X_j), Y_seg)
    cf[is.na(cf)] <- 0
    cf
  }, error = function(e) matrix(0, nrow(X_j), d))

  E_j <- Y_seg - X_j %*% B_j
  sum(E_j^2)
}

# ------------------------------------------------------------------
# Internal: profile likelihood loss for one segment.
# Returns (n_j/2) * log(det(Sigma_j)),  Sigma_j = E_j'E_j / n_j.
# Returns Inf if underdetermined or det <= 0.
# ------------------------------------------------------------------

seg_profile_loss_sb <- function(Y_seg, p) {
  n_j <- nrow(Y_seg)
  d   <- ncol(Y_seg)
  if (n_j <= p) return(Inf)

  Y_ext <- rbind(matrix(0, p, d), Y_seg)
  X_j   <- matrix(0.0, n_j, d * p)
  for (lag in seq_len(p)) {
    cols <- ((lag - 1L) * d + 1L):(lag * d)
    X_j[, cols] <- Y_ext[(p - lag + 1L):(n_j + p - lag), , drop = FALSE]
  }

  B_j <- tryCatch({
    cf <- qr.coef(qr(X_j), Y_seg)
    cf[is.na(cf)] <- 0
    cf
  }, error = function(e) matrix(0, nrow(X_j), d))

  E_j     <- Y_seg - X_j %*% B_j
  Sigma_j <- crossprod(E_j) / n_j
  Sigma_j <- (Sigma_j + t(Sigma_j)) * 0.5

  ld <- tryCatch(determinant(Sigma_j, logarithm = TRUE)$modulus,
                 error = function(e) -Inf)
  if (!is.finite(ld)) return(Inf)
  (n_j / 2) * ld
}

# ------------------------------------------------------------------
# Internal: sum of OLS Frobenius RSS over all segments.
# ------------------------------------------------------------------

compute_g_sb_ols <- function(Y, p, cps, n) {
  breaks <- c(1L, cps, n + 1L)
  total  <- 0.0
  for (b in seq_len(length(breaks) - 1L)) {
    rows  <- breaks[b]:(breaks[b + 1L] - 1L)
    total <- total + seg_ols_loss_sb(Y[rows, , drop = FALSE], p)
  }
  total
}

# ------------------------------------------------------------------
# Internal: sum of profile likelihood losses over all segments.
# ------------------------------------------------------------------

compute_g_sb_profile <- function(Y, p, cps, n) {
  breaks <- c(1L, cps, n + 1L)
  total  <- 0.0
  for (b in seq_len(length(breaks) - 1L)) {
    rows  <- breaks[b]:(breaks[b + 1L] - 1L)
    total <- total + seg_profile_loss_sb(Y[rows, , drop = FALSE], p)
  }
  total
}

# ------------------------------------------------------------------
# Internal: global sigma2 estimate for sigma_scaled variant.
# Fits no-break OLS on full Y; returns tr(Sigma_global) / d.
# ------------------------------------------------------------------

global_sigma2_sb <- function(Y, p, n, d) {
  rss <- seg_ols_loss_sb(Y, p)
  rss / (max(n - p, 1L) * d)
}

# ------------------------------------------------------------------
# Internal: LASSO objective on one segment.
# Returns the pred.error from var.lasso.brk, or Inf if underdetermined.
# ------------------------------------------------------------------

seg_loss_sb <- function(Y_seg, p, lambda) {
  n_j <- nrow(Y_seg)
  if (n_j <= p) {
    return(Inf)
  }

  tryCatch(
    var.lasso.brk(
      data          = Y_seg,
      weight        = NULL,
      lambda        = lambda,
      p             = p,
      max.iteration = 1000L,
      tol           = 1e-4
    )$pred.error[[1L]],
    error = function(e) Inf
  )
}

# ------------------------------------------------------------------
# Internal: total loss = sum of seg_loss_sb over all segments defined
# by the changepoint vector cps.
# ------------------------------------------------------------------

compute_g_sb <- function(Y, p, lambda, cps, n) {
  breaks <- c(1L, cps, n + 1L)
  total <- 0.0
  for (b in seq_len(length(breaks) - 1L)) {
    rows <- breaks[b]:(breaks[b + 1L] - 1L)
    total <- total + seg_loss_sb(Y[rows, , drop = FALSE], p, lambda)
  }
  total
}

# ------------------------------------------------------------------
# Internal: segment-wise MLE B and Sigma for the pruned model (OLS).
# Returns B (n x d*dp) and Sigma (n x d^2) in the same vec'd layout
# as B / Sigma in hsbvar-bea.R: column-major vec, constant per segment.
# Uses zero pre-sample values for the first p observations, identical
# to the Y_lag convention in hsbvar-bea.R.
# ------------------------------------------------------------------

refit_segments_sb <- function(Y, p, cps, n, d, dp) {
  Y_ext <- rbind(matrix(0, p, d), Y) # (n+p) x d
  Y_lag <- matrix(0.0, n, dp)
  for (k in seq_len(p)) {
    cols <- ((k - 1L) * d + 1L):(k * d)
    Y_lag[, cols] <- Y_ext[(p - k + 1L):(n + p - k), , drop = FALSE]
  }

  breaks <- c(1L, cps, n + 1L)
  B_out <- matrix(NA_real_, n, d * dp)
  Sig_out <- matrix(NA_real_, n, d * d)

  for (b in seq_len(length(breaks) - 1L)) {
    rows <- breaks[b]:(breaks[b + 1L] - 1L)
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
    Sigma_j <- crossprod(E_j) / length(rows)
    Sigma_j <- (Sigma_j + t(Sigma_j)) * 0.5

    B_out[rows, ] <- matrix(as.vector(B_j), nrow = length(rows), ncol = d * dp, byrow = TRUE)
    Sig_out[rows, ] <- matrix(as.vector(Sigma_j), nrow = length(rows), ncol = d * d, byrow = TRUE)
  }

  list(B = B_out, Sigma = Sig_out)
}

# ------------------------------------------------------------------
# Internal: data-driven omega_n selection (largest-gap rule).
# Runs BEA to empty, records the loss increment Delta_i at each
# removal, splits {Delta_i} at the largest gap, returns the minimum
# jump in the "large" (true-break) cluster.
# Identical logic to select_omega_hsbvar in hsbvar-bea.R.
# ------------------------------------------------------------------

select_omega_sb <- function(Y, p, lambda, cps_init, n) {
  cps <- sort(cps_init)
  deltas <- numeric(length(cps))
  idx <- 0L

  while (length(cps) > 0L) {
    m <- length(cps)
    g_cur <- compute_g_sb(Y, p, lambda, cps, n)

    g_rem <- vapply(seq_len(m), function(i) {
      compute_g_sb(Y, p, lambda, cps[-i], n)
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
# Public: BEA pruning of Stage 1 changepoints for SB-VAR
# ------------------------------------------------------------------

#' Backward Elimination Algorithm for SB-VAR (safikhani-supplemental.md, Appendix C)
#'
#' Prunes the over-selected candidate changepoints produced by Stage 1.
#' Three IC variants are available via \code{ic_type}:
#' \describe{
#'   \item{"rss"}{LASSO objective: \eqn{IC = \sum_j (RSS_j + \lambda\|\Phi_j\|_1) + m\,\omega_n}.
#'     Original Safikhani & Shojaie criterion. Not scale-invariant in \eqn{\sigma}.}
#'   \item{"sigma_scaled"}{OLS RSS with scaled penalty:
#'     \eqn{IC = \sum_j \mathrm{tr}(E_j'E_j) + m\,\hat\sigma^2\,\omega_n},
#'     where \eqn{\hat\sigma^2 = \mathrm{tr}(\hat\Sigma_{global})/d}.  Restores natural
#'     scaling when \eqn{\sigma \ne 1}.}
#'   \item{"profile_lik"}{Profile likelihood:
#'     \eqn{IC = \sum_j (n_j/2)\log\det(\hat\Sigma_j) + m\,\omega_n},
#'     \eqn{\hat\Sigma_j = E_j'E_j/n_j}.  Fully scale-invariant; multivariate
#'     generalisation of \eqn{n_j \log(RSS_j/n_j)}.}
#' }
#'
#' @param fit      Output of \code{sbvar()} or \code{first.step.cv()}.
#'                 Must contain \code{$cp} or \code{$brk.points}.
#'                 \code{$cv.final} is used to infer \code{lambda} when not supplied
#'                 (only used for \code{ic_type = "rss"}).
#' @param Y        n x d numeric matrix (same data passed to Stage 1).
#' @param p        VAR lag order.  Inferred from \code{fit$p} when not supplied.
#' @param lambda   LASSO penalty for segment fits.  Only used when
#'                 \code{ic_type = "rss"}.  Inferred from \code{fit$cv.final}.
#' @param omega_n  Penalty per break.  Default: BIC \eqn{d \cdot dp \cdot \log(n)/2}.
#'                 Pass \code{"data"} to use largest-gap selection (rss only).
#' @param ic_type  IC variant: \code{"rss"} (default), \code{"sigma_scaled"},
#'                 or \code{"profile_lik"}.
#'
#' @return A list with:
#' \describe{
#'   \item{cp}{Integer vector of refined changepoint locations.}
#'   \item{B}{n x (d*dp) matrix of segment-wise OLS VAR coefficients.}
#'   \item{Sigma}{n x d^2 matrix of segment-wise OLS covariances.}
#'   \item{ic}{IC value at the final changepoint set.}
#'   \item{omega_n}{Effective penalty per break (after any sigma scaling).}
#'   \item{ic_type}{The IC variant used.}
#' }
sbvar_bea <- function(fit, Y, p = NULL, lambda = NULL, omega_n = NULL,
                      ic_type = c("rss", "sigma_scaled", "profile_lik")) {
  ic_type <- match.arg(ic_type)
  n <- nrow(Y)
  d <- ncol(Y)

  if (is.null(p)) {
    if (!is.null(fit$p)) {
      p <- as.integer(fit$p)
    } else if (!is.null(fit$phi.hat)) {
      p <- ncol(fit$phi.hat) / d
      if (p != round(p) || p < 1L) stop("Cannot infer p from fit; supply p explicitly.")
      p <- as.integer(round(p))
    } else {
      stop("Cannot infer p from fit; supply p explicitly.")
    }
  }

  dp <- d * p

  # ---- lambda: only needed for ic_type = "rss" --------------------------------
  if (ic_type == "rss") {
    if (is.null(lambda)) {
      if (is.null(fit$cv.final)) {
        stop("Cannot infer lambda from fit; supply lambda explicitly.")
      }
      lambda <- fit$cv.final
    }
  }

  # ---- initial changepoints: prefer fit$cp (sbvar), fall back to fit$brk.points
  cps_init_raw <- if (!is.null(fit$cp)) fit$cp else fit$brk.points

  # ---- resolve omega_n -------------------------------------------------------
  if (identical(omega_n, "data")) {
    if (ic_type != "rss") stop('"data" omega selection requires ic_type = "rss".')
    omega_n <- select_omega_sb(Y, p, lambda, sort(cps_init_raw), n)
  } else if (is.null(omega_n)) {
    omega_n <- d * dp * log(n) / 2
  }

  # ---- sigma_scaled: scale omega_n by global noise variance ------------------
  omega_n_eff <- if (ic_type == "sigma_scaled") {
    omega_n * global_sigma2_sb(Y, p, n, d)
  } else {
    omega_n
  }

  # ---- dispatch IC function --------------------------------------------------
  ic_fn <- switch(ic_type,
    rss          = function(cps) compute_g_sb(Y, p, lambda, cps, n) + length(cps) * omega_n_eff,
    sigma_scaled = function(cps) compute_g_sb_ols(Y, p, cps, n)     + length(cps) * omega_n_eff,
    profile_lik  = function(cps) compute_g_sb_profile(Y, p, cps, n) + length(cps) * omega_n_eff
  )

  # ---- initialise ------------------------------------------------------------
  cps <- sort(cps_init_raw)

  if (length(cps) == 0L) {
    seg <- refit_segments_sb(Y, p, cps, n, d, dp)
    return(list(
      cp      = integer(0L),
      B       = seg$B,
      Sigma   = seg$Sigma,
      ic      = ic_fn(cps),
      omega_n = omega_n_eff,
      ic_type = ic_type
    ))
  }

  current_ic <- ic_fn(cps)

  # ---- iterative backward elimination ----------------------------------------
  repeat {
    m <- length(cps)
    if (m == 0L) break

    ic_without <- vapply(seq_len(m), function(i) ic_fn(cps[-i]), numeric(1L))
    best_i     <- which.min(ic_without)

    if (ic_without[best_i] <= current_ic) {
      cps        <- cps[-best_i]
      current_ic <- ic_without[best_i]
    } else {
      break
    }
  }

  # ---- refit pruned segments -------------------------------------------------
  seg <- refit_segments_sb(Y, p, cps, n, d, dp)

  list(
    cp      = cps,
    B       = seg$B,
    Sigma   = seg$Sigma,
    ic      = current_ic,
    omega_n = omega_n_eff,
    ic_type = ic_type
  )
}
