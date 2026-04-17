# sb_stage1_clean.R
# Stage 1: Joint group-LASSO change-point detection for SB-VAR.
#           (Safikhani & Shojaie 2020, Appendix C algorithm.)
#
# Minimises the TV-penalised squared loss
#
#   (1/n) ||Y - Z Theta||_F^2
#     + lambda_1 sum_i ||theta_i||_1          (sparsity)
#     + lambda_2 sum_k ||cumsum(theta)_k||_1  (total variation / break detection)
#
# via Gauss-Seidel block coordinate descent (var.break.fit).  Both penalties
# are collapsed to a single tuning parameter `lambda` as in the existing
# implementation (the "LASSO" method in var.break.fit merges lambda_1 and
# lambda_2 into one soft-threshold).
#
# Interface matches hsbvar-nll.R:
#   - Input/output names align with hsbvar() for drop-in comparison
#   - Theta: n x (d*dp) increment matrix (row t = vec(theta_t), first p rows zero)
#   - B_arr: n x (d*dp) cumulative VAR coefficients (= cumsum of Theta per column)
#   - Sigma_arr: n x d^2 segment-wise OLS covariances (constant within segments)
#   - cp: detected changepoints, same boundary/spacing filter as hsbvar
#
# Key difference from hsbvar: no Psi / Phi_arr (covariance not modelled).
# Sigma_arr is the sample covariance of OLS residuals within each segment,
# not a running precision inverse.

source(file.path("ss", "sb_estimation.R"))

#' Fit SB-VAR via TV-penalised Gauss-Seidel (Safikhani & Shojaie 2020)
#'
#' @param Y        n x d numeric matrix (rows = observations, cols = series)
#' @param p        VAR lag order (default 1)
#' @param lambda   Penalty strength; collapses lambda_1 and lambda_2 from the
#'                 paper into a single tuning parameter (default 0.1)
#' @param max_iter Maximum Gauss-Seidel iterations (default 1000)
#' @param tol      Convergence tolerance: stop when
#'                 max|phi_new - phi_old| < tol (default 1e-4)
#' @param thr      Frobenius-norm threshold for changepoint detection; any
#'                 theta_i with ||theta_i||_F <= thr is treated as zero
#'                 (default 1e-6)
#' @param verbose  Print progress? (default FALSE; note var.break.fit prints
#'                 internally regardless of this flag)
#'
#' @return A list with components:
#'   \describe{
#'     \item{Theta}{n x (d*dp) matrix of coefficient increments.
#'       Row t holds vec(theta_t) column-major; rows 1:p are zero (pre-sample).}
#'     \item{B_arr}{n x (d*dp) matrix of cumulative VAR coefficients.
#'       B_arr[t,] = vec(B_t) where B_t = sum_{i=1}^t theta_i.}
#'     \item{Sigma_arr}{n x d^2 matrix of segment-wise OLS covariances.
#'       Constant within each segment defined by cp.}
#'     \item{cp}{Integer vector of detected changepoints (indices into Y[1:n,]).}
#'     \item{brk.points}{Alias of cp; for compatibility with sbvar_bea.}
#'     \item{phi_hat}{Raw d x (dp*(n-p)) matrix from var.break.fit.}
#'     \item{obj_val}{Mean squared prediction error under the fitted model.}
#'     \item{n_iter}{NA -- var.break.fit does not expose iteration count.}
#'     \item{stop_crit}{"not_tracked".}
#'     \item{cv.final}{The lambda used; stored so sbvar_bea can infer lambda.}
#'     \item{p}{The lag order used; stored so sbvar_bea can infer p.}
#'   }
sbvar <- function(Y,
                  p        = 1,
                  lambda   = 0.1,
                  max_iter = 1000,
                  tol      = 1e-4,
                  thr      = 1e-6,
                  verbose  = FALSE) {
  n     <- nrow(Y)
  d     <- ncol(Y)
  dp    <- d * p
  n_eff <- n - p   # effective observations (responses at times p+1 ... n)

  if (n_eff < 1L) stop("n - p must be >= 1; reduce p or supply more data.")
  if (verbose) message("sbvar: n=", n, "  d=", d, "  p=", p,
                       "  lambda=", lambda, "  n_eff=", n_eff)

  # -----------------------------------------------------------------------
  # 1. Lag matrix  Y_lag[t, ] = c(Y[t-1,], ..., Y[t-p,])
  #    Pre-sample rows (t <= p) padded with zeros -- identical to hsbvar-nll.R
  # -----------------------------------------------------------------------
  Y_ext <- rbind(matrix(0, p, d), Y)   # (n+p) x d
  Y_lag <- matrix(0.0, n, dp)
  for (k in seq_len(p)) {
    cols <- ((k - 1L) * d + 1L):(k * d)
    Y_lag[, cols] <- Y_ext[(p - k + 1L):(n + p - k), , drop = FALSE]
  }

  # -----------------------------------------------------------------------
  # 2. Gauss-Seidel block coordinate descent (var.break.fit, "LASSO" method)
  #    phi_hat: d x (dp * n_eff)  -- block i = d x dp matrix theta_i
  #             for effective time i (original time p + i).
  # -----------------------------------------------------------------------
  raw     <- var.break.fit("LASSO", Y,
                           weight        = NULL,
                           lambda        = lambda,
                           p             = p,
                           max.iteration = max_iter,
                           tol           = tol)
  phi_hat <- raw$phi.hat   # d x (dp * n_eff)

  # -----------------------------------------------------------------------
  # 3. Reshape phi_hat into Theta: n x (d*dp)
  #    Theta[p + i, ] = vec(theta_i)  for i = 1, ..., n_eff
  #    Theta[1:p, ]   = 0             (no increment defined for pre-sample)
  #    Mirrors the n x (d*dp) layout of hsbvar-nll.R's Theta.
  # -----------------------------------------------------------------------
  Theta <- matrix(0.0, n, d * dp)
  for (i in seq_len(n_eff)) {
    Theta[p + i, ] <- as.vector(phi_hat[, ((i - 1L) * dp + 1L):(i * dp)])
  }

  # -----------------------------------------------------------------------
  # 4. Cumulative VAR coefficients  B_arr = cumsum(Theta) per column
  #    B_arr[t, ] = vec(B_t) where B_t = sum_{i=1}^t theta_i
  # -----------------------------------------------------------------------
  B_arr <- apply(Theta, 2L, cumsum)   # n x (d*dp)

  # -----------------------------------------------------------------------
  # 5. Changepoint detection
  #    cp at row t means ||theta_t||_F > thr (nonzero increment at time t).
  #    Since Theta[1:p, ] = 0, only rows p+2 ... n can be nonzero.
  #    Same boundary/spacing filter as hsbvar-nll.R.
  # -----------------------------------------------------------------------
  Theta_tail <- Theta[-1L, , drop = FALSE]          # (n-1) x (d*dp)
  norms      <- sqrt(rowSums(Theta_tail^2))          # length n-1
  cp_raw     <- which(norms > thr) + 1L             # indices into [2, n]

  filter_cp <- function(candidates) {
    x <- candidates[candidates > p + 3L & candidates < n]
    if (length(x) > 1L) {
      too_close <- which(diff(x) <= p + 1L)
      if (length(too_close) > 0L) x <- x[-too_close]
    }
    x
  }
  cp <- filter_cp(cp_raw)

  # -----------------------------------------------------------------------
  # 6. Segment-wise OLS Sigma_arr
  #    OLS residuals within each segment defined by cp give the sample
  #    covariance, stored constant within the segment.
  #    Comparable to hsbvar's Sigma_arr = Phi_t^{-1}, but not time-varying
  #    (sb does not model covariance changes).
  # -----------------------------------------------------------------------
  breaks  <- c(1L, cp, n + 1L)
  Sig_out <- matrix(NA_real_, n, d * d)

  for (b in seq_len(length(breaks) - 1L)) {
    rows <- breaks[b]:(breaks[b + 1L] - 1L)
    Y_j  <- Y[rows, , drop = FALSE]
    X_j  <- Y_lag[rows, , drop = FALSE]

    B_j <- tryCatch(
      t(solve(crossprod(X_j), crossprod(X_j, Y_j))),   # d x dp
      error = function(e) {
        cf <- qr.coef(qr(X_j), Y_j)
        cf[is.na(cf)] <- 0
        t(cf)
      }
    )

    E_j     <- Y_j - X_j %*% t(B_j)
    Sigma_j <- crossprod(E_j) / length(rows)
    Sigma_j <- (Sigma_j + t(Sigma_j)) * 0.5
    Sig_out[rows, ] <- matrix(as.vector(Sigma_j),
                              nrow = length(rows), ncol = d * d, byrow = TRUE)
  }

  # -----------------------------------------------------------------------
  # 7. Prediction RSS as obj_val (mean squared error over effective obs)
  # -----------------------------------------------------------------------
  resid <- matrix(NA_real_, n_eff, d)
  for (i in seq_len(n_eff)) {
    t         <- p + i
    resid[i, ] <- Y[t, ] - matrix(B_arr[t, ], d, dp) %*% Y_lag[t, ]
  }
  obj_val <- sum(resid^2) / n_eff

  if (verbose) message("sbvar: done  |cp|=", length(cp),
                       "  obj_val=", round(obj_val, 6L))

  list(
    Theta      = Theta,
    B_arr      = B_arr,
    Sigma_arr  = Sig_out,
    cp         = cp,
    brk.points = cp,      # alias for sbvar_bea compatibility
    phi_hat    = phi_hat,  # raw d x (dp * n_eff) output from var.break.fit
    obj_val    = obj_val,
    n_iter     = NA_integer_,
    stop_crit  = "not_tracked",
    cv.final   = lambda,   # sbvar_bea uses this to infer lambda
    p          = p         # sbvar_bea uses this to infer p
  )
}
