# hsbvar-nll.R
# H-SBVAR fitted by FISTA (accelerated proximal gradient) with backtracking.
#
# Minimises Q = f + g where:
#   f = (1 / (2 * n_tr)) * L_nat    smooth NLL loss (eq:var-loss of h-sbvar.tex)
#   g = lambda * sum_{i>=2} sqrt(||Theta_i||_F^2 + c * ||Psi_i||_F^2)
#
# Natural parametrisation (eq:var-nat):
#   Gamma_t = cumsum(Theta, dim=1)[t,,]   d x dp   precision-scaled VAR coefficients
#   Phi_t   = cumsum(Psi,   dim=1)[t,,]   d x d    precision matrix  (must be PD)
#
# Theta and Psi are stored as n x (d*dp) and n x d^2 matrices respectively,
# where row i holds vec(Theta_i) / vec(Psi_i) in column-major order.
# This mirrors the n x p / n x 1 layout of hsbar.R.
#
# See h-sbvar.tex Sec. 1 for the full mathematical development.

#' Fit H-SBVAR via FISTA with backtracking line search
#'
#' @param Y         n x d numeric matrix (rows = observations, cols = series)
#' @param p         VAR lag order
#' @param lambda    Joint penalty strength
#' @param c_scale   Scale parameter c > 0 balancing coefficient / precision breaks
#'                  (default 1)
#' @param keep_rows Integer vector of row indices to include in the likelihood.
#'                  Defaults to all rows.  Rows t <= p have zero regressors and
#'                  contribute the marginal (pre-sample) log-likelihood.
#' @param alpha0    Initial step size for backtracking (default 1)
#' @param beta      Backtracking shrinkage factor in (0, 1) (default 0.5)
#' @param max_iter  Maximum FISTA iterations (default 1000)
#' @param tol       Convergence tolerance on proximal gradient mapping norm
#'                  (default 1e-6)
#' @param thr       Frobenius-norm threshold for changepoint detection (default 1e-6)
#' @param restart   Apply gradient-based momentum restart? (default TRUE)
#' @param eps_tol   Stop if |dQ| < eps_tol * (1 + |Q|) (default 1e-10)
#' @param verbose   Print iteration log? (default FALSE)
#' @param init_Theta  n x (d*d*p) warm-start for Theta (optional)
#' @param init_Psi    n x d^2    warm-start for Psi   (optional)
#'
#' @return List with components:
#'   \describe{
#'     \item{Theta}{n x (d*d*p) matrix of natural coefficient increments}
#'     \item{Psi}{n x d^2 matrix of precision increments}
#'     \item{Phi_arr}{n x d^2 matrix of cumulative precisions (vec'd, col-major)}
#'     \item{Sigma_arr}{n x d^2 matrix of cumulative covariances}
#'     \item{B_arr}{n x (d*d*p) matrix of cumulative VAR coefficients}
#'     \item{cp, cp_theta, cp_psi}{detected changepoint indices}
#'     \item{obj_val, n_iter, stop_crit}{}
#'   }
hsbvar <- function(Y,
                   p = 1,
                   lambda = 0.1,
                   c_scale = 1,
                   keep_rows = NULL,
                   alpha0 = 1,
                   beta = 0.5,
                   max_iter = 1000,
                   tol = 1e-6,
                   thr = 1e-6,
                   restart = TRUE,
                   eps_tol = 1e-10,
                   verbose = FALSE,
                   init_Theta = NULL,
                   init_Psi = NULL) {
  n <- nrow(Y)
  d <- ncol(Y)
  dp <- d * p
  d2 <- d * d
  ddp <- d * dp

  if (is.null(keep_rows)) keep_rows <- seq_len(n)
  n_tr <- length(keep_rows)

  # -----------------------------------------------------------------------
  # 1. Lagged regressor matrix  Y_lag[t, ] = c(Y[t-1,], ..., Y[t-p,])
  #    Pre-sample rows (t <= p) are padded with zeros
  # -----------------------------------------------------------------------
  Y_ext <- rbind(matrix(0, p, d), Y) # (n+p) x d
  Y_lag <- matrix(0.0, n, dp)
  for (k in seq_len(p)) {
    cols <- ((k - 1L) * d + 1L):(k * d)
    Y_lag[, cols] <- Y_ext[(p - k + 1L):(n + p - k), , drop = FALSE]
  }
  Y_lag_tr <- Y_lag[keep_rows, , drop = FALSE] # n_tr x dp
  Y_tr <- Y[keep_rows, , drop = FALSE] # n_tr x d

  # -----------------------------------------------------------------------
  # 2. Cumulative-sum helpers
  #    fwd_cumsum_mat: cumsum of each column of an n x k matrix, return keep_rows
  #    bwd_cumsum_mat: scatter n_tr x k to n x k, then reverse cumsum per column
  #      This replaces  t(L_sub) %*% C_tr  from the dense formulation.
  # -----------------------------------------------------------------------
  fwd_cumsum_mat <- function(A) {
    apply(A, 2L, cumsum)[keep_rows, , drop = FALSE]
  }

  bwd_cumsum_mat <- function(C_tr) {
    out <- matrix(0.0, n, ncol(C_tr))
    out[keep_rows, ] <- C_tr
    apply(out, 2L, function(x) rev(cumsum(rev(x))))
  }

  # -----------------------------------------------------------------------
  # 3. Forward pass: cumulative Gamma_t, Phi_t and derived quantities
  #    For each training time ti:
  #      Phi_t   = Psi[1,,] + ... + Psi[t,,]   (symmetrised, PD checked via Cholesky)
  #      Phi_inv = Phi_t^{-1}                   (via chol2inv; NULL if not PD)
  #      v_t     = G_t X_t  in R^d
  # -----------------------------------------------------------------------
  compute_gp <- function(Theta, Psi) {
    Gamma_tr <- fwd_cumsum_mat(Theta) # n_tr x (d*dp)
    Phi_tr <- fwd_cumsum_mat(Psi) # n_tr x d^2

    Phi_inv <- vector("list", n_tr)
    logdet_v <- numeric(n_tr)
    v_mat <- matrix(0.0, n_tr, d)

    for (ti in seq_len(n_tr)) {
      Phi_t <- matrix(Phi_tr[ti, ], d, d)
      Phi_t <- (Phi_t + t(Phi_t)) * 0.5 # enforce symmetry
      ch <- tryCatch(chol(Phi_t), error = function(e) NULL)
      if (is.null(ch)) {
        logdet_v[ti] <- -Inf
        Phi_inv[[ti]] <- NULL
      } else {
        logdet_v[ti] <- 2.0 * sum(log(diag(ch)))
        Phi_inv[[ti]] <- chol2inv(ch)
      }
      G_t <- matrix(Gamma_tr[ti, ], d, dp)
      v_mat[ti, ] <- G_t %*% Y_lag_tr[ti, ] # d-vector  G_t X_t
    }

    list(
      Gamma_tr = Gamma_tr, Phi_tr = Phi_tr,
      Phi_inv = Phi_inv, logdet_v = logdet_v, v_mat = v_mat
    )
  }

  # -----------------------------------------------------------------------
  # 4. Smooth loss  f = (1 / (2 * n_tr)) * L_nat
  #
  #    Per-time term (eq:var-loss):
  #      -log det(Phi_t) + Y_t^T Phi_t Y_t - 2 Y_t^T v_t + v_t^T Phi_t^{-1} v_t
  #    where v_t = G_t X_t.
  # -----------------------------------------------------------------------
  compute_f <- function(gp) {
    if (!all(is.finite(gp$logdet_v))) {
      return(Inf)
    }
    total <- 0.0
    for (ti in seq_len(n_tr)) {
      Phi_t <- matrix(gp$Phi_tr[ti, ], d, d)
      Phi_t <- (Phi_t + t(Phi_t)) * 0.5
      Y_t <- Y_tr[ti, ]
      v_t <- gp$v_mat[ti, ]
      w_t <- drop(gp$Phi_inv[[ti]] %*% v_t) # Phi_t^{-1} v_t  (d-vector)
      total <- total + (
        -gp$logdet_v[ti]
        + sum(Y_t * (Phi_t %*% Y_t))
          - 2.0 * sum(Y_t * v_t)
          + sum(v_t * w_t)
      )
    }
    total / (2.0 * n_tr)
  }

  # -----------------------------------------------------------------------
  # 5. Penalty  g = lambda * sum_{i>=2} sqrt(||Theta_i||_F^2 + c*||Psi_i||_F^2)
  # -----------------------------------------------------------------------
  compute_g_pen <- function(Theta, Psi) {
    if (n == 1L) {
      return(0.0)
    }
    norms <- sqrt(rowSums(Theta[-1L, , drop = FALSE]^2) +
      c_scale * rowSums(Psi[-1L, , drop = FALSE]^2))
    lambda * sum(norms)
  }

  # -----------------------------------------------------------------------
  # 6. Gradients of f  (eq:grad-theta / eq:grad-psi of h-sbvar.tex)
  #
  #    Per-time contributions (reverse-cumsummed to get gradient per increment):
  #      grad_Theta: (w_t - Y_t) X_t^T   d x dp
  #      grad_Psi:   -1/2 Phi_t^{-1} + 1/2 Y_t Y_t^T - 1/2 w_t w_t^T   d x d
  #    where w_t = Phi_t^{-1} G_t X_t.
  # -----------------------------------------------------------------------
  compute_grad <- function(gp) {
    cth <- matrix(0.0, n_tr, ddp) # per-time d x dp contributions (vec'd)
    cps <- matrix(0.0, n_tr, d2) # per-time d x d  contributions (vec'd)

    for (ti in seq_len(n_tr)) {
      Y_t <- Y_tr[ti, ]
      X_t <- Y_lag_tr[ti, ]
      v_t <- gp$v_mat[ti, ]
      w_t <- drop(gp$Phi_inv[[ti]] %*% v_t) # Phi_t^{-1} G_t X_t  (d-vector)
      cth[ti, ] <- as.vector(outer(w_t - Y_t, X_t))
      cps[ti, ] <- as.vector(
        -0.5 * gp$Phi_inv[[ti]]
          + 0.5 * outer(Y_t, Y_t)
          - 0.5 * outer(w_t, w_t)
      )
    }

    list(
      grad_Theta = bwd_cumsum_mat(cth) / n_tr, # n x (d*dp)
      grad_Psi   = bwd_cumsum_mat(cps) / n_tr # n x d^2
    )
  }

  # -----------------------------------------------------------------------
  # 7. Proximal operator: joint Frobenius-norm soft-threshold for i >= 2
  #    prox_{alpha*g}(A, B)_i = (1 - alpha*lambda / ||[A_i; sqrt(c)*B_i]||_F)_+ [A_i; B_i]
  # -----------------------------------------------------------------------
  prox_g <- function(a_Th, a_Ps, alpha) {
    if (n == 1L) {
      return(list(Th = a_Th, Ps = a_Ps))
    }
    tail_Th <- a_Th[-1L, , drop = FALSE]
    tail_Ps <- a_Ps[-1L, , drop = FALSE]
    norms <- sqrt(rowSums(tail_Th^2) + c_scale * rowSums(tail_Ps^2))
    shrink <- ifelse(norms > 0, pmax(0.0, 1.0 - alpha * lambda / norms), 0.0)
    out_Th <- a_Th
    out_Ps <- a_Ps
    out_Th[-1L, ] <- shrink * tail_Th
    out_Ps[-1L, ] <- shrink * tail_Ps
    list(Th = out_Th, Ps = out_Ps)
  }

  # -----------------------------------------------------------------------
  # 8. Initialise
  #    w = true iterate; m = extrapolated point (start equal); s = momentum
  #    Psi[1,] = vec(Omega_0) where Omega_0 = sample precision of Y_tr.
  # -----------------------------------------------------------------------
  if (!is.null(init_Theta) && !is.null(init_Psi)) {
    Theta <- init_Theta
    Psi <- init_Psi
  } else {
    Theta <- matrix(0.0, n, ddp)
    S0 <- cov(Y_tr)
    if (rcond(S0) < 1e-10) S0 <- S0 + 1e-6 * diag(d)
    Psi <- matrix(0.0, n, d2)
    Psi[1L, ] <- as.vector(solve(S0))
  }

  m_Th <- Theta
  m_Ps <- Psi
  s <- 1.0
  n_iter <- max_iter
  obj_prev <- Inf
  stop_crit <- "max_iter"

  # -----------------------------------------------------------------------
  # 9. FISTA main loop
  # -----------------------------------------------------------------------
  for (iter in seq_len(max_iter)) {
    gp_m <- compute_gp(m_Th, m_Ps)
    f_m <- compute_f(gp_m)
    grad_m <- compute_grad(gp_m)

    # --- Backtracking from extrapolated point m ---
    alpha <- alpha0
    repeat {
      Th_tent <- m_Th - alpha * grad_m$grad_Theta
      Ps_tent <- m_Ps - alpha * grad_m$grad_Psi
      prx <- prox_g(Th_tent, Ps_tent, alpha)
      Th_new <- prx$Th
      Ps_new <- prx$Ps

      gp_new <- compute_gp(Th_new, Ps_new)
      if (all(is.finite(gp_new$logdet_v))) {
        f_new <- compute_f(gp_new)
        d_Th <- Th_new - m_Th
        d_Ps <- Ps_new - m_Ps
        inner <- sum(grad_m$grad_Theta * d_Th) + sum(grad_m$grad_Psi * d_Ps)
        sq_nrm <- sum(d_Th^2) + sum(d_Ps^2)
        if (f_new <= f_m + inner + sq_nrm / (2.0 * alpha)) break
      }

      alpha <- beta * alpha
      if (alpha < 1e-16) {
        warning("FISTA: backtracking step size collapsed at iteration ", iter)
        break
      }
    }

    # --- Gradient-based momentum restart (O'Donoghue-Candes) ---
    # Reset when the new iterate overshoots the extrapolation point.
    do_restart <- restart && (
      sum((Th_new - Theta) * (m_Th - Th_new)) +
        sum((Ps_new - Psi) * (m_Ps - Ps_new)) > 0
    )

    # --- Nesterov momentum update ---
    s_new <- (1.0 + sqrt(1.0 + 4.0 * s^2)) / 2.0

    if (do_restart) {
      m_Th <- Th_new
      m_Ps <- Ps_new
      s_new <- 1.0
    } else {
      coef <- (s - 1.0) / s_new
      m_Th <- Th_new + coef * (Th_new - Theta)
      m_Ps <- Ps_new + coef * (Ps_new - Psi)
    }

    Theta <- Th_new
    Psi <- Ps_new
    s <- s_new

    # --- Proximal gradient norm at the true iterate ---
    # G_alpha(w) = (w - prox_{alpha*g}(w - alpha * grad_f(w))) / alpha
    # Zero at the optimum; used instead of ||Th_new - Theta_old|| / alpha
    # to avoid momentum contamination.
    gp_t <- gp_new
    grad_t <- compute_grad(gp_t)
    prx_t <- prox_g(
      Theta - alpha * grad_t$grad_Theta,
      Psi - alpha * grad_t$grad_Psi,
      alpha
    )
    pg_norm <- sqrt(
      sum((prx_t$Th - Theta)^2) + sum((prx_t$Ps - Psi)^2)
    ) / alpha

    f_cur <- compute_f(gp_t)
    g_pen <- compute_g_pen(Theta, Psi)
    obj_cur <- f_cur + g_pen

    if (verbose) {
      rst <- if (do_restart) " R" else "  "
      message(sprintf(
        "iter %4d%s| f=%.6f | g=%.6f | Q=%.6f | pg=%.2e | a=%.2e",
        iter, rst, f_cur, g_pen, obj_cur, pg_norm, alpha
      ))
    }

    if (pg_norm < tol) {
      stop_crit <- "pg_norm"
      n_iter <- iter
      break
    }
    if (abs(obj_cur - obj_prev) < eps_tol * (1.0 + abs(obj_cur))) {
      stop_crit <- "obj_change"
      n_iter <- iter
      break
    }
    obj_prev <- obj_cur
  }

  # -----------------------------------------------------------------------
  # 10. Recover original parameters from natural parametrisation
  #     B_t   = Phi_t^{-1} Gamma_t   (d x dp VAR coefficients)
  #     Sigma_t = Phi_t^{-1}          (d x d covariance)
  # -----------------------------------------------------------------------
  Gamma_full <- apply(Theta, 2L, cumsum) # n x (d*dp)
  Phi_full <- apply(Psi, 2L, cumsum) # n x d^2
  B_arr <- matrix(NA_real_, n, ddp)
  Sigma_arr <- matrix(NA_real_, n, d2)

  for (t in seq_len(n)) {
    Phi_t <- matrix(Phi_full[t, ], d, d)
    Phi_t <- (Phi_t + t(Phi_t)) * 0.5
    ch <- tryCatch(chol(Phi_t), error = function(e) NULL)
    if (!is.null(ch)) {
      Phi_inv_t <- chol2inv(ch)
      Sigma_arr[t, ] <- as.vector(Phi_inv_t)
      G_t <- matrix(Gamma_full[t, ], d, dp)
      B_arr[t, ] <- as.vector(Phi_inv_t %*% G_t)
    }
  }

  # -----------------------------------------------------------------------
  # 11. Detect changepoints
  # -----------------------------------------------------------------------
  Theta_tail <- Theta[-1L, , drop = FALSE] # (n-1) x (d*dp)
  Psi_tail <- Psi[-1L, , drop = FALSE] # (n-1) x d^2

  norm_joint <- sqrt(rowSums(Theta_tail^2) + c_scale * rowSums(Psi_tail^2))
  norm_theta <- sqrt(rowSums(Theta_tail^2))
  norm_psi <- sqrt(rowSums(Psi_tail^2))

  cp <- which(norm_joint > thr) + 1L
  cp_theta <- which(norm_theta > thr) + 1L
  cp_psi <- which(norm_psi > thr) + 1L

  # -----------------------------------------------------------------------
  # 12. Post-processing: boundary trim and minimum-spacing filter
  #     - Drop candidates within p+3 of the start or at the final index n
  #     - Of any two adjacent candidates with gap <= p+1, drop the earlier
  # -----------------------------------------------------------------------
  filter_cp <- function(candidates) {
    x <- candidates[candidates > p + 3L & candidates < n]
    if (length(x) > 1L) {
      too_close <- which(diff(x) <= p + 1L)
      if (length(too_close) > 0L) x <- x[-too_close]
    }
    x
  }
  cp <- filter_cp(cp)
  cp_theta <- filter_cp(cp_theta)
  cp_psi <- filter_cp(cp_psi)

  list(
    Theta     = Theta,
    Psi       = Psi,
    Phi_arr   = Phi_full,
    Sigma_arr = Sigma_arr,
    B_arr     = B_arr,
    cp        = cp,
    cp_theta  = cp_theta,
    cp_psi    = cp_psi,
    obj_val   = f_cur + compute_g_pen(Theta, Psi),
    n_iter    = n_iter,
    stop_crit = stop_crit
  )
}
