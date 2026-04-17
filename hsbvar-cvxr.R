# hsbvar-cvxr.R
# H-SBVAR: Heterogeneous Structural Break VAR, solved exactly via CVXR.
#
# Minimises:
#   Q = (1/(2*n_tr)) * L_nat + lambda * sum_{i>=2} sqrt(||Theta_i||_F^2 + c*||Psi_i||_F^2)
#
# where L_nat = sum_t [-logdet(Phi_t) + Y_t^T Phi_t Y_t - 2 Y_t^T v_t + v_t^T Phi_t^{-1} v_t]
# and v_t = G_t X_t,  G_t = cumsum(Theta)[t,,],  Phi_t = cumsum(Psi)[t,,].
#
# This is the multivariate analogue of hsbar-cvxr.R.  Intended as a reference
# solver to validate hsbvar-nll.R (FISTA) on small problems.

library(CVXR)

#' Fit H-SBVAR via CVXR (reference solver)
#'
#' @param Y         n x d numeric matrix
#' @param p         VAR lag order
#' @param lambda    Joint penalty strength
#' @param c_scale   Scale c > 0 weighting precision vs coefficient breaks (default 1)
#' @param keep_rows Integer vector of row indices to include (defaults to all)
#' @param solver    CVXR solver (default "CLARABEL")
#' @param thr       Frobenius-norm threshold for changepoint detection (default 1e-3)
#' @param verbose   Print solver output?
#'
#' @return List: Theta, Psi, Phi_arr, Sigma_arr, B_arr,
#'               cp, cp_theta, cp_psi, status, obj_val
hsbvar_cvxr <- function(Y,
                        p         = 1,
                        lambda    = 0.1,
                        c_scale   = 1,
                        keep_rows = NULL,
                        solver    = "CLARABEL",
                        thr       = 1e-3,
                        verbose   = FALSE) {
  n   <- nrow(Y)
  d   <- ncol(Y)
  dp  <- d * p
  d2  <- d * d
  ddp <- d * dp

  if (is.null(keep_rows)) keep_rows <- seq_len(n)
  n_tr <- length(keep_rows)

  # -----------------------------------------------------------------------
  # 1. Lagged regressor matrix  Y_lag[t, ] = c(Y[t-1,], ..., Y[t-p,])
  # -----------------------------------------------------------------------
  Y_ext <- rbind(matrix(0, p, d), Y)
  Y_lag <- matrix(0.0, n, dp)
  for (k in seq_len(p)) {
    cols <- ((k - 1L) * d + 1L):(k * d)
    Y_lag[, cols] <- Y_ext[(p - k + 1L):(n + p - k), , drop = FALSE]
  }
  Y_lag_tr <- Y_lag[keep_rows, , drop = FALSE]  # n_tr x dp
  Y_tr     <- Y[keep_rows, , drop = FALSE]       # n_tr x d

  # -----------------------------------------------------------------------
  # 2. Cumsum selector: l_sub[r, i] = 1 if i <= keep_rows[r]
  # -----------------------------------------------------------------------
  l_sub <- base::outer(keep_rows, seq_len(n), ">=") + 0L  # n_tr x n

  # -----------------------------------------------------------------------
  # 3. CVXR decision variables
  # -----------------------------------------------------------------------
  Theta_var <- Variable(c(n, ddp))  # row i = vec(Theta_i), column-major
  Psi_var   <- Variable(c(n, d2))   # row i = vec(Psi_i),   column-major

  # -----------------------------------------------------------------------
  # 4. Affine cumulative-sum expressions at training rows
  # -----------------------------------------------------------------------
  Gamma_tr <- l_sub %*% Theta_var  # n_tr x ddp  (row t = vec(G_t))
  Phi_tr   <- l_sub %*% Psi_var    # n_tr x d2   (row t = vec(Phi_t))

  # Helper: build a d x d CVXR matrix from a 1 x d^2 row expression.
  # Entries stored column-major: Phi_t[r, cc] = row_expr[1, (cc-1)*d + r].
  make_mat_expr <- function(row_expr) {
    bmat(lapply(seq_len(d), function(r) {
      lapply(seq_len(d), function(cc) {
        row_expr[1L, (cc - 1L) * d + r, drop = FALSE]
      })
    }))
  }

  # -----------------------------------------------------------------------
  # 5. NLL contributions: one term per training time step
  # -----------------------------------------------------------------------
  loss_terms <- lapply(seq_len(n_tr), function(ti) {
    Y_t   <- Y_tr[ti, ]                           # d-vector (fixed)
    X_t   <- Y_lag_tr[ti, ]                       # dp-vector (fixed)
    yyvec <- as.vector(outer(Y_t, Y_t))            # d^2-vector, column-major

    # Phi_t as d x d CVXR matrix expression
    Phi_t_expr <- make_mat_expr(Phi_tr[ti, , drop = FALSE])

    # v_t = G_t X_t via the Kronecker identity:
    #   v_t = (X_t^T ⊗ I_d) vec(G_t)
    M_t <- kronecker(t(X_t), diag(d))              # d x ddp (fixed)
    v_t <- M_t %*% t(Gamma_tr[ti, , drop = FALSE]) # d x 1 CVXR expr

    # -log det(Phi_t)  [convex in Phi_t]
    neg_logdet <- -log_det(Phi_t_expr)

    # Y_t^T Phi_t Y_t = vec(Phi_t)^T (Y_t ⊗ Y_t)  [linear in vec(Phi_t)]
    quad_phi <- sum_entries(
      Phi_tr[ti, , drop = FALSE] * matrix(yyvec, nrow = 1L, ncol = d2)
    )

    # -2 Y_t^T v_t  [linear in Theta]
    cross <- -2.0 * sum_entries(matrix(Y_t, ncol = 1L) * v_t)

    # v_t^T Phi_t^{-1} v_t  [matrix_frac: convex in (v_t, Phi_t)]
    mf <- matrix_frac(v_t, Phi_t_expr)

    neg_logdet + quad_phi + cross + mf
  })

  total_loss <- Reduce("+", loss_terms) / (2.0 * n_tr)

  # -----------------------------------------------------------------------
  # 6. Joint co-location penalty
  #    sqrt(||Theta_i||_F^2 + c * ||Psi_i||_F^2)  for i = 2, ..., n
  # -----------------------------------------------------------------------
  sq_c      <- sqrt(c_scale)
  pen_joint <- lambda * Reduce("+", lapply(seq(2L, n), function(i) {
    cvxr_norm(
      hstack(
        Theta_var[i, , drop = FALSE],
        sq_c * Psi_var[i, , drop = FALSE]
      ), "F"
    )
  }))

  # -----------------------------------------------------------------------
  # 7. Solve
  # -----------------------------------------------------------------------
  problem <- Problem(Minimize(total_loss + pen_joint))
  psolve(problem, solver = solver, verbose = verbose)

  sol_status <- status(problem)
  sol_obj    <- value(problem)
  if (!sol_status %in% c("optimal", "optimal_inaccurate")) {
    warning("H-SBVAR CVXR solver status: ", sol_status)
  }

  # -----------------------------------------------------------------------
  # 8. Extract results
  # -----------------------------------------------------------------------
  Theta_hat  <- value(Theta_var)
  Psi_hat    <- value(Psi_var)

  # -----------------------------------------------------------------------
  # 9. Recover original parameters
  #    B_t = Phi_t^{-1} G_t,   Sigma_t = Phi_t^{-1}
  # -----------------------------------------------------------------------
  Gamma_full <- apply(Theta_hat, 2L, cumsum)  # n x ddp
  Phi_full   <- apply(Psi_hat,   2L, cumsum)  # n x d2
  B_arr      <- matrix(NA_real_, n, ddp)
  Sigma_arr  <- matrix(NA_real_, n, d2)

  for (t in seq_len(n)) {
    Phi_t <- matrix(Phi_full[t, ], d, d)
    Phi_t <- (Phi_t + t(Phi_t)) * 0.5
    ch    <- tryCatch(chol(Phi_t), error = function(e) NULL)
    if (!is.null(ch)) {
      Pi             <- chol2inv(ch)
      Sigma_arr[t, ] <- as.vector(Pi)
      G_t            <- matrix(Gamma_full[t, ], d, dp)
      B_arr[t, ]     <- as.vector(Pi %*% G_t)
    }
  }

  # -----------------------------------------------------------------------
  # 10. Detect changepoints
  # -----------------------------------------------------------------------
  Theta_tail <- Theta_hat[-1L, , drop = FALSE]
  Psi_tail   <- Psi_hat[-1L,   , drop = FALSE]
  norm_joint <- sqrt(rowSums(Theta_tail^2) + c_scale * rowSums(Psi_tail^2))
  norm_theta <- sqrt(rowSums(Theta_tail^2))
  norm_psi   <- sqrt(rowSums(Psi_tail^2))

  cp       <- which(norm_joint > thr) + 1L
  cp_theta <- which(norm_theta > thr) + 1L
  cp_psi   <- which(norm_psi   > thr) + 1L

  filter_cp <- function(candidates) {
    x <- candidates[candidates > p + 3L & candidates < n]
    if (length(x) > 1L) {
      too_close <- which(diff(x) <= p + 1L)
      if (length(too_close) > 0L) x <- x[-too_close]
    }
    x
  }
  cp       <- filter_cp(cp)
  cp_theta <- filter_cp(cp_theta)
  cp_psi   <- filter_cp(cp_psi)

  list(
    Theta     = Theta_hat,
    Psi       = Psi_hat,
    Phi_arr   = Phi_full,
    Sigma_arr = Sigma_arr,
    B_arr     = B_arr,
    cp        = cp,
    cp_theta  = cp_theta,
    cp_psi    = cp_psi,
    status    = sol_status,
    obj_val   = sol_obj
  )
}
