# cv_sbdetect_stage1.R
# Interpolation-based cross-validation for Safikhani-Shojaie Stage 1 over a
# lambda path.  Interface mirrors cv_hsbvar_nll() in cv-hsbvar-nll.R.
#
# Methodology follows first.step.cv.new() in sb_stage1.R:
#   1. Choose equally spaced validation points with spacing val_spacing > p.
#   2. Remove validation rows (and optionally the p downstream rows) to form
#      the training set Y_train.
#   3. For each lambda fit var.break.fit() on Y_train (warm-started from the
#      previous lambda on the decreasing path).
#   4. Reconstruct cumulative coefficient phi_t at each validation point from
#      the increment matrix, using the original Y for lag computation.
#   5. MSPE(lambda) = mean over val points of mean((Y_t - phi_t %*% X_t)^2).
#
# @param Y                 n x d numeric matrix.
# @param p                 VAR lag order.
# @param lambda            Numeric vector of lambda candidates.
#                          Default: 10-point log grid on [1e-3, 1].
# @param val_spacing       Spacing between validation points. Must be > p.
#                          Default 10.
# @param drop_poisoned_rows If TRUE, also exclude the p rows following each
#                          validation point from training. Default FALSE.
# @param lambda_rule       "min" or "1se". Default "min".
# @param verbose           Print per-lambda progress? Default TRUE.
# @param method            Estimation method forwarded to var.break.fit.
#                          Default "LASSO".
# @param max_iter          Maximum coordinate-descent iterations. Default 1000.
# @param tol               Convergence tolerance. Default 1e-4.
#
# @return List with the same structure as cv_hsbvar_nll():
#   cv_table   — data.frame(lambda, mspe, mspe_se), one row per lambda
#   best       — single-row data.frame for the selected lambda
#   val_points — integer vector of validation row indices into Y
#   keep_rows  — integer vector of training row indices into Y

source(file.path("ss", "sb_stage1.R"))

cv_sbdetect_stage1 <- function(Y,
                               p = 1L,
                               lambda = NULL,
                               val_spacing = 10L,
                               drop_poisoned_rows = FALSE,
                               lambda_rule = c("min", "1se"),
                               verbose = TRUE,
                               method = "LASSO",
                               max_iter = 1000,
                               tol = 1e-4) {
  n <- nrow(Y)
  d <- ncol(Y)

  # ---- lambda path (decreasing for warm starts) ----------------------------
  lambda_vec <- if (is.null(lambda)) 10^seq(-3, 0, length.out = 10) else lambda
  lambda_vec <- sort(lambda_vec, decreasing = TRUE)
  n_lambda <- length(lambda_vec)

  # ---- validation points ---------------------------------------------------
  if (val_spacing <= p) {
    stop(sprintf("val_spacing (%d) must be > p (%d).", val_spacing, p))
  }
  val_points <- seq(p + val_spacing, n - p, by = val_spacing)
  k_val <- length(val_points)
  if (k_val < 2L) {
    stop("Fewer than 2 validation points. Reduce val_spacing or increase n.")
  }

  # ---- excluded / training rows --------------------------------------------
  excl_rows <- sort(unique(as.integer(unlist(lapply(val_points, function(t) {
    if (drop_poisoned_rows) t:(t + p) else t
  })))))
  keep_rows <- setdiff(seq_len(n), excl_rows)

  if (length(keep_rows) < 2L * p) {
    stop("Too few training rows after exclusion. Reduce val_spacing.")
  }

  Y_train <- Y[keep_rows, , drop = FALSE]
  T_tr <- nrow(Y_train)

  # For each val point t, the phi_full_all index is the number of training rows
  # strictly before t (i.e. length of keep_rows < t), minus p (because
  # phi_full_all[[i]] holds the cumulative phi through training response i,
  # indexed from 1 for the first response row at training row p+1).
  # This generalises first.step.cv.new's formula (cv.index[j]-1-p-j+1) to
  # arbitrary exclusion sets including drop_poisoned_rows = TRUE.
  phi_idx <- vapply(val_points, function(t) sum(keep_rows < t) - p, integer(1L))

  lambda_rule <- match.arg(lambda_rule, c("min", "1se"))

  # ---- main CV loop --------------------------------------------------------
  sq_err <- matrix(NA_real_, n_lambda, k_val)
  mspe <- numeric(n_lambda)
  cp_list <- vector("list", n_lambda)
  prev_phi <- NULL # warm-start increments (initial.phi in var.break.fit)

  for (j in seq_len(n_lambda)) {
    ln <- lambda_vec[j]

    fit <- tryCatch(
      var.break.fit(method, Y_train,
        weight        = NULL,
        lambda        = ln,
        p             = p,
        max.iteration = max_iter,
        tol           = tol,
        initial.phi   = prev_phi
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      mspe[j] <- NA_real_
      if (verbose) message(sprintf("[%d/%d]  lambda=%.4g  FAILED", j, n_lambda, ln))
      next
    }

    prev_phi <- fit$phi.hat
    phi_hat_full <- fit$phi.hat # d x (d*p*(T_tr-p))

    # -- reconstruct cumulative phi at every training response step ----------
    phi_full_all <- vector("list", T_tr - p)
    phi_cum <- matrix(0.0, d, d * p)
    for (jj in seq(p + 1L, T_tr)) {
      s <- (jj - p - 1L) * d * p + 1L
      e <- (jj - p) * d * p
      phi_cum <- phi_cum + phi_hat_full[, s:e, drop = FALSE]
      phi_full_all[[jj - p]] <- phi_cum
    }

    # -- predict at each validation point using original Y for lags ----------
    for (jv in seq_len(k_val)) {
      idx <- phi_idx[jv]
      if (idx < 1L || idx > length(phi_full_all)) next
      phi_t <- phi_full_all[[idx]]
      if (is.null(phi_t)) next

      # pred(Y_transposed, phi, p, T, k, h) returns d-vector prediction
      Y_hat <- as.vector(pred(t(Y), phi_t, p, val_points[jv] - 1L, d, 1L))
      sq_err[j, jv] <- mean((Y[val_points[jv], ] - Y_hat)^2)
    }
    mspe[j] <- mean(sq_err[j, ], na.rm = TRUE)

    # -- extract break points (training response space -> original Y rows) ---
    n_tr <- T_tr - p
    m_hat <- 0L
    brk <- integer(n_tr)
    for (i in seq(2L, n_tr)) {
      s <- (i - 1L) * d * p + 1L
      e <- i * d * p
      if (sum(phi_hat_full[, s:e]^2) != 0) {
        m_hat <- m_hat + 1L
        brk[m_hat] <- i
      }
    }
    brk <- brk[seq_len(m_hat)]
    brk <- brk[brk > (p + 3L) & brk < n_tr]
    if (length(brk) > 1L) {
      too_close <- which(diff(brk) <= (p + 1L))
      if (length(too_close)) brk <- brk[-too_close]
    }
    # translate training response indices to original Y rows
    cp_list[[j]] <- if (length(brk)) keep_rows[p + brk] else integer(0L)

    if (verbose) {
      message(sprintf(
        "[%d/%d]  lambda=%.4g  MSPE=%.5g  ncp=%d",
        j, n_lambda, ln, mspe[j], length(cp_list[[j]])
      ))
    }
  }

  # ---- assemble output -----------------------------------------------------
  mspe_se <- apply(sq_err, 1L, sd, na.rm = TRUE) / sqrt(k_val)
  cv_table <- data.frame(lambda = lambda_vec, mspe = mspe, mspe_se = mspe_se)

  min_idx <- which.min(mspe)

  best_idx <- if (lambda_rule == "min") {
    min_idx
  } else {
    # segment-based variance rule (mirrors cv_hsbvar_nll "1se" branch)
    cp_min <- cp_list[[min_idx]]
    seg_id <- vapply(val_points, function(t) sum(t > cp_min) + 1L, integer(1L))
    seg_vars <- tapply(sq_err[min_idx, ], seg_id, stats::var)
    pt_var <- as.numeric(seg_vars[as.character(seg_id)])
    pt_var[is.na(pt_var)] <- stats::var(sq_err[min_idx, ], na.rm = TRUE)
    cv_var <- (1.0 / (d * k_val)) * sqrt(sum(pt_var, na.rm = TRUE))
    threshold <- mspe[min_idx] + cv_var
    candidates <- which(!is.na(mspe) & seq_len(n_lambda) <= min_idx &
      mspe <= threshold)
    if (length(candidates) == 0L) min_idx else min(candidates)
  }

  list(
    cv_table   = cv_table,
    best       = cv_table[best_idx, , drop = FALSE],
    val_points = val_points,
    keep_rows  = keep_rows
  )
}

#' Plot MSPE along the lambda path
#'
#' @param cv_result Output of cv_sbdetect_stage1().
#' @param log_x     Log-scale the lambda axis? Default TRUE.
#' @param ...       Extra arguments passed to plot().
plot_cv_sbdetect_stage1 <- function(cv_result, log_x = TRUE, ...) {
  tbl <- cv_result$cv_table
  best <- cv_result$best

  x <- tbl$lambda
  bx <- best$lambda
  if (log_x) {
    x <- log10(x)
    bx <- log10(bx)
  }

  xlab <- if (log_x) expression(log[10](lambda)) else expression(lambda)
  plot(x, tbl$mspe,
    type = "b", pch = 19, col = "steelblue",
    xlab = xlab, ylab = "CV MSPE",
    main = "Safikhani Stage 1 \u2013 cross-validation MSPE",
    ...
  )
  abline(v = bx, lty = 2, col = "firebrick")
  legend("topright",
    legend = sprintf("Best = %.4g\nMSPE = %.5g", best$lambda, best$mspe),
    bty = "n", text.col = "firebrick"
  )
  invisible(cv_result)
}
