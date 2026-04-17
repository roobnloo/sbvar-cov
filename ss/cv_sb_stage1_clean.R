# cv_sb_stage1_clean.R
# Interpolation-based cross-validation for SB-VAR over a lambda path.
#
# Mirrors cv-hsbvar-nll.R in structure and parameters.  Key methodological
# difference: hsbvar supports keep_rows, so the full Y is passed and only
# training rows contribute to the loss (lag structure stays intact).  sbvar
# (var.break.fit) has no keep_rows mechanism, so we pass the thinned training
# data Y[train_rows, ] directly.  This introduces a small lag misalignment at
# the gaps left by removed validation points; the bias is negligible when
# val_spacing >> p.
#
# Prediction at validation point t uses the carry-forward coefficient:
#   Theta_orig is an n x (d*dp) zero matrix with Theta_train rows scattered
#   into their original time indices.  B_cumsum = cumsum(Theta_orig) gives
#   the piecewise coefficient at every original row, including val points, via
#   the same cumsum-selector logic as cv-hsbvar-nll.R.

source(file.path("ss", "sb_stage1_clean.R"))

#' Interpolation CV for SB-VAR over a lambda path
#'
#' @param Y                 n x d numeric matrix.
#' @param p                 VAR lag order.
#' @param lambda            Numeric vector of lambda candidates.
#'                          Default: 10-point log grid on [1e-3, 1].
#' @param val_spacing       Spacing between validation points. Must be > p.
#'                          Default 10.
#' @param drop_poisoned_rows If TRUE, also exclude the p rows following each
#'                          validation point from the training set. Default FALSE.
#' @param lambda_rule       \code{"min"} or \code{"1se"}. Default \code{"min"}.
#' @param verbose           Print per-lambda progress? Default TRUE.
#' @param ...               Additional arguments forwarded to \code{sbvar()}.
#'
#' @return List with components:
#'   \describe{
#'     \item{cv_table}{data.frame with columns lambda, mspe, mspe_se.}
#'     \item{best}{Single-row data.frame for the selected lambda.}
#'     \item{val_points}{Integer vector of validation row indices.}
#'     \item{keep_rows}{Integer vector of training row indices.}
#'   }
cv_sbvar <- function(Y,
                     p                  = 1L,
                     lambda             = NULL,
                     val_spacing        = 10L,
                     drop_poisoned_rows = FALSE,
                     lambda_rule        = c("min", "1se"),
                     verbose            = TRUE,
                     ...) {
  n  <- nrow(Y)
  d  <- ncol(Y)
  dp <- d * p

  # ---- lambda path (decreasing for warm restarts) --------------------------
  lambda_vec <- if (is.null(lambda)) 10^seq(-3, 0, length.out = 10L) else lambda
  lambda_vec <- sort(lambda_vec, decreasing = TRUE)
  n_lambda   <- length(lambda_vec)

  # ---- validation set -------------------------------------------------------
  if (val_spacing <= p) {
    stop(sprintf("val_spacing (%d) must be > p (%d).", val_spacing, p))
  }
  first_val  <- p + val_spacing
  last_val   <- n - p
  val_points <- seq(first_val, last_val, by = val_spacing)
  k          <- length(val_points)

  if (k < 2L) {
    stop("Fewer than 2 validation points. Reduce val_spacing or increase n.")
  }

  # ---- excluded rows -------------------------------------------------------
  excl_rows <- sort(unique(as.integer(unlist(lapply(val_points, function(t) {
    if (drop_poisoned_rows) t:(t + p) else t
  })))))
  keep_rows  <- setdiff(seq_len(n), excl_rows)
  n_train    <- length(keep_rows)
  n_eff_tr   <- n_train - p

  if (n_eff_tr < 1L) {
    stop("Too few training rows after exclusion. Reduce val_spacing or p.")
  }

  # ---- predictor matrix at validation points  (k x dp) ---------------------
  # x_val[jv, ] = c(Y[t-1,], ..., Y[t-p,]) for t = val_points[jv]
  Y_ext  <- rbind(matrix(0, p, d), Y)      # (n+p) x d, zero-padded
  x_val  <- matrix(0.0, k, dp)
  for (lag in seq_len(p)) {
    cols <- ((lag - 1L) * d + 1L):(lag * d)
    x_val[, cols] <- Y_ext[val_points - lag + p, , drop = FALSE]
  }
  Y_val <- Y[val_points, , drop = FALSE]   # k x d

  lambda_rule <- match.arg(lambda_rule, c("min", "1se"))

  # ---- row mapping: thinned → original -------------------------------------
  # After fitting sbvar on Y[keep_rows, ], its Theta has n_train rows:
  #   rows 1:p are zero (pre-sample)
  #   rows (p+1):n_train correspond to keep_rows[(p+1):n_train] in original Y
  # We scatter these increments back into an n x (d*dp) sparse Theta_orig,
  # then cumsum to get B_cumsum -- the same carry-forward logic as
  # l_val %*% fit$Theta in cv-hsbvar-nll.R.
  orig_eff_rows <- keep_rows[(p + 1L):n_train]  # original indices of effective obs

  # ---- main CV loop --------------------------------------------------------
  mspe    <- numeric(n_lambda)
  sq_err  <- matrix(NA_real_, n_lambda, k)
  cp_list <- vector("list", n_lambda)
  warm_phi <- NULL

  Y_train <- Y[keep_rows, , drop = FALSE]  # thinned data (n_train x d)

  for (j in seq_len(n_lambda)) {
    ln  <- lambda_vec[j]
    fit <- tryCatch(
      sbvar(Y_train,
            p        = p,
            lambda   = ln,
            ...),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      mspe[j] <- NA_real_
    } else {
      warm_phi <- fit$phi_hat

      # Scatter thinned Theta back to original time axis
      Theta_orig <- matrix(0.0, n, d * dp)
      Theta_orig[orig_eff_rows, ] <- fit$Theta[(p + 1L):n_train, ]

      # Cumulative coefficient at every original row (carry-forward at gaps)
      B_cumsum <- apply(Theta_orig, 2L, cumsum)   # n x (d*dp)

      for (jv in seq_len(k)) {
        t   <- val_points[jv]
        B_t <- matrix(B_cumsum[t, ], d, dp)
        Y_hat <- B_t %*% x_val[jv, ]
        sq_err[j, jv] <- mean((Y_val[jv, ] - Y_hat)^2)
      }
      mspe[j]    <- mean(sq_err[j, ], na.rm = TRUE)
      cp_list[[j]] <- fit$cp
    }

    if (verbose) {
      n_cp <- if (is.null(fit)) NA_integer_ else length(fit$cp)
      message(sprintf(
        "[%d/%d]  lambda=%.4g  MSPE=%.5g  ncp=%s",
        j, n_lambda, ln, mspe[j], n_cp
      ))
    }
  }

  # ---- assemble output -----------------------------------------------------
  mspe_se  <- apply(sq_err, 1L, sd, na.rm = TRUE) / sqrt(k)
  cv_table <- data.frame(lambda = lambda_vec, mspe = mspe, mspe_se = mspe_se)

  min_idx  <- which.min(mspe)
  best_idx <- if (lambda_rule == "min") {
    min_idx
  } else {
    # 1se rule: segment-aware variance estimate, same as cv-hsbvar-nll.R
    cp_min  <- cp_list[[min_idx]]
    seg_id  <- vapply(val_points, function(t) sum(t > cp_min) + 1L, integer(1L))
    seg_vars <- tapply(sq_err[min_idx, ], seg_id, stats::var)
    pt_var  <- as.numeric(seg_vars[as.character(seg_id)])
    pt_var[is.na(pt_var)] <- stats::var(sq_err[min_idx, ], na.rm = TRUE)
    cv_var  <- (1 / k) * sqrt(sum(pt_var, na.rm = TRUE))
    threshold  <- mspe[min_idx] + cv_var
    candidates <- which(!is.na(mspe) &
                          seq_len(n_lambda) <= min_idx &
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
#' @param cv_result Output of \code{cv_sbvar()}.
#' @param log_x     Log-scale the lambda axis? Default TRUE.
#' @param ...       Extra arguments passed to \code{plot()}.
plot_cv_sbvar <- function(cv_result, log_x = TRUE, ...) {
  tbl  <- cv_result$cv_table
  best <- cv_result$best

  x  <- tbl$lambda
  bx <- best$lambda
  if (log_x) {
    x  <- log10(x)
    bx <- log10(bx)
  }

  xlab <- if (log_x) expression(log[10](lambda)) else expression(lambda)
  plot(x, tbl$mspe,
       type = "b", pch = 19L, col = "steelblue",
       xlab = xlab, ylab = "CV MSPE",
       main = "SB-VAR \u2013 cross-validation MSPE",
       ...)
  abline(v = bx, lty = 2L, col = "firebrick")
  legend("topright",
         legend = sprintf("Best = %.4g\nMSPE = %.5g", best$lambda, best$mspe),
         bty = "n", text.col = "firebrick")
  invisible(cv_result)
}
