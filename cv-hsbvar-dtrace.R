# cv-hsbvar-dtrace.R
# Interpolation-based cross-validation for H-SBVAR with D-trace loss.
#
# Methodology mirrors cv-hsbar.R (Safikhani & Shojaie 2022):
#   1. Choose k equally spaced validation points with spacing s > p.
#   2. For each t in T remove the validation row and (optionally) the p
#      downstream rows that would contain Y_t as a lagged predictor.
#   3. Fit H-SBVAR on the thinned training set via hsbvar_dtrace(keep_rows=...).
#   4. Recover B_t = Phi_t^{-1} G_t at each validation point and predict
#        Y_hat_t = B_t %*% X_t
#   5. MSPE(lambda) = mean over val points of mean((Y_t - Y_hat_t)^2).
#
# Prediction always uses the NLL-natural B_t = Phi_t^{-1} G_t, regardless of
# which loss was used for training.  Lambda must be tuned separately from the
# NLL version because the D-trace loss has different curvature.
#
# Varies lambda over a log-linear path; c_scale is fixed.

if (!exists("hsbvar_dtrace", mode = "function")) source("hsbvar-dtrace.R")

#' Interpolation CV for H-SBVAR (D-trace loss) over a lambda path
#'
#' @param Y                 n x d numeric matrix.
#' @param p                 VAR lag order.
#' @param lambda            Numeric vector of lambda candidates.
#'                          Default: 10-point log grid on [1e-3, 1].
#' @param c_scale           Fixed scale c (not cross-validated). Default 1.
#' @param val_spacing       Spacing between validation points. Must be > p.
#'                          Default 10.
#' @param drop_poisoned_rows If TRUE, also exclude the p rows following each
#'                          validation point. Default FALSE.
#' @param lambda_rule       \code{"min"} or \code{"1se"}. Default \code{"min"}.
#' @param verbose           Print per-lambda progress? Default TRUE.
#' @param ...               Additional arguments forwarded to \code{hsbvar_dtrace()}.
#'
#' @return List: cv_table (lambda, mspe, mspe_se), best, val_points, keep_rows.
cv_hsbvar_dtrace <- function(Y,
                             p                  = 1L,
                             lambda             = NULL,
                             c_scale            = 1,
                             val_spacing        = 10L,
                             drop_poisoned_rows = FALSE,
                             lambda_rule        = c("min", "1se"),
                             verbose            = TRUE,
                             ...) {
  n  <- nrow(Y)
  d  <- ncol(Y)
  dp <- d * p

  # ---- lambda path (decreasing for warm restarts) ------------------------
  lambda_vec <- if (is.null(lambda)) 10^seq(-3, 0, length.out = 10) else lambda
  lambda_vec <- sort(lambda_vec, decreasing = TRUE)
  n_lambda   <- length(lambda_vec)

  # ---- validation set -----------------------------------------------------
  if (val_spacing <= p) {
    stop(sprintf("val_spacing (%d) must be > p (%d).", val_spacing, p))
  }
  first_valid <- p + val_spacing
  last_valid  <- n - p
  val_points  <- seq(first_valid, last_valid, by = val_spacing)
  k           <- length(val_points)

  if (k < 2L) {
    stop("Fewer than 2 validation points. Reduce val_spacing or increase n.")
  }

  # ---- excluded rows ------------------------------------------------------
  excl_rows <- sort(unique(as.integer(unlist(lapply(val_points, function(t) {
    if (drop_poisoned_rows) t:(t + p) else t
  })))))
  keep_rows <- setdiff(seq_len(n), excl_rows)

  if (length(keep_rows) < 2L * p) {
    stop("Too few training rows after exclusion. Reduce val_spacing or n.")
  }

  # ---- predictor matrix for validation points  (k x dp) ------------------
  Y_ext <- rbind(matrix(0, p, d), Y)
  x_val <- matrix(0.0, k, dp)
  for (lag in seq_len(p)) {
    cols          <- ((lag - 1L) * d + 1L):(lag * d)
    x_val[, cols] <- Y_ext[val_points - lag + p, , drop = FALSE]
  }
  Y_val <- Y[val_points, , drop = FALSE]  # k x d

  # ---- cumsum selector for validation rows --------------------------------
  l_val <- base::outer(val_points, seq_len(n), ">=") + 0L  # k x n

  lambda_rule <- match.arg(lambda_rule, c("min", "1se"))

  # ---- main CV loop -------------------------------------------------------
  mspe       <- numeric(n_lambda)
  sq_err     <- matrix(NA_real_, n_lambda, k)
  cp_list    <- vector("list", n_lambda)
  warm_Theta <- NULL
  warm_Psi   <- NULL

  for (j in seq_len(n_lambda)) {
    ln  <- lambda_vec[j]
    fit <- tryCatch(
      hsbvar_dtrace( # nolint: object_usage_linter
        Y, p,
        lambda     = ln,
        c_scale    = c_scale,
        keep_rows  = keep_rows,
        init_Theta = warm_Theta,
        init_Psi   = warm_Psi,
        ...
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      mspe[j] <- NA_real_
    } else {
      warm_Theta <- fit$Theta
      warm_Psi   <- fit$Psi

      Gamma_val <- l_val %*% fit$Theta  # k x ddp
      Phi_val   <- l_val %*% fit$Psi    # k x d^2

      for (jv in seq_len(k)) {
        Phi_t <- matrix(Phi_val[jv, ], d, d)
        Phi_t <- (Phi_t + t(Phi_t)) * 0.5
        ch    <- tryCatch(chol(Phi_t), error = function(e) NULL)
        if (is.null(ch)) next
        B_t           <- chol2inv(ch) %*% matrix(Gamma_val[jv, ], d, dp)
        Y_hat         <- B_t %*% x_val[jv, ]
        sq_err[j, jv] <- mean((Y_val[jv, ] - Y_hat)^2)
      }
      mspe[j]      <- mean(sq_err[j, ], na.rm = TRUE)
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

  # ---- assemble output ----------------------------------------------------
  mspe_se  <- apply(sq_err, 1, sd) / sqrt(k)
  cv_table <- data.frame(lambda = lambda_vec, mspe = mspe, mspe_se = mspe_se)

  min_idx  <- which.min(mspe)
  best_idx <- if (lambda_rule == "min") {
    min_idx
  } else {
    cp_min   <- cp_list[[min_idx]]
    seg_id   <- vapply(val_points, function(t) sum(t > cp_min) + 1L, integer(1L))
    seg_vars <- tapply(sq_err[min_idx, ], seg_id, stats::var)
    pt_var   <- as.numeric(seg_vars[as.character(seg_id)])
    pt_var[is.na(pt_var)] <- stats::var(sq_err[min_idx, ], na.rm = TRUE)
    cv_var    <- (1 / k) * sqrt(sum(pt_var, na.rm = TRUE))
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
#' @param cv_result Output of \code{cv_hsbvar_dtrace()}.
#' @param log_x     Log-scale the lambda axis? Default TRUE.
#' @param ...       Extra arguments passed to \code{plot()}.
plot_cv_hsbvar_dtrace <- function(cv_result, log_x = TRUE, ...) {
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
    type = "b", pch = 19, col = "steelblue",
    xlab = xlab, ylab = "CV MSPE",
    main = "H-SBVAR D-trace \u2013 cross-validation MSPE",
    ...
  )
  abline(v = bx, lty = 2, col = "firebrick")
  legend("topright",
    legend = sprintf("Best = %.4g\nMSPE = %.5g", best$lambda, best$mspe),
    bty = "n", text.col = "firebrick"
  )
  invisible(cv_result)
}
