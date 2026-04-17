# sbdetect_stage1.R
# Thin wrapper around the Safikhani-Shojaie Stage 1 initial change point
# detection (before BEA) with an interface matching hsbvar() in hsbvar-nll.R.
#
# @param Y         n x d numeric matrix (rows = observations, cols = series)
# @param p         VAR lag order (default 1)
# @param lambda    Penalty strength, single numeric (default 0.1)
# @param method    Estimation method passed to var.break.fit (default "LASSO")
# @param max_iter  Maximum coordinate-descent iterations (default 1000)
# @param tol       Convergence tolerance (default 1e-4)
#
# @return List with:
#   cp      — integer vector of changepoint row indices into Y, same scale as
#             hsbvar()'s $cp (i.e. 1-indexed rows of Y)
#   phi_hat — k x (k*p*(T-p)) coefficient increment matrix from var.break.fit

source(file.path("ss", "sb_stage1.R"))

sbdetect_stage1 <- function(Y,
                            p = 1,
                            lambda = 0.1,
                            method = "LASSO",
                            max_iter = 1000,
                            tol = 1e-4) {
  # first.step() reads T and k as free globals from its enclosing (global) env.
  # T = nrow(Y) and k = ncol(Y) — both are properties of the data.
  # Inject them via a local closure
  local_first_step <- first.step
  e <- new.env(parent = environment(first.step))
  e$T <- nrow(Y)
  e$k <- ncol(Y)
  environment(local_first_step) <- e

  result <- local_first_step(
    method        = method,
    data.temp     = Y,
    weight        = NULL,
    lambda        = lambda,
    p             = p,
    max.iteration = max_iter,
    tol           = tol
  )

  # first.step returns brk.points as a single-element list (ll = c(0)).
  # Indices are in n = T-p response space; add p to get Y row indices.
  cp_n <- result$brk.points[[1]]
  if (is.null(cp_n) || length(cp_n) == 0) cp_n <- integer(0)
  cp <- cp_n + p # translate to Y row indices (matches hsbvar $cp convention)

  list(
    cp      = cp,
    phi_hat = result$phi.hat
  )
}
