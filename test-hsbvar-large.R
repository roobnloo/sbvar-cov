# hsbvar-test-large.R
# Timing study for H-SBVAR (NLL and D-trace FISTA solvers) as d increases.
# Fixed: n=300, p=1, lambda=0.3, one changepoint at t=150.

source("hsbvar-nll.R")
source("hsbvar-dtrace.R")
source("generate-var-data.R")

n <- 300
p <- 1
lambda <- 0.3
cp_true <- 150

d_seq <- c(1, 2, 3, 5, 8, 10, 15)

hdr <- sprintf("%-4s  %-12s  %-6s  %-25s  %-12s  %-6s  %s",
               "d", "NLL time(s)", "iter", "NLL cp",
               "Dtrace time(s)", "iter", "Dtrace cp")
cat(hdr, "\n")
cat(strrep("-", nchar(hdr)), "\n")

for (d in d_seq) {
  # Diagonal-dominant coefficient matrices; spectral radius well inside unit disk
  a1 <- matrix(0.04, d, d)
  diag(a1) <- 0.40
  a2 <- matrix(0.04, d, d)
  diag(a2) <- -0.35

  sim <- generate_var_piecewise(
    n            = n,
    break_points = c(cp_true),
    a_list       = list(list(a1), list(a2)),
    sigma_list   = list(diag(d), diag(d)),
    seed         = 42L
  )

  t_nll <- system.time(
    fit_nll <- hsbvar(sim$y, p = p, lambda = lambda,
                      max_iter = 2000, verbose = FALSE)
  )[["elapsed"]]

  t_dtrace <- system.time(
    fit_dtrace <- hsbvar_dtrace(sim$y, p = p, lambda = lambda,
                                max_iter = 2000, verbose = FALSE)
  )[["elapsed"]]

  cp_nll    <- if (length(fit_nll$cp)    == 0) "none" else paste(fit_nll$cp,    collapse = " ")
  cp_dtrace <- if (length(fit_dtrace$cp) == 0) "none" else paste(fit_dtrace$cp, collapse = " ")

  cat(sprintf("%-4d  %-12.3f  %-6d  %-25s  %-12.3f  %-6d  %s\n",
              d,
              t_nll,    fit_nll$n_iter,    cp_nll,
              t_dtrace, fit_dtrace$n_iter, cp_dtrace))
}
