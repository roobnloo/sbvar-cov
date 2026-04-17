# sim-scenario1.R
# Runs CV + BEA for H-SBVAR (NLL and D-trace) and SB-VAR on Scenario 1.
#
# Scenario 1 DGP: d=2, p=1, n=200, one true break at t=100.
#   Regime 1: A = [[0.60, 0.10], [0.10, 0.55]],
#             Sigma = [[1.0, 0.2], [0.2, 1.0]]
#   Regime 2: A = [[-0.50, 0.10], [0.10, -0.45]],
#             Sigma = [[1.5, -0.3], [-0.3, 1.5]]

setwd("/Users/boethius/code/sbvar-cov")

source("generate-var-data.R")
source("cv-hsbvar-nll.R")
source("cv-hsbvar-dtrace.R")
source("hsbvar-bea.R")
source(file.path("ss", "cv_sb_stage1_clean.R"))
source(file.path("ss", "sb_stage2_bea_clean.R"))

fmt_cp <- function(cp) {
  if (length(cp) == 0) "none" else paste(cp, collapse = ", ")
}

# ---- data --------------------------------------------------------------------
sim <- generate_scenario1()
y <- sim$y
n <- sim$n
d <- sim$d
p <- sim$p
true_cp <- sim$break_points # 100

cat(sprintf(
  "\n=== Scenario 1: n=%d  d=%d  p=%d  true break at t=%s ===\n\n",
  n, d, p, fmt_cp(true_cp)
))

# ---- NLL: CV -----------------------------------------------------------------
cat("--- NLL method: cross-validation ---\n")
cv_nll <- cv_hsbvar_nll(y, p = p, verbose = TRUE)

cat(sprintf(
  "\n  Best lambda (NLL) = %.4g  (CV MSPE = %.5g)\n",
  cv_nll$best$lambda, cv_nll$best$mspe
))

# Re-fit at best lambda to extract Stage 1 changepoints
fit_nll <- hsbvar(y, p = p, lambda = cv_nll$best$lambda)
cp_nll_stage1 <- fit_nll$cp

cat(sprintf(
  "  Stage 1 changepoints (NLL) : %s\n",
  fmt_cp(cp_nll_stage1)
))

# BEA pruning
bea_nll <- hsbvar_bea(fit_nll, y, p = p)
cp_nll_bea <- bea_nll$cp

cat(sprintf(
  "  After BEA             (NLL) : %s\n",
  fmt_cp(cp_nll_bea)
))

# ---- D-trace: CV -------------------------------------------------------------
cat("\n--- D-trace method: cross-validation ---\n")
cv_dtrace <- cv_hsbvar_dtrace(y, p = p, verbose = TRUE)

cat(sprintf(
  "\n  Best lambda (D-trace) = %.4g  (CV MSPE = %.5g)\n",
  cv_dtrace$best$lambda, cv_dtrace$best$mspe
))

# Re-fit at best lambda
fit_dtrace <- hsbvar_dtrace(y, p = p, lambda = cv_dtrace$best$lambda)
cp_dtrace_stage1 <- fit_dtrace$cp

cat(sprintf(
  "  Stage 1 changepoints (D-trace) : %s\n",
  fmt_cp(cp_dtrace_stage1)
))

# BEA pruning
bea_dtrace <- hsbvar_bea(fit_dtrace, y, p = p)
cp_dtrace_bea <- bea_dtrace$cp

cat(sprintf(
  "  After BEA             (D-trace) : %s\n",
  fmt_cp(cp_dtrace_bea)
))

# ---- SB-VAR: CV --------------------------------------------------------------
cat("\n--- SB-VAR method: cross-validation ---\n")
cv_sb <- cv_sbvar(y, p = p, verbose = TRUE)

cat(sprintf(
  "\n  Best lambda (SB-VAR) = %.4g  (CV MSPE = %.5g)\n",
  cv_sb$best$lambda, cv_sb$best$mspe
))

# Re-fit at best lambda to extract Stage 1 changepoints
fit_sb <- sbvar(y, p = p, lambda = cv_sb$best$lambda)
cp_sb_stage1 <- fit_sb$cp

cat(sprintf(
  "  Stage 1 changepoints (SB-VAR) : %s\n",
  fmt_cp(cp_sb_stage1)
))

# BEA pruning — three IC variants
n_eff        <- n - p
omega_author <- (1 / 1.75) * log(n_eff) * log(d)
lambda2      <- log(n_eff) * log(d) / n_eff

bea_sb_rss     <- sbvar_bea(fit_sb, y, p = p, lambda = lambda2, omega_n = omega_author,
                            ic_type = "rss")
bea_sb_scaled  <- sbvar_bea(fit_sb, y, p = p, omega_n = omega_author,
                            ic_type = "sigma_scaled")
bea_sb_profile <- sbvar_bea(fit_sb, y, p = p, omega_n = omega_author,
                            ic_type = "profile_lik")

cat(sprintf("  After BEA rss         (SB-VAR) : %s\n", fmt_cp(bea_sb_rss$cp)))
cat(sprintf("  After BEA sigma_scaled(SB-VAR) : %s\n", fmt_cp(bea_sb_scaled$cp)))
cat(sprintf("  After BEA profile_lik (SB-VAR) : %s\n", fmt_cp(bea_sb_profile$cp)))

# ---- Summary -----------------------------------------------------------------
cat(sprintf(
  "\n=== Summary (true break: t=%s) ===\n",
  fmt_cp(true_cp)
))
cat(sprintf(
  "  NLL          Stage 1: %-30s  BEA: %s\n",
  fmt_cp(cp_nll_stage1), fmt_cp(cp_nll_bea)
))
cat(sprintf(
  "  Dtrace       Stage 1: %-30s  BEA: %s\n",
  fmt_cp(cp_dtrace_stage1), fmt_cp(cp_dtrace_bea)
))
cat(sprintf(
  "  SB-VAR (rss) Stage 1: %-30s  BEA: %s\n",
  fmt_cp(cp_sb_stage1), fmt_cp(bea_sb_rss$cp)
))
cat(sprintf(
  "  SB-VAR (sig) Stage 1: %-30s  BEA: %s\n",
  fmt_cp(cp_sb_stage1), fmt_cp(bea_sb_scaled$cp)
))
cat(sprintf(
  "  SB-VAR (pro) Stage 1: %-30s  BEA: %s\n",
  fmt_cp(cp_sb_stage1), fmt_cp(bea_sb_profile$cp)
))
