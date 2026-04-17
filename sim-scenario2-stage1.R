# sim-scenario2-stage1.R
# Runs Safikhani-Shojaie Stage 1 (cv_sbdetect_stage1) on scenario 2 data.
#
# Scenario 2 DGP mirrors SBDetection_simulation_1.R:
#   k=20 series, p=1, n=300, two breaks at t=101 and t=201.
#   Regime coefficient matrices have a sparse superdiagonal:
#     Regime 1: A[j, j+1] = -0.6   (j = 1..19)
#     Regime 2: A[j, j+1] =  0.75
#     Regime 3: A[j, j+1] = -0.8
#   Innovation covariance: 0.01 * I_20.
#
# CV settings are chosen to match the Stage 1 step of SBDetection_simulation_1.R:
#   lambda grid : seq(0.10, 0.50, 0.05)  (= lambda.1.cv)
#   val_spacing : 10                      (= cv.length)
#   max_iter    : 100                     (= max.iteration)
#   tol         : 4e-2                    (= tol)
#   lambda_rule : "1se"                   (= lll + ind.new selection)

# ---- working directory -------------------------------------------------------
# All source() paths assume the sbvar-cov root.

# ---- sources -----------------------------------------------------------------
source("generate-var-data.R")
source(file.path("ss", "cv_sbdetect_stage1.R"))

# ---- data generation ---------------------------------------------------------
sim <- generate_scenario2()
Y <- sim$y
n <- sim$n
d <- sim$d
p <- sim$p
true_cp <- sim$break_points # c(101, 201)

cat(sprintf(
  "\n=== Scenario 2: n=%d  d=%d  p=%d  true breaks at t=%s ===\n\n",
  n, d, p, paste(true_cp, collapse = ", ")
))

# ---- Stage 1: CV break point detection ---------------------------------------
# Parameters match SBDetection_simulation_1.R Stage 1 (first.step.cv.new).
lambda_grid <- seq(0.10, 0.50, 0.05)

cv_fit <- cv_sbdetect_stage1(
  Y,
  p           = p,
  lambda      = lambda_grid,
  val_spacing = 10L,
  max_iter    = 100L,
  tol         = 4e-2,
  lambda_rule = "1se",
  verbose     = TRUE
)

# Re-run at best lambda to extract break points.
source(file.path("ss", "sbdetect_stage1.R"))
fit <- sbdetect_stage1(Y, p = p, lambda = cv_fit$best$lambda)
cp <- fit$cp

cat(sprintf("\n  Best lambda = %.4g  (MSPE = %.5g)\n", cv_fit$best$lambda, cv_fit$best$mspe))
cat(sprintf(
  "  Detected changepoints : %s\n",
  if (length(cp) == 0) "none" else paste(cp, collapse = ", ")
))
cat(sprintf(
  "  True    changepoints  : %s\n\n",
  paste(true_cp, collapse = ", ")
))

# ---- CV curve ----------------------------------------------------------------
plot_cv_sbdetect_stage1(cv_fit)
