# sb_stage2_bea.R
# Stage 2: Backward Elimination Algorithm (BEA) for change point screening.
#   break.var    — evaluates segmented loss L_n for a given set of candidate break points
#   backward.sel — one BEA pass: computes loss with each candidate removed
#   second.step  — iterates BEA with information criterion until convergence

source(file.path("ss", "sb_estimation.R"))

break.var <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts) {
  k <- length(data[1, ])
  T <- length(data[, 1])
  m <- length(pts)
  L.n <- rep(0, m + 1)
  if (m == 0) {
    pts.temp <- c(1, T + 1)
  }
  if (m > 0) {
    pts.temp <- rep(0, m + 2)
    pts.temp[1] <- 1
    pts.temp[m + 2] <- T + 1
    pts.temp[(2):(m + 1)] <- pts
  }
  for (mm in 1:(m + 1)) {
    # print("mmmmm"); print(mm)
    data.temp <- data[(pts.temp[mm]):(pts.temp[(mm + 1)] - 1), ]
    try <- var.lasso.brk(data = data.temp, weight = NULL, lambda, p, max.iteration = 1000, tol = 10^(-4))
    L.n[mm] <- try$pred.error
  }
  return(list(L.n = sum(L.n)))
}

backward.sel <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts) {
  k <- length(data[1, ])
  T <- length(data[, 1])
  m <- length(pts)
  L.n <- rep(0, m)
  L.n.curr <- rep(0, 1)

  ###### SHVAR FUN FULL PTS ##################
  try <- break.var(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts)
  L.n.curr <- try$L.n

  for (mm in 1:m) {
    pts.temp <- pts[-mm]

    ### SHVAR FIT FUNCTION ###################
    try <- break.var(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts.temp)
    L.n[mm] <- try$L.n
  }
  return(list(L.n = L.n, L.n.curr = L.n.curr))
}

second.step <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts, omega) {
  m <- length(pts)
  if (m == 0) {
    break
  }
  mm <- 0
  while (mm < m) {
    mm <- mm + 1
    try <- backward.sel(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts)
    L.n <- try$L.n
    L.n.curr <- try$L.n.curr
    if (min(L.n) + (m - 1) * omega >= L.n.curr + m * omega) {
      ic <- L.n.curr + (m - mm + 1) * omega
      break
    }
    if (min(L.n) + (m - 1) * omega < L.n.curr + m * omega) {
      pts <- pts[-which(L.n == min(L.n))]
      print(pts)
    }
  }


  return(list(pts = pts, ic = ic))
}
