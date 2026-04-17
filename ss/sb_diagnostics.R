# sb_diagnostics.R
# Model evaluation and visualisation utilities:
#   estimation.check      — L2 error, FPR/TPR, lag detection vs. ground truth
#   estimation.check.new  — enhanced version for segmented models
#   mspe.plot             — MSPE vs. tuning parameter plot
#   plot.new              — multivariate time series plot
#   MTSplot.new           — MTS plot (legacy)
#   plot.matrix           — heatmap of VAR coefficient matrices

source(file.path("ss", "sb_utils.R"))

estimation.check <- function(phi, phi.hat) {
  k <- length(phi[, 1])
  p <- (length(phi[1, ])) / k
  N <- (length(phi.hat[1, ])) / (k * p)
  L <- matrix(0, k, k)
  l2.error <- rep(0, N)
  true.zero <- rep(0, N)
  true.non.zero <- rep(0, N)
  L.hat <- matrix(0, k, k * N)
  L.hat.final <- rep(0, N)
  false.zero <- rep(0, N)
  false.non.zero <- rep(0, N)

  count.non <- 0
  count.zero <- 0
  for (i in 1:k) {
    for (j in 1:(k * p)) {
      if (phi[i, j] != 0) {
        count.non <- count.non + 1
      }
      if (phi[i, j] == 0) {
        count.zero <- count.zero + 1
      }
    }
  }
  for (i in 1:k) {
    for (j in 1:k) {
      temp <- 0
      for (l in 1:p) {
        if (phi[i, ((l - 1) * k + j)] != 0) {
          temp <- temp + 1
        }
      }
      L[i, j] <- temp
    }
  }
  for (jj in 1:N) {
    phi.temp <- phi.hat[, ((jj - 1) * k * p + 1):(jj * k * p)]
    count.false.zero <- 0
    count.false.non.zero <- 0
    count.true.non.zero <- 0
    count.true.zero <- 0
    for (i in 1:k) {
      for (j in 1:(k * p)) {
        if (phi[i, j] != 0 && phi.hat[i, ((jj - 1) * k * p + j)] == 0) {
          count.false.zero <- count.false.zero + 1
        }
        if (phi[i, j] == 0 && phi.hat[i, ((jj - 1) * k * p + j)] != 0) {
          count.false.non.zero <- count.false.non.zero + 1
        }
        if (phi[i, j] == 0 && phi.hat[i, ((jj - 1) * k * p + j)] == 0) {
          count.true.zero <- count.true.zero + 1
        }
        if (phi[i, j] != 0 && phi.hat[i, ((jj - 1) * k * p + j)] != 0) {
          count.true.non.zero <- count.true.non.zero + 1
        }
        if (phi[i, j] == 0) {
          phi.temp[i, j] <- 0
        }
      }
    }
    l2.error[jj] <- sum((phi.temp - phi)^2)
    true.zero[jj] <- count.true.zero
    true.non.zero[jj] <- count.true.non.zero
    false.zero[jj] <- count.false.zero
    false.non.zero[jj] <- count.false.non.zero
  }
  for (jj in 1:N) {
    for (i in 1:k) {
      for (j in 1:k) {
        temp <- 0
        for (l in 1:p) {
          if (phi.hat[i, (((jj - 1) * k * p) + ((l - 1) * k) + j)] != 0) {
            temp <- temp + 1
          }
        }
        L.hat[i, (((jj - 1) * k) + j)] <- temp
      }
    }
  }
  L.hat.final <- sapply(c(1:N), function(jj) (sum((L.hat[, ((jj - 1) * k + 1):(jj * k)] - L)^2)) / (sum(L)))

  return(list(
    l2.error.mean = mean(l2.error), l2.error.sd = sd(l2.error), true.zero.median = median(true.zero),
    true.non.zero.median = median(true.non.zero), false.zero.median = median(false.zero),
    false.non.zero.median = median(false.non.zero), lag.error.mean = mean(L.hat.final),
    lag.error.sd = sd(L.hat.final)
  ))
}

estimation.check.new <- function(phi, phi.final, k, p, m.hat, pts.final.full) {
  N <- length(phi.final)
  l2.error <- rep(0, N)
  true.zero <- rep(0, N)
  true.non.zero <- rep(0, N)
  false.zero <- rep(0, N)
  false.non.zero <- rep(0, N)

  count.non <- 0
  count.zero <- 0
  for (i in 1:k) {
    for (j in 1:(k * m.hat * p)) {
      if (phi[i, j] != 0) {
        count.non <- count.non + 1
      }
      if (phi[i, j] == 0) {
        count.zero <- count.zero + 1
      }
    }
  }

  for (jj in 1:N) {
    phi.temp <- phi.final[[jj]]
    pt.temp <- pts.final.full[[jj]]
    if (length(phi.temp[1, ]) < k * m.hat * p) {
      phi.temp.new <- matrix(0, k, k * m.hat * p)
      phi.temp.new[, (1):(k * p)] <- phi.temp[, (1):(k * p)]
      ind <- 0
      for (i in 2:m.hat) {
        if (pt.temp[i - 1] == 0) {
          ind <- ind + 1
          phi.temp.new[, ((i - 1) * k * p + 1):(i * k * p)] <- phi.temp.new[, ((i - 2) * k * p + 1):((i - 1) * k * p)]
        }
        if (pt.temp[i - 1] != 0) {
          phi.temp.new[, ((i - 1) * k * p + 1):(i * k * p)] <- phi.temp[, ((i - 1 - ind) * k * p + 1):((i - ind) * k * p)]
        }
      }
      phi.temp <- phi.temp.new
    }

    count.false.zero <- 0
    count.false.non.zero <- 0
    count.true.non.zero <- 0
    count.true.zero <- 0
    for (i in 1:k) {
      for (j in 1:(k * m.hat * p)) {
        if (phi[i, j] != 0 && phi.temp[i, j] == 0) {
          count.false.zero <- count.false.zero + 1
        }
        if (phi[i, j] == 0 && phi.temp[i, j] != 0) {
          count.false.non.zero <- count.false.non.zero + 1
        }
        if (phi[i, j] == 0 && phi.temp[i, j] == 0) {
          count.true.zero <- count.true.zero + 1
        }
        if (phi[i, j] != 0 && phi.temp[i, j] != 0) {
          count.true.non.zero <- count.true.non.zero + 1
        }
        if (phi[i, j] == 0) {
          phi.temp[i, j] <- 0
        }
      }
    }
    l2.error[jj] <- sqrt(sum((phi.temp - phi)^2) / sum(phi^2))
    true.zero[jj] <- count.true.zero
    true.non.zero[jj] <- count.true.non.zero
    false.zero[jj] <- count.false.zero
    false.non.zero[jj] <- count.false.non.zero
  }


  return(list(
    l2.error.mean = mean(l2.error), l2.error.sd = sd(l2.error), true.zero.median = median(true.zero),
    true.non.zero.median = median(true.non.zero), false.zero.median = median(false.zero),
    false.non.zero.median = median(false.non.zero)
  ))
}

mspe.plot <- function(method, pred.error, lambda, tune.final, jj) {
  plot(lambda, pred.error, type = "o", col = "blue", main = c(jj, method))
  abline(v = tune.final)
}

plot.new <- function(data, caltime = NULL) {
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  if (is.ts(data)) {
    plot(data)
  } else {
    nT <- dim(data)[1]
    tdx <- c(1:nT)
    if (length(caltime) > 1) {
      tdx <- caltime
    }
    k <- dim(data)[2]
    # if (k < 4) {
    #   par(mfcol = c(k, 1))
    #   for (j in 1:k) {
    #     plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j],
    #          type = "l")
    #   }
    # }

    if (k >= 1) {
      par(mfcol = c(1, 1))
      yl <- range(data) * 1.05
      # plot(tdx, data[, 1], xlab = "time", ylab = " ", type = "l",
      #      ylim = yl)
      plot(tdx, data[, 1],
        xlab = "seconds", ylab = " ", type = "l",
        ylim = yl, xaxt = "n"
      )
      axis(1, c(1, 500, 1000, 1500, 2000), c(0, 50, 100, 150, 200))
      axis(2)
      for (j in 2:k) {
        lines(tdx, data[, j], lty = j, col = j)
      }
    }
  }
  par(mfcol = c(1, 1))
}

MTSplot.new <- function(data, caltime = NULL) {
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  if (is.ts(data)) {
    plot(data)
  } else {
    nT <- dim(data)[1]
    tdx <- c(1:nT)
    if (length(caltime) > 1) {
      tdx <- caltime
    }
    k <- dim(data)[2]
    if (k < 0) {
      par(mfcol = c(k, 1))
      for (j in 1:k) {
        plot(tdx, data[, j],
          xlab = "time", ylab = colnames(data)[j],
          type = "l"
        )
      }
    }
    if (k == 0) {
      par(mfcol = c(2, 2))
      for (j in 1:k) {
        plot(tdx, data[, j],
          xlab = "time", ylab = colnames(data)[j],
          type = "l"
        )
      }
    }
    if ((k > 0) && (k < 1)) {
      par(mfcol = c(3, 2), mai = c(0.3, 0.3, 0.3, 0.3))
      k1 <- 6
      jcnt <- 0
      for (j in 1:k) {
        plot(tdx, data[, j],
          xlab = "time", ylab = colnames(data)[j],
          type = "l", cex.axis = 0.8
        )
        jcnt <- jcnt + 1
        if ((jcnt == k1) && (k > 6)) {
          jcnt <- 0
          cat("Hit return for more plots: ", "\n")
          readline()
        }
      }
    }
    if (k > 0) {
      par(mfcol = c(1, 1))
      yl <- range(data) * 1.05
      plot(tdx, data[, 1],
        xlab = "time", ylab = " ", type = "l",
        ylim = yl
      )
      for (j in 2:k) {
        lines(tdx, data[, j], lty = j, col = j)
      }
    }
  }
  par(mfcol = c(1, 1))
}

plot.matrix <- function(phi, p, name = NULL) {
  B <- phi
  if (nrow(B) == 1) {
    B <- matrix(B[, 1:ncol(B)], nrow = 1)
  } else {
    B <- B[, 1:ncol(B)]
  }
  k <- nrow(B)
  s1 <- 0
  m <- 0
  s <- 0
  s <- s + s1
  text <- c()
  for (i in 1:p) {
    text1 <- as.expression(bquote(bold(Phi)^(.(i))))
    text <- append(text, text1)
  }
  if (m > 0) {
    for (i in (p + 1):(p + s + 1)) {
      text1 <- as.expression(bquote(bold(beta)^(.(i - p -
        s1))))
      text <- append(text, text1)
    }
  }
  f <- function(m) t(m)[, nrow(m):1]
  rgb.palette <- colorRampPalette(c("white", "blue"), space = "Lab")
  at <- seq(k / 2 + 0.5, p * (k) + 0.5, by = k)
  if (m > 0) {
    at2 <- seq(p * k + m / 2 + 0.5, p * k + s * m + 0.5, by = m)
  } else {
    at2 <- c()
  }
  at <- c(at, at2)
  se2 <- seq(1.75, by = k, length = k)
  L2 <- levelplot(as.matrix(f(B)),
    col.regions = rgb.palette,
    colorkey = NULL, xlab = NULL, ylab = NULL, main = list(
      label = name,
      cex = 1
    ), panel = function(...) {
      panel.levelplot(...)
      panel.abline(a = NULL, b = 1, h = seq(1.5, m * s +
        p * k + 0.5, by = 1), v = seq(1.5, by = 1, length = p *
        k + m * s), lwd = 0.5)
      bl1 <- seq(k + 0.5, p * k + 0.5, by = k)
      b23 <- seq(p * k + 0.5, p * k + 0.5 + s * m, by = m)
      b1 <- c(bl1, b23)
      panel.abline(a = NULL, b = 1, v = p * k + 0.5, lwd = 3)
      panel.abline(a = NULL, b = 1, v = b1, lwd = 2)
    }, scales = list(x = list(
      alternating = 1, labels = text,
      cex = 1, at = at, tck = c(0, 0)
    ), y = list(
      alternating = 0,
      tck = c(0, 0)
    ))
  )
  return(L2)
}
