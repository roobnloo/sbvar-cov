# sb_stage3_refit.R
# Stage 3: Post-break parameter estimation (refitting).
#   ar.est — estimates VAR coefficients per segment after removing observations
#             within r.n of each detected break point

source(file.path("ss", "sb_estimation.R"))

ar.est <- function(method, data, weight = NULL, lambda, p, break.pts, r.n, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3)) {
  method.full <- c("LASSO")
  if (!(method %in% method.full)) {
    print("ERROR")
    break
  }

  k <- length(data[1, ])
  T <- length(data[, 1])
  T.1 <- T
  m.hat <- length(break.pts) + 1
  ind.remain <- rep(0, (2 + 2 * length(break.pts)))
  ind.remain[1] <- p
  ind.remain[(2 + 2 * length(break.pts))] <- T.1
  for (i in 1:length(break.pts)) {
    ind.remain[(2 * i)] <- break.pts[i] - r.n - 1
    ind.remain[(2 * i + 1)] <- break.pts[i] + r.n + 1
  }
  iter <- matrix(0, k, length(lambda))
  phi.hat <- matrix(0, k, k * m.hat * p)
  phi.hat.fista <- matrix(0, max.iteration, k * m.hat * p)
  pred.error <- rep(0, length(lambda))
  phi.hat.temp <- matrix(0, k, k * m.hat * p * length(lambda))
  std.res <- rep(0, length(lambda))

  Y <- as.matrix(t(data))
  Y <- Y[, -seq(1, p, 1)]
  Z <- matrix(0, k * p, T.1 - p)
  # Z <- Z%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  for (i in 1:(T.1 - p)) {
    for (j in 1:p) {
      Z[((j - 1) * k + 1):(j * k), i] <- t(data[i + p - j, ])
    }
  }

  for (i in 1:length(break.pts)) {
    Y <- Y[, -seq(break.pts[i] - r.n, break.pts[i] + r.n, 1)]
    # Z <- Z[,-seq(break.pts[i]-r.n,break.pts[i]+r.n,1)];
  }

  n <- length(Y[1, ])
  Z.new <- matrix(0, T.1 - p, k * m.hat * p)
  Z.new[(1:(break.pts[1] - r.n - 1)), 1:(k * p)] <- t(Z[1:(k * p), (1:(break.pts[1] - r.n - 1))])
  if (m.hat > 2) {
    for (i in 1:(m.hat - 2)) {
      # ind <- break.pts[i]-r.n-1;
      Z.new[((break.pts[i] + r.n + 1):(break.pts[i + 1] - r.n - 1)), ((i) * k * p + 1):((i + 1) * k * p)] <- t(Z[1:(k * p), ((break.pts[i] + r.n + 1):(break.pts[i + 1] - r.n - 1))])
    }
  }

  Z.new[((break.pts[m.hat - 1] + r.n + 1):(T.1 - p)), ((m.hat - 1) * k * p + 1):((m.hat) * k * p)] <- t(Z[1:(k * p), ((break.pts[m.hat - 1] + r.n + 1):(T.1 - p))])
  del.ind <- c()
  for (i in 1:(T.1 - p)) {
    if (sum(Z.new[i, ]^2) == 0) {
      del.ind <- c(del.ind, i)
    }
  }
  Z.new <- Z.new[-del.ind, ]
  Z <- t(Z.new)

  step.size <- 1 / (max(svd(Z)$d))^2
  # print(step.size)
  # step.size <- 10^(-3);


  for (ll in 1:length(lambda)) {
    for (ii in 1:k) {
      l <- 2
      while (l < max.iteration) {
        l <- l + 1
        phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
        phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
        phi.new <- soft(phi.new, rep(0, k * m.hat * p), lambda[ll])
        if (max(abs(phi.new - phi.temp)) < tol) {
          break
        }
        if (max(abs(phi.new - phi.temp)) > tol) {
          phi.hat.fista[l, ] <- phi.new
          # print(l);
        }
      }

      # print("l="); print(l)
      iter[ii, ll] <- l
      phi.hat.temp[ii, ((ll - 1) * k * m.hat * p + 1):(ll * k * m.hat * p)] <- phi.new
    }

    forecast <- matrix(0, k, T.1)
    for (i in 1:m.hat) {
      len <- ind.remain[(2 * i)] - ind.remain[(2 * i - 1)]
      delay <- floor(2 * len / 3)
      delay <- 0
      l.b <- ind.remain[(2 * i - 1)] + delay + 1
      u.b <- ind.remain[(2 * i)]
      forecast[, (l.b):(u.b)] <- sapply(
        c((l.b - 1):(u.b - 1)),
        function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * m.hat * p + (i - 1) * k * p + 1):((ll - 1) * k * m.hat * p + (i) * k * p)], p, jjj, k, 1)
      )
    }

    if (ll == 1) {
      del.ind <- c()
      for (i in 1:(T.1)) {
        if (sum(forecast[, i]^2) == 0) {
          del.ind <- c(del.ind, i)
        }
      }
    }

    # foo <- t(forecast);
    # plot(c(1:300),foo[,1],type="l",col="red")
    # lines(c(1:300),data[,1],col="green")

    residual <- t(data) - forecast
    BIC.temp <- 0
    for (i in 1:m.hat) {
      l.b <- ind.remain[(2 * i - 1)] + 0 + 1
      u.b <- ind.remain[(2 * i)]
      temp <- AIC.BIC(residual[, (l.b):(u.b)], phi.hat.temp[, ((ll - 1) * k * m.hat * p + (i - 1) * k * p + 1):((ll - 1) * k * m.hat * p + (i) * k * p)])
      BIC.temp <- BIC.temp + temp$BIC
    }


    residual <- residual[, -del.ind]
    nn <- length(residual[1, ])
    # check <- AIC.BIC(residual,phi.hat.temp[,((ll-1)*k*m.hat*p+1):(ll*k*m.hat*p)]);
    # pred.error[ll] <- sum(residual^2)/(k*nn); std.res[ll] <- sd(residual[,]^2);
    # pred.error[ll] <- check$BIC;
    pred.error[ll] <- BIC.temp
    # pred.error[ll] <- check$AIC;
    # print(ll)
  }
  ll.final <- which(pred.error == min(pred.error))
  ll.final <- min(ll.final)
  # sd.final <- 1*median(std.res)/(sqrt(k*nn));
  # # sd.final <- 1*max(std.res)/(sqrt(k*nn));
  # # sd.final <- 3*sd(pred.error);
  # # sd.final <- quantile(pred.error[1:ll.final],0.975) + pred.error[ll.final];
  # # sd.final <- 1*median(std.res)/(sqrt(nn));
  # p.error <- abs(pred.error[ll.final:length(lambda)] - (pred.error[ll.final] + sd.final));
  # ll.final <- which(p.error==min(p.error)); ll.final <- min(ll.final);
  phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * m.hat * p + 1):(ll.final * k * m.hat * p)]

  # for(i in 1:k){
  #   for(j in 1:(k*m.hat*p)){
  #     if(abs(phi.hat[i,j]) < lambda[ll.final]  ){phi.hat[i,j] <- 0;}
  #   }
  # }

  return(list(phi.hat = phi.hat, iter = iter, pred.error = pred.error, tune.final = lambda[ll.final]))
}
