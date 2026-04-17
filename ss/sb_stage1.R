# sb_stage1.R
# Stage 1: Initial change point detection via TV-penalised VAR regression.
#   first.step        — single lambda; returns candidate break points
#   first.step.cv     — CV over lambda grid; selects best lambda by prediction error
#   first.step.cv.new — CV with held-out index set; variance-adjusted lambda selection
#
# Note: first.step() reads T and k as free variables from its enclosing environment.
# Use sbdetect_stage1() (in sbdetect_stage1.R) for a clean interface that handles
# this automatically without touching .GlobalEnv.

source(file.path("ss", "sb_estimation.R"))

first.step <- function(method, data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1 / 2) * 10^(-4)) {
  test <- var.break.fit(method, data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  phi.hat.full <- test$phi.hat
  ll <- seq(0.05, max(abs(phi.hat.full)), 0.05)
  temp.lam <- ll
  # temp.lam <- quantile(phi.hat.full,ll)
  ll <- c(0)
  brk.points.list <- vector("list", length(ll))

  for (j in 1:length(ll)) {
    phi.hat <- phi.hat.full
    # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);

    n <- T - p
    m.hat <- 0
    brk.points <- rep(0, n)

    for (i in 2:n)
    {
      if (sum((phi.hat[, ((i - 1) * k * p + 1):(i * k * p)])^2) != 0) {
        m.hat <- m.hat + 1
        brk.points[m.hat] <- i
      }
    }

    loc <- rep(0, m.hat)
    brk.points <- brk.points[1:m.hat]
    brk.points <- brk.points[which(brk.points > (p + 3))]
    brk.points <- brk.points[which(brk.points < (n))]
    m.hat <- length(brk.points)
    if (m.hat > 1) {
      for (mm in 1:(m.hat - 1)) {
        if (abs(brk.points[mm + 1] - brk.points[mm]) <= (p + 1)) {
          loc[mm] <- mm
        }
      }
    }
    loc <- loc[which(loc != 0)]
    if (length(loc) > 0) {
      brk.points <- brk.points[-loc]
    }
    brk.points.list[[j]] <- brk.points
  }


  return(list(brk.points = brk.points.list, phi.hat = phi.hat.full))
}


first.step.cv <- function(method, data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1 / 2) * 10^(-4)) {
  kk <- length(lambda)
  cv <- rep(0, kk)
  phi.final <- vector("list", kk)
  T <- length(data.temp[, 1])
  k <- length(data.temp[1, ])
  brk.points.final <- vector("list", kk)

  for (i in 1:kk) {
    test <- var.break.fit(method, data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol)
    phi.hat.full <- test$phi.hat
    phi.final[[i]] <- phi.hat.full

    # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
    # temp.lam <- ll;
    # temp.lam <- quantile(phi.hat.full,ll)
    ll <- c(0)
    brk.points.list <- vector("list", length(ll))

    for (j in 1:length(ll)) {
      phi.hat <- phi.hat.full
      # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);

      n <- T - p
      m.hat <- 0
      brk.points <- rep(0, n)

      for (iii in 2:n)
      {
        if (sum((phi.hat[, ((iii - 1) * k * p + 1):(iii * k * p)])^2) != 0) {
          m.hat <- m.hat + 1
          brk.points[m.hat] <- iii
        }
      }

      loc <- rep(0, m.hat)
      brk.points <- brk.points[1:m.hat]
      brk.points <- brk.points[which(brk.points > (p + 3))]
      brk.points <- brk.points[which(brk.points < (n))]
      m.hat <- length(brk.points)
      if (m.hat > 1) {
        for (mm in 1:(m.hat - 1)) {
          if (abs(brk.points[mm + 1] - brk.points[mm]) <= (p + 1)) {
            loc[mm] <- mm
          }
        }
      }
      loc <- loc[which(loc != 0)]
      if (length(loc) > 0) {
        brk.points <- brk.points[-loc]
      }
      brk.points.list[[j]] <- brk.points
    }

    brk.points.final[[i]] <- brk.points
    m.hat <- length(brk.points)



    phi.full.all <- vector("list", T - p)
    phi.temp.cv <- matrix(0, k, k * p)
    forecast <- matrix(0, k, T - p)
    for (j in (p + 1):T) {
      phi.temp.cv <- phi.temp.cv + phi.hat.full[, ((j - p - 1) * k * p + 1):((j - p) * k * p)]
      phi.full.all[[(j - p)]] <- phi.temp.cv
      forecast[, (j - p)] <- pred(t(data), phi.temp.cv, p, j - 1, k, 1)
    }

    phi.hat.temp <- matrix(0, k, (k * p * (m.hat + 1)))
    if (m.hat > 1) {
      for (jj in 1:m.hat) {
        phi.hat.temp[, ((jj - 1) * k * p + 1):((jj) * k * p)] <- phi.full.all[[(brk.points[jj] - p - 1)]]
      }
    }
    phi.hat.temp[, ((m.hat) * k * p + 1):((m.hat + 1) * k * p)] <- phi.full.all[[(T - p)]]

    residual <- t(data[(p + 1):T, ]) - forecast
    BIC.temp <- 0
    if (m.hat == 0) {
      temp <- AIC.BIC.CV(residual[, (1):(T - p)], phi.hat.temp[, ((1 - 1) * k * m.hat * p + (1 - 1) * k * p + 1):((1 - 1) * k * m.hat * p + (1) * k * p)])
      BIC.temp <- BIC.temp + temp$BIC
    }
    if (m.hat >= 1) {
      temp <- AIC.BIC.CV(residual[, (1):(brk.points[1] - 1)], phi.hat.temp[, ((1 - 1) * k * m.hat * p + (1 - 1) * k * p + 1):((1 - 1) * k * m.hat * p + (1) * k * p)])
      BIC.temp <- BIC.temp + temp$BIC
      if (m.hat >= 2) {
        for (ii in 2:m.hat) {
          l.b <- brk.points[(ii - 1)]
          u.b <- brk.points[ii] - 1
          temp <- AIC.BIC.CV(residual[, (l.b):(u.b)], phi.hat.temp[, ((1 - 1) * k * m.hat * p + (ii - 1) * k * p + 1):((1 - 1) * k * m.hat * p + (ii) * k * p)])
          BIC.temp <- BIC.temp + temp$BIC
        }
      }
      temp <- AIC.BIC.CV(residual[, (brk.points[m.hat]):(T - p)], phi.hat.temp[, ((1 - 1) * k * m.hat * p + (m.hat) * k * p + 1):((1 - 1) * k * m.hat * p + (m.hat + 1) * k * p)])
      BIC.temp <- BIC.temp + temp$BIC
    }


    cv[i] <- sum((forecast - t(data[(p + 1):T, ]))^2)
    # cv[i] <- BIC.temp;
    # temp <- AIC.BIC.CV(residual,phi.hat.temp);
    # cv[i] <- temp$BIC;
    print("====================================")
  }

  lll <- min(which(cv == min(cv)))
  phi.hat.full <- phi.final[[lll]]


  # test <- var.break.fit(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  # phi.hat.full <- test$phi.hat;
  # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
  # temp.lam <- ll;
  # # temp.lam <- quantile(phi.hat.full,ll)
  # ll <- c(0);
  # brk.points.list <- vector("list",length(ll));
  #
  # for(j in 1:length(ll)){
  #
  #   phi.hat <- phi.hat.full;
  #   # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
  #
  #   n <- T - p;
  #   m.hat <- 0; brk.points <- rep(0,n);
  #
  #   for (i in 2:n)
  #   {
  #     if ( sum((phi.hat[,((i-1)*k*p+1):(i*k*p)] )^2 ) != 0   ){
  #       m.hat <- m.hat + 1; brk.points[m.hat] <- i;
  #     }
  #   }
  #
  #   loc <- rep(0,m.hat);
  #   brk.points <- brk.points[1:m.hat];
  #   brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
  #   m.hat <- length(brk.points);
  #   if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
  #   loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
  #   brk.points.list[[j]] <- brk.points;
  # }


  return(list(brk.points = brk.points.final[[lll]], cv = cv, cv.final = lambda[lll]))
}



first.step.cv.new <- function(method, data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1 / 2) * 10^(-4), cv.index) {
  cv.l <- length(cv.index)
  data.org <- data.temp
  T.org <- length(data.temp[, 1])
  k.org <- length(data.temp[1, ])
  data.temp <- data.temp[-cv.index, ]
  kk <- length(lambda)
  cv <- rep(0, kk)
  cv.var <- rep(0, kk)
  phi.final <- vector("list", kk)
  T <- length(data.temp[, 1])
  k <- length(data.temp[1, ])
  brk.points.final <- vector("list", kk)

  for (i in 1:kk) {
    if (i == 1) {
      test <- var.break.fit(method, data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol)
    }
    if (i > 1) {
      initial.phi <- phi.final[[(i - 1)]]
      test <- var.break.fit(method, data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol, initial.phi = initial.phi)
    }
    phi.hat.full <- test$phi.hat
    phi.final[[i]] <- phi.hat.full

    # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
    # temp.lam <- ll;
    # temp.lam <- quantile(phi.hat.full,ll)
    ll <- c(0)
    brk.points.list <- vector("list", length(ll))

    for (j in 1:length(ll)) {
      phi.hat <- phi.hat.full
      # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);

      n <- T - p
      m.hat <- 0
      brk.points <- rep(0, n)

      for (iii in 2:n)
      {
        if (sum((phi.hat[, ((iii - 1) * k * p + 1):(iii * k * p)])^2) != 0) {
          m.hat <- m.hat + 1
          brk.points[m.hat] <- iii
        }
      }

      loc <- rep(0, m.hat)
      brk.points <- brk.points[1:m.hat]
      brk.points <- brk.points[which(brk.points > (p + 3))]
      brk.points <- brk.points[which(brk.points < (n))]
      m.hat <- length(brk.points)
      if (m.hat > 1) {
        for (mm in 1:(m.hat - 1)) {
          if (abs(brk.points[mm + 1] - brk.points[mm]) <= (p + 1)) {
            loc[mm] <- mm
          }
        }
      }
      loc <- loc[which(loc != 0)]
      if (length(loc) > 0) {
        brk.points <- brk.points[-loc]
      }
      brk.points.list[[j]] <- brk.points
    }

    brk.points.final[[i]] <- brk.points
    m.hat <- length(brk.points)



    phi.full.all <- vector("list", T - p)
    phi.temp.cv <- matrix(0, k, k * p)
    forecast <- matrix(0, k, T - p)
    forecast.new <- matrix(0, k, cv.l)
    for (j in (p + 1):T) {
      phi.temp.cv <- phi.temp.cv + phi.hat.full[, ((j - p - 1) * k * p + 1):((j - p) * k * p)]
      phi.full.all[[(j - p)]] <- phi.temp.cv
      forecast[, (j - p)] <- pred(t(data.temp), phi.temp.cv, p, j - 1, k, 1)
    }

    for (j in (1):cv.l) {
      forecast.new[, j] <- pred(t(data.org), phi.full.all[[(cv.index[j] - 1 - p - j + 1)]], p, cv.index[j] - 1, k, 1)
    }

    forecast.all <- matrix(0, k, T.org - p)
    forecast.all[, cv.index] <- forecast.new
    forecast.all[, -cv.index] <- forecast
    residual <- (forecast.all - t(data.org[((p + 1):T.org), ]))^2
    var.matrix <- matrix(0, k, m.hat + 1)
    if (m.hat == 0) {
      var.matrix <- sapply(c(1:k), function(jjj) var(residual[jjj, ]))
    }
    if (m.hat >= 1) {
      var.matrix[, 1] <- t(as.vector(sapply(c(1:k), function(jjj) var(residual[jjj, (1):(brk.points[1] - 1)]))))
      var.matrix[, (m.hat + 1)] <- sapply(c(1:k), function(jjj) var(residual[jjj, (brk.points[(m.hat)]):(T.org - p)]))
    }
    if (m.hat >= 2) {
      for (mm in 2:m.hat) {
        var.matrix[, mm] <- sapply(c(1:k), function(jjj) var(residual[jjj, (brk.points[(mm - 1)]):(brk.points[mm] - 1)]))
      }
    }


    if (m.hat == 0) {
      cv.var[i] <- sum(var.matrix)
    }
    if (m.hat >= 1) {
      for (i.1 in 1:cv.l) {
        ind <- 0
        for (i.2 in 1:m.hat) {
          if (cv.index[i.1] > brk.points[i.2]) {
            ind <- ind + 1
          }
        }
        cv.var[i] <- cv.var[i] + sum(var.matrix[, (ind + 1)])
      }
    }

    cv.var[i] <- (1 / (k * cv.l)) * sqrt(cv.var[i])

    # phi.hat.temp <- matrix(0,k,(k*p*(m.hat+1)));
    # if( m.hat > 1){
    #   for(jj in 1:m.hat){
    #     phi.hat.temp[,((jj-1)*k*p+1):((jj)*k*p)] <- phi.full.all[[(brk.points[jj]-p-1)]];
    #   }
    # }
    # phi.hat.temp[,((m.hat)*k*p+1):((m.hat+1)*k*p)] <- phi.full.all[[(T-p)]];
    #
    # residual <- t(data.temp[(p+1):T,]) - forecast;
    # BIC.temp <- 0;
    # if (m.hat == 0){
    #   temp <- AIC.BIC.CV(residual[,(1):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
    #   BIC.temp <- BIC.temp + temp$BIC;
    # }
    # if(m.hat >=1){
    #   temp <- AIC.BIC.CV(residual[,(1):(brk.points[1]-1)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
    #   BIC.temp <- BIC.temp + temp$BIC;
    #   if ( m.hat >= 2){
    #     for(ii in 2:m.hat){
    #       l.b <-  brk.points[(ii-1)]; u.b <- brk.points[ii]-1;
    #       temp <- AIC.BIC.CV(residual[,(l.b):(u.b)],phi.hat.temp[,((1-1)*k*m.hat*p+(ii-1)*k*p+1):((1-1)*k*m.hat*p+(ii)*k*p)]);
    #       BIC.temp <- BIC.temp + temp$BIC;
    #     }
    #   }
    #   temp <- AIC.BIC.CV(residual[,(brk.points[m.hat]):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(m.hat)*k*p+1):((1-1)*k*m.hat*p+(m.hat+1)*k*p)]);
    #   BIC.temp <- BIC.temp + temp$BIC;
    # }


    cv[i] <- (1 / (k * cv.l)) * sum((forecast.new - t(data.org[cv.index, ]))^2)
    # cv[i] <- BIC.temp;
    # temp <- AIC.BIC.CV(residual,phi.hat.temp);
    # cv[i] <- temp$BIC;
    print("====================================")
  }

  lll <- min(which(cv == min(cv)))
  ind.new <- 0
  if (lll < kk) {
    for (i.3 in (lll + 1):(kk)) {
      if (cv[i.3] < (cv[lll] + cv.var[lll])) {
        ind.new <- ind.new + 1
      }
    }
  }
  lll <- lll + ind.new
  phi.hat.full <- phi.final[[lll]]
  print("CV")
  print(cv)
  print("CV.VAR")
  print(cv.var)

  # test <- var.break.fit(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  # phi.hat.full <- test$phi.hat;
  # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
  # temp.lam <- ll;
  # # temp.lam <- quantile(phi.hat.full,ll)
  # ll <- c(0);
  # brk.points.list <- vector("list",length(ll));
  #
  # for(j in 1:length(ll)){
  #
  #   phi.hat <- phi.hat.full;
  #   # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
  #
  #   n <- T - p;
  #   m.hat <- 0; brk.points <- rep(0,n);
  #
  #   for (i in 2:n)
  #   {
  #     if ( sum((phi.hat[,((i-1)*k*p+1):(i*k*p)] )^2 ) != 0   ){
  #       m.hat <- m.hat + 1; brk.points[m.hat] <- i;
  #     }
  #   }
  #
  #   loc <- rep(0,m.hat);
  #   brk.points <- brk.points[1:m.hat];
  #   brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
  #   m.hat <- length(brk.points);
  #   if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
  #   loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
  #   brk.points.list[[j]] <- brk.points;
  # }


  return(list(brk.points = brk.points.final[[lll]], cv = cv, cv.final = lambda[lll]))
}
