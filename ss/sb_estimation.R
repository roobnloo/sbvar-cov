# sb_estimation.R
# Core VAR estimation routines for the Safikhani-Shojaie method:
#   shvar.fit       — single-segment VAR fitting (many penalty methods)
#   var.lasso.brk   — LASSO VAR for a single segment (used in BEA)
#   var.break.fit   — full TV-penalised VAR (Stage 1 workhorse)
#   AIC.BIC         — information criteria
#   AIC.BIC.CV      — cross-validated information criteria

source(file.path("ss", "sb_utils.R"))

shvar.fit <- function(method, data, weight = NULL, lambda, p, T.1, T.2, max.iteration = 1000, tol = 10^(-4)) {
  method.full <- c(
    "VAR", "LASSO", "HVARC", "HVAROO", "HVARELEM", "SLASSO", "SHVARC", "SHVAROO", "SHVARELEM", "SSHVARC",
    "SSHVAROO", "SSHVARELEM", "SSLASSO", "DHVAR", "DHVARC", "SSDHVAR", "SSDHVARC"
  )
  if (!(method %in% method.full)) {
    print("ERROR")
    break
  }

  k <- length(data[1, ])
  T <- length(data[, 1])
  iter <- matrix(0, k, length(lambda))
  phi.hat <- matrix(0, k, k * p)
  phi.hat.fista <- matrix(0, max.iteration, k * p)
  pred.error <- rep(0, length(lambda))
  phi.hat.temp <- matrix(0, k, k * p * length(lambda))
  Y <- as.matrix(t(data))
  Y <- as.matrix(Y[, (p + 1):T.1])
  # Y <- Y%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  Z <- matrix(0, k * p, T.1 - p)
  # Z <- Z%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  for (i in 1:(T.1 - p)) {
    for (j in 1:p) {
      Z[((j - 1) * k + 1):(j * k), i] <- t(data[i + p - j, ])
    }
  }
  step.size <- 1 / (max(svd(Z)$d))^2
  # Z.tilde <- vector("list",k);
  # for ( i in 1:k){
  #   Z.new <- matrix(0,k*p,T.1-p);
  #   for (j in 1:(T.1-p)){Z.new[,j] <- Z[,j]/(rep(1,k*p)+rep(weight[i,],p))}
  #   Z.tilde[[i]] <- Z.new;
  # }
  # step.size.new <- rep(0,k)
  # for (i in 1:k){step.size.new[i] <- 1/( max(svd(Z.tilde[[i]])$d)  )^2; }
  #

  if (method == "VAR") {
    ll <- 1
    for (ii in 1:k) {
      l <- 2
      while (l < max.iteration) {
        l <- l + 1
        phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
        phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
        if (max(abs(phi.new - phi.temp)) < tol) {
          break
        }
        if (max(abs(phi.new - phi.temp)) > tol) {
          phi.hat.fista[l, ] <- phi.new
        }
      }
      iter[ii, ll] <- l
      phi.hat[ii, ] <- phi.new
    }
    ll.final <- 1
  }

  # if ( method == 'VAR'){
  #   ll <- 1;
  #   for ( ii in 1:k){
  #     phi.temp <- matrix(0,1,k*p); l <- 0;
  #     while( l < max.iteration){
  #       l <- l+1;
  #       phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
  #       if ( max(abs(phi.new - phi.temp)) < tol) {break;}
  #       if (max(abs(phi.new - phi.temp)) > tol ) {
  #         phi.temp <- phi.new; }
  #     }
  #     iter[ii,ll] <- l;
  #     phi.hat[ii,] <- phi.new;
  #   }
  #   ll.final <- 1;
  # }

  # else if (method == 'LASSO'){
  #   for (ll in 1:length(lambda)){
  #     for ( ii in 1:k){
  #       phi.temp <- matrix(0,1,k*p); l <- 0;
  #       while( l < max.iteration){
  #         l <- l+1;
  #         # print(k)
  #         phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
  #         phi.new <- soft(phi.new,rep(0,k*p),lambda[ll]);
  #         if ( max(abs(phi.new - phi.temp)) < tol) {break;}
  #         if (max(abs(phi.new - phi.temp)) > tol ) {
  #           # print("ERROR"); print(max(abs(phi.new - phi.temp)));
  #           phi.temp <- phi.new; }
  #       }
  #       # print("l="); print(l)
  #       iter[ii,ll] <- l;
  #       phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
  #     }
  #     #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
  #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
  #     pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
  #   }
  #   ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
  #   phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
  # }

  else if (method == "LASSO") {
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        l <- 2
        while (l < max.iteration) {
          l <- l + 1
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
          phi.new <- soft(phi.new, rep(0, k * p), lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.1 - p)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, p + jjj, k, 1))
      res <- t(data[(p + 1):(T.1), ]) - forecast
      temp <- AIC.BIC(res, phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)])
      pred.error[ll] <- temp$BIC
      # forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1+jjj-1,k,1)  )
      # pred.error[ll] <- (1/(k*(T.2-T.1)))*sum((t(data[(T.1+1):(T.2),])-forecast)^2);
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "SLASSO") {
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        l <- 2
        while (l < max.iteration) {
          l <- l + 1
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
          # phi.new <- soft(phi.new,rep(weight[ii,],p)-rep(1,k*p),lambda[ll]);
          phi.new <- soft(phi.new, rep(weight[ii, ], p), lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "SSLASSO") {
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        l <- 2
        Z.new <- Z.tilde[[ii]]
        while (l < max.iteration) {
          l <- l + 1
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size.new[ii] * t(Z.new %*% (t(Y[ii, ] - phi.temp %*% Z.new)))
          # phi.new <- soft(phi.new,rep(weight[ii,],p)-rep(1,k*p),lambda[ll]);
          phi.new <- soft(phi.new, rep(0, k * p), lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new / (rep(1, k * p) + rep(weight[ii, ], p))
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "HVARC") {
    group.lag <- matrix(0, p, p * k)
    for (l in 1:p) {
      group.lag[l, (((p - l) * k + 1):(p * k))] <- c((((p - l) * k + 1):(p * k)))
    }
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        l <- 2
        while (l < max.iteration) {
          l <- l + 1
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
          phi.new <- soft.group(phi.new, rep(0, k * p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "SHVARC") {
    group.lag <- matrix(0, p, p * k)
    for (l in 1:p) {
      group.lag[l, (((p - l) * k + 1):(p * k))] <- c((((p - l) * k + 1):(p * k)))
    }
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        l <- 2
        while (l < max.iteration) {
          l <- l + 1
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
          phi.new <- soft.group(phi.new, rep(weight[ii, ], p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "SSHVARC") {
    group.lag <- matrix(0, p, p * k)
    for (l in 1:p) {
      group.lag[l, (((p - l) * k + 1):(p * k))] <- c((((p - l) * k + 1):(p * k)))
    }
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        l <- 2
        Z.new <- Z.tilde[[ii]]
        while (l < max.iteration) {
          l <- l + 1
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size.new[ii] * t(Z.new %*% (t(Y[ii, ] - phi.temp %*% Z.new)))
          phi.new <- soft.group(phi.new, rep(0, k * p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new / (rep(1, k * p) + rep(weight[ii, ], p))
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "HVAROO") {
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        ## OO ###
        group.lag <- matrix(0, 2 * p, p * k)
        for (l.new in 1:(2 * p)) {
          if ((l.new %% 2) == 1) {
            l.temp <- floor(l.new / 2) + 1
            lag.temp <- c(((p - l.temp) * k + 1):((p - l.temp + 1) * k))
            lag.temp <- lag.temp[-c(ii)]
            group.lag[l.new, 1:(k - 1)] <- lag.temp
            l.temp <- floor(l.new / 2)
            if (l.temp != 0) {
              group.lag[l.new, (((p - l.temp) * k + 1):(p * k))] <- c((((p - l.temp) * k + 1):(p * k)))
            }
          }
          if ((l.new %% 2) == 0) {
            l.temp <- l.new / 2
            group.lag[l.new, (((p - l.temp) * k + 1):(p * k))] <- c((((p - l.temp) * k + 1):(p * k)))
          }
        }
        l <- 2
        while (l < max.iteration) {
          l <- l + 1
          # print(k)
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
          phi.new <- soft.group(phi.new, rep(0, k * p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "SHVAROO") {
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        ## OO ###
        group.lag <- matrix(0, 2 * p, p * k)
        for (l.new in 1:(2 * p)) {
          if ((l.new %% 2) == 1) {
            l.temp <- floor(l.new / 2) + 1
            lag.temp <- c(((p - l.temp) * k + 1):((p - l.temp + 1) * k))
            lag.temp <- lag.temp[-c(ii)]
            group.lag[l.new, 1:(k - 1)] <- lag.temp
            l.temp <- floor(l.new / 2)
            if (l.temp != 0) {
              group.lag[l.new, (((p - l.temp) * k + 1):(p * k))] <- c((((p - l.temp) * k + 1):(p * k)))
            }
          }
          if ((l.new %% 2) == 0) {
            l.temp <- l.new / 2
            group.lag[l.new, (((p - l.temp) * k + 1):(p * k))] <- c((((p - l.temp) * k + 1):(p * k)))
          }
        }
        l <- 2
        while (l < max.iteration) {
          l <- l + 1
          # print(k)
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
          phi.new <- soft.group(phi.new, rep(weight[ii, ], p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "SSHVAROO") {
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        ## OO ###
        group.lag <- matrix(0, 2 * p, p * k)
        for (l.new in 1:(2 * p)) {
          if ((l.new %% 2) == 1) {
            l.temp <- floor(l.new / 2) + 1
            lag.temp <- c(((p - l.temp) * k + 1):((p - l.temp + 1) * k))
            lag.temp <- lag.temp[-c(ii)]
            group.lag[l.new, 1:(k - 1)] <- lag.temp
            l.temp <- floor(l.new / 2)
            if (l.temp != 0) {
              group.lag[l.new, (((p - l.temp) * k + 1):(p * k))] <- c((((p - l.temp) * k + 1):(p * k)))
            }
          }
          if ((l.new %% 2) == 0) {
            l.temp <- l.new / 2
            group.lag[l.new, (((p - l.temp) * k + 1):(p * k))] <- c((((p - l.temp) * k + 1):(p * k)))
          }
        }
        l <- 2
        Z.new <- Z.tilde[[ii]]
        while (l < max.iteration) {
          l <- l + 1
          # print(k)
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size.new[ii] * t(Z.new %*% (t(Y[ii, ] - phi.temp %*% Z.new)))
          phi.new <- soft.group(phi.new, rep(0, k * p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new / (rep(1, k * p) + rep(weight[ii, ], p))
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "HVARELEM") {
    ## ELEM ###
    group.lag <- matrix(0, k * p, p)
    for (j in 1:k) {
      for (l in 1:p) {
        temp <- j + seq(0, ((p - 1) * k), k)
        group.lag[((j - 1) * p + l), (1:l)] <- temp[(p - l + 1):(p)]
      }
    }
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        l <- 2
        while (l < max.iteration) {
          l <- l + 1
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
          phi.new <- soft.group(phi.new, rep(0, k * p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "SHVARELEM") {
    ## ELEM ###
    group.lag <- matrix(0, k * p, p)
    for (j in 1:k) {
      for (l in 1:p) {
        temp <- j + seq(0, ((p - 1) * k), k)
        group.lag[((j - 1) * p + l), (1:l)] <- temp[(p - l + 1):(p)]
      }
    }
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        l <- 2
        while (l < max.iteration) {
          l <- l + 1
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
          phi.new <- soft.group(phi.new, rep(weight[ii, ], p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "SSHVARELEM") {
    ## ELEM ###
    group.lag <- matrix(0, k * p, p)
    for (j in 1:k) {
      for (l in 1:p) {
        temp <- j + seq(0, ((p - 1) * k), k)
        group.lag[((j - 1) * p + l), (1:l)] <- temp[(p - l + 1):(p)]
      }
    }
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        l <- 2
        Z.new <- Z.tilde[[ii]]
        while (l < max.iteration) {
          l <- l + 1
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size.new[ii] * t(Z.new %*% (t(Y[ii, ] - phi.temp %*% Z.new)))
          phi.new <- soft.group(phi.new, rep(0, k * p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new / (rep(1, k * p) + rep(weight[ii, ], p))
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "DHVAR") {
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        ## OO ###
        group.lag <- matrix(0, p * k, k)
        for (l.new in 1:(p)) {
          for (j1 in 1:k) {
            lag.temp <- rank(weight[ii, ], ties.method = "first")
            group.lag[((l.new - 1) * k + j1), (1):(j1)] <- rep((l.new - 1) * k, j1) + lag.temp[(k - j1 + 1):(k)]
          }
        }

        l <- 2
        while (l < max.iteration) {
          l <- l + 1
          # print(k)
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
          phi.new <- soft.group(phi.new, rep(0, k * p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "DHVARC") {
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        ## OO ###
        lag.temp <- rank(weight[ii, ], ties.method = "first")
        group.lag <- matrix(0, k, p * k)
        for (j1 in 1:k) {
          for (l.new in 1:(p)) {
            group.lag[(j1), ((l.new - 1) * k + 1):((l.new - 1) * k + j1)] <- rep((l.new - 1) * k, j1) + lag.temp[(k - j1 + 1):(k)]
          }
        }

        l <- 2
        while (l < max.iteration) {
          l <- l + 1
          # print(k)
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
          phi.new <- soft.group(phi.new, rep(0, k * p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "SSDHVAR") {
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        ## OO ###
        group.lag <- matrix(0, p * k, k)
        for (l.new in 1:(p)) {
          for (j1 in 1:k) {
            lag.temp <- rank(weight[ii, ], ties.method = "first")
            group.lag[((l.new - 1) * k + j1), (1):(j1)] <- rep((l.new - 1) * k, j1) + lag.temp[(k - j1 + 1):(k)]
          }
        }

        l <- 2
        Z.new <- Z.tilde[[ii]]
        while (l < max.iteration) {
          l <- l + 1
          # print(k)
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size.new[ii] * t(Z.new %*% (t(Y[ii, ] - phi.temp %*% Z.new)))
          phi.new <- soft.group(phi.new, rep(0, k * p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new / (rep(1, k * p) + rep(weight[ii, ], p))
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  } else if (method == "SSDHVARC") {
    for (ll in 1:length(lambda)) {
      for (ii in 1:k) {
        ## OO ###
        lag.temp <- rank(weight[ii, ], ties.method = "first")
        group.lag <- matrix(0, k, p * k)
        for (j1 in 1:k) {
          for (l.new in 1:(p)) {
            group.lag[(j1), ((l.new - 1) * k + 1):((l.new - 1) * k + j1)] <- rep((l.new - 1) * k, j1) + lag.temp[(k - j1 + 1):(k)]
          }
        }

        l <- 2
        Z.new <- Z.tilde[[ii]]
        while (l < max.iteration) {
          l <- l + 1
          # print(k)
          phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
          phi.new <- phi.temp + step.size.new[ii] * t(Z.new %*% (t(Y[ii, ] - phi.temp %*% Z.new)))
          phi.new <- soft.group(phi.new, rep(0, k * p), group.lag, lambda[ll])
          if (max(abs(phi.new - phi.temp)) < tol) {
            break
          }
          if (max(abs(phi.new - phi.temp)) > tol) {
            # print("ERROR"); print(max(abs(phi.new - phi.temp)));
            phi.hat.fista[l, ] <- phi.new
          }
        }
        # print("l="); print(l)
        iter[ii, ll] <- l
        phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new / (rep(1, k * p) + rep(weight[ii, ], p))
      }
      #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- sapply(c(1:(T.2 - T.1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, T.1 + jjj - 1, k, 1))
      pred.error[ll] <- (1 / (k * (T.2 - T.1))) * sum((t(data[(T.1 + 1):(T.2), ]) - forecast)^2)
    }
    ll.final <- which(pred.error == min(pred.error))
    ll.final <- min(ll.final)
    phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]
  }


  return(list(phi.hat = phi.hat, iter = iter, pred.error = pred.error, tune.final = lambda[ll.final]))
}

var.lasso.brk <- function(data = data.temp, weight = NULL, lambda, p, max.iteration = 1000, tol = 10^(-4)) {
  k <- length(data[1, ])
  T <- length(data[, 1])
  T.1 <- T
  iter <- matrix(0, k, length(lambda))
  phi.hat <- matrix(0, k, k * p)
  phi.hat.fista <- matrix(0, max.iteration, k * p)
  pred.error <- rep(0, length(lambda))
  phi.hat.temp <- matrix(0, k, k * p * length(lambda))
  Y <- as.matrix(t(data))
  Y <- as.matrix(Y[, (p + 1):T.1])
  # Y <- Y%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  Z <- matrix(0, k * p, T.1 - p)
  # Z <- Z%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  for (i in 1:(T.1 - p)) {
    for (j in 1:p) {
      Z[((j - 1) * k + 1):(j * k), i] <- t(data[i + p - j, ])
    }
  }
  step.size <- 1 / (max(svd(Z)$d))^2
  # Z.tilde <- vector("list",k);
  # for ( i in 1:k){
  #   Z.new <- matrix(0,k*p,T.1-p);
  #   for (j in 1:(T.1-p)){Z.new[,j] <- Z[,j]/(rep(1,k*p)+rep(weight[i,],p))}
  #   Z.tilde[[i]] <- Z.new;
  # }
  # step.size.new <- rep(0,k)
  # for (i in 1:k){step.size.new[i] <- 1/( max(svd(Z.tilde[[i]])$d)  )^2; }





  for (ll in 1:length(lambda)) {
    for (ii in 1:k) {
      l <- 2
      while (l < max.iteration) {
        l <- l + 1
        phi.temp <- phi.hat.fista[l - 1, ] + ((l - 2) / (l + 1)) * (phi.hat.fista[l - 1, ] - phi.hat.fista[l - 2, ])
        phi.new <- phi.temp + step.size * t(Z %*% (t(Y[ii, ] - phi.temp %*% Z)))
        phi.new <- soft(phi.new, rep(0, k * p), lambda[ll])
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
      phi.hat.temp[ii, ((ll - 1) * k * p + 1):(ll * k * p)] <- phi.new
    }
    #     forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
    forecast <- matrix(0, k, T.1 - p)
    forecast <- sapply(c((p):(T.1 - 1)), function(jjj) pred(t(data), phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)], p, jjj, k, 1))
    pred.error[ll] <- sum((t(data[(p + 1):(T.1), ]) - forecast)^2) + lambda[ll] * sum(abs(phi.hat.temp[, ((ll - 1) * k * p + 1):(ll * k * p)]))
  }
  ll.final <- 1
  phi.hat <- phi.hat.temp[, ((ll.final - 1) * k * p + 1):(ll.final * k * p)]







  return(list(phi.hat = phi.hat, iter = iter, pred.error = pred.error, tune.final = lambda[ll.final]))
}

var.break.fit <- function(method, data, weight = NULL, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), initial.phi = NULL) {
  method.full <- c("LASSO")
  if (!(method %in% method.full)) {
    print("ERROR")
    break
  }

  k <- length(data[1, ])
  T <- length(data[, 1])
  iter <- matrix(0, k, length(lambda))
  n <- T - p
  Y <- matrix(0, k * p, n)
  for (i in p:(T - 1)) {
    Y[, (i - p + 1)] <- sapply(c(1:p), function(jjj) data[i - jjj + 1, ])
  }

  C <- vector("list", n)
  # for(jjj in 1:n){C[[jjj]] <- as.matrix(Y[,jjj],k*p,1)%*%t(as.matrix(data[jjj+p,],1,k)); }
  C <- lapply(c(1:n), function(jjj) as.matrix(Y[, jjj], k * p, 1) %*% t(as.matrix(data[jjj + p, ], 1, k)))
  C.sum <- matrix(0, k * p * n, k)
  C.sum[1:(k * p), ] <- C[[1]]
  for (i in 2:n) {
    C.sum[((i - 1) * k * p + 1):(i * k * p), ] <- C.sum[((i - 2) * k * p + 1):((i - 1) * k * p), ] + C[[i]]
  }
  C.sum.new <- matrix(0, k * p * n, k)
  C.sum.new[1:(k * p), ] <- C.sum[((n - 1) * k * p + 1):(n * k * p), ]
  for (i in 2:n) {
    C.sum.new[((i - 1) * k * p + 1):(i * k * p), ] <- C.sum[((n - 1) * k * p + 1):(n * k * p), ] - C.sum[((i - 2) * k * p + 1):((i - 1) * k * p), ]
  }

  D <- vector("list", n)
  D <- lapply(c(1:n), function(jjj) as.matrix(Y[, jjj], k * p, 1) %*% t(as.matrix(Y[, jjj], k * p, 1)))
  D.sum <- matrix(0, k * p * n, k * p)
  D.sum[1:(k * p), ] <- D[[1]]
  for (i in 2:n) {
    D.sum[((i - 1) * k * p + 1):(i * k * p), ] <- D.sum[((i - 2) * k * p + 1):((i - 1) * k * p), ] + D[[i]]
  }
  D.sum.new <- matrix(0, k * p * n, k * p)
  D.sum.new[1:(k * p), ] <- D.sum[((n - 1) * k * p + 1):(n * k * p), ]
  # D.sum.new[(k*p+1):(n*k*p),] <- sapply(c(2:n), function(jjj) D.sum[((n-1)*k*p+1):(n*k*p),] - D.sum[((jjj-2)*k*p+1):((jjj-1)*k*p),]    )
  for (i in 2:n) {
    D.sum.new[((i - 1) * k * p + 1):(i * k * p), ] <- D.sum[((n - 1) * k * p + 1):(n * k * p), ] - D.sum[((i - 2) * k * p + 1):((i - 1) * k * p), ]
  }
  D.sum.new.inv <- matrix(0, k * p * n, k * p)
  for (jjj in 1:n) {
    D.sum.new.inv[((jjj - 1) * k * p + 1):(jjj * k * p), ] <- solve(D.sum.new[((jjj - 1) * k * p + 1):(jjj * k * p), ] + (tol) * diag(k * p))
  }
  # D.sum.new.inv <- sapply(c(1:n), function(jjj)  solve(D.sum.new[((jjj-1)*k*p+1):(jjj*k*p),])   )


  phi.hat <- matrix(0, k, k * p * n)
  if (!is.null(initial.phi)) {
    phi.hat <- initial.phi
  }
  active <- rep(0, n)
  active <- sapply(c(1:n), function(jjj) {
    if (sum((phi.hat[, ((jjj - 1) * k * p + 1):(jjj * k * p)])^2) != 0) {
      jjj
    } else {
      0
    }
  })
  active <- active[which(active != 0)]
  phi.new <- matrix(0, k, k * p * n)
  # phi.temp <- phi.hat;



  # X <- matrix(0,n,n*k*p);
  # for ( i in 1:n){
  #   for(j in 1:i){
  #     X[i,((j-1)*p*k+1):(j*p*k)] <-  t(as.matrix(Y[,i],k*p,1));
  #   }
  # }
  # step.size <- 1/( max(svd(X)$d)  )^2;
  # tol <- 0.5*step.size;
  # step.size <- (0.25)*tol;
  # step.size <- 0.01;

  # print(step.size)




  if (method == "LASSO") {
    for (ll in 1:length(lambda)) {
      l <- 2

      while (l < max.iteration) {
        l <- l + 1
        phi.compare <- phi.hat
        # Pre-compute D_phi_right = sum_jjj D_jjj %*% phi[jjj]^T once before the sweep.
        # phi_left_sum accumulates updated phi blocks as ii advances (Gauss-Seidel).
        # Together they let us compute E_net(ii) in O(k^2 p^2) instead of O(n k^2 p^2),
        # reducing the per-while-iteration cost from O(n^2) to O(n).
        phi_left_sum <- matrix(0, k * p, k)
        D_phi_right <- matrix(0, k * p, k)
        for (jjj in 1:n) {
          jdx <- ((jjj - 1) * k * p + 1):(jjj * k * p)
          D_phi_right <- D_phi_right + D.sum.new[jdx, ] %*% t(phi.hat[, jdx])
        }

        for (ii in 1:n) {
          idx <- ((ii - 1) * k * p + 1):(ii * k * p)
          D_ii <- D.sum.new[idx, ]
          phi_old_ii <- t(phi.hat[, idx]) # (k*p) x k

          E_net <- D_ii %*% phi_left_sum + D_phi_right - D_ii %*% phi_old_ii
          S <- C.sum.new[idx, ] - E_net
          S <- soft.full(S, lambda[ll], k, p, n)

          phi.temp <- t(D.sum.new.inv[idx, ] %*% S) # k x (k*p)
          phi.hat[, idx] <- phi.temp
          phi.new[, idx] <- phi.temp

          # Incremental updates for next ii
          phi_left_sum <- phi_left_sum + t(phi.temp)
          D_phi_right <- D_phi_right - D_ii %*% phi_old_ii
        }


        if (max(abs(phi.new - phi.compare)) < tol) {
          # print(max(abs(phi.new - phi.compare)))
          break
        }
        if (max(abs(phi.new - phi.compare)) > tol) {
          phi.hat <- phi.new
          # print(max(abs(phi.new - phi.compare)))
        }
      }
      # print("l="); print(l);
    }
  }

  # stopCluster(cl)

  # for(i in 1:(n*p*k)){
  #   for( j in 1:k){
  #     if ( abs (phi.hat[j,i] <= tol) ){phi.hat[j,i] <- 0;}
  #   }
  # }
  return(list(phi.hat = phi.hat, iter = iter))
}

AIC.BIC <- function(residual, phi) {
  k <- length(phi[, 1])
  k.lam <- length(phi[1, ])
  T.new <- length(residual[1, ])
  count <- 0
  for (i in 1:k) {
    for (j in 1:k.lam) {
      if (phi[i, j] != 0) {
        count <- count + 1
      }
    }
  }
  # for ( i in 1:T.new){residual[,i] <- residual[,i] - as.matrix(rowMeans(residual))}
  # sigma.hat <- (1/(T.new-1))*(residual%*%t(residual)) + 10^(-10)*diag(k);
  sigma.hat <- 0 * diag(k)
  for (i in 1:T.new) {
    sigma.hat <- sigma.hat + residual[, i] %*% t(residual[, i])
  }
  sigma.hat <- (1 / (T.new)) * sigma.hat
  ee.temp <- min(eigen(sigma.hat)$values)
  if (ee.temp <= 0) {
    sigma.hat <- sigma.hat + (2.0) * (abs(ee.temp) + 10^(-4)) * diag(k)
  }
  # sigma.hat <- (1/(T.new-1))*(t(residual)%*%(residual)) + 0*10^(-10)*diag(T.new);
  log.det <- log(det(sigma.hat) + 0 * 10^(-10))
  # print(det(sigma.hat));
  return(list(AIC = log.det + 2 * count / T.new, BIC = log.det + log(T.new) * count / T.new))
}

AIC.BIC.CV <- function(residual, phi) {
  k <- length(phi[, 1])
  k.lam <- length(phi[1, ])
  T.new <- length(residual[1, ])
  count <- 0
  for (i in 1:k) {
    for (j in 1:k.lam) {
      if (phi[i, j] != 0) {
        count <- count + 1
      }
    }
  }
  for (i in 1:T.new) {
    residual[, i] <- residual[, i] - as.matrix(rowMeans(residual))
  }
  sigma.hat <- (1 / (T.new - 1)) * (residual %*% t(residual))
  log.det <- log(det(sigma.hat))
  if (abs(det(sigma.hat)) < 10^(-8)) {
    log.det <- -8^(1)
  }
  return(list(AIC = log.det + 2 * count / T.new, BIC = log.det + log(T.new) * count / T.new))
}
