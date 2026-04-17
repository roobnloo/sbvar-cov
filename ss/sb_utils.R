# sb_utils.R
# Simulation helpers, proximal operators, and miscellaneous utilities
# extracted from functions_SBDetection.R.
#
# Requires: mvtnorm (for var.sim, var.sim.break)

var.sim <- function (nobs, arlags = NULL, malags = NULL, cnst = NULL, phi = NULL,theta = NULL, skip = 200, sigma) {
  if (!is.matrix(sigma))
    sigma = as.matrix(sigma)
  k = nrow(sigma)
  nT = nobs + skip
  at = rmvnorm(nT, rep(0, k), sigma)
  nar = length(arlags)
  p = 0
  if (nar > 0) {
    arlags = sort(arlags)
    p = arlags[nar]
  }
  q = 0
  nma = length(malags)
  if (nma > 0) {
    malags = sort(malags)
    q = malags[nma]
  }
  ist = max(p, q) + 1
  zt = matrix(0, nT, k)
  if (length(cnst) == 0)
    cnst = rep(0, k)
  for (it in ist:nT) {
    tmp = matrix(at[it, ], 1, k)
    if (nma > 0) {
      for (j in 1:nma) {
        jdx = (j - 1) * k
        thej = theta[, (jdx + 1):(jdx + k)]
        atm = matrix(at[it - malags[j], ], 1, k)
        tmp = tmp - atm %*% t(thej)
      }
    }
    if (nar > 0) {
      for (i in 1:nar) {
        idx = (i - 1) * k
        phj = phi[, (idx + 1):(idx + k)]
        ztm = matrix(zt[it - arlags[i], ], 1, k)
        tmp = tmp + ztm %*% t(phj)
      }
    }
    zt[it, ] = cnst + tmp
  }
  zt = zt[(1 + skip):nT, ]
  at = at[(1 + skip):nT, ]
  VARMAsim <- list(series = zt, noises = at)
}

var.sim.break <- function (nobs, arlags = NULL, malags = NULL, cnst = NULL, phi = NULL,theta = NULL, skip = 200, sigma, brk = nobs+1) {
  if (!is.matrix(sigma))
    sigma = as.matrix(sigma)
  k = nrow(sigma)
  m <- length(brk)
  nT = nobs + skip
  at = rmvnorm(nT, rep(0, k), sigma)
  nar = length(arlags)
  p = 0
  if (nar > 0) {
    arlags = sort(arlags)
    p = arlags[nar]
  }
  q = 0
  nma = length(malags)
  if (nma > 0) {
    malags = sort(malags)
    q = malags[nma]
  }
  ist = max(p, q) + 1
  zt = matrix(0, nT, k)
  if (length(cnst) == 0)
    cnst = rep(0, k)
  if (m == 1){
    for (it in ist:nT) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
  }

  if (m > 1){

    for (it in ist:(skip+brk[1]-1)) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
    for ( mm in 1:(m-1)){
      for (it in (skip+brk[mm]):(skip+brk[mm+1]-1) ) {
        tmp = matrix(at[it, ], 1, k)
        if (nma > 0) {
          for (j in 1:nma) {
            jdx = (j - 1) * k
            thej = theta[, (jdx + 1):(jdx + k)]
            atm = matrix(at[it - malags[j], ], 1, k)
            tmp = tmp - atm %*% t(thej)
          }
        }
        if (nar > 0) {
          for (i in 1:nar) {
            idx = (i - 1) * k
            phj = phi[, ((mm)*p*k+idx + 1):((mm)*p*k+idx + k)]
            ztm = matrix(zt[it - arlags[i], ], 1, k)
            tmp = tmp + ztm %*% t(phj)
          }
        }
        zt[it, ] = cnst + tmp
      }
    }
  }

  zt = zt[(1 + skip):nT, ]
  at = at[(1 + skip):nT, ]
  VARMAsim <- list(series = zt, noises = at)
}

pred <- function(Y,phi,p,T,k,h){
  concat.Y <- matrix(0,k,T+h); concat.Y[,1:T] <- Y[,1:T];
  for ( j in 1:h){
    temp <- matrix(0,k,1);
    for (i in 1:p){temp <- temp +  phi[,((i-1)*k+1):(i*k)]%*%concat.Y[,T+j-i];}
    concat.Y[,T+j] <- temp;
  }
  return(as.matrix(concat.Y[,T+h]))
}

soft <- function(L,weight,lambda){
  for (i in 1:length(L[1,])){
    lambda <- lambda*(1+weight[i])
    if ( L[i] > lambda){L[i] <- L[i] - lambda}
    if ( L[i] < -lambda){L[i] <- L[i] + lambda}
    if ( abs(L[i]) <= lambda){L[i] <- 0}
  }
  return(L)
}

soft.full <- function(L,lambda,k,p,n){

  # for(kk in 1:n){
  #   temp <- L[,((kk-1)*k*p+1):(kk*k*p)];
  #   nrm <- sum(abs(temp))
  #   if ( nrm <= lambda){ L[,((kk-1)*k*p+1):(kk*k*p)] <- matrix(0,k,k*p)  }
  #   if ( nrm > lambda) { L[,((kk-1)*k*p+1):(kk*k*p)] <- L[,((kk-1)*k*p+1):(kk*k*p)] - matrix((lambda/(p*k^2)),k,k*p); }
  # }
  # nrm <- sum(abs(L))
  # if ( nrm <= lambda){ L <- matrix(0,k,k*p)  }
  # if ( nrm > lambda) { L <- L - matrix((lambda/(p*k^2)),k,k*p);  }


  sign(L) * pmax(abs(L) - lambda, 0)
}

soft.group <- function(L,weight,group.lag,lambda){
  L <- L/(1+weight);
  for (i in 1:length(group.lag[,1])){
    temp <- group.lag[i,]; temp <- temp[which(temp!=0)];
    L[temp] <- (max(0,1-lambda/(sqrt(sum((L[temp])^2)))))*L[temp]
  }
  return(L*(1+weight))
}

myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }

  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

  # # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  # ColorRamp <- rgb( seq(0,1,length=256),  # Red
  #                   seq(0,1,length=256),  # Green
  #                   seq(1,0,length=256))  # Blue
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- grey(seq(1, 0, length = 256))
  ColorLevels <- seq(min, max, length=length(ColorRamp))

  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]

  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)

  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")

  layout(1)
}

generate.phi.new <- function(k,p.t,l.bound,max.iter,seed){

  for ( try in 1:max.iter){
    phi <- matrix(0,k,k*p.t)
    set.seed(try+seed)
    base <- matrix(2*(runif((k^2)*p.t)-1/2),k,k*p.t);
    for(i in 1:k){
      for(j in 1:(k*p.t)){
        d <- (abs(abs(i-j)-1)+1);
        if ( abs(base[i,j]/d) > l.bound   ){phi[i,j] <- base[i,j]/d; }
        if ( abs(base[i,j]/d) <= l.bound   ){phi[i,j] <- 0; }
      }
    }

    companion.phi <- matrix(0,k*p.t,k*p.t);
    companion.phi[1:k,] <- phi;
    if( p.t > 1){
      for(i in 1:(p.t-1)){companion.phi[(i*k+1):((i+1)*k),((i-1)*k+1):((i)*k)] <- diag(k)       }
    }
    aaa <- eigen(companion.phi)$values; aaa <- Mod(aaa);
    # print("TRY=="); print(try)
    if(max(aaa) < 1){break;}
  }

  return(list(phi = phi, try = try))
}
