#===============================================================================
# This is an R function named 'twosremrFun'
# which computes a two-stage random-effects meta-regression estimator for
# the two-stage random-effects meta-analysis estimators
# described in Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. (2024).
#===============================================================================

twosremrFun <- function(Y,X,SE,ID,rho){

  N   <- length(Y)
  IDs <- unique(ID)
  I   <- length(IDs)
  if(length(rho)!=I){rho <- matrix(rho[1],I,1)}

  # Observations per group:
  J <- matrix(0,I,1)
  for(i in 1:I){
    J[i] <- sum(ID==IDs[i])
  }

  # Dimension of X:
  K <- length(X[1,])

  # Y and SE as a IxJmax matrices (y and se):
  {
    y  <- matrix(0,I,max(J))
    se <- matrix(0,I,max(J))
    for(i in 1:I){
      y[i,1:J[i]]  <- Y[which(ID==IDs[i])]
      se[i,1:J[i]] <- SE[which(ID==IDs[i])]
    }
  }

  # Use OLS to initialize bhat:
  XpXi     <- solve(t(X)%*%X)
  bhat.OLS <- XpXi %*% t(X) %*% Y
  ehat     <- X %*% bhat.OLS - Y

  # Cluster- and heteroskedastic-robust VCV matrix
  {
    # (Cameron and Miller 2013)
    BB <- matrix(0,K,K)
    for(i in 1:I){
      rows  <- which(ID==IDs[i])
      nrows <- length(rows)
      Xi    <- matrix(X[rows,],nrows,K)
      ei    <- ehat[rows,]
      BB    <- BB + t(Xi) %*% (ei %*% t(ei)) %*% Xi
    }
    VCVr      <- XpXi %*% BB %*% XpXi
    se.ro.OLS <- as.matrix(sqrt(diag(VCVr)),K,1)
  }

  # UNCONSTRAINED SIG.MUHAT:
  done <- 0
  bhat <- bhat.OLS
  while(done==0){

    # Estimate sig.mu.UNC for each group:
    {
      sig.muhat.UNC <- matrix(0,I,1)
      for(i in 1:I){

        if(J[i]>1){ # Estimate sig.mu only for groups with >1 observation

          CV <- se[i,1:J[i]] %*% t(se[i,1:J[i]])

          xbeta <- X[which(ID==IDs[i]),] %*% bhat

          sig.mu2hat <- var(y[i,1:J[i]]) -
            1/J[i]*(sum(diag(CV)-rho[i]/(J[i]^2)*
                          (sum(CV)-sum(diag(CV))))) -
            var(xbeta)

          # Only estimate sig.mu if estimated sig.mu^2 is positive:
          if(sig.mu2hat>0){sig.muhat.UNC[i] <- sqrt(sig.mu2hat)}

        }
        sig.muhat.UNC[i] <- sig.muhat.UNC[i] * (J[i]/(J[i]-1))
      }

      # For all groups with one observation, impute sig.mu as sample size
      # weighted average of all other group's sig.muhat':
      sig.muhat.UNC[which(J==1)] <- sum( sig.muhat.UNC[which(J>1)] * J[which(J>1)]
                                         / sum(J[which(J>1)]) )
    }

    # Calculate ghat:
    {
      outs <- ghFun(se,0,sig.muhat.UNC,rho)
      ghat <- outs[[1]]
    }

    # Estimate sig.eta.UNC:
    {
      # Calculate yhati and Vi for each group w/ more than one observation:
      yhati <- matrix(0,I,1)
      Vi    <- matrix(0,I,1)
      for(i in 1:I){
        if(J[i]>1){
          yhati[i] <- sum(y[i,1:J[i]]*ghat[i,1:J[i]])
          for(j in 1:J[i]){
            Vi[i] <- Vi[i] + ghat[i,j]^2*(sig.muhat.UNC[i]^2+se[i,j]^2)
            for(k in 1:J[i]){
              if(k!=j){
                Vi[i] <- Vi[i] + rho[i]*ghat[i,j]*ghat[i,k]*se[i,j]*se[i,k]
              }
            }
          }
        }
        if(J[i]==1){yhati[i] <- y[i,1]}
      }

      # Use yhati's and Vi's to calculate sig.etahat [Equation 49]:
      part1 <- 1/(I-1)*sum((yhati-mean(yhati))^2)
      part2 <- 1/sum(J>1)*sum(Vi[J>1])
      if(part1>part2){sig.etahat.UNC <- sqrt(part1-part2)}
      else{sig.etahat.UNC <- 0}
    }

    # Calculate ghat and hhat and weights for WLS:
    {
      outs <- ghFun(se,sig.etahat.UNC,sig.muhat.UNC,rho)
      ghat <- outs[[1]]
      hhat <- outs[[2]]
      W    <- matrix(0,N,1)
      ij   <- 0
      for(i in 1:I){
        for(j in 1:J[i]){
          ij <- ij + 1
          W[ij] <- ghat[i,j]*hhat[i]
        }
      }
    }

    # Estimate beta using WLS:
    {
      YW    <- Y * sqrt(W)
      XW    <- X * matrix(sqrt(W),N,length(X[1,]))
      XpXi  <- solve(t(XW)%*%XW)
      bhat2 <- XpXi %*% t(XW) %*% YW
      ewhat <- XW %*% bhat2 - YW
    }

    # Check for convergence:
    {
      if( max(abs((bhat2[1:K]-bhat[1:K])/bhat[1:K])) < 0.00001 ){
        done <- 1
      }
      bhat <- bhat + .1*(bhat2-bhat)
    }
  }

  bhat.UNC   <- bhat
  ehat       <- X %*% bhat.UNC - Y
  R2.UNC     <- 1-sum(ewhat^2)/sum((YW-mean(YW))^2)
  R2.adj.UNC <- 1-sum(ewhat^2)/(N-K)/(sum((YW-mean(YW))^2)/(N-1))

  # UNC Leave-one-out cross-validation errors:
  {
    XWpXWi    <- solve(t(XW)%*%XW)
    H         <- XW %*% XWpXWi %*% t(XW)
    h         <- diag(H)
    e.cv.UNC  <- ewhat/(1-h)
    R2.cv.UNC <- 1-sum(e.cv.UNC^2)/sum((YW-mean(YW))^2)
  }

  # Cluster- and heteroskedastic-robust VCV matrix
  {
    # (Cameron and Miller 2013)
    BB <- matrix(0,K,K)
    for(i in 1:I){
      rows  <- which(ID==IDs[i])
      nrows <- length(rows)
      Xi    <- matrix(XW[rows,],nrows,K)
      ei    <- ewhat[rows,]
      BB    <- BB + t(Xi) %*% (ei %*% t(ei)) %*% Xi
    }
    VCVr      <- XpXi %*% BB %*% XpXi
    se.ro.UNC <- as.matrix(sqrt(diag(VCVr)),K,1)
  }

  # CONSTRAINED SIG.MUHAT:
  done <- 0
  bhat <- bhat.OLS
  while(done==0){

    # Estimate sig.mu.CON (same for all groups)
    {
      sig.muhat.CON <- matrix(mean(sig.muhat.UNC[which(J>1)]),I,1)
    }

    # Calculate ghat:
    {
      outs<- ghFun(se,0,sig.muhat.CON,rho)
      ghat <- outs[[1]]
    }

    # Estimate sig.eta.CON:
    {
      # Calculate yhati and Vi for each group w/ more than one observation:
      yhati <- matrix(0,I,1)
      Vi    <- matrix(0,I,1)
      for(i in 1:I){
        if(J[i]>1){
          yhati[i] <- sum(y[i,1:J[i]]*ghat[i,1:J[i]])
          for(j in 1:J[i]){
            Vi[i] <- Vi[i] + ghat[i,j]^2*(sig.muhat.CON[i]^2+se[i,j]^2)
            for(k in 1:J[i]){
              if(k!=j){
                Vi[i] <- Vi[i] + rho[i]*ghat[i,j]*ghat[i,k]*se[i,j]*se[i,k]
              }
            }
          }
        }
        if(J[i]==1){yhati[i] <- y[i,1]}
      }

      # Use yhati's and Vi's to calculate sig.etahat [Equation 49]:
      part1 <- 1/(I-1)*sum((yhati-mean(yhati))^2)
      part2 <- 1/sum(J>1)*sum(Vi[J>1])
      if(part1>part2){sig.etahat.CON <- sqrt(part1-part2)}
      else{sig.etahat.CON <- 0}
    }

    # Calculate ghat and hhat and weights for WLS:
    {
      outs <- ghFun(se,sig.etahat.CON,sig.muhat.CON,rho)
      ghat <- outs[[1]]
      hhat <- outs[[2]]
      W    <- matrix(0,N,1)
      ij   <- 0
      for(i in 1:I){
        for(j in 1:J[i]){
          ij <- ij + 1
          W[ij] <- ghat[i,j]*hhat[i]
        }
      }
    }

    # Estimate beta using WLS:
    {
      YW    <- Y * sqrt(W)
      XW    <- X * matrix(sqrt(W),N,length(X[1,]))
      XpXi  <- solve(t(XW)%*%XW)
      bhat2 <- XpXi %*% t(XW) %*% YW
      ewhat <- XW %*% bhat2 - YW
    }

    # Check for convergence:
    {
      if(max(abs((bhat2[1:K]-bhat[1:K])/bhat[1:K]))<0.00001){done <- 1}
      bhat <- bhat2
    }
  }

  bhat.CON   <- bhat
  ehat       <- X%*%bhat.CON-Y
  R2.CON     <- 1-sum(ewhat^2)/(sum((YW-mean(YW))^2))
  R2.adj.CON <- 1-sum(ewhat^2)/(N-K)/(sum((YW-mean(YW))^2)/(N-1))

  # log likelihood:
  lnL <- -N/2*log(2*pi)-sum(1/2*log(as.vector(sig.etahat.CON^2)+SE^2))-
    1/2*sum((Y-X%*%bhat.CON)^2/(as.vector(sig.etahat.CON^2)+SE^2))

  # CON Leave-one-out cross-validation errors:
  {
    XWpXWi    <- solve(t(XW)%*%XW)
    H         <- XW %*% XWpXWi %*% t(XW)
    h         <- diag(H)
    e.cv.CON  <- ewhat/(1-h)
    R2.cv.CON <- 1-sum(e.cv.CON^2)/sum((YW-mean(YW))^2)
    }

  # Cluster- and heteroskedastic-robust VCV matrix
  {
    # (Cameron and Miller 2013)
    BB <- matrix(0,K,K)
    for(i in 1:I){
      rows  <- which(ID==IDs[i])
      nrows <- length(rows)
      Xi    <- matrix(XW[rows,],nrows,K)
      ei    <- ewhat[rows,]
      BB    <- BB + t(Xi) %*% (ei %*% t(ei)) %*% Xi
    }
    VCVr      <- XpXi %*% BB %*% XpXi
    se.ro.CON <- as.matrix(sqrt(diag(VCVr)),K,1)
  }

  return(list( bhat.OLS       = bhat.OLS,
               se.ro.OLS      = se.ro.OLS,
               bhat.UNC       = bhat.UNC,
               se.ro.UNC      = se.ro.UNC,
               sig.muhat.UNC  = sig.muhat.UNC,
               sig.etahat.UNC = sig.etahat.UNC,
               R2.UNC         = R2.UNC,
               R2.adj.UNC     = R2.adj.UNC,
               e.cv.UNC       = e.cv.UNC,
               R2.cv.UNC      = R2.cv.UNC,
               bhat.CON       = bhat.CON,
               se.ro.CON      = se.ro.CON,
               sig.muhat.CON  = sig.muhat.CON,
               sig.etahat.CON = sig.etahat.CON,
               R2.CON         = R2.CON,
               R2.adj.CON     = R2.adj.CON,
               e.cv.CON       = e.cv.CON,
               R2.cv.CON      = R2.cv.CON))

}
