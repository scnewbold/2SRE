yhatFun <- function(ID,Y,SE,rho){

  IDs <- unique(ID)
  I   <- length(IDs)
  J   <- matrix(0,I,1); for(i in 1:I){J[i] <- sum(ID==IDs[i])}

  if(length(rho)!=I){rho <- matrix(rho[1],I,1)}

  # Y and SE as a IxJmax matrices (y and se):
  {
    y  <- matrix(0,I,max(J))
    se <- matrix(0,I,max(J))
    J  <- matrix(0,I,1)
    for(i in 1:I){
      zz <- which(ID==IDs[i])
      J[i] <- length(zz)
      y[i,1:J[i]]  <- Y[zz]
      se[i,1:J[i]] <- SE[zz]
    }
  }

  # Eq 2.1 [sig2.mu.hat]
  {
    sig2.mu <- matrix(0,I,1)
    for(i in 1:I){
      if(J[i]>1){
        sig2.mu[i] <- var(y[i,1:J[i]])
        for(j in 1:J[i]){
          sig2.mu[i] <- sig2.mu[i] - 1/J[i]*se[i,j]^2
          for(k in 1:J[i]){
            sig2.mu[i] <- sig2.mu[i] + 1/J[i]^2*rho[i]*se[i,j]*se[i,k]*(k!=j)
          }
        }
        sig2.mu[i] <- max(sig2.mu[i],0)
        sig2.mu[i] <- sig2.mu[i] * J[i]/(J[i]-1)
      }
    }
    # impute for J=1 using weighted average of sig.mu
    # sig2.mu[J==1] <- sum( sqrt(sig2.mu[J>1]) * J[J>1] / sum(J[J>1]) )^2
    sig2.mu[J==1] <- sum(sig2.mu[J>1]*J[J>1]/sum(J[J>1]))
  }

  # Eq 2.2 [A.hat]
  {
    A <- matrix(0,I,max(J))
    for(i in 1:I){
      for(j in 1:J[i]){
        for(k in 1:J[i]){
          A[i,j] <- A[i,j] +
             (sig2.mu[i] + (1-rho[i])*se[i,j]^2)/(sig2.mu[i]+(1-rho[i])*se[i,k]^2)
        }
        A[i,j] <- 1/A[i,j]
      }
    }
  }

  # Eq 2.3 [B.hat]
  {
    B <- array(0,dim=c(max(J),max(J),I))
    for(i in 1:I){
      for(j in 1:J[i]){
        num <- sum(se[i,1:J[i]]/(sig2.mu[i]+(1-rho[i])*se[i,1:J[i]]^2))
        den <- sum((sig2.mu[i]+(1-rho[i])*se[i,j]^2)/(sig2.mu[i]+(1-rho[i])*se[i,1:J[i]]^2))
        for(k in 1:J[i]){
          B[j,k,i] <- rho[i]*(num/den-se[i,j]/(sig2.mu[i]+(1-rho[i])*se[i,j]^2))*se[i,k]
        }
      }
    }
  }

  # Eq 2.4 [g.hat]
  {
    g <- matrix(0,I,max(J))
    for(i in 1:I){
      g[i,1:J[i]] <- solve(diag(1,J[i])-B[1:J[i],1:J[i],i]) %*% A[i,1:J[i]]
      if(sum(g[i,])!=0){
        g[i,g[i,]<0] <- 0
        g[i,] <- g[i,]/sum(g[i,])
      }
    }
  }

  # Eq 2.5 [yi.hat]
  {
    yihat <- matrix(0,I,1)
    for(i in 1:I){
      yihat[i] <- sum(g[i,]*y[i,])
    }
  }

  # Eq 2.6 [sig2eta.hat]
  {
    V <- matrix(0,I,1)
    for(i in 1:I){
      for(j in 1:J[i]){
        V[i] <- V[i] + g[i,j]^2*(sig2.mu[i]+se[i,j]^2)
        for(k in 1:J[i]){
          V[i] <- V[i] + rho[i]*g[i,j]*g[i,k]*se[i,j]*se[i,k]*(k!=j)
        }
      }
    }
    sig2.eta <- var(yihat) - 1/sum(J>1) * sum(V[J>1])
    sig2.eta <- max(sig2.eta,0)
  }

  # Eq 2.7 [vi.hat]
  {
    v <- matrix(0,I,1)
    for(i in 1:I){
      v[i] <- sig2.eta + V[i]
    }
  }

  # Eq 2.8 [hi.hat]
  {
    h <- (1/v)/sum(1/v)
  }

  # Eq 2.9 [wij.hat]
  {
    w <- matrix(0,I,max(J))
    for(i in 1:I){
      for(j in 1:J[i]){
        w[i,j] <- g[i,j]*h[i]
      }
    }
  }

  # Eq 2.10 [yhat]
  {
    yhat <- sum(h*yihat)
    yhat <- sum(w*y)
  }

  return(list(sig2.mu=sig2.mu,sig2.eta=sig2.eta,yhat=yhat))

}


