#===============================================================================
# ghFun.R:
# This function will compute level 1 and level 2 observation weights for
# the two-stage random-effects meta-analysis estimators
# described in Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. (2023)
#===============================================================================

ghFun <- function(se,sig.eta,sig.mu,rho){

  # Replace any negative g's with 0 and re-scale remaining g's (1=yes,0=no):
  g.scale <- 1

  I    <- nrow(se)
  J.hi <- ncol(se)
  if(J.hi> 1){J <- matrix(rowSums(se>0),I,1)}
  if(J.hi==1){J <- matrix(1,I,1)}


  ## Stage 1: compute within group weights (g) using equation (26) or (30)
  g <- matrix(0,I,J.hi)
  A <- matrix(0,I,J.hi)
  B <- array(0,dim=c(J.hi,J.hi,I))

  for(i in 1:I){
    for(j in 1:J[i]){

      num <- 0; den <- 0
      for(kk in 1:J[i]){
        num <- num + se[i,kk]/(sig.mu[i]^2+(1-rho[i])*se[i,kk]^2)
        den <- den + (sig.mu[i]^2+(1-rho[i])*se[i,j]^2)/
                     (sig.mu[i]^2+(1-rho[i])*se[i,kk]^2)
      }

      for(k in 1:J[i]){
        A[i,j] <- A[i,j] + (sig.mu[i]^2+(1-rho[i])*se[i,j]^2)/
                           (sig.mu[i]^2+(1-rho[i])*se[i,k]^2)

        B[j,k,i] <- rho[i]*(num/den-se[i,j]/(sig.mu[i]^2+(1-rho[i])*se[i,j]^2))*se[i,k]
      }
      A[i,j] <- 1/A[i,j]
    }
    g[i,1:J[i]] <- solve(diag(J[i])-B[1:J[i],1:J[i],i])%*%as.matrix(A[i,1:J[i]])

    if(g.scale==1){
      # If any elements of g<0 replace with 0 and re-scale remaining elements
      if(sum(g[i,1:J[i]]<0)>0){
        g[i,g[i,1:J[i]]<0] <- 0
        g[i,] <- g[i,]/sum(g[i,])
      }
    }

  }

  ## stage 2: compute across group weights (h)
  # find variance of studies first
  V <- matrix(0,I,1)   ## 2nd part of equation (7)
  for(i in 1:I){
    for(j in 1:J[i]){
      V[i] <- V[i] + ( g[i,j]^2*(sig.mu[i]^2+se[i,j]^2) )
      rho.term <- 0
      for(k in 1:J[i]){
        if(k!=j){rho.term <- rho.term + rho[i]*g[i,j]*g[i,k]*se[i,j]*se[i,k]}
      }
      V[i] <- V[i] + rho.term
    }
  }
  varyhati <- sig.eta^2 + V  ## equation(7)

  h <- matrix(0,I,1)   ## equation(43)
  for(i in 1:I){
    h[i] <- 1/varyhati[i] / sum( 1/varyhati )
  }

  return(list(g,h))

}
