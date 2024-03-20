#===============================================================================
# This is an R function named 'twosremaFun'
# which computes the two-stage random-effects meta-analysis estimators
# described in Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. (2024).
#===============================================================================
twosremaFun <- function(Y,SE,ID,rho){

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

  # UNCONSTRAINED SIGMU.HAT:

  # Estimate sig.mu for each group:
  {
    sig.muhat <- matrix(0,I,1)
    for(i in 1:I){

      if(J[i]>1){ # Estimate sig.mu only for groups with >1 observation

        CV <- se[i,1:J[i]] %*% t(se[i,1:J[i]])

        sig.mu2hat <- var(y[i,1:J[i]])-
        1/J[i]*(sum(diag(CV)-rho[i]/(J[i]^2)*(sum(CV)-sum(diag(CV)))))

        # Only estimate sig.mu if estimated sig.mu^2 is positive:
        if(sig.mu2hat>0){sig.muhat[i] <- sqrt(sig.mu2hat)}

      }
      # Small sample correction:
      sig.muhat[i] <- sig.muhat[i] * (J[i]/(J[i]-1))
    }

    # For all groups with one observation, impute sig.mu as sample size
    # weighted average of all other group's sig.muhat':
    sig.muhat[J==1] <- sum( sig.muhat[J>1] * J[J>1] / sum(J[J>1]) )
  }

  # Calculate ghat:
  {
    # g's don't use sig.eta
    outs <- ghFun(se,0,sig.muhat,rho)
    ghat <- outs[[1]]
  }

  # Estimate sig.eta:
  {
    # Calculate yhati and Vi for each group w/ more than one observation:
    yhati <- matrix(0,I,1)
    Vi    <- matrix(0,I,1)
    for(i in 1:I){
      if(J[i]>1){
        yhati[i] <- sum(y[i,1:J[i]]*ghat[i,1:J[i]])
        for(j in 1:J[i]){
          Vi[i] <- Vi[i] + ghat[i,j]^2*(sig.muhat[i]^2+se[i,j]^2)
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
    sig.etahat <- 0
    if(I>1 & sum(J>1)>0){
      part1 <- 1/(I-1)*sum((yhati-mean(yhati))^2)
      part2 <- 1/sum(J>1)*sum(Vi[J>1])
      if(part1>part2){sig.etahat <- sqrt(part1-part2)}
    }
  }

  # Calculate ghat and hhat:
  {
    outs <- ghFun(se,sig.etahat,sig.muhat,rho)
    ghat <- outs[[1]]
    hhat <- outs[[2]]
  }

  # CONSTRAINED SIGMU.HAT:

  # Estimate sig.muC (same for all groups)
  {
    sig.muhatC <- matrix(mean(sig.muhat[which(J>1)]),I,1)
  }

  # Calculate ghatC:
  {
    outs  <- ghFun(se,0,sig.muhatC,rho)
    ghatC <- outs[[1]]
  }

  # Estimate sig.etaC:
  {
    # Calculate yhati and Vi for each group w/ more than one observation:
    yhati <- matrix(0,I,1)
    Vi    <- matrix(0,I,1)
    for(i in 1:I){
      if(J[i]>1){
        yhati[i] <- sum(y[i,1:J[i]]*ghatC[i,1:J[i]])
        for(j in 1:J[i]){
          Vi[i] <- Vi[i] + ghatC[i,j]^2*(sig.muhatC[i]^2+se[i,j]^2)
          for(k in 1:J[i]){
            if(k!=j){
              Vi[i] <- Vi[i] + rho[i]*ghatC[i,j]*ghatC[i,k]*se[i,j]*se[i,k]
            }
          }
        }
      }
      if(J[i]==1){yhati[i] <- y[i,1]}
    }

    # Use yhati's and Vi's to calculate sig.etahat [Equation 49]:
    part1 <- 1/(I-1)*sum((yhati-mean(yhati))^2)
    part2 <- 1/sum(J>1)*sum(Vi[J>1])
    if(part1>part2){sig.etahatC <- sqrt(part1-part2)}else{sig.etahatC <- 0}
  }

  # Calculate ghatC and hhatC:
  {
    outs  <- ghFun(se,sig.etahatC,sig.muhatC,rho)
    ghatC <- outs[[1]]
    hhatC <- outs[[2]]
  }

  # COMPUTE ESTIMATES:
  {
    yhat.sm <- mean(y[which(y!=0)])
    yhat.mm <- mean((1/J)*rowSums(y))
    yhat.ru <- sum(rowSums(ghat*y)*hhat)
    yhat.rc <- sum(rowSums(ghatC*y)*hhatC)
    se.ru   <- seroFun(I,J,y,yhat.ru,ghat,hhat)
    se.rc   <- seroFun(I,J,y,yhat.rc,ghatC,hhatC)
  }

  return(list(yhat.sm     = yhat.sm,
              yhat.mm     = yhat.mm,
              yhat.ru     = yhat.ru,
              yhat.rc     = yhat.rc,
              se.ru       = se.ru,
              se.rc       = se.rc,
              sig.muhat   = sig.muhat,
              sig.etahat  = sig.etahat,
              sig.muhatC  = sig.muhatC,
              sig.etahatC = sig.etahatC))
}
