# This is a function named 'trimfillFun'

# which computes the trim-and-fill publication bias correction estimator for
# the two-stage random-effects meta-analysis estimators
# described in Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. 2023.

# 1. Apply 2SRE meta-analysis to full dataset to estimate VSL.i
# 2. Compute Wilcoxen rank test statistic, Tn = sum of ranks of individual
#    observations larger than VSL.i
# 3. Compute k = (4*Tn-n*(n+1))/(2*n-1) rounded to integer and truncated at 0
# 4. Trim right-most [largest] k observations
# 5. Return to step 1 and repeat until k does not change.

trimfillFun <- function(Y,SE,ID){


  # contrained (2SRE-equal)
  {
    
    N       <- length(Y)
    Yhati.c <- -99
    fail.c  <- 0
    outs    <- twosremaFun(Y,SE,ID,rho)
    Yhat    <- outs$yhat.rc
    Tn      <- sum(rank(Y)*(Y>Yhat))
    k       <- max(0,floor((4*Tn-N*(N+1))/(2*N-1)))
    done  <- 0
    while(done==0){
      keep.c <- which(rank(Y)<(N-k+1))
      Yi     <- Y[keep.c]
      SEi    <- SE[keep.c]
      IDi    <- ID[keep.c]
      if(length(unique(IDi))==length(IDi)){done <- 1; fail.c <- 1}else{
        outs    <- twosremaFun(Yi,SEi,IDi,rho)
        Yhati.c <- outs$yhat.rc
        Tni     <- sum(rank(Y)*(Y>Yhati.c))
        ki      <- max(0,floor((4*Tni-N*(N+1))/(2*N-1)))
        if(ki<=k){done <- 1}else{k <- ki}
      }
    }
    
  }
  
  # uncontrained (2SRE-free)
  {
    
    N       <- length(Y)
    Yhati.u <- -99
    fail.u  <- 0
    outs    <- twosremaFun(Y,SE,ID,rho)
    Yhat    <- outs$yhat.ru
    Tn      <- sum(rank(Y)*(Y>Yhat))
    k       <- max(0,floor((4*Tn-N*(N+1))/(2*N-1)))
    done    <- 0
    while(done==0){
      keep.u <- which(rank(Y)<(N-k+1))
      Yi     <- Y[keep.u]
      SEi    <- SE[keep.u]
      IDi    <- ID[keep.u]
      if(length(unique(IDi))==length(IDi)){done <- 1; fail.u <- 1}else{
        outs    <- twosremaFun(Yi,SEi,IDi,rho)
        Yhati.u <- outs$yhat.ru
        Tni     <- sum(rank(Y)*(Y>Yhati.u))
        ki      <- max(0,floor((4*Tni-N*(N+1))/(2*N-1)))
        if(ki<=k){done <- 1}else{k <- ki}
      }
    }
    
  }

  return(list(Yhati.c = Yhati.c,
              fail.c  = fail.c,
              keep.c  = keep.c,
              Yhati.u = Yhati.u,
              fail.u  = fail.u,
              keep.u  = keep.u ))

}
