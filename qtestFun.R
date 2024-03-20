#===============================================================================
# This is an R function named 'qtestFUN'
# which computes the Q statistic to test for excess heterogeneity among
# the primary estimates for meta-anlaysis, as used in
# Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. (2024)
#===============================================================================

qtestFun <- function(ID,y,se){
  IDs  <- unique(ID)
  I    <- length(IDs)
  Q    <- matrix(0,I,1)
  crit <- matrix(0,I,1)
  for(i in 1:I){
    ii      <- which(ID==IDs[i])
    Q[i]    <- sum((y[ii]-mean(y[ii]))/(se[ii]^2))
    crit[i] <- qchisq(.95,sum(ID==IDs[i]))
  }
  return(list(Q,crit))
}
