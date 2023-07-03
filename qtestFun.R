# This is a function named 'qtestFun'

# which [...does what...]
# described in Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. 2023.

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
