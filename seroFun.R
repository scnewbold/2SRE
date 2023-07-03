# This is a function named 'seroFun'

# which computes robust standard errors (following Hedges et al 2010) for the
# the two-stage random-effects meta-analysis estimators
# described in Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. 2023.

seroFun <- function(I,J,y,yhat,g,h){
  outer <- 0
  middle <- 0
  for(i in 1:I){
    Xi <- matrix(1,J[i],1)
    ei <- matrix(yhat-y[i,1:J[i]],J[i],1)
    if(J[i]>1){Wi <- diag(g[i,1:J[i]]*h[i])}
    else      {Wi <- as.matrix(g[i,1]*h[i])}
    outer  <- outer  + t(Xi) %*% Wi %*% Xi
    middle <- middle + t(Xi) %*% Wi %*% ei %*% t(ei) %*% Wi %*% Xi
  }
  outer.inv <- solve(outer)
  sero      <- sqrt(outer.inv %*% middle %*% outer.inv)*I/(I-1)

  return(sero)
}
