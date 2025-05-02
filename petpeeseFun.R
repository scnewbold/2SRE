#===============================================================================
# This is an R function named 'petpeeseFun'
# which computes the PET-PEESE publication bias correction estimator for
# the two-stage random-effects meta-analysis estimators
# described in Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. (2025)
#===============================================================================

petpeeseFun <- function(Y,SE){

  fail <- 0

  SE2 <- SE^2
  # PET (Egger's regression):
  outs  <- lm(Y ~ SE,weights=1/SE2)
  tstat <- outs$coefficients[2]/summary(outs)$coefficients[2,2]
  if(tstat>1.96){
    # PEESE
    outs       <- lm(Y ~ SE2,weights=1/SE2)
    yhat.PP    <- outs$coefficients[1]
    se.PP      <- summary(outs)$coefficients[1,2]
    pub.bias   <- 1
  }else{
    outs       <- lm(Y ~ 0 + matrix(1,length(Y),1),weights=1/SE2)
    yhat.PP    <- outs$coefficients[1]
    se.PP      <- summary(outs)$coefficients[1,2]
    pub.bias   <- 0
  }

  return(list(yhat.PP=yhat.PP,se.PP=se.PP,fail=fail,pub.bias=pub.bias))

}
