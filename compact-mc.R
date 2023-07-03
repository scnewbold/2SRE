#===============================================================================
# compact-mc.R
# This script runs a streamlined (compact) Monte Carlo experiment to verify
# that our estimators for the non-sampling error variances are valid.
# Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. (2023)
#===============================================================================

source("getdataFun.R")
source("yhatFun.R")
source("pauseFun.R")

set.seed(1234)

MC       <- 1000
sig2.mu5 <- matrix(0,MC,1)
sig2.eta <- matrix(0,MC,1)
yhat     <- matrix(0,MC,1)

par(mfrow=c(1,3))

for(mc in 1:MC){

  outs   <- getdataFun('sim',floor(runif(1,100,10000)))
  ID     <- outs$ID
  Y      <- outs$Y
  SE     <- outs$SE
  y.true <- outs$y.true

  outs         <- yhatFun(ID,Y,SE,rho)
  yhat[mc]     <- outs$yhat
  sig2.eta[mc] <- outs$sig2.eta
  sig2.mu5[mc] <- outs$sig2.mu[5]

  if(mc>10){

    cat('\014Working on mc =',sprintf('%-.0f',mc),'of',sprintf('%-.0f',MC))

    hist(yhat[1:mc])
    lines(c(y.true,y.true),c(0,MC),col='red')
    lines(c(mean(yhat[1:mc]),mean(yhat[1:mc])),c(0,MC),col='blue')

    hist(sig2.eta[1:mc])
    lines(c(4,4),c(0,MC),col='red')
    lines(c(mean(sig2.eta[1:mc]),mean(sig2.eta[1:mc])),c(0,MC),col='blue')

    hist(sig2.mu5[1:mc])
    lines(c(2,2),c(0,MC),col='red')
    lines(c(mean(sig2.mu5[1:mc]),mean(sig2.mu5[1:mc])),c(0,MC),col='blue')

    pauseFun(.05)
  }

}

