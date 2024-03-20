#===============================================================================
# This is an R function named 'getadataFun'
# which will construct a small simulated dataset or will import the
# VSL demonstration dataset used in
# Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. (2024)
#===============================================================================

getdataFun <- function(sim.or.app,seed){

  if(sim.or.app == 'sim'){

    set.seed(seed)

    y.true   <- 10
    sig2.mu  <- c(2,3,2,3,2)
    sig2.eta <- 4
    ID       <- c(1,2,2,3,3,3,4,4,4,4,5,5,5,5,5)
    SE       <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)
    I        <- 5
    J        <- c(1,2,3,4,5)
    rho      <- c(.5,.5,.5,.5,.5)

    y <- matrix(0,I,max(J))
    for(i in 1:I){

      temp <- matrix(1,J[i],J[i])*rho[i]
      for(j in 1:J[i]){temp[j,j] <- 1}
      CF  <- chol(temp)
      sei <- SE[which(ID==i)]
      CV  <- sei %*% t(sei)

      eta.i <- rnorm(1)*sqrt(sig2.eta)

      yi <- matrix(rnorm(J[i]),1,J[i]) %*% CF
      yi <- matrix(yi,J[i],1)
      for(j in 1:J[i]){
        y[i,j] <- yi[j]*sei[j] + y.true + rnorm(1)*sqrt(sig2.mu[i]) + eta.i
      }
    }

    Y <- c()
    for(i in 1:I){
      Y <- c(Y,y[i,1:J[i]])
    }

  }

  if(sim.or.app == 'app'){

    # IMPORT DATA USED IN NEWBOLD ET AL (2024):
    if(TRUE){
      gs4_deauth()

      # Input data are in the Google Sheet file titled "EPA_VSL_metadata" at the
      # following link:

      data <- read_sheet('https://docs.google.com/spreadsheets/d/1wxWgCSZKYWBuX55i4vCCw-ZyldcUP3SfbdvkOnAqdvA/edit?usp=sharing',sheet='Sheet1')

      ID <- data$GroupID             # unique group id
      Y  <- data$VSLestimate.2013    # VSL estimate [2013$US]
      SE <- data$StandardError.2013  # standard error [2013$US]
    }

  }

  return(list(ID=ID,Y=Y,SE=SE,y.true=y.true))

}
