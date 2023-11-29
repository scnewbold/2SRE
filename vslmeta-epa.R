#===============================================================================
# vslmeta-epa.R:
# This script will replicate the results of our application of the 2SRE meta-
# analysis estimator to the EPA VSL demonstration dataset as reported in
# Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. (2023)
#===============================================================================

#-------------------------------------------------------------------------------
# PRELIMINARIES:
#-------------------------------------------------------------------------------
{
  # Clear environment to start fresh
  rm(list=ls()) 
  
  # Clear console:
  cat('\014');
  
  # Clear all plots:
  try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
  try(dev.off(),silent=TRUE)
  
  # Grab script name for file handling:
  script.name <- basename(rstudioapi::getSourceEditorContext()$path) 
  script.name <- gsub(".R","",script.name)
  
  # Packages
  list.of.packages <-
    c('googlesheets4',
      'NlcOptim',
      'MASS',
      'metafor',
      'robumeta',
      'MAd',
      'tikzDevice',
      'readxl',
      'stringr')

  new.packages <- list.of.packages[!(list.of.packages %in%
                  installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages,function(x){library(x,character.only=TRUE)})

  # Define paths for file handling:
  this.dir <- dirname(parent.frame(2)$ofile) # source file dir
  setwd(this.dir)                            # set wd to source file dir
  code.path <- getwd()                       # define code path
  output.path <- getwd()                     # define output path

  # Create output file in working directory:
  date.time     <- gsub(" ","_",Sys.time())
  date.time     <- gsub("-","_",date.time)
  date.time     <- gsub(":","_",date.time)
  out.file.name <- paste(output.path,'/',script.name,'-',date.time,'.out',sep='')
  outfile       <- file.create(out.file.name)

  # WRITE SOURCE FILE TO OUTPUT FILE:
  {
    source.file.name <- paste(code.path,'/',script.name,'.R',sep='')

    # Read lines of source file:
    Rscript <- readLines(source.file.name)

    # Write lines of source file to output file:
    for(i in 1:length(Rscript)){cat('\n',Rscript[i],file=out.file.name,append=TRUE)}
    cat('\n\n|---------------------------------------------------------------------------|',file=out.file.name,append=TRUE)
    cat('\n| R SCRIPT ABOVE                                                            |',file=out.file.name,append=TRUE)
    cat('\n|---------------------------------------------------------------------------|',file=out.file.name,append=TRUE)
    cat('\n| R OUTPUT BELOW                                                            |',file=out.file.name,append=TRUE)
    cat('\n|---------------------------------------------------------------------------|\n',file=out.file.name,append=TRUE)

  }

  s2 <- function(x){return(sprintf('%-.2f',x))}
  s3 <- function(x){return(sprintf('%-.3f',x))}
  s4 <- function(x){return(sprintf('%-.4f',x))}

}

#-------------------------------------------------------------------------------
# FUNCTIONS:
#-------------------------------------------------------------------------------
{
  source("pauseFun.R")
  source("ghFun.R")
  source("seroFun.R")
  source("qtestFun.R")
  source("trimfillFun.R")
  source("petpeeseFun.R")
  source("twosremaFun.R")
  source("twosremrFun.R")
}

#-------------------------------------------------------------------------------
# MAIN PROGRAM:
#-------------------------------------------------------------------------------

# Assumed correlation among observations within studies (rho):
# Our estimation approach does not allow estimation of rho, so we must
# assume a value for this parameter. We can test the robustness of our results
# by varying this assumption and reporting the associated VSL estimates.
rho.all <- 0.5

# Replace any negative g's with 0 and re-scale remaining g's (1=yes,0=no):
g.scale <- 1

# In the paper, meta-analysis results are reported in Table 8.
TABLE8 <- matrix(0,14,11)
for(case in 1:8){

  set.seed(1234)
  B <- 1000 # Number of bootstrap reps:

  if(case==1){include.means<-1;include.medians<-1;include.HW<-1;include.SP<-0;balanced<-0} # HW mm
  if(case==2){include.means<-1;include.medians<-0;include.HW<-1;include.SP<-0;balanced<-0} # HW m
  if(case==3){include.means<-1;include.medians<-1;include.HW<-0;include.SP<-1;balanced<-0} # SP mm
  if(case==4){include.means<-1;include.medians<-0;include.HW<-0;include.SP<-1;balanced<-0} # SP m
  if(case==5){include.means<-1;include.medians<-1;include.HW<-1;include.SP<-1;balanced<-0} # pooled mm
  if(case==6){include.means<-1;include.medians<-0;include.HW<-1;include.SP<-1;balanced<-0} # pooled m
  if(case==7){include.means<-1;include.medians<-1;include.HW<-1;include.SP<-1;balanced<-1} # balanced mm
  if(case==8){include.means<-1;include.medians<-0;include.HW<-1;include.SP<-1;balanced<-1} # balanced m

  # IMPORT DATA FOR THIS CASE:
  {
    gs4_deauth()

    # Input data are in the Google Sheet file titled "EPA_VSL_metadata" at the
    # following link:
    data <- read_sheet('https://docs.google.com/spreadsheets/d/1wxWgCSZKYWBuX55i4vCCw-ZyldcUP3SfbdvkOnAqdvA/edit?usp=sharing',
    sheet='meta-data')

    if(include.means==0)  {data <- data[data$MeanDummy==0,]}
    if(include.medians==0){data <- data[data$MeanDummy==1,]}
    if(include.HW==0)     {data <- data[data$SPDummy==1,]}
    if(include.SP==0)     {data <- data[data$SPDummy==0,]}

    # Extract variables:
    ID              <- data$GroupID          # unique group id
    pubyear         <- data$PubYear          # year of publication
    datayear        <- data$DataYear         # year of data collection
    dollaryear      <- data$DollarYear       # US$ year of reported VSL estimates
    mnD             <- data$MeanDummy        # mean dummy variable
    mdD             <- matrix(1,length(mnD),1)-mnD # median dummy variable
    spD             <- data$SPDummy          # stated preference dummy variable
    samplesize      <- data$SampleSize       # sample size of original study
    Y               <- data$VSL              # VSL estimate [DataYear $US]
    SE              <- data$SE               # standard error [DataYear $US]
    income          <- data$Income / 10000   # [2013 $US]

    # Group ids:
    IDs <- unique(ID)

    # Number of groups:
    I <- length(IDs)

    # Assumed correlation among observations within studies (rho):
    rho <- matrix(rho.all,I,1)

    # Observations per group:
    J <- matrix(0,I,1)
    for(i in 1:I){
      J[i] <- sum(ID==IDs[i])
    }

    # Total number of observations:
    N <- length(Y)

    # Convert all estimates from DataYear $US to 2020$US using the CPI index:
    # Source: https://data.bls.gov/pdq/SurveyOutputServlet
    CPI <- read_sheet('https://docs.google.com/spreadsheets/d/1wxWgCSZKYWBuX55i4vCCw-ZyldcUP3SfbdvkOnAqdvA/edit?usp=sharing',
                       sheet='CPI',range='A13:P122',col_names=FALSE)
    CPI      <- as.matrix(CPI)
    CPIyr    <- CPI[,1]
    CPIindex <- CPI[,14]
    
    Y2 <- Y
    SE2 <- SE
    
    Y2013 <- Y
    SE2013 <- SE
    
    income2 <- income
    for(i in 1:length(Y)){
      Y2[i]  <- Y[i]  * CPIindex[which(CPIyr==2020)] / CPIindex[which(CPIyr==dollaryear[i])]
      SE2[i] <- SE[i] * CPIindex[which(CPIyr==2020)] / CPIindex[which(CPIyr==dollaryear[i])]
      income2[i] <- income[i] * CPIindex[which(CPIyr==2020)] / CPIindex[which(CPIyr==2013)]
      Y2013[i]   <- Y[i]  * CPIindex[which(CPIyr==2013)] / CPIindex[which(CPIyr==dollaryear[i])]
      SE2013[i]  <- SE[i] * CPIindex[which(CPIyr==2013)] / CPIindex[which(CPIyr==dollaryear[i])]
    }
    Y <- Y2
    SE <- SE2
    income <- income2

    # Number of observations from parent group of each observation:
    m <- matrix(0,N,1)
    for(ij in 1:N){
      m[ij] <- sum(ID==ID[ij])
    }

    #if(case==5){data1<-cbind(Y,SE)}

  }

  # TWO-STAGE RANDOM-EFFECTS META-ANALYSIS (no covariates) (2SREMA):
  if(TRUE){

    Y0  <- Y
    SE0 <- SE
    ID0 <- ID

    if(balanced==0){

      outs     <- twosremaFun(Y,SE,ID,rho)
      yhat.sm  <- outs$yhat.sm
      yhat.mm  <- outs$yhat.mm
      yhat.ru  <- outs$yhat.ru
      yhat.rc  <- outs$yhat.rc
      sero.ru  <- outs$se.ru
      sero.rc  <- outs$se.rc

      outs     <- twosremrFun(Y,matrix(1,length(Y),1),SE,ID,rho)

      # Bootstrap standard errors:
      yhat.sm.BS <- matrix(0,B,1)
      yhat.mm.BS <- matrix(0,B,1)
      yhat.ru.BS <- matrix(0,B,1)
      yhat.rc.BS <- matrix(0,B,1)

      for(b in 1:B){

        if(floor(b/10)==b/10){
          cat('\014')
          cat('Case ',sprintf('%-.0f',case),'-- meta-analysis\n')
          cat('Working on bootstrap rep',sprintf('%-.0f',b),'of',sprintf('%-.0f',B))
        }

        # Resample data with replacement:
        done <- 0
        while(done==0){
          BSgroups <- sample(IDs,I,replace=TRUE)
          z   <- which(ID==BSgroups[1])
          Yb  <- Y[z]
          SEb <- SE[z]
          IDb <- matrix(1,length(Yb),1)
          for(i in 2:I){
            z   <- which(ID==BSgroups[i])
            Yb  <- c(Yb,Y[z])
            SEb <- c(SEb,SE[z])
            IDb <- rbind(IDb,matrix(i,length(z),1))
          }
          done <- 1
        }
        outs <- twosremaFun(Yb,SEb,IDb,rho)
        yhat.sm.BS[b]  <- outs[[1]]
        yhat.mm.BS[b]  <- outs[[2]]
        yhat.ru.BS[b]  <- outs[[3]]
        yhat.rc.BS[b]  <- outs[[4]]

      }

      sebs.sm <- sqrt(var(yhat.sm.BS))
      sebs.mm <- sqrt(var(yhat.mm.BS))
      sebs.ru <- sqrt(var(yhat.ru.BS))
      sebs.rc <- sqrt(var(yhat.rc.BS))

    }

    if(balanced==1){

      # SP:
      Y   <- Y0[which(spD==1)]
      SE  <- SE0[which(spD==1)]
      ID  <- ID0[which(spD==1)]
      IDs <- unique(ID)
      I   <- length(IDs)

      outs       <- twosremaFun(Y,SE,ID,rho)
      yhat.sm.SP <- outs[[1]]
      yhat.mm.SP <- outs[[2]]
      yhat.ru.SP <- outs[[3]]
      yhat.rc.SP <- outs[[4]]
      sero.ru.SP <- outs[[5]]
      sero.rc.SP <- outs[[6]]

      # Bootstrap standard errors:
      yhat.sm.BS.SP <- matrix(0,B,1)
      yhat.mm.BS.SP <- matrix(0,B,1)
      yhat.ru.BS.SP <- matrix(0,B,1)
      yhat.rc.BS.SP <- matrix(0,B,1)

      for(b in 1:B){

        if(floor(b/10)==b/10){
          cat('\014')
          cat('Case ',sprintf('%-.0f',case),'-- meta-analysis\n')
          cat('Working on bootstrap rep',sprintf('%-.0f',b),'of',sprintf('%-.0f',B))
        }

        # Resample data with replacement:
        done <- 0
        while(done==0){
          BSgroups <- sample(IDs,I,replace=TRUE)
          z   <- which(ID==BSgroups[1])
          Yb  <- Y[z]
          SEb <- SE[z]
          IDb <- matrix(1,length(Yb),1)
          for(i in 2:I){
            z   <- which(ID==BSgroups[i])
            Yb  <- c(Yb,Y[z])
            SEb <- c(SEb,SE[z])
            IDb <- rbind(IDb,matrix(i,length(z),1))
          }
          done <- 1
        }
        outs <- twosremaFun(Yb,SEb,IDb,rho)
        yhat.sm.BS.SP[b]  <- outs[[1]]
        yhat.mm.BS.SP[b]  <- outs[[2]]
        yhat.ru.BS.SP[b]  <- outs[[3]]
        yhat.rc.BS.SP[b]  <- outs[[4]]

      }

      # HW:
      Y   <- Y0[which(spD==0)]
      SE  <- SE0[which(spD==0)]
      ID  <- ID0[which(spD==0)]
      IDs <- unique(ID)
      I   <- length(IDs)

      outs       <- twosremaFun(Y,SE,ID,rho)
      yhat.sm.HW <- outs[[1]]
      yhat.mm.HW <- outs[[2]]
      yhat.ru.HW <- outs[[3]]
      yhat.rc.HW <- outs[[4]]
      sero.ru.HW <- outs[[5]]
      sero.rc.HW <- outs[[6]]

      # Bootstrap standard errors:
      yhat.sm.BS.HW <- matrix(0,B,1)
      yhat.mm.BS.HW <- matrix(0,B,1)
      yhat.ru.BS.HW <- matrix(0,B,1)
      yhat.rc.BS.HW <- matrix(0,B,1)

      for(b in 1:B){

        if(floor(b/10)==b/10){
          cat('\014')
          cat('Working on bootstrap rep',sprintf('%-.0f',b),'of',sprintf('%-.0f',B))
        }

        # Resample data with replacement:
        done <- 0
        while(done==0){
          BSgroups <- sample(IDs,I,replace=TRUE)
          z   <- which(ID==BSgroups[1])
          Yb  <- Y[z]
          SEb <- SE[z]
          IDb <- matrix(1,length(Yb),1)
          for(i in 2:I){
            z   <- which(ID==BSgroups[i])
            Yb  <- c(Yb,Y[z])
            SEb <- c(SEb,SE[z])
            IDb <- rbind(IDb,matrix(i,length(z),1))
          }
          done <- 1
        }
        outs <- twosremaFun(Yb,SEb,IDb,rho)
        yhat.sm.BS.HW[b]  <- outs[[1]]
        yhat.mm.BS.HW[b]  <- outs[[2]]
        yhat.ru.BS.HW[b]  <- outs[[3]]
        yhat.rc.BS.HW[b]  <- outs[[4]]

      }

      yhat.sm = .5*(yhat.sm.SP+yhat.sm.HW)
      yhat.mm = .5*(yhat.mm.SP+yhat.mm.HW)
      yhat.ru = .5*(yhat.ru.SP+yhat.ru.HW)
      yhat.rc = .5*(yhat.rc.SP+yhat.rc.HW)

      yhat.sm.BS = .5*(yhat.sm.BS.SP+yhat.sm.BS.HW)
      yhat.mm.BS = .5*(yhat.mm.BS.SP+yhat.mm.BS.HW)
      yhat.ru.BS = .5*(yhat.ru.BS.SP+yhat.ru.BS.HW)
      yhat.rc.BS = .5*(yhat.rc.BS.SP+yhat.rc.BS.HW)

      sebs.sm <- sqrt(var(yhat.sm.BS))
      sebs.mm <- sqrt(var(yhat.mm.BS))
      sebs.ru <- sqrt(var(yhat.ru.BS))
      sebs.rc <- sqrt(var(yhat.rc.BS))

    }

  }

  # PUBLICATION BIAS (Trim-and-Fill estimator):
  if(TRUE){

    # yhat.TF <- -99; sebs.TF <- -99

    Y  <- Y0
    SE <- SE0
    ID <- ID0

    outs     <- trimfillFun(Y,SE,ID)
    yhat.TFu <- outs$Yhati.u
    yhat.TFc <- outs$Yhati.c

    # Bootstrap standard errors:
    yhat.TF.BSu <- matrix(0,B,1)
    yhat.TF.BSc <- matrix(0,B,1)
    for(b in 1:B){

      if(floor(b/10)==b/10){
        cat('\014')
        cat('Case ',sprintf('%-.0f',case),'-- trim and fill\n')
        cat('Working on bootstrap rep',sprintf('%-.0f',b),'of',sprintf('%-.0f',B))
      }

      # Resample data with replacement:
      done <- 0
      while(done==0){
        BSgroups <- sample(IDs,I,replace=TRUE)
        z   <- which(ID==BSgroups[1])
        Yb  <- Y[z]
        SEb <- SE[z]
        IDb <- matrix(1,length(Yb),1)
        for(i in 2:I){
          z   <- which(ID==BSgroups[i])
          Yb  <- c(Yb,Y[z])
          SEb <- c(SEb,SE[z])
          IDb <- rbind(IDb,matrix(i,length(z),1))
        }
        outs <- trimfillFun(Yb,SEb,IDb)
        if((outs$fail.c+outs$fail.u)==0){done <- 1} # Discards re-sampled data sets that won't estimate.
      }
      yhat.TF.BSu[b] <- outs$Yhati.u
      yhat.TF.BSc[b] <- outs$Yhati.c

    }
    sebs.TFu <- sd(yhat.TF.BSu)
    sebs.TFc <- sd(yhat.TF.BSc)

  }

  # PUBLICATION BIAS (PET-PEESE estimator):
  if(TRUE){

    outs    <- petpeeseFun(Y,SE)
    yhat.PP <- outs$yhat.PP
    se.PP   <- outs$se.PP

    # Bootstrap standard errors:
    yhat.PP.BS <- matrix(0,B,1)
    for(b in 1:B){

      if(floor(b/10)==b/10){
        cat('\014')
        cat('Case ',sprintf('%-.0f',case),'-- PET-PEESE\n')
        cat('Working on bootstrap rep',sprintf('%-.0f',b),'of',sprintf('%-.0f',B))
      }

      # Resample data with replacement:
      done <- 0
      while(done==0){
        BSgroups <- sample(IDs,I,replace=TRUE)
        z   <- which(ID==BSgroups[1])
        Yb  <- Y[z]
        SEb <- SE[z]
        IDb <- matrix(1,length(Yb),1)
        for(i in 2:I){
          z   <- which(ID==BSgroups[i])
          Yb  <- c(Yb,Y[z])
          SEb <- c(SEb,SE[z])
          IDb <- rbind(IDb,matrix(i,length(z),1))
        }
        outs <- petpeeseFun(Yb,SEb)
        if(outs$fail==0){done <- 1} # Discards re-sampled data sets that won't estimate.
      }
      yhat.PP.BS[b] <- outs$yhat.PP

    }
    sebs.PP <- sd(yhat.PP.BS)

  }

  # SAVE RESULTS FOR THIS CASE TO TABLE8 MATRIX
  {
    if(case==1){
      TABLE8[1,1] <-yhat.sm; TABLE8[1,2] <-sebs.sm
      TABLE8[3,1] <-yhat.mm; TABLE8[3,2] <-sebs.mm
      TABLE8[5,1] <-yhat.ru; TABLE8[5,2] <-sebs.ru
      TABLE8[7,1] <-yhat.rc; TABLE8[7,2] <-sebs.rc
      TABLE8[9,1] <-yhat.TFu;TABLE8[9,2] <-sebs.TFu
      TABLE8[11,1]<-yhat.TFc;TABLE8[11,2]<-sebs.TFc
      TABLE8[13,1]<-yhat.PP; TABLE8[13,2]<-sebs.PP
    }
    if(case==2){
      TABLE8[2,1] <-yhat.sm; TABLE8[2,2] <-sebs.sm
      TABLE8[4,1] <-yhat.mm; TABLE8[4,2] <-sebs.mm
      TABLE8[6,1] <-yhat.ru; TABLE8[6,2] <-sebs.ru
      TABLE8[8,1] <-yhat.rc; TABLE8[8,2] <-sebs.rc
      TABLE8[10,1]<-yhat.TFu;TABLE8[10,2]<-sebs.TFu
      TABLE8[12,1]<-yhat.TFc;TABLE8[12,2]<-sebs.TFc
      TABLE8[14,1]<-yhat.PP; TABLE8[14,2]<-sebs.PP
    }
    if(case==3){
      TABLE8[1,3] <-yhat.sm; TABLE8[1,4] <-sebs.sm
      TABLE8[3,3] <-yhat.mm; TABLE8[3,4] <-sebs.mm
      TABLE8[5,3] <-yhat.ru; TABLE8[5,4] <-sebs.ru
      TABLE8[7,3] <-yhat.rc; TABLE8[7,4] <-sebs.rc
      TABLE8[9,3] <-yhat.TFu;TABLE8[9,4] <-sebs.TFu
      TABLE8[11,3]<-yhat.TFc;TABLE8[11,4]<-sebs.TFc
      TABLE8[13,3]<-yhat.PP; TABLE8[13,4]<-sebs.PP
    }
    if(case==4){
      TABLE8[2,3] <-yhat.sm; TABLE8[2,4] <-sebs.sm
      TABLE8[4,3] <-yhat.mm; TABLE8[4,4] <-sebs.mm
      TABLE8[6,3] <-yhat.ru; TABLE8[6,4] <-sebs.ru
      TABLE8[8,3] <-yhat.rc; TABLE8[8,4] <-sebs.rc
      TABLE8[10,3]<-yhat.TFu;TABLE8[10,4]<-sebs.TFu
      TABLE8[12,3]<-yhat.TFc;TABLE8[12,4]<-sebs.TFc
      TABLE8[14,3]<-yhat.PP; TABLE8[14,4]<-sebs.PP
    }
    if(case==5){
      TABLE8[1,6] <-yhat.sm; TABLE8[1,7] <-sebs.sm
      TABLE8[3,6] <-yhat.mm; TABLE8[3,7] <-sebs.mm
      TABLE8[5,6] <-yhat.ru; TABLE8[5,7] <-sebs.ru
      TABLE8[7,6] <-yhat.rc; TABLE8[7,7] <-sebs.rc
      TABLE8[9,6] <-yhat.TFu;TABLE8[9,7] <-sebs.TFu
      TABLE8[11,6]<-yhat.TFc;TABLE8[11,7]<-sebs.TFc
      TABLE8[13,6]<-yhat.PP; TABLE8[13,7]<-sebs.PP
    }
    if(case==6){
      TABLE8[2,6] <-yhat.sm; TABLE8[2,7] <-sebs.sm
      TABLE8[4,6] <-yhat.mm; TABLE8[4,7] <-sebs.mm
      TABLE8[6,6] <-yhat.ru; TABLE8[6,7] <-sebs.ru
      TABLE8[8,6] <-yhat.rc; TABLE8[8,7] <-sebs.rc
      TABLE8[10,6]<-yhat.TFu;TABLE8[10,7]<-sebs.TFu
      TABLE8[12,6]<-yhat.TFc;TABLE8[12,7]<-sebs.TFc
      TABLE8[14,6]<-yhat.PP; TABLE8[14,7]<-sebs.PP
    }
    # Balanced cases:
    if(case==7){
      TABLE8[1,9]  <- mean(c(TABLE8[1,1],TABLE8[1,3]))
      TABLE8[1,10] <- .5^2*TABLE8[1,2]+.5^2*TABLE8[1,4]

      TABLE8[3,9]  <- mean(c(TABLE8[3,1],TABLE8[3,3]))
      TABLE8[3,10] <- .5^2*TABLE8[3,2]+.5^2*TABLE8[3,4]

      TABLE8[5,9]  <- mean(c(TABLE8[5,1],TABLE8[5,3]))
      TABLE8[5,10] <- .5^2*TABLE8[5,2]+.5^2*TABLE8[5,4]

      TABLE8[7,9]  <- mean(c(TABLE8[7,1],TABLE8[7,3]))
      TABLE8[7,10] <- .5^2*TABLE8[7,2]+.5^2*TABLE8[7,4]

      TABLE8[9,9]  <- mean(c(TABLE8[9,1],TABLE8[9,3]))
      TABLE8[9,10] <- .5^2*TABLE8[9,2]+.5^2*TABLE8[9,4]

      TABLE8[11,9] <- mean(c(TABLE8[11,1],TABLE8[11,3]))
      TABLE8[11,10]<- .5^2*TABLE8[11,2]+.5^2*TABLE8[11,4]
      
      TABLE8[13,9] <- mean(c(TABLE8[13,1],TABLE8[13,3]))
      TABLE8[13,10]<- .5^2*TABLE8[13,2]+.5^2*TABLE8[13,4]
    }
    if(case==8){
      TABLE8[2,9]  <- mean(c(TABLE8[2,1],TABLE8[2,3]))
      TABLE8[2,10] <- .5^2*TABLE8[2,2]+.5^2*TABLE8[2,4]

      TABLE8[4,9]  <- mean(c(TABLE8[4,1],TABLE8[4,3]))
      TABLE8[4,10] <- .5^2*TABLE8[4,2]+.5^2*TABLE8[4,4]

      TABLE8[6,9]  <- mean(c(TABLE8[6,1],TABLE8[6,3]))
      TABLE8[6,10] <- .5^2*TABLE8[6,2]+.5^2*TABLE8[6,4]

      TABLE8[8,9]  <- mean(c(TABLE8[8,1],TABLE8[8,3]))
      TABLE8[8,10] <- .5^2*TABLE8[8,2]+.5^2*TABLE8[8,4]

      TABLE8[10,9]  <- mean(c(TABLE8[10,1],TABLE8[10,3]))
      TABLE8[10,10] <- .5^2*TABLE8[10,2]+.5^2*TABLE8[10,4]

      TABLE8[12,9] <- mean(c(TABLE8[12,1],TABLE8[12,3]))
      TABLE8[12,10]<- .5^2*TABLE8[12,2]+.5^2*TABLE8[12,4]
      
      TABLE8[14,9] <- mean(c(TABLE8[14,1],TABLE8[14,3]))
      TABLE8[14,10]<- .5^2*TABLE8[14,2]+.5^2*TABLE8[14,4]
    }
    
    # RMSEs:
    {
      TABLE8[1,5] <- sqrt(TABLE8[1,4]^2 +(TABLE8[1,3] -TABLE8[2,3])^2) 
      TABLE8[3,5] <- sqrt(TABLE8[3,4]^2 +(TABLE8[3,3] -TABLE8[4,3])^2) 
      TABLE8[5,5] <- sqrt(TABLE8[5,4]^2 +(TABLE8[5,3] -TABLE8[6,3])^2)
      TABLE8[7,5] <- sqrt(TABLE8[7,4]^2 +(TABLE8[7,3] -TABLE8[8,3])^2) 
      TABLE8[9,5] <- sqrt(TABLE8[9,4]^2 +(TABLE8[9,3] -TABLE8[10,3])^2) 
      TABLE8[11,5]<- sqrt(TABLE8[11,4]^2+(TABLE8[11,3]-TABLE8[12,3])^2)
      TABLE8[13,5]<- sqrt(TABLE8[13,4]^2+(TABLE8[13,3]-TABLE8[14,3])^2) 
      
      TABLE8[1,8] <- sqrt(TABLE8[1,7]^2 +(TABLE8[1,6] -TABLE8[2,6])^2) 
      TABLE8[3,8] <- sqrt(TABLE8[3,7]^2 +(TABLE8[3,6] -TABLE8[4,6])^2) 
      TABLE8[5,8] <- sqrt(TABLE8[5,7]^2 +(TABLE8[5,6] -TABLE8[6,6])^2)
      TABLE8[7,8] <- sqrt(TABLE8[7,7]^2 +(TABLE8[7,6] -TABLE8[8,6])^2) 
      TABLE8[9,8] <- sqrt(TABLE8[9,7]^2 +(TABLE8[9,6] -TABLE8[10,6])^2) 
      TABLE8[11,8]<- sqrt(TABLE8[11,7]^2+(TABLE8[11,6]-TABLE8[12,6])^2)
      TABLE8[13,8]<- sqrt(TABLE8[13,7]^2+(TABLE8[13,6]-TABLE8[14,6])^2) 
      
      TABLE8[1,11] <- sqrt(TABLE8[1,10]^2 +(TABLE8[1,9] -TABLE8[2,9])^2) 
      TABLE8[3,11] <- sqrt(TABLE8[3,10]^2 +(TABLE8[3,9] -TABLE8[4,9])^2) 
      TABLE8[5,11] <- sqrt(TABLE8[5,10]^2 +(TABLE8[5,9] -TABLE8[6,9])^2)
      TABLE8[7,11] <- sqrt(TABLE8[7,10]^2 +(TABLE8[7,9] -TABLE8[8,9])^2) 
      TABLE8[9,11] <- sqrt(TABLE8[9,10]^2 +(TABLE8[9,9] -TABLE8[10,9])^2) 
      TABLE8[11,11]<- sqrt(TABLE8[11,10]^2+(TABLE8[11,9]-TABLE8[12,9])^2)
      TABLE8[13,11]<- sqrt(TABLE8[13,10]^2+(TABLE8[13,9]-TABLE8[14,9])^2) 
    }
  }

}

# WRITE TABLE8 TO OUTPUT FILE:
if(TRUE){
  cat('\nTable 8\n',file=out.file.name,append=TRUE)
  cat('\\hline\\hline\n',file=out.file.name,append=TRUE)
  cat('Estimator & mm/m & HW & SP & pooled & balanced \\\\ \n',file=out.file.name,append=TRUE)
  cat('\\hline\n',file=out.file.name,append=TRUE)
  cat('simple mean    & mm & ',
      s2(TABLE8[1,1]),' & (',s2(TABLE8[1,2]),') & ',
      s2(TABLE8[1,3]),' & (',s2(TABLE8[1,4]),') & [',s2(TABLE8[1,5]),'] & ',
      s2(TABLE8[1,6]),' & (',s2(TABLE8[1,7]),') & [',s2(TABLE8[1,8]),'] & ',
      s2(TABLE8[1,9]),' & (',s2(TABLE8[1,10]),') & [',s2(TABLE8[1,11]),'] \\\\ \n',sep='',
      file=out.file.name,append=TRUE)

  cat('               & m  & ',
      s2(TABLE8[2,1]),' & (',s2(TABLE8[2,2]),') & ',
      s2(TABLE8[2,3]),' & (',s2(TABLE8[2,4]),') & & ',
      s2(TABLE8[2,6]),' & (',s2(TABLE8[2,7]),') & & ',
      s2(TABLE8[2,9]),' & (',s2(TABLE8[2,10]),') & \\\\ \n',sep='',
      file=out.file.name,append=TRUE)

  cat('group means    & mm & ',
      s2(TABLE8[3,1]),' & (',s2(TABLE8[3,2]),') & ',
      s2(TABLE8[3,3]),' & (',s2(TABLE8[3,4]),') & [',s2(TABLE8[3,5]),'] & ',
      s2(TABLE8[3,6]),' & (',s2(TABLE8[3,7]),') & [',s2(TABLE8[3,8]),'] & ',
      s2(TABLE8[3,9]),' & (',s2(TABLE8[3,10]),') & [',s2(TABLE8[3,11]),'] \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('               & m  & ',
      s2(TABLE8[4,1]),' & (',s2(TABLE8[4,2]),') & ',
      s2(TABLE8[4,3]),' & (',s2(TABLE8[4,4]),') & & ',
      s2(TABLE8[4,6]),' & (',s2(TABLE8[4,7]),') & & ',
      s2(TABLE8[4,9]),' & (',s2(TABLE8[4,10]),') & \\\\ \n',sep='',
      file=out.file.name,append=TRUE)

  cat('2SRE--free     & mm & ',
      s2(TABLE8[5,1]),' & (',s2(TABLE8[5,2]),') & ',
      s2(TABLE8[5,3]),' & (',s2(TABLE8[5,4]),') & [',s2(TABLE8[5,5]),'] & ',
      s2(TABLE8[5,6]),' & (',s2(TABLE8[5,7]),') & [',s2(TABLE8[5,8]),'] & ',
      s2(TABLE8[5,9]),' & (',s2(TABLE8[5,10]),') & [',s2(TABLE8[5,11]),'] \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('               & m  & ',
      s2(TABLE8[6,1]),' & (',s2(TABLE8[6,2]),') & ',
      s2(TABLE8[6,3]),' & (',s2(TABLE8[6,4]),') & & ',
      s2(TABLE8[6,6]),' & (',s2(TABLE8[6,7]),') & & ',
      s2(TABLE8[6,9]),' & (',s2(TABLE8[6,10]),') & \\\\ \n',sep='',
      file=out.file.name,append=TRUE)

  cat('\\,\\,--equal  & mm & ',
      s2(TABLE8[7,1]),' & (',s2(TABLE8[7,2]),') & ',
      s2(TABLE8[7,3]),' & (',s2(TABLE8[7,4]),') & [',s2(TABLE8[7,5]),'] & ',
      s2(TABLE8[7,6]),' & (',s2(TABLE8[7,7]),') & [',s2(TABLE8[7,8]),'] & ',
      s2(TABLE8[7,9]),' & (',s2(TABLE8[7,10]),') & [',s2(TABLE8[7,11]),'] \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('               & m  & ',
      s2(TABLE8[8,1]),' & (',s2(TABLE8[8,2]),') & ',
      s2(TABLE8[8,3]),' & (',s2(TABLE8[8,4]),') & & ',
      s2(TABLE8[8,6]),' & (',s2(TABLE8[8,7]),') & & ',
      s2(TABLE8[8,9]),' & (',s2(TABLE8[8,10]),') & \\\\ \n',sep='',
      file=out.file.name,append=TRUE)

  cat('\\,\\,--free T\\&F  & mm & ',
      s2(TABLE8[9,1]),' & (',s2(TABLE8[9,2]),') & ',
      s2(TABLE8[9,3]),' & (',s2(TABLE8[9,4]),') & [',s2(TABLE8[9,5]),'] & ',
      s2(TABLE8[9,6]),' & (',s2(TABLE8[9,7]),') & [',s2(TABLE8[9,8]),'] & ',
      s2(TABLE8[9,9]),' & (',s2(TABLE8[9,10]),') & [',s2(TABLE8[9,11]),'] \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('               & m  & ',
      s2(TABLE8[10,1]),' & (',s2(TABLE8[10,2]),') & ',
      s2(TABLE8[10,3]),' & (',s2(TABLE8[10,4]),') & & ',
      s2(TABLE8[10,6]),' & (',s2(TABLE8[10,7]),') & & ',
      s2(TABLE8[10,9]),' & (',s2(TABLE8[10,10]),') & \\\\ \n',sep='',
      file=out.file.name,append=TRUE)

  cat('\\,\\,--equal T\\&F  & mm & ',
      s2(TABLE8[11,1]),' & (',s2(TABLE8[11,2]),') & ',
      s2(TABLE8[11,3]),' & (',s2(TABLE8[11,4]),') & [',s2(TABLE8[11,5]),'] & ',
      s2(TABLE8[11,6]),' & (',s2(TABLE8[11,7]),') & [',s2(TABLE8[11,8]),'] & ',
      s2(TABLE8[11,9]),' & (',s2(TABLE8[11,10]),') & [',s2(TABLE8[11,11]),'] \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('               & m  & ',
      s2(TABLE8[12,1]),' & (',s2(TABLE8[12,2]),') & ',
      s2(TABLE8[12,3]),' & (',s2(TABLE8[12,4]),') & & ',
      s2(TABLE8[12,6]),' & (',s2(TABLE8[12,7]),') & & ',
      s2(TABLE8[12,9]),' & (',s2(TABLE8[12,10]),') & \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('\\,\\,--P-P  & mm & ',
      s2(TABLE8[13,1]),' & (',s2(TABLE8[13,2]),') & ',
      s2(TABLE8[13,3]),' & (',s2(TABLE8[13,4]),') & [',s2(TABLE8[13,5]),'] & ',
      s2(TABLE8[13,6]),' & (',s2(TABLE8[13,7]),') & [',s2(TABLE8[13,8]),'] & ',
      s2(TABLE8[13,9]),' & (',s2(TABLE8[13,10]),') & [',s2(TABLE8[13,11]),'] \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('               & m  & ',
      s2(TABLE8[14,1]),' & (',s2(TABLE8[14,2]),') & ',
      s2(TABLE8[14,3]),' & (',s2(TABLE8[14,4]),') & & ',
      s2(TABLE8[14,6]),' & (',s2(TABLE8[14,7]),') & & ',
      s2(TABLE8[14,9]),' & (',s2(TABLE8[14,10]),') & \\\\ \n',sep='',
      file=out.file.name,append=TRUE)

  cat('\\hline\\hline\n',file=out.file.name,append=TRUE)
  mean8 <- (sum(TABLE8[,1])/2 + sum(TABLE8[,3]) + sum(TABLE8[,6]) + sum(TABLE8[,9]))/49
  cat('\nMean of all non-repeated estimates in Table 8 =',sprintf('%-.2f\n',mean8),file=out.file.name,append=TRUE)
}

# TWO-STAGE RANDOM-EFFECTS META-REGRESSION (2SREMR):
if(TRUE){

  RSS.m  <- matrix(0,28,1)
  yhat.m <- matrix(0,28,1)
  
  # IMPORT DATA:
  {
    gs4_deauth()
    
    # Input data are in the Google Sheet file titled "EPA_VSL_metadata" at the
    # following link:
    data <- read_sheet('https://docs.google.com/spreadsheets/d/1wxWgCSZKYWBuX55i4vCCw-ZyldcUP3SfbdvkOnAqdvA/edit?usp=sharing',
                       sheet='meta-data')

    # Extract variables:
    ID              <- data$GroupID          # unique group id
    pubyear         <- data$PubYear          # year of publication
    datayear        <- data$DataYear         # year of data collection
    dollaryear      <- data$DollarYear       # US$ year of reported VSL estimates
    mnD             <- data$MeanDummy        # mean dummy variable
    mdD             <- matrix(1,length(mnD),1)-mnD # median dummy variable
    spD             <- data$SPDummy          # stated preference dummy variable
    samplesize      <- data$SampleSize       # sample size of original study
    Y               <- data$VSL              # VSL estimate [DataYear $US]
    SE              <- data$SE               # standard error [DataYear $US]
    income          <- data$Income / 10000   # [2013 $US]
    
    # Group ids:
    IDs <- unique(ID)
    
    # Number of groups:
    I <- length(IDs)
    
    # Assumed correlation among observations within studies (rho):
    rho <- matrix(rho.all,I,1)
    
    # Observations per group:
    J <- matrix(0,I,1)
    for(i in 1:I){
      J[i] <- sum(ID==IDs[i])
    }
    
    # Total number of observations:
    N <- length(Y)
    
    # Convert all estimates from DataYear $US to 2020$US using the CPI index:
    # Source: https://data.bls.gov/pdq/SurveyOutputServlet
    CPI <- read_sheet('https://docs.google.com/spreadsheets/d/1wxWgCSZKYWBuX55i4vCCw-ZyldcUP3SfbdvkOnAqdvA/edit?usp=sharing',
           sheet='CPI',range='A13:P122',col_names=FALSE)
    CPI      <- as.matrix(CPI)
    CPIyr    <- CPI[,1]
    CPIindex <- CPI[,14]
    
    Y2 <- Y
    SE2 <- SE
    
    Y2013 <- Y
    SE2013 <- SE
    
    income2 <- income
    for(i in 1:length(Y)){
      Y2[i]  <- Y[i]  * CPIindex[which(CPIyr==2020)] / CPIindex[which(CPIyr==dollaryear[i])]
      SE2[i] <- SE[i] * CPIindex[which(CPIyr==2020)] / CPIindex[which(CPIyr==dollaryear[i])]
      income2[i] <- income[i] * CPIindex[which(CPIyr==2020)] / CPIindex[which(CPIyr==2013)]
      Y2013[i]   <- Y[i]  * CPIindex[which(CPIyr==2013)] / CPIindex[which(CPIyr==dollaryear[i])]
      SE2013[i]  <- SE[i] * CPIindex[which(CPIyr==2013)] / CPIindex[which(CPIyr==dollaryear[i])]
    }
    Y <- Y2
    SE <- SE2
    income <- income2
    
    # Number of observations from parent group of each observation:
    m <- matrix(0,N,1)
    for(ij in 1:N){
      m[ij] <- sum(ID==ID[ij])
    }
    
    #if(case==5){data1<-cbind(Y,SE)}
    
  }

  for(version in c(9,10,11)){ # For Tables 9/12, 10/13, 11/14

    if(version==9 ){PET <- FALSE; PEESE <- FALSE}
    if(version==10){PET <- TRUE;  PEESE <- FALSE}
    if(version==11){PET <- FALSE; PEESE <- TRUE}

    # Models [specifications] to estimate:
    #
    # 0 constant only (i.e., meta-analysis not meta-regression model)
    # 1 SP dummy (to estimate a fixed effect for SP vs RP)
    # 2 SP dummy + median dummy + data year
    # 3 SP dummy + median dummy + income
    # 4 SP dummy + median dummy + data year + income
    # 5 SP dummy + median dummy + data year + SP dummy * data year
    # 6 SP dummy + median dummy + income + SP dummy * income
    #
    # Note: Set up each model so that the intercept corresponds to the average VSL
    # from hedonic wage studies in the most recent data year

    TABLEcon <- matrix(NA,24,7)
    TABLEunc <- matrix(NA,24,7)

    for(model in 0:6){

      if(model==0){
        if(version==8){data2 <- cbind(Y,SE)}
        X <- matrix(1,N,1);                    var.names <- 'constant'
        if(PET)  {X <- cbind(X,SE);            var.names <- c(var.names,'s.e.') }
        if(PEESE){X <- cbind(X,SE^2);          var.names <- c(var.names,'s.e.^2')}
      }
      if(model==1){
        X <- matrix(1,N,1);                    var.names <- 'constant'
        X <- cbind(X,spD);                     var.names <- c(var.names,'sp dummy')
        X <- cbind(X,mdD);                     var.names <- c(var.names,'md dummy')
        if(PET)  {X <- cbind(X,SE);            var.names <- c(var.names,'s.e.') }
        if(PEESE){X <- cbind(X,SE^2);          var.names <- c(var.names,'s.e.^2')}
      }
      if(model==2){
        X <- matrix(1,N,1);                    var.names <- 'constant'
        X <- cbind(X,spD);                     var.names <- c(var.names,'sp dummy')
        X <- cbind(X,mdD);                     var.names <- c(var.names,'md dummy')
        X <- cbind(X,datayear-mean(datayear));  var.names <- c(var.names,'year')
        if(PET)  {X <- cbind(X,SE);            var.names <- c(var.names,'s.e.') }
        if(PEESE){X <- cbind(X,SE^2);          var.names <- c(var.names,'s.e.^2')}
      }
      if(model==3){
        X <- matrix(1,N,1);                    var.names <- 'constant'
        X <- cbind(X,spD);                     var.names <- c(var.names,'sp dummy')
        X <- cbind(X,mdD);                     var.names <- c(var.names,'md dummy')
        X <- cbind(X,income-mean(income));     var.names <- c(var.names,'income')
        if(PET)  {X <- cbind(X,SE);            var.names <- c(var.names,'s.e.') }
        if(PEESE){X <- cbind(X,SE^2);          var.names <- c(var.names,'s.e.^2')}
      }
      if(model==4){
        X <- matrix(1,N,1);                    var.names <- 'constant'
        X <- cbind(X,spD);                     var.names <- c(var.names,'sp dummy')
        X <- cbind(X,mdD);                     var.names <- c(var.names,'md dummy')
        X <- cbind(X,datayear-mean(datayear));  var.names <- c(var.names,'year')
        X <- cbind(X,income-mean(income));     var.names <- c(var.names,'income')
        if(PET)  {X <- cbind(X,SE);            var.names <- c(var.names,'s.e.') }
        if(PEESE){X <- cbind(X,SE^2);          var.names <- c(var.names,'s.e.^2')}
      }
      if(model==5){
        X <- matrix(1,N,1);                    var.names <- 'constant'
        X <- cbind(X,spD);                     var.names <- c(var.names,'sp dummy')
        X <- cbind(X,mdD);                     var.names <- c(var.names,'md dummy')
        X <- cbind(X,datayear-mean(datayear)); var.names <- c(var.names,'year')
        X <- cbind(X,spD*(datayear-mean(datayear))); var.names <- c(var.names,'spD*year')
        if(PET)  {X <- cbind(X,SE);            var.names <- c(var.names,'s.e.') }
        if(PEESE){X <- cbind(X,SE^2);          var.names <- c(var.names,'s.e.^2')}
      }
      if(model==6){
        X <- matrix(1,N,1);                    var.names <- 'constant'
        X <- cbind(X,spD);                     var.names <- c(var.names,'sp dummy')
        X <- cbind(X,mdD);                     var.names <- c(var.names,'md dummy')
        X <- cbind(X,income-mean(income));     var.names <- c(var.names,'income')
        X <- cbind(X,spD*(income-mean(income))); var.names <- c(var.names,'spD*income')
        if(PET)  {X <- cbind(X,SE);            var.names <- c(var.names,'s.e.') }
        if(PEESE){X <- cbind(X,SE^2);          var.names <- c(var.names,'s.e.^2')}
      }

      K <- length(X[1,])

      outs        <- twosremrFun(Y,X,SE,ID,rho)
      bhat.OLS    <- outs$bhat.OLS
      se.ro.OLS   <- outs$se.ro.OLS
      bhat.UNC    <- outs$bhat.UNC
      se.ro.UNC   <- outs$se.ro.UNC
      sig.mu.UNC  <- mean(outs$sig.muhat.UNC)
      sig.eta.UNC <- outs$sig.etahat.UNC
      bhat.CON    <- outs$bhat.CON
      se.ro.CON   <- outs$se.ro.CON
      sig.mu.CON  <- mean(outs$sig.muhat.CON)
      sig.eta.CON <- outs$sig.etahat.CON
      R2.CON      <- outs$R2.CON
      R2.adj.CON  <- outs$R2.adj.CON
      R2.cv.CON   <- outs$R2.cv.CON
      R2.UNC      <- outs$R2.UNC
      R2.adj.UNC  <- outs$R2.adj.UNC
      R2.cv.UNC   <- outs$R2.cv.UNC
      
      # Save cross validation errors for model averaging:
      {
        if(model==0 & PET==0 & PEESE==0){RSS.m[1] <- sum(outs$e.cv.CON^2); yhat.m[1] <- bhat.CON[1]}
        if(model==0 & PET==0 & PEESE==0){RSS.m[2] <- sum(outs$e.cv.UNC^2); yhat.m[2] <- bhat.UNC[1]}
        if(model==0 & PET==0 & PEESE==1){RSS.m[3] <- sum(outs$e.cv.CON^2); yhat.m[3] <- bhat.CON[1]}
        if(model==0 & PET==0 & PEESE==1){RSS.m[4] <- sum(outs$e.cv.UNC^2); yhat.m[4] <- bhat.UNC[1]}
        
        if(model==1 & PET==0 & PEESE==0){RSS.m[5] <- sum(outs$e.cv.CON^2); yhat.m[5] <- bhat.CON[1]+0.5*bhat.CON[2]}
        if(model==1 & PET==0 & PEESE==0){RSS.m[6] <- sum(outs$e.cv.UNC^2); yhat.m[6] <- bhat.UNC[1]+0.5*bhat.UNC[2]}
        if(model==1 & PET==0 & PEESE==1){RSS.m[7] <- sum(outs$e.cv.CON^2); yhat.m[7] <- bhat.CON[1]+0.5*bhat.CON[2]}
        if(model==1 & PET==0 & PEESE==1){RSS.m[8] <- sum(outs$e.cv.UNC^2); yhat.m[8] <- bhat.UNC[1]+0.5*bhat.UNC[2]}
        
        if(model==2 & PET==0 & PEESE==0){RSS.m[9] <- sum(outs$e.cv.CON^2); yhat.m[9] <- bhat.CON[1]+0.5*bhat.CON[2]}
        if(model==2 & PET==0 & PEESE==0){RSS.m[10]<- sum(outs$e.cv.UNC^2); yhat.m[10]<- bhat.UNC[1]+0.5*bhat.UNC[2]}
        if(model==2 & PET==0 & PEESE==1){RSS.m[11]<- sum(outs$e.cv.CON^2); yhat.m[11]<- bhat.CON[1]+0.5*bhat.CON[2]}
        if(model==2 & PET==0 & PEESE==1){RSS.m[12]<- sum(outs$e.cv.UNC^2); yhat.m[12]<- bhat.UNC[1]+0.5*bhat.UNC[2]}
        
        if(model==3 & PET==0 & PEESE==0){RSS.m[13]<- sum(outs$e.cv.CON^2); yhat.m[13]<- bhat.CON[1]+0.5*bhat.CON[2]}
        if(model==3 & PET==0 & PEESE==0){RSS.m[14]<- sum(outs$e.cv.UNC^2); yhat.m[14]<- bhat.UNC[1]+0.5*bhat.UNC[2]}
        if(model==3 & PET==0 & PEESE==1){RSS.m[15]<- sum(outs$e.cv.CON^2); yhat.m[15]<- bhat.CON[1]+0.5*bhat.CON[2]}
        if(model==3 & PET==0 & PEESE==1){RSS.m[16]<- sum(outs$e.cv.UNC^2); yhat.m[16]<- bhat.UNC[1]+0.5*bhat.UNC[2]}
        
        if(model==4 & PET==0 & PEESE==0){RSS.m[17]<- sum(outs$e.cv.CON^2); yhat.m[17]<- bhat.CON[1]+0.5*bhat.CON[2]}
        if(model==4 & PET==0 & PEESE==0){RSS.m[18]<- sum(outs$e.cv.UNC^2); yhat.m[18]<- bhat.UNC[1]+0.5*bhat.UNC[2]}
        if(model==4 & PET==0 & PEESE==1){RSS.m[19]<- sum(outs$e.cv.CON^2); yhat.m[19]<- bhat.CON[1]+0.5*bhat.CON[2]}
        if(model==4 & PET==0 & PEESE==1){RSS.m[20]<- sum(outs$e.cv.UNC^2); yhat.m[20]<- bhat.UNC[1]+0.5*bhat.UNC[2]}
        
        if(model==5 & PET==0 & PEESE==0){RSS.m[21]<- sum(outs$e.cv.CON^2); yhat.m[21]<- bhat.CON[1]+0.5*bhat.CON[2]}
        if(model==5 & PET==0 & PEESE==0){RSS.m[22]<- sum(outs$e.cv.UNC^2); yhat.m[22]<- bhat.UNC[1]+0.5*bhat.UNC[2]}
        if(model==5 & PET==0 & PEESE==1){RSS.m[23]<- sum(outs$e.cv.CON^2); yhat.m[23]<- bhat.CON[1]+0.5*bhat.CON[2]}
        if(model==5 & PET==0 & PEESE==1){RSS.m[24]<- sum(outs$e.cv.UNC^2); yhat.m[24]<- bhat.UNC[1]+0.5*bhat.UNC[2]}
        
        if(model==6 & PET==0 & PEESE==0){RSS.m[25]<- sum(outs$e.cv.CON^2); yhat.m[25]<- bhat.CON[1]+0.5*bhat.CON[2]}
        if(model==6 & PET==0 & PEESE==0){RSS.m[26]<- sum(outs$e.cv.UNC^2); yhat.m[26]<- bhat.UNC[1]+0.5*bhat.UNC[2]}
        if(model==6 & PET==0 & PEESE==1){RSS.m[27]<- sum(outs$e.cv.CON^2); yhat.m[27]<- bhat.CON[1]+0.5*bhat.CON[2]}
        if(model==6 & PET==0 & PEESE==1){RSS.m[28]<- sum(outs$e.cv.UNC^2); yhat.m[28]<- bhat.UNC[1]+0.5*bhat.UNC[2]}
      }

      kk <- 0
      if(sum(var.names=='income')>0){
        kk <- which(var.names=='income')
        IE.OLS <- bhat.OLS[kk] * mean(income) / mean(Y)
        IE.UNC <- bhat.UNC[kk] * mean(income) / mean(Y)
        IE.CON <- bhat.CON[kk] * mean(income) / mean(Y)
      }else{
        IE.OLS <- -99
        IE.UNC <- -99
        IE.CON <- -99
      }

      # Bootstrap standard errors:
      bhat.OLS.BS <- matrix(0,B,K)
      bhat.UNC.BS <- matrix(0,B,K)
      bhat.CON.BS <- matrix(0,B,K)
      IE.OLS.BS   <- matrix(0,B,1)
      IE.UNC.BS   <- matrix(0,B,1)
      IE.CON.BS   <- matrix(0,B,1)
      if(FALSE)for(b in 1:B){

        cat('\014')
        cat('Working on bootstrap rep',sprintf('%-.0f',b),'of',sprintf('%-.0f',B))

        # Re-sample data with replacement:
        done <- 0
        while(done==0){
          BSgroups <- sample(IDs,I,replace=TRUE)
          z   <- which(ID==BSgroups[1])
          Yb  <- Y[z]
          SEb <- SE[z]
          Xb  <- X[z,]
          IDb <- matrix(1,length(Yb),1)
          for(i in 2:I){
            z   <- which(ID==BSgroups[i])
            Yb  <- c(Yb,Y[z])
            SEb <- c(SEb,SE[z])
            Xb  <- rbind(Xb,X[z,])
            IDb <- rbind(IDb,matrix(i,length(z),1))
          }
          if (rcond(Xb)>10^-7){done <- 1}
        }
        outs <- twosremrFun(Yb,Xb,SEb,IDb,0)
        bhat.OLS.BS[b,] <- outs$bhat.OLS
        bhat.UNC.BS[b,] <- outs$bhat.UNC
        bhat.CON.BS[b,] <- outs$bhat.CON

        # IE.OLS.BS[b]    <- bhat.OLS.BS[b,2] * mean(Xb[,2]) / mean(Yb)
        # IE.UNC.BS[b]    <- bhat.UNC.BS[b,2] * mean(Xb[,2]) / mean(Yb)
        # IE.CON.BS[b]    <- bhat.CON.BS[b,2] * mean(Xb[,2]) / mean(Yb)

      }

      se.bs.OLS <- apply(bhat.OLS.BS,2,sd)
      se.bs.UNC <- apply(bhat.UNC.BS,2,sd)
      se.bs.CON <- apply(bhat.CON.BS,2,sd)

      if(model==0){

        TABLEcon[1,1]  <- bhat.CON[1]
        TABLEcon[2,1]  <- se.ro.CON[1]
        TABLEcon[19,1] <- sig.mu.CON
        TABLEcon[20,1] <- sig.eta.CON
        TABLEcon[23,1] <- R2.CON
        TABLEcon[24,1] <- R2.cv.CON
        if(PET){
          TABLEcon[15,1] <- bhat.CON[2]
          TABLEcon[16,1] <- se.ro.CON[2]
        }
        if(PEESE){
          TABLEcon[17,1] <- bhat.CON[2]
          TABLEcon[18,1] <- se.ro.CON[2]
        }

        TABLEunc[1,1]  <- bhat.UNC[1]
        TABLEunc[2,1]  <- se.ro.UNC[1]
        TABLEunc[19,1] <- sig.mu.UNC
        TABLEunc[20,1] <- sig.eta.UNC
        TABLEunc[23,1] <- R2.UNC
        TABLEunc[24,1] <- R2.cv.UNC
        if(PET){
          TABLEunc[15,1] <- bhat.UNC[2]
          TABLEunc[16,1] <- se.ro.UNC[2]
        }
        if(PEESE){
          TABLEunc[17,1] <- bhat.UNC[2]
          TABLEunc[18,1] <- se.ro.UNC[2]
        }

      }
      if(model==1){

        TABLEcon[1,2]  <- bhat.CON[1]
        TABLEcon[2,2]  <- se.ro.CON[1]
        TABLEcon[3,2]  <- bhat.CON[2]
        TABLEcon[4,2]  <- se.ro.CON[2]
        TABLEcon[5,2]  <- bhat.CON[3]
        TABLEcon[6,2]  <- se.ro.CON[3]
        TABLEcon[19,2] <- sig.mu.CON
        TABLEcon[20,2] <- sig.eta.CON
        TABLEcon[23,2] <- R2.CON
        TABLEcon[24,2] <- R2.cv.CON
        if(PET){
          TABLEcon[15,2] <- bhat.CON[4]
          TABLEcon[16,2] <- se.ro.CON[4]
        }
        if(PEESE){
          TABLEcon[17,2] <- bhat.CON[4]
          TABLEcon[18,2] <- se.ro.CON[4]
        }

        TABLEunc[1,2]  <- bhat.UNC[1]
        TABLEunc[2,2]  <- se.ro.UNC[1]
        TABLEunc[3,2]  <- bhat.UNC[2]
        TABLEunc[4,2]  <- se.ro.UNC[2]
        TABLEunc[5,2]  <- bhat.UNC[3]
        TABLEunc[6,2]  <- se.ro.UNC[3]
        TABLEunc[19,2] <- sig.mu.UNC
        TABLEunc[20,2] <- sig.eta.UNC
        TABLEunc[23,2] <- R2.UNC
        TABLEunc[24,2] <- R2.cv.UNC
        if(PET){
          TABLEunc[15,2] <- bhat.UNC[4]
          TABLEunc[16,2] <- se.ro.UNC[4]
        }
        if(PEESE){
          TABLEunc[17,2] <- bhat.UNC[4]
          TABLEunc[18,2] <- se.ro.UNC[4]
        }

      }
      if(model==2){

        TABLEcon[1,3]  <- bhat.CON[1]
        TABLEcon[2,3]  <- se.ro.CON[1]
        TABLEcon[3,3]  <- bhat.CON[2]
        TABLEcon[4,3]  <- se.ro.CON[2]
        TABLEcon[5,3]  <- bhat.CON[3]
        TABLEcon[6,3]  <- se.ro.CON[3]
        TABLEcon[7,3]  <- bhat.CON[4]
        TABLEcon[8,3]  <- se.ro.CON[4]
        TABLEcon[19,3] <- sig.mu.CON
        TABLEcon[20,3] <- sig.eta.CON
        TABLEcon[23,3] <- R2.CON
        TABLEcon[24,3] <- R2.cv.CON
        if(PET){
          TABLEcon[15,3] <- bhat.CON[5]
          TABLEcon[16,3] <- se.ro.CON[5]
        }
        if(PEESE){
          TABLEcon[17,3] <- bhat.UNC[5]
          TABLEcon[18,3] <- se.ro.UNC[5]
        }

        TABLEunc[1,3]  <- bhat.UNC[1]
        TABLEunc[2,3]  <- se.ro.UNC[1]
        TABLEunc[3,3]  <- bhat.UNC[2]
        TABLEunc[4,3]  <- se.ro.UNC[2]
        TABLEunc[5,3]  <- bhat.UNC[3]
        TABLEunc[6,3]  <- se.ro.UNC[3]
        TABLEunc[7,3]  <- bhat.UNC[4]
        TABLEunc[8,3]  <- se.ro.UNC[4]
        TABLEunc[19,3] <- sig.mu.UNC
        TABLEunc[20,3] <- sig.eta.UNC
        TABLEunc[23,3] <- R2.UNC
        TABLEunc[24,3] <- R2.cv.UNC
        if(PET){
          TABLEunc[15,3] <- bhat.UNC[5]
          TABLEunc[16,3] <- se.ro.UNC[5]
        }
        if(PEESE){
          TABLEunc[17,3] <- bhat.UNC[5]
          TABLEunc[18,3] <- se.ro.UNC[5]
        }

      }
      if(model==3){

        TABLEcon[1,4]  <- bhat.CON[1]
        TABLEcon[2,4]  <- se.ro.CON[1]
        TABLEcon[3,4]  <- bhat.CON[2]
        TABLEcon[4,4]  <- se.ro.CON[2]
        TABLEcon[5,4]  <- bhat.CON[3]
        TABLEcon[6,4]  <- se.ro.CON[3]
        TABLEcon[9,4]  <- bhat.CON[4]
        TABLEcon[10,4] <- se.ro.CON[4]
        TABLEcon[19,4] <- sig.mu.CON
        TABLEcon[20,4] <- sig.eta.CON
        TABLEcon[21,4] <- IE.CON
        TABLEcon[22,4] <- se.ro.CON[kk]*mean(income)/mean(Y)
        TABLEcon[23,4] <- R2.CON
        TABLEcon[24,4] <- R2.cv.CON
        if(PET){
          TABLEcon[15,4] <- bhat.CON[5]
          TABLEcon[16,4] <- se.ro.CON[5]
        }
        if(PEESE){
          TABLEcon[17,4] <- bhat.CON[5]
          TABLEcon[18,4] <- se.ro.CON[5]
        }

        TABLEunc[1,4]  <- bhat.UNC[1]
        TABLEunc[2,4]  <- se.ro.UNC[1]
        TABLEunc[3,4]  <- bhat.UNC[2]
        TABLEunc[4,4]  <- se.ro.UNC[2]
        TABLEunc[5,4]  <- bhat.UNC[3]
        TABLEunc[6,4]  <- se.ro.UNC[3]
        TABLEunc[9,4]  <- bhat.UNC[4]
        TABLEunc[10,4] <- se.ro.UNC[4]
        TABLEunc[19,4] <- sig.mu.UNC
        TABLEunc[20,4] <- sig.eta.UNC
        TABLEunc[21,4] <- IE.UNC
        TABLEunc[22,4] <- se.ro.UNC[kk]*mean(income)/mean(Y)
        TABLEunc[23,4] <- R2.UNC
        TABLEunc[24,4] <- R2.cv.UNC
        if(PET){
          TABLEunc[15,4] <- bhat.UNC[5]
          TABLEunc[16,4] <- se.ro.UNC[5]
        }
        if(PEESE){
          TABLEunc[17,4] <- bhat.UNC[5]
          TABLEunc[18,4] <- se.ro.UNC[5]
        }

      }
      if(model==4){

        TABLEcon[1,5]  <- bhat.CON[1]
        TABLEcon[2,5]  <- se.ro.CON[1]
        TABLEcon[3,5]  <- bhat.CON[2]
        TABLEcon[4,5]  <- se.ro.CON[2]
        TABLEcon[5,5]  <- bhat.CON[3]
        TABLEcon[6,5]  <- se.ro.CON[3]
        TABLEcon[7,5]  <- bhat.CON[4]
        TABLEcon[8,5]  <- se.ro.CON[4]
        TABLEcon[9,5]  <- bhat.CON[5]
        TABLEcon[10,5]  <- se.ro.CON[5]
        TABLEcon[19,5] <- sig.mu.CON
        TABLEcon[20,5] <- sig.eta.CON
        TABLEcon[21,5] <- IE.CON
        TABLEcon[22,5] <- se.ro.CON[kk]*mean(income)/mean(Y)
        TABLEcon[23,5] <- R2.CON
        TABLEcon[24,5] <- R2.cv.CON
        if(PET){
          TABLEcon[15,5] <- bhat.CON[6]
          TABLEcon[16,5] <- se.ro.CON[6]
        }
        if(PEESE){
          TABLEcon[17,5] <- bhat.CON[6]
          TABLEcon[18,5] <- se.ro.CON[6]
        }

        TABLEunc[1,5]  <- bhat.UNC[1]
        TABLEunc[2,5]  <- se.ro.UNC[1]
        TABLEunc[3,5]  <- bhat.UNC[2]
        TABLEunc[4,5]  <- se.ro.UNC[2]
        TABLEunc[5,5]  <- bhat.UNC[3]
        TABLEunc[6,5]  <- se.ro.UNC[3]
        TABLEunc[7,5]  <- bhat.UNC[4]
        TABLEunc[8,5]  <- se.ro.UNC[4]
        TABLEunc[9,5]  <- bhat.UNC[5]
        TABLEunc[10,5]  <- se.ro.UNC[5]
        TABLEunc[19,5] <- sig.mu.UNC
        TABLEunc[20,5] <- sig.eta.UNC
        TABLEunc[21,5] <- IE.UNC
        TABLEunc[22,5] <- se.ro.UNC[kk]*mean(income)/mean(Y)
        TABLEunc[23,5] <- R2.UNC
        TABLEunc[24,5] <- R2.cv.UNC
        if(PET){
          TABLEunc[15,5] <- bhat.UNC[6]
          TABLEunc[16,5] <- se.ro.UNC[6]
        }
        if(PEESE){
          TABLEunc[17,5] <- bhat.UNC[6]
          TABLEunc[18,5] <- se.ro.UNC[6]
        }

      }
      if(model==5){

        TABLEcon[1,6]  <- bhat.CON[1]
        TABLEcon[2,6]  <- se.ro.CON[1]
        TABLEcon[3,6]  <- bhat.CON[2]
        TABLEcon[4,6]  <- se.ro.CON[2]
        TABLEcon[5,6]  <- bhat.CON[3]
        TABLEcon[6,6]  <- se.ro.CON[3]
        TABLEcon[7,6]  <- bhat.CON[4]
        TABLEcon[8,6]  <- se.ro.CON[4]
        TABLEcon[11,6] <- bhat.CON[5]
        TABLEcon[12,6] <- se.ro.CON[5]
        TABLEcon[19,6] <- sig.mu.CON
        TABLEcon[20,6] <- sig.eta.CON
        TABLEcon[23,6] <- R2.CON
        TABLEcon[24,6] <- R2.cv.CON
        if(PET){
          TABLEcon[15,6] <- bhat.CON[6]
          TABLEcon[16,6] <- se.ro.CON[6]
        }
        if(PEESE){
          TABLEcon[17,6] <- bhat.CON[6]
          TABLEcon[18,6] <- se.ro.CON[6]
        }

        TABLEunc[1,6]  <- bhat.UNC[1]
        TABLEunc[2,6]  <- se.ro.UNC[1]
        TABLEunc[3,6]  <- bhat.UNC[2]
        TABLEunc[4,6]  <- se.ro.UNC[2]
        TABLEunc[5,6]  <- bhat.UNC[3]
        TABLEunc[6,6]  <- se.ro.UNC[3]
        TABLEunc[7,6]  <- bhat.UNC[4]
        TABLEunc[8,6]  <- se.ro.UNC[4]
        TABLEunc[11,6] <- bhat.UNC[5]
        TABLEunc[12,6] <- se.ro.UNC[5]
        TABLEunc[19,6] <- sig.mu.UNC
        TABLEunc[20,6] <- sig.eta.UNC
        TABLEunc[23,6] <- R2.UNC
        TABLEunc[24,6] <- R2.cv.UNC
        if(PET){
          TABLEunc[15,6] <- bhat.UNC[6]
          TABLEunc[16,6] <- se.ro.UNC[6]
        }
        if(PEESE){
          TABLEunc[17,6] <- bhat.UNC[6]
          TABLEunc[18,6] <- se.ro.UNC[6]
        }

      }
      if(model==6){

        TABLEcon[1,7]  <- bhat.CON[1]
        TABLEcon[2,7]  <- se.ro.CON[1]
        TABLEcon[3,7]  <- bhat.CON[2]
        TABLEcon[4,7]  <- se.ro.CON[2]
        TABLEcon[5,7]  <- bhat.CON[3]
        TABLEcon[6,7]  <- se.ro.CON[3]
        TABLEcon[9,7]  <- bhat.CON[4]
        TABLEcon[10,7]  <- se.ro.CON[4]
        TABLEcon[13,7] <- bhat.CON[5]
        TABLEcon[14,7] <- se.ro.CON[5]
        TABLEcon[19,7] <- sig.mu.CON
        TABLEcon[20,7] <- sig.eta.CON
        TABLEcon[21,7] <- IE.CON
        TABLEcon[22,7] <- se.ro.CON[kk]*mean(income)/mean(Y)
        TABLEcon[23,7] <- R2.CON
        TABLEcon[24,7] <- R2.cv.CON
        if(PET){
          TABLEcon[15,7] <- bhat.CON[6]
          TABLEcon[16,7] <- se.ro.CON[6]
        }
        if(PEESE){
          TABLEcon[17,7] <- bhat.CON[6]
          TABLEcon[18,7] <- se.ro.CON[6]
        }

        TABLEunc[1,7]  <- bhat.UNC[1]
        TABLEunc[2,7]  <- se.ro.UNC[1]
        TABLEunc[3,7]  <- bhat.UNC[2]
        TABLEunc[4,7]  <- se.ro.UNC[2]
        TABLEunc[5,7]  <- bhat.UNC[3]
        TABLEunc[6,7]  <- se.ro.UNC[3]
        TABLEunc[9,7]  <- bhat.UNC[4]
        TABLEunc[10,7]  <- se.ro.UNC[4]
        TABLEunc[13,7] <- bhat.UNC[5]
        TABLEunc[14,7] <- se.ro.UNC[5]
        TABLEunc[19,7] <- sig.mu.UNC
        TABLEunc[20,7] <- sig.eta.UNC
        TABLEunc[21,7] <- IE.UNC
        TABLEunc[22,7] <- se.ro.UNC[kk]*mean(income)/mean(Y)
        TABLEunc[23,7] <- R2.UNC
        TABLEunc[24,7] <- R2.cv.UNC
        if(PET){
          TABLEunc[15,7] <- bhat.UNC[6]
          TABLEunc[16,7] <- se.ro.UNC[6]
        }
        if(PEESE){
          TABLEunc[17,7] <- bhat.UNC[6]
          TABLEunc[18,7] <- se.ro.UNC[6]
        }

      }

    }

    # WRITE TABLEcon [9/10/11] TO OUTPUT FILE:
    if(TRUE){
      cat('\nTable ',sprintf('%-.0f',version),'\n',sep='',file=out.file.name,append=TRUE)
      cat('\\hline\\hline\n',file=out.file.name,append=TRUE)
      cat(' & S0 & S1 & S2 & S3 & S4 & S5 & S6 \\\\ \n',file=out.file.name,append=TRUE)
      cat('\\hline\n',file=out.file.name,append=TRUE)

      cat('constant ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[1,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEcon[2,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('SP ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[3,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEcon[4,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('median ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[5,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEcon[6,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('year ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[7,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEcon[8,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('income ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[9,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEcon[10,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('SP$\\times$year ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[11,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEcon[12,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('SP$\\times$income ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[13,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEcon[14,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('\\textit{se} ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[15,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEcon[16,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('\\textit{se}$^2$ ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[17,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEcon[18,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('$\\sigma_{\\mu}$ ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[19,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('$\\sigma_{\\eta}$',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[20,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('\\textit{IEVSL} ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[21,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEcon[22,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('$R^2$ ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[23,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('$R_{CV}^2$ ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEcon[24,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('\\hline\\hline\n',file=out.file.name,append=TRUE)
    }

    # WRITE TABLEunc [12/13/14] TO OUTPUT FILE:
    if(TRUE){
      cat('\nTable ',sprintf('%-.0f',version+3),'\n',sep='',file=out.file.name,append=TRUE)
      cat('\\hline\\hline\n',file=out.file.name,append=TRUE)
      cat(' & 0 & 1 & 2 & 3 & 4 & 5 & 6 \\\\ \n',file=out.file.name,append=TRUE)
      cat('\\hline\n',file=out.file.name,append=TRUE)

      cat('constant ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[1,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEunc[2,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('SP ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[3,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEunc[4,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('median ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[5,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEunc[6,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('year ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[7,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEunc[8,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('income ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[9,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEunc[10,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('SP$\\times$year ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[11,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEunc[12,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('SP$\\times$income ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[13,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEunc[14,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('\\textit{se} ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[15,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEunc[16,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('\\textit{se}$^2$ ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[17,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEunc[18,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('$\\sigma_{\\mu}$ ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[19,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('$\\sigma_{\\eta}$',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[20,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('\\textit{IEVSL} ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[21,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)
      cat('            ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & (',s3(TABLEunc[22,k]),')',sep='',file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('$R^2$ ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[23,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('$R_{CV}^2$ ',file=out.file.name,append=TRUE)
      for(k in 1:7){cat(' & ',s3(TABLEunc[24,k]),file=out.file.name,append=TRUE)}
      cat('\\\\ \n',file=out.file.name,append=TRUE)

      cat('\\hline\\hline\n',file=out.file.name,append=TRUE)
    }

  }

  # JACKKNIFE MODEL AVERAGING:
  {
    M <- 28
    w.m <- matrix(1/M,M,1)
    # w.m <- w.m*0; w.m[5] <- 1
    done <- 0
    while(done==0){
      w.mm <- (1/M) * sum(w.m*RSS.m) / RSS.m
      w.mm <- w.mm/sum(w.mm)
      if(max(abs(w.mm-w.m))<1e-6){done <- 1}else{w.m <- w.mm}
    }
    yhat.ma <- sum(yhat.m*w.m)
    cat('yhat.ma =',sprintf('%-.3f\n',yhat.ma))
  }
  
}


