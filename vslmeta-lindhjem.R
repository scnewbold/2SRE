#===============================================================================
# This is an R script named 'vslmeta-lindhjem'
# which will replicate the results of our application of the 2SRE meta-
# analysis estimator to the Lindhjem et al. (2011) dataset as reported in
# Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. (2025)
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
  # out.file.name <- paste(output.path,'/',script.name,'-',date.time,'.out',sep='')
  out.file.name <- paste(output.path,'/',script.name,'.out',sep='')
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
  source("tictocFun.R")
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

tictocFun('tic')

# Assumed correlation among observations within studies (rho):
# Our estimation approach does not allow estimation of rho, so we must
# assume a value for this parameter. We can test the robustness of our results
# by varying this assumption and reporting the associated VSL estimates.
rho.all <- 0.5

# Replace any negative g's with 0 and re-scale remaining g's (1=yes,0=no):
g.scale <- 1

# In the paper, meta-analysis results are reported in Table 8.
TABLE8b <- matrix(0,7,2)
for(case in 1:1){ # only case 1 applies to Lindhjem data bc only SP and only means
  
  set.seed(1234)
  B <- 10 # Number of bootstrap reps [1000]

  # IMPORT DATA:
  {
    gs4_deauth()
    
    # Input data are in the Google Sheet file titled "EPA_VSL_metadata" at the
    # following link:
    data <- read_sheet('https://docs.google.com/spreadsheets/d/1wxWgCSZKYWBuX55i4vCCw-ZyldcUP3SfbdvkOnAqdvA/edit?usp=sharing',
                       sheet='lindhjem-meta-data',range="A1:M1172")

    for(j in 1:dim(data)[2]){
      data[which(is.na(data[,j])),j] <- -99
    }
    
    # Extract variables:
    ID              <- data$surveyid     # unique group id
    datayear        <- data$collyear     # year of data collection
    Y               <- data$vslaic       # VSL estimate [$US]
    SE              <- data$sd_vslaic    # standard error [$US]
    income          <- data$income_aic   # [$US]
    include         <- data$include_2025 # Our best guess of included estimated based on close reading of Lindhjem et al. (2011)
    
    # Replicate Table 1 in Lindhjem et al 2011:
    {
      Yfull <- Y[which(include == 1)]
      Ytrim <- Yfull[which(Yfull > quantile(Yfull,.025) & Yfull < quantile(Yfull,.975))]
      
      IDfull <- ID[which(include == 1)]
      IDtrim <- IDfull[which(Yfull > quantile(Yfull,.025) & Yfull < quantile(Yfull,.975))]
      
      IDs.full <- unique(IDfull)
      Ystudy.full <- matrix(0,length(IDs.full),1)
      for(j in 1:length(IDs.full)){
        Ystudy.full[j] <- mean(Yfull[which(IDfull==IDs.full[j])])
      }
      
      IDs.trim <- unique(IDtrim)
      Ystudy.trim <- matrix(0,length(IDs.trim),1)
      for(j in 1:length(IDs.trim)){
        Ystudy.trim[j] <- mean(Ytrim[which(IDtrim==IDs.trim[j])])
      }
      
      cat("\nReplication of Table 1 from Lindhjem et al. (2011):\n",file=out.file.name,append=TRUE)
      cat("==============================================\n",file=out.file.name,append=TRUE)
      cat("                   Full Sample  Trimmed Sample\n",file=out.file.name,append=TRUE)
      cat("----------------------------------------------\n",file=out.file.name,append=TRUE)
      cat(sprintf("Mean VSL            %10.0f      %10.0f\n",mean(Yfull),mean(Ytrim)),file=out.file.name,append=TRUE)
      cat(sprintf("Weighted mean VSL   %10.0f      %10.0f\n",mean(Ystudy.full),mean(Ystudy.trim)),file=out.file.name,append=TRUE)
      cat(sprintf("Median              %10.0f      %10.0f\n",median(Yfull),median(Ytrim)),file=out.file.name,append=TRUE)
      cat(sprintf("Minimum value       %10.0f      %10.0f\n",min(Yfull),min(Ytrim)),file=out.file.name,append=TRUE)
      cat(sprintf("Maximum value       %10.0f      %10.0f\n",max(Yfull),max(Ytrim)),file=out.file.name,append=TRUE)
      cat(sprintf("Number of estimates %10.0f      %10.0f\n",length(Yfull),length(Ytrim)),file=out.file.name,append=TRUE)
      cat("==============================================\n",file=out.file.name,append=TRUE)
    }

    # Filter:
    keep     <- which(include==1 & Y!=-99 & SE!=-99)
    ID       <- ID[keep]
    datayear <- datayear[keep]
    Y        <- Y[keep]/1e6
    SE       <- SE[keep]/1e6
    income   <- income[keep]/1e4
    
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
    
    # Number of observations from parent group of each observation:
    m <- matrix(0,N,1)
    for(ij in 1:N){
      m[ij] <- sum(ID==ID[ij])
    }

  }
  
  # TWO-STAGE RANDOM-EFFECTS META-ANALYSIS (no covariates) (2SREMA):
  if(TRUE){
    
    Y0  <- Y
    SE0 <- SE
    ID0 <- ID
    
    # estimation and bootstrapped standard errors:
    {
      
      outs     <- twosremaFun(Y,SE,ID,rho)
      yhat.sm  <- outs$yhat.sm
      yhat.mm  <- outs$yhat.mm
      yhat.ru  <- outs$yhat.ru
      yhat.rc  <- outs$yhat.rc
      sero.ru  <- outs$se.ru
      sero.rc  <- outs$se.rc
      
      # outs     <- twosremrFun(Y,matrix(1,length(Y),1),SE,ID,rho)
      
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
        done <- FALSE
        while(done==FALSE){
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
          done <- TRUE
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
      done <- FALSE
      while(done==FALSE){
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
        if((outs$fail.c+outs$fail.u)==0){done <- TRUE} # Discards re-sampled data sets that won't estimate.
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
    bias.PP <- outs$pub.bias
    
    # Bootstrap standard errors:
    yhat.PP.BS <- matrix(0,B,1)
    bias.PP.BS <- matrix(0,B,1)
    for(b in 1:B){
      
      if(floor(b/10)==b/10){
        cat('\014')
        cat('Case ',sprintf('%-.0f',case),'-- PET-PEESE\n')
        cat('Working on bootstrap rep',sprintf('%-.0f',b),'of',sprintf('%-.0f',B))
      }
      
      # Resample data with replacement:
      done <- FALSE
      while(done==FALSE){
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
        if(outs$fail==0){done <- TRUE} # Discards re-sampled data sets that won't estimate.
      }
      yhat.PP.BS[b] <- outs$yhat.PP
      bias.PP.BS[b] <- outs$pub.bias
      
    }
    sebs.PP <- sd(yhat.PP.BS)
    bias.PP <- mean(bias.PP.BS)
    
  }
  
  # SAVE RESULTS TO TABLE8b MATRIX
  {
    TABLE8b[1,1] <-yhat.sm; TABLE8b[1,2] <-sebs.sm
    TABLE8b[2,1] <-yhat.mm; TABLE8b[2,2] <-sebs.mm
    TABLE8b[3,1] <-yhat.ru; TABLE8b[3,2] <-sebs.ru
    TABLE8b[4,1] <-yhat.rc; TABLE8b[4,2] <-sebs.rc
    TABLE8b[5,1] <-yhat.TFu;TABLE8b[5,2] <-sebs.TFu
    TABLE8b[6,1] <-yhat.TFc;TABLE8b[6,2] <-sebs.TFc
    TABLE8b[7,1] <-yhat.PP; TABLE8b[7,2] <-sebs.PP
  }
  
}

# WRITE TABLE8b TO OUTPUT FILE:
if(TRUE){
  cat('\nTable 8b\n',file=out.file.name,append=TRUE)
  cat('\\hline\\hline\n',file=out.file.name,append=TRUE)
  cat('Estimator & estimate & s.e. \\\\ \n',file=out.file.name,append=TRUE)
  cat('\\hline\n',file=out.file.name,append=TRUE)
  cat('simple mean    & ',
      s2(TABLE8b[1,1]),' & (',s2(TABLE8b[1,2]),') \\\\ \n',sep='',
      file=out.file.name,append=TRUE)

  cat('group means    & ',
      s2(TABLE8b[2,1]),' & (',s2(TABLE8b[2,2]),') \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('2SRE--free     & ',
      s2(TABLE8b[3,1]),' & (',s2(TABLE8b[3,2]),') \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('\\,\\,--equal  & ',
      s2(TABLE8b[4,1]),' & (',s2(TABLE8b[4,2]),') \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('\\,\\,--free T\\&F  & ',
      s2(TABLE8b[5,1]),' & (',s2(TABLE8b[5,2]),') \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('\\,\\,--equal T\\&F  & ',
      s2(TABLE8b[6,1]),' & (',s2(TABLE8b[6,2]),') \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('\\,\\,--P-P  & ',
      s2(TABLE8b[7,1]),' & (',s2(TABLE8b[7,2]),') \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('freq bias & ',s2(bias.PP),' & \\\\ \n',sep='',
      file=out.file.name,append=TRUE)
  
  cat('\\hline\\hline\n',file=out.file.name,append=TRUE)
}

tictocFun('toc')


