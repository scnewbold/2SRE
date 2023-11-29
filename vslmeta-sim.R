#===============================================================================
# vslmeta-sim.R:
# This script will replicate the results of our application of the 2SRE meta-
# analysis estimator to the constructed datasets as reported in
# Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. (2023)
#
# [calls vslmeta-sim-iteration.R]
#===============================================================================

#===============================================================================
# PRELIMINARIES:
#===============================================================================
if(TRUE){
  
  # Clear environment to start fresh
  rm(list=ls()) 
  
  # Grab script name for file handling:
  script.name <- basename(rstudioapi::getSourceEditorContext()$path) 
  script.name <- gsub(".R","",script.name)
  
  # Packages
  list.of.packages <-
    c('metafor',
      'robumeta',
      'MAd')

  new.packages <- list.of.packages[!(list.of.packages %in%
                                       installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages,function(x){library(x,character.only=TRUE)})

  # Clear console:
  cat('\014');

  # Clear all plots:
  try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
  try(dev.off(),silent=TRUE)

  # Define paths for file handling:
  this.dir <- dirname(parent.frame(2)$ofile) # source file dir
  setwd(this.dir)                            # set wd to source file dir
  code.path   <- getwd()                     # define code path
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

  }

}

#===============================================================================
# FUNCTIONS:
#===============================================================================

source("ghFun.R")
source("seroFun.R")
source("tictocFun.R")

#===============================================================================
# MAIN PROGRAM:
#===============================================================================

tictocFun('tic')

set.seed(12345)

# Set simulation parameters:
VSL       <- 10
J.lo      <- 1
se.lo     <- 0.5 
se.hi     <- 5.0 
sig.eta   <- 0   # set in cases
sig.mu.lo <- 0.5
sig.mu.hi <- 0   # set in cases
K         <- 0   # Number of moderator variables for meta-regression

# Set Monte Carlo and bootstrap reps:
MC <- 1000  # 1000
BS <- 0     # 100 (to turn off use BS = 0)

eps <- 1e-15 # epsilon used later to avoid divide by zero

# Define cases for simulation:
#  I  J.hi   sig.eta  sig.mu.hi
cases <- matrix(c(
  20,    5,        1,         1,
  20,    5,        1,         3,
  20,    5,        3,         1,
  20,    5,        3,         3,
  20,   15,        1,         1,
  20,   15,        1,         3,
  20,   15,        3,         1,
  20,   15,        3,         3,
  60,    5,        1,         1,
  60,    5,        1,         3,
  60,    5,        3,         1,
  60,    5,        3,         3,
  60,   15,        1,         1,
  60,   15,        1,         3,
  60,   15,        3,         1,
  60,   15,        3,         3),16,4,byrow=TRUE)

for(iteration in 1:1){
  
  if(iteration==1){rho.lo <- 0.0; rho.hi <- 0.0; rho.hat <- 0.0}
  if(iteration==2){rho.lo <- 0.5; rho.hi <- 0.5; rho.hat <- 0.0}
  if(iteration==3){rho.lo <- 0.5; rho.hi <- 0.5; rho.hat <- 0.5}
  if(iteration==4){rho.lo <- 0.0; rho.hi <- 0.0; rho.hat <- 0.5}
  
  source("vslmeta-sim-iteration.R")
  
}

tictocFun('toc')

