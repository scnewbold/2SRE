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
      'MAd',
      'MASS',
      'metafor',
      'NlcOptim',
      'readxl',
      'robumeta',
      'stringr',
      'tikzDevice',
      'tinytex')
  
  new.packages <- list.of.packages[!(list.of.packages %in%
                                       installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages,function(x){library(x,character.only=TRUE)})
  
  # Define paths for file handling:
  this.dir <- dirname(parent.frame(2)$ofile) # source file dir
  setwd(this.dir)                            # set wd to source file dir
  code.path <- getwd()                       # define code path
  output.path <- getwd()                     # define output path
  
}

# The quantities plotted in the figures below were computed manually as follows:
# The center, bottom, and top horizontal lines of each box are the average, 
# minimum, and maximum relative precision measures for each estimator across the
# 16 experimental design settings tested. The relative precision was computed as
# the difference between the standard error of the estimator and the minimum 
# theoretical standard error divided by the minimum theoretical standard error,
# i.e., (se.hat−se.min)/se.min. See Tables S2.1-S2.4 in the supplemental
# information file S4 for the raw results used to compute these summary
# statistics.

# Figure 1a (rho=0,rho.hat=0)
{
  center <- c(0.155 , 0.207 , 0.251 , 0.172 , 0.066 , 0.092 , 0.188 , 0.022)
  plus   <- c(0.153 , 0.243 , 0.232 , 0.167 , 0.060 , 0.107 , 0.515 , 0.045)
  minus  <- c(0.079 , 0.175 , 0.144 , 0.130 , 0.042 , 0.074 , 0.174 , 0.021)
  
  hw <- 0.4 # half-width of bars
  
  fig.name <- paste(output.path,'/fig1a',sep='')
  tikz(paste(fig.name,'.tex',sep=''),
       width      = 3.2,
       height     = 3.2,
       pointsize  = 12,
       standAlone = TRUE)
  
  par(mfrow=c(1,1))
  par(mar=c(3,1.5,1.5,0.5))
  
  plot(1:8,center,xlim=c(0.5,8.5),cex=0.3,ylim=c(0,1.2),xlab='',ylab='',axes=FALSE)
  box(lwd=1)
  for(k in 1:8){
    lines(c(k-hw,k+hw),c(center[k],center[k]))
    lines(c(k-hw,k+hw),c(center[k]+plus[k],center[k]+plus[k]))
    lines(c(k-hw,k+hw),c(center[k]-minus[k],center[k]-minus[k]))
    lines(c(k-hw,k-hw),c(center[k]-minus[k],center[k]+plus[k]))
    lines(c(k+hw,k+hw),c(center[k]-minus[k],center[k]+plus[k]))
    
  }
  mtext(c('simple','group','metafor','robum.','robum.','MAd','2SRE','2SRE'), 
        side = 1, at = seq(1,8,1), las = 1, line=-.3, cex=.5) # x-axis labels
  mtext(c('','','','CORR','HIER','','free','equal'), 
        side = 1, at = seq(1,8,1), las = 1, line=.2, cex=.5) # x-axis labels
  mtext(seq(0,1.2,.2), side = 2, at = seq(0,1.2,.2) , las = 1 , cex = .5, line = .3 )
  mtext('true $\\rho=0$, assumed $\\rho=0$', side = 3, at = 4.5, las = 1 , cex = .6, line = .3 )
  
  dev.off()
  tinytex::latexmk(paste(fig.name,".tex",sep=""))
  file.remove(paste(fig.name,".tex",sep=""))
}

# Figure 1b (rho=0.5,rho.hat=0)
{
  center <- c(0.361 , 0.298 , 0.220 , 0.241 , 0.170 , 0.229 , 0.086 , 0.068)
  plus   <- c(0.422 , 0.352 , 0.264 , 0.365 , 0.136 , 0.305 , 0.116 , 0.096)
  minus  <- c(0.219 , 0.223 , 0.120 , 0.171 , 0.068 , 0.179 , 0.067 , 0.058)
  
  hw <- 0.4 # half-width of bars
  
  fig.name <- paste(output.path,'/fig1b.tex',sep='')
  tikz(paste(fig.name,'.tex',sep=''),
       width      = 3.2,
       height     = 3.2,
       pointsize  = 12,
       standAlone = TRUE)
  
  par(mfrow=c(1,1))
  par(mar=c(3,1.5,1.5,0.5))
  
  plot(1:8,center,xlim=c(0.5,8.5),cex=0.3,ylim=c(0,1.2),xlab='',ylab='',axes=FALSE)
  box(lwd=1)
  for(k in 1:8){
    lines(c(k-hw,k+hw),c(center[k],center[k]))
    lines(c(k-hw,k+hw),c(center[k]+plus[k],center[k]+plus[k]))
    lines(c(k-hw,k+hw),c(center[k]-minus[k],center[k]-minus[k]))
    lines(c(k-hw,k-hw),c(center[k]-minus[k],center[k]+plus[k]))
    lines(c(k+hw,k+hw),c(center[k]-minus[k],center[k]+plus[k]))
    
  }
  mtext(c('simple','group','metafor','robum.','robum.','MAd','2SRE','2SRE'), 
        side = 1, at = seq(1,8,1), las = 1, line=-.3, cex=.5) # x-axis labels
  mtext(c('','','','CORR','HIER','','free','equal'), 
        side = 1, at = seq(1,8,1), las = 1, line=0.2, cex=.5) # x-axis labels
  mtext(seq(0,1.2,.2), side = 2, at = seq(0,1.2,.2) , las = 1 , cex = .5, line = .3 )
  mtext('true $\\rho=0.5$, assumed $\\rho=0$', side = 3, at = 4.5, las = 1 , cex = .6, line = .3 )
  
  dev.off()
  tinytex::latexmk(paste(fig.name,".tex",sep=""))
  file.remove(paste(fig.name,".tex",sep=""))
}

# Figure 1c (rho=0.5,rho.hat=0.5)
{
  center <- c(0.346 , 0.279 , 0.328 , 0.230 , 0.164 , 0.213 , 0.108 , 0.024)
  plus   <- c(0.343 , 0.352 , 0.191 , 0.364 , 0.109 , 0.375 , 0.310 , 0.050)
  minus  <- c(0.213 , 0.202 , 0.238 , 0.169 , 0.094 , 0.163 , 0.101 , 0.025)
  
  hw <- 0.4 # half-width of bars
  
  fig.name <- paste(output.path,'/fig1c.tex',sep='')
  tikz(paste(fig.name,'.tex',sep=''),
       width      = 2.6,
       height     = 2.6,
       pointsize  = 12,
       standAlone = TRUE)
  
  par(mfrow=c(1,1))
  par(mar=c(3,1.5,1.5,0.5))
  
  plot(1:8,center,xlim=c(0.5,8.5),cex=0.3,ylim=c(0,1.2),xlab='',ylab='',axes=FALSE)
  box(lwd=1)
  for(k in 1:8){
    lines(c(k-hw,k+hw),c(center[k],center[k]))
    lines(c(k-hw,k+hw),c(center[k]+plus[k],center[k]+plus[k]))
    lines(c(k-hw,k+hw),c(center[k]-minus[k],center[k]-minus[k]))
    lines(c(k-hw,k-hw),c(center[k]-minus[k],center[k]+plus[k]))
    lines(c(k+hw,k+hw),c(center[k]-minus[k],center[k]+plus[k]))
    
  }
  mtext(c('simple','group','metafor','robum.','robum.','MAd','2SRE','2SRE'), 
        side = 1, at = seq(1,8,1), las = 1, line=-.3, cex=.5) # x-axis labels
  mtext(c('','','','CORR','HIER','','free','equal'), 
        side = 1, at = seq(1,8,1), las = 1, line=0.2, cex=.5) # x-axis labels
  mtext(seq(0,1.2,.2), side = 2, at = seq(0,1.2,.2) , las = 1 , cex = .5, line = .3 )
  mtext('true $\\rho=0.5$, assumed $\\rho=0.5$', side = 3, at = 4.5, las = 1 , cex = .6, line = .3 )
  
  dev.off()
  tinytex::latexmk(paste(fig.name,".tex",sep=""))
  file.remove(paste(fig.name,".tex",sep=""))
}

# Figure 1d (rho=0,rho.hat=0.5)
{
  center <- c(0.168 , 0.205 , 0.494 , 0.182 , 0.068 , 0.149 , 0.367 , 0.035)
  plus   <- c(0.126 , 0.275 , 0.370 , 0.157 , 0.067 , 0.185 , 0.361 , 0.072)
  minus  <- c(0.119 , 0.181 , 0.192 , 0.161 , 0.039 , 0.130 , 0.336 , 0.031)
  
  hw <- 0.4 # half-width of bars
  
  fig.name <- paste(output.path,'/fig1d.tex',sep='')
  tikz(paste(fig.name,'.tex',sep=''),
       width      = 2.6,
       height     = 2.6,
       pointsize  = 12,
       standAlone = TRUE)
  
  par(mfrow=c(1,1))
  par(mar=c(3,1.5,1.5,0.5))
  
  plot(1:8,center,xlim=c(0.5,8.5),cex=0.3,ylim=c(0,1.2),xlab='',ylab='',axes=FALSE)
  box(lwd=1)
  for(k in 1:8){
    lines(c(k-hw,k+hw),c(center[k],center[k]))
    lines(c(k-hw,k+hw),c(center[k]+plus[k],center[k]+plus[k]))
    lines(c(k-hw,k+hw),c(center[k]-minus[k],center[k]-minus[k]))
    lines(c(k-hw,k-hw),c(center[k]-minus[k],center[k]+plus[k]))
    lines(c(k+hw,k+hw),c(center[k]-minus[k],center[k]+plus[k]))
    
  }
  mtext(c('simple','group','metafor','robum.','robum.','MAd','2SRE','2SRE'), 
        side = 1, at = seq(1,8,1), las = 1, line=-.3, cex=.5) # x-axis labels
  mtext(c('','','','CORR','HIER','','free','equal'), 
        side = 1, at = seq(1,8,1), las = 1, line=0.2, cex=.5) # x-axis labels
  mtext(seq(0,1.2,.2), side = 2, at = seq(0,1.2,.2) , las = 1 , cex = .5, line = .3 )
  mtext('true $\\rho=0$, assumed $\\rho=0.5$', side = 3, at = 4.5, las = 1 , cex = .6, line = .3 )
  
  dev.off()
  tinytex::latexmk(paste(fig.name,".tex",sep=""))
  file.remove(paste(fig.name,".tex",sep=""))
}

# COMBINED FIGURE:
{
  # Figure 1a (rho=0,rho.hat=0)
  center <- c(0.155 , 0.207 , 0.251 , 0.172 , 0.066 , 0.092 , 0.188 , 0.022)
  plus   <- c(0.153 , 0.243 , 0.232 , 0.167 , 0.060 , 0.107 , 0.515 , 0.045)
  minus  <- c(0.079 , 0.175 , 0.144 , 0.130 , 0.042 , 0.074 , 0.174 , 0.021)
  
  hw <- 0.4 # half-width of bars
  
  fig.name <- paste(output.path,'/fig1.tex',sep='')
  tikz(paste(fig.name,'.tex',sep=''),
       width      = 5,
       height     = 5,
       pointsize  = 12,
       standAlone = TRUE)
  
  par(mfrow=c(2,2))
  par(mar=c(3,1.5,1.5,0.5))
  
  plot(1:8,center,xlim=c(0.5,8.5),cex=0.3,ylim=c(0,1.2),xlab='',ylab='',axes=FALSE)
  box(lwd=1)
  for(k in 1:8){
    lines(c(k-hw,k+hw),c(center[k],center[k]))
    lines(c(k-hw,k+hw),c(center[k]+plus[k],center[k]+plus[k]))
    lines(c(k-hw,k+hw),c(center[k]-minus[k],center[k]-minus[k]))
    lines(c(k-hw,k-hw),c(center[k]-minus[k],center[k]+plus[k]))
    lines(c(k+hw,k+hw),c(center[k]-minus[k],center[k]+plus[k]))
    
  }
  mtext(c('simple','group','metafor','robum.','robum.','MAd','2SRE','2SRE'), 
        side = 1, at = seq(1,8,1), las = 1, line=-.3, cex=.5) # x-axis labels
  mtext(c('','','','CORR','HIER','','free','equal'), 
        side = 1, at = seq(1,8,1), las = 1, line=.2, cex=.5) # x-axis labels
  mtext(seq(0,1.2,.2), side = 2, at = seq(0,1.2,.2) , las = 1 , cex = .5, line = .3 )
  mtext('true $\\rho=0$, assumed $\\rho=0$', side = 3, at = 4.5, las = 1 , cex = .6, line = .3 )
  
  # Figure 1b (rho=0.5,rho.hat=0)
  center <- c(0.361 , 0.298 , 0.220 , 0.241 , 0.170 , 0.229 , 0.086 , 0.068)
  plus   <- c(0.422 , 0.352 , 0.264 , 0.365 , 0.136 , 0.305 , 0.116 , 0.096)
  minus  <- c(0.219 , 0.223 , 0.120 , 0.171 , 0.068 , 0.179 , 0.067 , 0.058)
  
  plot(1:8,center,xlim=c(0.5,8.5),cex=0.3,ylim=c(0,1.2),xlab='',ylab='',axes=FALSE)
  box(lwd=1)
  for(k in 1:8){
    lines(c(k-hw,k+hw),c(center[k],center[k]))
    lines(c(k-hw,k+hw),c(center[k]+plus[k],center[k]+plus[k]))
    lines(c(k-hw,k+hw),c(center[k]-minus[k],center[k]-minus[k]))
    lines(c(k-hw,k-hw),c(center[k]-minus[k],center[k]+plus[k]))
    lines(c(k+hw,k+hw),c(center[k]-minus[k],center[k]+plus[k]))
    
  }
  mtext(c('simple','group','metafor','robum.','robum.','MAd','2SRE','2SRE'), 
        side = 1, at = seq(1,8,1), las = 1, line=-.3, cex=.5) # x-axis labels
  mtext(c('','','','CORR','HIER','','free','equal'), 
        side = 1, at = seq(1,8,1), las = 1, line=0.2, cex=.5) # x-axis labels
  mtext(seq(0,1.2,.2), side = 2, at = seq(0,1.2,.2) , las = 1 , cex = .5, line = .3 )
  mtext('true $\\rho=0.5$, assumed $\\rho=0$', side = 3, at = 4.5, las = 1 , cex = .6, line = .3 )

  # Figure 1c (rho=0.5,rho.hat=0.5)
  center <- c(0.346 , 0.279 , 0.328 , 0.230 , 0.164 , 0.213 , 0.108 , 0.024)
  plus   <- c(0.343 , 0.352 , 0.191 , 0.364 , 0.109 , 0.375 , 0.310 , 0.050)
  minus  <- c(0.213 , 0.202 , 0.238 , 0.169 , 0.094 , 0.163 , 0.101 , 0.025)
  
  plot(1:8,center,xlim=c(0.5,8.5),cex=0.3,ylim=c(0,1.2),xlab='',ylab='',axes=FALSE)
  box(lwd=1)
  for(k in 1:8){
    lines(c(k-hw,k+hw),c(center[k],center[k]))
    lines(c(k-hw,k+hw),c(center[k]+plus[k],center[k]+plus[k]))
    lines(c(k-hw,k+hw),c(center[k]-minus[k],center[k]-minus[k]))
    lines(c(k-hw,k-hw),c(center[k]-minus[k],center[k]+plus[k]))
    lines(c(k+hw,k+hw),c(center[k]-minus[k],center[k]+plus[k]))
    
  }
  mtext(c('simple','group','metafor','robum.','robum.','MAd','2SRE','2SRE'), 
        side = 1, at = seq(1,8,1), las = 1, line=-.3, cex=.5) # x-axis labels
  mtext(c('','','','CORR','HIER','','free','equal'), 
        side = 1, at = seq(1,8,1), las = 1, line=0.2, cex=.5) # x-axis labels
  mtext(seq(0,1.2,.2), side = 2, at = seq(0,1.2,.2) , las = 1 , cex = .5, line = .3 )
  mtext('true $\\rho=0.5$, assumed $\\rho=0.5$', side = 3, at = 4.5, las = 1 , cex = .6, line = .3 )

  # Figure 1d (rho=0,rho.hat=0.5)
  center <- c(0.168 , 0.205 , 0.494 , 0.182 , 0.068 , 0.149 , 0.367 , 0.035)
  plus   <- c(0.126 , 0.275 , 0.370 , 0.157 , 0.067 , 0.185 , 0.361 , 0.072)
  minus  <- c(0.119 , 0.181 , 0.192 , 0.161 , 0.039 , 0.130 , 0.336 , 0.031)

  plot(1:8,center,xlim=c(0.5,8.5),cex=0.3,ylim=c(0,1.2),xlab='',ylab='',axes=FALSE)
  box(lwd=1)
  for(k in 1:8){
    lines(c(k-hw,k+hw),c(center[k],center[k]))
    lines(c(k-hw,k+hw),c(center[k]+plus[k],center[k]+plus[k]))
    lines(c(k-hw,k+hw),c(center[k]-minus[k],center[k]-minus[k]))
    lines(c(k-hw,k-hw),c(center[k]-minus[k],center[k]+plus[k]))
    lines(c(k+hw,k+hw),c(center[k]-minus[k],center[k]+plus[k]))
    
  }
  mtext(c('simple','group','metafor','robum.','robum.','MAd','2SRE','2SRE'), 
        side = 1, at = seq(1,8,1), las = 1, line=-.3, cex=.5) # x-axis labels
  mtext(c('','','','CORR','HIER','','free','equal'), 
        side = 1, at = seq(1,8,1), las = 1, line=0.2, cex=.5) # x-axis labels
  mtext(seq(0,1.2,.2), side = 2, at = seq(0,1.2,.2) , las = 1 , cex = .5, line = .3 )
  mtext('true $\\rho=0$, assumed $\\rho=0.5$', side = 3, at = 4.5, las = 1 , cex = .6, line = .3 )
  
  dev.off()
  tinytex::latexmk(paste(fig.name,".tex",sep=""))
  file.remove(paste(fig.name,".tex",sep=""))
}