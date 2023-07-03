#===============================================================================
# vslmeta-sim-iteration.R:
# This script will replicate the results of our application of the 2SRE meta-
# analysis estimator to the constructed datasets as reported in
# Newbold SC, Dockins C, Simon N, Maguire K, Sakib A. (2023)
#
# [called by vslmeta-sim.R]
#===============================================================================

if(iteration==1){
  
  # WRITE SOURCE FILE TO OUTPUT FILE:
  {
    script.name <- 'vslmeta-sim-iteration'
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
  
}

for (qq in 1:length(cases[,1])){
  
  # Number of groups:
  I <- cases[qq,1]
  
  # Max number of observations per group:
  J.hi <- cases[qq,2]
  
  # Group level non-sampling error's sd (stays same for a group):
  sig.eta <- cases[qq,3]
  
  # Max sig.mu:
  sig.mu.hi <- cases[qq,4]
  
  # Number of observations in each group:
  J <- matrix(sample(J.lo:J.hi,I,replace=TRUE),I,1)
  if(J.lo==J.hi){J <- matrix(J.lo,I,1)}
  
  # Observation level sampling error:
  se <- matrix(runif(I*J.hi),I,J.hi)*(se.hi-se.lo) + se.lo
  for(i in 1:I){
    if(J[i]<J.hi){se[i,(J[i]+1):J.hi] <- 0}
  }
  
  # Correlation among sampling errors in each group:
  rho <- matrix(runif(I)*(rho.hi-rho.lo)+rho.lo, I, 1) # to simulated data
  
  # Observation level non-sampling error's sd
  sig.mu <- matrix(runif(I)*(sig.mu.hi-sig.mu.lo)+sig.mu.lo, I, 1)
  
  outs <- ghFun(se,sig.eta,sig.mu,rho)
  g    <- outs[[1]]
  h    <- outs[[2]]
  
  N    <- sum(J)
  XX   <- matrix(0,MC,1+BS)
  XXX  <- matrix(0,MC,1)
  y.rk <- XX; sero.rk <- XXX; sebs.rk <- XXX # 2 stage RE known variances
  y.ru <- XX; sero.ru <- XXX; sebs.ru <- XXX # 2 stage RE unconstrained sig.mu
  y.rc <- XX; sero.rc <- XXX; sebs.rc <- XXX # 2 stage RE constrained sig.mu
  y.sm <- XX; sero.sm <- XXX; sebs.sm <- XXX # simple mean
  y.mm <- XX; sero.mm <- XXX; sebs.mm <- XXX # mean of means
  y.mf <- XX # metafor package
  y.rr <- XX # robumeta package CORR
  y.rh <- XX # robumeta package HIER
  y.ma <- XX # MAd package
  
  varyhat     <- matrix(0,MC,1)
  sig.muhat   <- matrix(0,MC,I)
  sig.etahat  <- matrix(0,MC,1)
  sig.muChat  <- matrix(0,MC,I)
  sig.etaChat <- matrix(0,MC,1)
  
  start.time <- proc.time()
  start.time <- start.time[3]
  for(mc in 1:MC){
    
    # Report progress to console:
    {
      if(floor(mc/10)==mc/10){
        now.time <- proc.time()
        now.time <- now.time[3]
        rate <- (now.time-start.time)/(mc-1)
        cat('\014')
        cat('Working on iteration',sprintf('%-.0f',iteration),'of 4.\n')
        cat('Working on case',sprintf('%-.0f',qq),'of',sprintf('%-.0f\n',length(cases[,1])))
        cat('  Working on MC rep',sprintf('%-.0f',mc),'of',sprintf('%-.0f\n',MC))
        if(mc>1){
          cat('  [Time to completion for this case =',sprintf('%-.1f',(MC-mc)*rate/60),'min]')
        }
      }
    }
    
    # Simulate data:
    {
      done <- 0
      while(done==0){
        y <- matrix(0,I,J.hi)
        for(i in 1:I){
          
          temp <- matrix(1,J[i],J[i])*rho[i]
          for(j in 1:J[i]){temp[j,j] <- 1}
          CF <- chol(temp)
          CV <- se[i,1:J[i]] %*% t(se[i,1:J[i]])
          
          eta.i <- rnorm(1)*sig.eta
          
          yi <- matrix(rnorm(J[i]),1,J[i]) %*% CF
          yi <- matrix(yi,J[i],1)
          for(j in 1:J[i]){
            y[i,j] <- yi[j]*se[i,j] + VSL + rnorm(1)*sig.mu[i] + eta.i
          }
        }
        #if(sum(y<0)==0){done <- 1}
        done <- 1
      }
      
      # Reformat data for other R meta-analysis packages:
      N <- sum(J)
      XX <- matrix(0,N,1)
      id <- XX
      yy <- XX
      vv <- XX
      n  <- 0
      for(i in 1:I){
        for(j in 1:J[i]){
          n     <- n + 1
          id[n] <- i
          yy[n] <- y[i,j]
          vv[n] <- se[i,j]^2
        }
      }
      # Convert reformated data to dataframe:
      DM.df <- as.data.frame(cbind(id, yy, vv))
    }
    
    # Calculate group means and variances:
    {
      m  <- matrix(0,I,1)
      v  <- matrix(0,I,1)
      va <- matrix(0,I,1)
      for(i in 1:I){
        m[i]  <- mean(y[i,1:J[i]])
        v[i]  <- var(y[i,1:J[i]])
        va[i] <- mean(se[i,1:J[i]]^2) # for use in metafor package
      }
      v[which(J==1)] <- sum(v[which(J>1)]*J[which(J>1)])/sum(J[which(J>1)])
    }
    
    # Calculate true g, h, and var:
    {
      outs <- ghFun(se,sig.eta,sig.mu,rho)
      g    <- outs[[1]]
      h    <- outs[[2]]
      
      varyhati <- matrix(0,I,1)
      for(i in 1:I){
        varyhati[i] <- sig.eta^2
        for(j in 1:J[i]){
          varyhati[i] <- varyhati[i] + g[i,j]^2*(sig.mu[i]^2+se[i,j]^2)
          for(k in 1:J[i]){
            if(j!=k){
              varyhati[i] <- varyhati[i] + rho[i] * g[i,j]*g[i,k]*se[i,j]*se[i,k]
            }
          }
        }
      }
      
      varyhat[mc] <- t(h^2) %*% varyhati
      
    }
    
    # UNCONSTRAINED SIG.MUHAT:
    
    # Estimate sig.mu for each group:
    {
      for(i in 1:I){
        
        if(J[i]>1){ # Estimate sig.mu only for groups with >1 observation
          
          CV <- se[i,1:J[i]] %*% t(se[i,1:J[i]])
          
          sig.mu2hat <- var(y[i,1:J[i]])-1/J[i]*(sum(diag(CV)-rho.hat/(J[i]^2)*
                                                       (sum(CV)-sum(diag(CV)))))
          
          # Only estimate sig.mu if estimated sig.mu^2 is positive:
          if(sig.mu2hat>0){sig.muhat[mc,i] <- sqrt(sig.mu2hat)}
          
        }
        sig.muhat[mc,i] <- sig.muhat[mc,i] * (J[i]/(J[i]-1))
      }
      
      # For all groups with one observation, impute sig.mu as sample size
      # weighted average of all other group's sig.muhat':
      sig.muhat[mc,which(J==1)] <- sum( sig.muhat[mc,which(J>1)] * J[which(J>1)]
                                        / sum(J[which(J>1)]) )
    }
    
    # Calculate ghat:
    {
      outs <- ghFun(se,0,sig.muhat[mc,],matrix(rho.hat,I,1))
      ghat <- outs[[1]]
    }
    
    # Estimate sig.eta:
    {
      # Calculate yhati and Vi for each group w/ more than one observation:
      yhati <- matrix(0,I,1)
      Vi    <- matrix(0,I,1)
      for(i in 1:I){
        if(J[i]>1){
          yhati[i] <- sum(y[i,1:J[i]]*ghat[i,1:J[i]])
          for(j in 1:J[i]){
            Vi[i] <- Vi[i] + ghat[i,j]^2*(sig.muhat[i]^2+se[i,j]^2)
            for(k in 1:J[i]){
              if(k!=j){
                Vi[i] <- Vi[i] + rho.hat*ghat[i,j]*ghat[i,k]*se[i,j]*se[i,k]
              }
            }
          }
        }
        if(J[i]==1){yhati[i] <- y[i,1]}
      }
      
      # Use yhati's and Vi's to calculate sig.etahat [Equation 49]:
      part1 <- 1/(I-1)*sum((yhati-mean(yhati))^2)
      part2 <- 1/sum(J>1)*sum(Vi[J>1])
      if(part1>part2){sig.etahat[mc] <- sqrt(part1-part2)}
    }
    
    # Calculate ghat and hhat:
    {
      outs <- ghFun(se,sig.etahat[mc],sig.muhat[mc,],matrix(rho.hat,I,1))
      ghat <- outs[[1]]
      hhat <- outs[[2]]
    }
    
    # CONSTRAINED SIG.MUHAT:
    
    # Estimate sig.muC (constrained)
    {
      sig.muChat[mc,] <- matrix(mean(sig.muhat[mc,which(J>1)]),I,1)
    }
    
    # Calculate gChat:
    {
      outs  <- ghFun(se,0,sig.muChat[mc,],matrix(rho.hat,I,1))
      gChat <- outs[[1]]
    }
    
    # Estimate sig.etaC:
    {
      # Calculate yhati and Vi for each group w/ more than one observation:
      yhati <- matrix(0,I,1)
      Vi    <- matrix(0,I,1)
      for(i in 1:I){
        if(J[i]>1){
          yhati[i] <- sum(y[i,1:J[i]]*gChat[i,1:J[i]])
          for(j in 1:J[i]){
            Vi[i] <- Vi[i] + gChat[i,j]^2*(sig.muChat[i]^2+se[i,j]^2)
            for(k in 1:J[i]){
              if(k!=j){
                Vi[i] <- Vi[i] + rho.hat*gChat[i,j]*gChat[i,k]*se[i,j]*se[i,k]
              }
            }
          }
        }
        if(J[i]==1){yhati[i] <- y[i,1]}
      }
      
      # Use yhati's and Vi's to calculate sig.etahat [Equation 49]:
      part1 <- 1/(I-1)*sum((yhati-mean(yhati))^2)
      part2 <- 1/sum(J>1)*sum(Vi[J>1])
      if(part1>part2){sig.etaChat[mc] <- sqrt(part1-part2)}
    }
    
    # Calculate gChat and hChat:
    {
      outs  <- ghFun(se,sig.etaChat[mc],sig.muChat[mc,],matrix(rho.hat,I,1))
      gChat <- outs[[1]]
      hChat <- outs[[2]]
    }
    
    # CALCULATE ALTERNATIVE ESTIMATES:
    if(TRUE){
      y.rk[mc,1]  <- sum(rowSums(g*y)*h)
      sero.rk[mc] <- seroFun(I,J,y,y.rk[mc,1],g,h)
      
      y.ru[mc,1]  <- sum(rowSums(ghat*y)*hhat)
      sero.ru[mc] <- seroFun(I,J,y,y.ru[mc,1],ghat,hhat)
      
      y.rc[mc,1]  <- sum(rowSums(gChat*y)*hChat)
      sero.rc[mc] <- seroFun(I,J,y,y.rc[mc,1],gChat,hChat)
      
      y.sm[mc,1]  <- mean(y[which(y!=0)])
      sero.sm[mc] <- seroFun(I,J,y,y.sm[mc,1],(g>0)/N,matrix(1,I,1))
      
      y.mm[mc,1]  <- mean((1/J)*rowSums(y))
      sero.mm[mc] <- seroFun(I,J,y,y.mm[mc,1],
                             matrix(1/J,I,J.hi)*(g>0),matrix(1/I,I,1))
      
      # metafor
      # run.1 <- rma(yi      = m,    # effect size
      #              vi      = va,   # variance
      #              method  = 'ML', # options: [DL], ML
      #              control = list(stepadj=0.1, maxit=50000000))
      
      # Create y vector for metafor:
      {
        y.rma <- {}
        for(i in 1:I){
          y.rma <- rbind(y.rma,matrix(y[i,1:J[i]],J[i],1))
        }
      }
      
      # Create variance-covariance matrix for metafor:
      {
        V1 <- matrix(0,J[1],J[1])
        for(j in 1:J[1]){
          for(k in 1:J[1]){
            if(j==k){V1[j,k] <- se[1,j]^2}else{V1[j,k] <- rho.hat*se[1,j]*se[1,k]}
          }
        }
        V.rma <- V1
        
        for(i in 2:I){
          
          V.i <- matrix(0,J[i],J[i])
          for(j in 1:J[i]){
            for(k in 1:J[i]){
              if(j==k){V.i[j,k] <- se[i,j]^2}else{V.i[j,k] <- rho.hat*se[i,j]*se[i,k]}
            }
          }
          
          temp <- matrix(0,nrow(V.rma)+nrow(V.i),ncol(V.rma)+ncol(V.i))
          temp[1:nrow(V.rma),1:ncol(V.rma)] <- V.rma
          temp[(nrow(V.rma)+1):(nrow(V.rma)+nrow(V.i)),(ncol(V.rma)+1):(ncol(V.rma)+ncol(V.i))] <- V.i
          V.rma <- temp
          
        }
      }
      
      outs  <- rma.mv(yi      = y.rma,  # effect size
                      V       = V.rma,  # variance
                      method  = 'REML', # options: [REML], ML
                      control = list(stepadj=0.01, maxit=50000000))
      
      y.mf[mc,1] <- coef.rma(outs)
      
      # robumeta CORR
      outs <- robu(formula      = yy ~ 1, # effect size
                   var.eff.size = vv,     # variance
                   data         = DM.df,  # dataset
                   studynum     = id,     # study-level ID
                   rho          = rho.hat,# (to match our 2SRE estimator)
                   modelweights = 'CORR', # options: [CORR], HIER
                   small        = TRUE)  # [TRUE], FALSE
      
      y.rr[mc,1] <- outs$b.r
      
      # robumeta HIER
      outs <- robu(formula      = yy ~ 1, # effect size
                   var.eff.size = vv,     # variance
                   data         = DM.df,  # dataset
                   studynum     = id,     # study-level ID
                   modelweights = 'HIER', # options: [CORR], HIER
                   small        = TRUE)  # [TRUE], FALSE
      
      y.rh[mc,1] <- outs$b.r
      
      # MAd
      outs <- agg(id     = id,
                  es     = yy,
                  var    = vv,
                  cor    = rho.hat,
                  method = "BHHR",
                  data   = DM.df)
      
      outs <- mareg(outs$es ~ 1,
                    var     = outs$var,
                    method  = 'REML',
                    control = list(stepadj=0.01, maxit=50000000),
                    data    = outs)
      
      y.ma[mc,1] <- coef(outs)
    }
    
    # Bootstrap:
    if(BS>0){
      
      for(b in 1:BS){
        
        # Re-sample data:
        {
          done <- 0
          while(done==0){
            draw       <- sample(1:I,I,replace=TRUE)
            yb         <- y[draw,]
            seb        <- se[draw,]
            sig.mub    <- sig.mu[draw]
            sig.muhatb <- sig.muhat[mc,draw]
            rhob       <- rho[draw]
            Jb         <- J[draw]
            if(sum(Jb>=2)>=2){done <- 1}
          }
          
          # Reformat data for other R meta-analysis packages:
          N   <- sum(Jb)
          XX  <- matrix(0,N,1)
          idb <- XX
          yyb <- XX
          vvb <- XX
          n   <- 0
          for(i in 1:I){
            for(j in 1:Jb[i]){
              n      <- n + 1
              idb[n] <- i
              yyb[n] <- yb[i,j]
              vvb[n] <- seb[i,j]^2
            }
          }
          # Convert reformated data to dataframe:
          DM.dfb <- as.data.frame(cbind(idb, yyb, vvb))
        }
        
        # Calculate group means and variances:
        {
          mb  <- matrix(0,I,1)
          vb  <- matrix(0,I,1)
          vab <- matrix(0,I,1)
          for(i in 1:I){
            mb[i]  <- mean(yb[i,1:Jb[i]])
            vb[i]  <- var(yb[i,1:Jb[i]])
            vab[i] <- mean(seb[i,1:Jb[i]]^2) # for use in metafor package
          }
          vb[which(Jb==1)] <- sum(vb[which(Jb>1)]*Jb[which(Jb>1)])/
            sum(Jb[which(Jb>1)])
        }
        
        # Calculate g and h:
        {
          outs <- ghFun(seb,sig.eta,sig.mub,rhob)
          gb   <- outs[[1]]
          hb   <- outs[[2]]
        }
        
        # UNCONSTRAINED SIG.MUHAT:
        
        # Estimate sig.mu for each group:
        {
          for(i in 1:I){
            if(Jb[i]>1){ # Estimate sig.mu only for groups with >1 observation
              
              CVb <- seb[i,1:Jb[i]] %*% t(seb[i,1:Jb[i]])
              
              sig.mu2hatb <- var(yb[i,1:Jb[i]])-1/Jb[i]*
                (sum(diag(CVb)-rhob[i]/(Jb[i]^2)*
                       (sum(CVb)-sum(diag(CVb)))))
              
              # Only estimate sig.mu if estimated sig.mu^2 is positive:
              if(sig.mu2hatb>0){sig.muhatb[i] <- sqrt(sig.mu2hatb)}
              
            }
            sig.muhatb[i] <- sig.muhatb[i] * (Jb[i]/(Jb[i]-1))
          }
          
          # For all groups with one observation, impute sig.mu as sample size
          # weighted average of all other group's sig.muhat':
          sig.muhatb[which(Jb==1)] <- sum( sig.muhatb[which(Jb>1)] *
                                             Jb[which(Jb>1)] / sum(Jb[which(Jb>1)]) )
        }
        
        # Calculate ghat:
        {
          outs  <- ghFun(seb,0,sig.muhatb,rhob)
          ghatb <- outs[[1]]
        }
        
        # Estimate sig.eta:
        {
          sig.etahatb <- 0
          # Calculate yhati and Vi for each group w/ more than one observation:
          yhatib <- matrix(0,I,1)
          Vib    <- matrix(0,I,1)
          for(i in 1:I){
            if(Jb[i]>1){
              yhatib[i] <- sum(yb[i,1:Jb[i]]*ghatb[i,1:Jb[i]])
              for(j in 1:Jb[i]){
                Vib[i] <- Vib[i] + ghatb[i,j]^2*(sig.muhatb[i]^2+seb[i,j]^2)
                for(k in 1:Jb[i]){
                  if(k!=j){
                    Vib[i] <- Vib[i] +
                      rhob[i]*ghatb[i,j]*ghatb[i,k]*seb[i,j]*seb[i,k]
                  }
                }
              }
            }
            if(Jb[i]==1){yhatib[i] <- yb[i,1]}
          }
          
          # Use yhati's and Vi's to calculate sig.etahat [Equation 49]:
          part1 <- 1/(I-1)*sum((yhatib-mean(yhatib))^2)
          part2 <- 1/sum(Jb>1)*sum(Vib[Jb>1])
          if(part1>part2){sig.etahatb <- sqrt(part1-part2)}
        }
        
        # Calculate ghat and hhat:
        {
          outs  <- ghFun(seb,sig.etahatb,sig.muhatb,rhob)
          ghatb <- outs[[1]]
          hhatb <- outs[[2]]
        }
        
        # CONSTRAINED SIG.MUHAT:
        
        # Estimate sig.muC (constrained)
        {
          sig.muChatb <- matrix(mean(sig.muhatb[which(Jb>1)]),I,1)
        }
        
        # Calculate gChat:
        {
          outs   <- ghFun(seb,0,sig.muChatb,rhob)
          gChatb <- outs[[1]]
        }
        
        # Estimate sig.etaC:
        {
          # Calculate yhati and Vi for each group w/ more than one observation:
          yhatib <- matrix(0,I,1)
          Vib    <- matrix(0,I,1)
          for(i in 1:I){
            if(Jb[i]>1){
              yhatib[i] <- sum(yb[i,1:Jb[i]]*gChatb[i,1:Jb[i]])
              for(j in 1:Jb[i]){
                Vib[i] <- Vib[i] + gChatb[i,j]^2*(sig.muChatb[i]^2+seb[i,j]^2)
                for(k in 1:Jb[i]){
                  if(k!=j){
                    Vib[i] <- Vib[i] +
                      rhob[i]*gChatb[i,j]*gChatb[i,k]*seb[i,j]*seb[i,k]
                  }
                }
              }
            }
            if(Jb[i]==1){yhatib[i] <- yb[i,1]}
          }
          
          # Use yhati's and Vi's to calculate sig.etahat [Equation 49]:
          part1 <- 1/(I-1)*sum((yhati-mean(yhati))^2)
          part2 <- 1/sum(J>1)*sum(Vi[J>1])
          if(part1>part2){sig.etaChatb <- sqrt(part1-part2)}
          else{sig.etaChatb <- 0}
        }
        
        # Calculate gChat and hChat:
        {
          outs   <- ghFun(seb,sig.etaChatb,sig.muChatb,rhob)
          gChatb <- outs[[1]]
          hChatb <- outs[[2]]
        }
        
        # CALCULATE ALTERNATIVE ESTIMATES:
        {
          y.rk[mc,1+b] <- sum(rowSums(gb*yb)*hb)
          y.ru[mc,1+b] <- sum(rowSums(ghatb*yb)*hhatb)
          y.rc[mc,1+b] <- sum(rowSums(gChatb*yb)*hChatb)
          y.sm[mc,1+b] <- mean(yb[which(yb!=0)])
          y.mm[mc,1+b] <- mean((1/Jb)*rowSums(yb))
        }
        
      }
      
      sebs.rk[mc] <- sqrt(var(y.rk[mc,2:(b+1)]))*I/(I-1)
      sebs.ru[mc] <- sqrt(var(y.ru[mc,2:(b+1)]))*I/(I-1)
      sebs.rc[mc] <- sqrt(var(y.rc[mc,2:(b+1)]))*I/(I-1)
      sebs.sm[mc] <- sqrt(var(y.sm[mc,2:(b+1)]))*I/(I-1)
      sebs.mm[mc] <- sqrt(var(y.mm[mc,2:(b+1)]))*I/(I-1)
      
    }
    
  }
  
  # Write results to output file:
  {
    if(qq==1){
      cat(sprintf('\nrho = %-.1f; rho.hat = %-.1f\n',rho.lo,rho.hat),file=out.file.name,append=TRUE)
    }
    # Monte Carlo estimates of standard errors:
    cat(sprintf('%4.0f & %1.0f,%2.0f & %4.1f & %3.1f,%3.1f & %5.3f & %5.3f & %5.3f & %5.3f & %5.3f & %5.3f & %5.3f & %5.3f & %5.3f & %5.3f\\\\ \n',
                I,
                J.lo,J.hi,
                sig.eta,
                sig.mu.lo,sig.mu.hi,
                sd(y.sm[,1]),   # simple mean
                sd(y.mm[,1]),   # mean of means
                sd(y.mf[,1]),   # metafor
                sd(y.rr[,1]),   # robumeta CORR
                sd(y.rh[,1]),   # robumeta HIER
                sd(y.ma[,1]),   # MAd
                sd(y.rk[,1]),   # 2SRE true
                sd(y.ru[,1]),   # 2SRE free
                sd(y.rc[,1]),   # 2SRE equal
                mean(sero.rc)), # robust s.e. for 2SRE constrained estimator
        file=out.file.name,append=TRUE)
  }
  
}

