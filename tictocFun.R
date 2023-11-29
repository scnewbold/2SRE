tictocFun <- function(tic.or.toc){
  
  if(tic.or.toc=='tic'){
    tic.time<<-proc.time()[3]
  }
  
  if(tic.or.toc=='toc'){
    toc.time <- proc.time()[3]
    elapsed <- toc.time - tic.time
    return(sprintf("elapsed time = %0.2f sec = %0.2f min = %0.2f hr",elapsed,elapsed/60,elapsed/3600))
  }
  
}
