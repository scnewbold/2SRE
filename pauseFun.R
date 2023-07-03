pauseFun <- function(sec){

  if(missing(sec)){
    readline(prompt='Paused. Press [enter] to continue.\n')
  }else{
    Sys.sleep(sec)
  }

}
