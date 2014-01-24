finneyGM      <- function(df,arg){
  # Initialize gm
  gm        <- numeric(length(arg));
  # Evaluate all args 
  for (i in 1:length(gm)){
    if (abs(arg[i])>50) {
      gm[i] <- NA
      cat('NA assigned to large arg[',i,']. Aborting...')
      # break
    }
    else {
      gm[i]   <- 1;
      bt      <- arg[i]*df^2/(2*(df+1));
      term    <- 1;
      # Proceed with the computation
      for (p in 1:1000){
        term  <- term * bt/((df/2+p-1)*p);
        gm[i] <- gm[i] + term;
        if (p>1 & abs(term)<1e-7) break
      }
    }
    # cat(format(i,digits=5),format(p,digits=5),format(gm[i],digits=5),"\n"); 
  }
  return(gm)
}
