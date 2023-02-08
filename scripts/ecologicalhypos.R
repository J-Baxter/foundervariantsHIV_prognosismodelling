## Functions outlining high-level ecological interactions between variants


exclusion <- function(transmitter, variants){
  
  recip <- max(variants) 
  
  return(recip)
  
}


additive <- function(transmitter, variants, scale = 0.7 ){
  
  delta <- sapply(variants, function(x) x - transmitter) 
  
  recip <- sum(delta)*scale + transmitter
  
  return(recip)
  
}


interaction <- function(transmitter, variants, scale = 0.7){
  
  delta <- sapply(variants, function(x) x - transmitter) 
  print(delta)
  
  recip <- prod(delta)**(1/length(delta))+ transmitter
  
  
  return(recip)
  
}

interaction <- function(transmitter, variants, scale = 0.7){
  
  
  recip <- prod(variants)**(1/length(variants))
  
  
  return(recip)
  
}
3.99929, c(3.402199,3.851453,4.062112)