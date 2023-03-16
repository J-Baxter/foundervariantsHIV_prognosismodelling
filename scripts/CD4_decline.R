SimCD4Decline <- function(delta_cd4, init_cd4, theshold = 200){
  z <- init_cd4
  i <-  1
  
  while(z[i]>threshold){
    
    z[i+1] = z[i] - delta_cd4 
    i = i + 1 
    
  }
  
  return(z)
  
}


          
          
          