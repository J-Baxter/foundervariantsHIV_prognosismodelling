SimCD4Decline <- function(delta_cd4, init_cd4, threshold = 200){
  z <- init_cd4
  i <-  1
  
  while(z[i]>threshold){
    z[i+1] = z[i] - abs(delta_cd4)*365.25
    i = i + 1 
    
  }
  
  return(i-1)
  
}


          
          
          