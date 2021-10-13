
###################################################################################################
###################################################################################################
# Miscellaneous functions for probabilstic models linking multiple founder variants to the associated
# increase in viral load.
###################################################################################################
###################################################################################################

###################################################################################################
# Execute in parallel
RunParallel <- function(func, v1, ...){
  options(warn = 1)
  require(parallel)
  
  # Set up cluster (fork)
  cl <- detectCores() %>% `-` (2) 
 
  start <- Sys.time()
  print(start)
    
  para <- mclapply(v1,
                   func, 
                   ...,
                   mc.cores = cl,
                   mc.set.seed = FALSE) #child process has the same initial random number generator (RNG) state as the current R session
    
    end <- Sys.time()
    elapsed <- end-start
    print(end)
    print(elapsed)
 
  return(para)
}


###################################################################################################