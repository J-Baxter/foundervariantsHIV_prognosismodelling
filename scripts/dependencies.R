
###################################################################################################
###################################################################################################
# Dependencies for probabilistic models linking multiple founder variants to the associated
# increase in viral load.
###################################################################################################
############################# Check necessary packages are available and install; import if not.
pkgs <- c('tidyverse', 'parallel', 'ggpmisc', 'ggpubr', 'ggprism', 'cowplot', 'scales',
          'sysfonts', 'showtext', 'showtextdb', 'ggExtra', 'ggdag', 'brms', 'magrittr',
          'RColorBrewer') 

for (pkg in pkgs){
  
  if (!(pkg %in% installed.packages())){
    
    install.packages(pkg)
    
    library(pkg, character.only = T)
    
  }else if (pkg %in% installed.packages()){
    
    library(pkg, character.only = T)
  }
}


###################################################################################################
# Create directories for results and figures
ddmonthyy <- format(Sys.Date(), '%d%b%y')
check_dirs <- paste(c('./results', './figures'), ddmonthyy, sep = '/')
dirs <- list.dirs()

for (check_dir in check_dirs){
  
  if (!(check_dir %in% dirs)){
    dir.create(check_dir, recursive = T)
  }
  
}

results_dir <- check_dirs[[1]]
figs_dir <- check_dirs[[2]]

###################################################################################################

# Execute in parallel
RunParallel <- function(func, v1, ...){
  options(warn = 1)
  require(parallel)
 
  # Set up cluster (fork)
  cl <- detectCores() %>% `-` (3) 
 
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