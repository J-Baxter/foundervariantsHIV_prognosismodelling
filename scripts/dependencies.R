
###################################################################################################
###################################################################################################
# Dependencies for probabilistic models linking multiple founder variants to the associated
# increase in viral load.
###################################################################################################
############################# Check necessary packages are available and install; import if not.
pkgs <- c('tidyverse', 'parallel', 'ggpmisc', 'ggpubr', 'ggprism', 'cowplot', 'scales',
          'sysfonts', 'showtext', 'showtextdb', 'ggExtra', 'ggdag', 'brms', 'magrittr',
          'RColorBrewer', 'tmvtnorm', 'stringi', 'matrixcalc', 'tidybayes', 'ggmcmc', 'epitools', 
          'performance', 'bayesplot', 'pdftools') 

# NB: package tmvtnorm will require gfortranand gcc to be installed and present in your system path
# /usr/local/lib/

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