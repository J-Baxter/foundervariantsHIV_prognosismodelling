# Initialise study population
# SPVL distributions from Hollingsworth et al., PLoS 2010 https://doi.org/10.1371/journal.ppat.1000876

# Initialisation function
InitPop <- function(N, H2, donor_vl, recipient_vl){
  require(MASS)
  require(tidyverse)
  
  matrix <- matrix(c(1, 1-(1/H2),
                     1-(1/H2), 1), ncol = 2)
  
  d <-diag(c(donor_vl[["sd"]], recipient_vl[["sd"]]))
  
  S <- d * matrix * d
  
  paired_spvl <- MASS::mvrnorm(n = N, 
                               mu = c('transmitter_log10SPVL' = donor_vl[["mean"]],
                                      'recipient_log10SPVL' = recipient_vl[["mean"]]), 
                               Sigma = S) %>% 
    data.frame() %>%
    mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SPVL')}"))
  
  return(paired_spvl)
}

## END ##
