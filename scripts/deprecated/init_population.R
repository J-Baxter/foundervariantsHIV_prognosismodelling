# Initialise study population
# SPVL distributions from Hollingsworth et al., PLoS 2010 https://doi.org/10.1371/journal.ppat.1000876

# Initialisation function
InitPop <- function(N, H2, donor_vl, recipient_vl){
  
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


WeightPDF <- function(viralload){
  alpha = -3.55
  sigma <- 0.78/(sqrt(1 - ((2*alpha**2)/(pi*(1 + alpha**2)))))
  mu <-  4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha**2)))
  weight <-  (2/sigma)*dnorm((log10(viralload) - mu)/sigma)*pnorm(alpha*(log10(viralload) - mu)/sigma)
  return(weight)
} 


InitSimTransmitter <- function(pop_size, transmitter_min, transmitter_max, sample_prob){
  sim_range <- seq(transmitter_min, transmitter_max, length.out = pop_size)
  sim_logspvl <- sample(sim_range, size = pop_size, prob = sample_prob, replace = T)
  
  return(sim_logspvl)
}


## END ##
