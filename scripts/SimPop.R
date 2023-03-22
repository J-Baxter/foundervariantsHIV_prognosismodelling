InitPop <- function(N, H2, donor_vl, recipient_vl){
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


# PDF ('g' from Thompson et al. 2019)
WeightPDF <- function(viralload){
  alpha = -3.55
  sigma <- 0.78/(sqrt(1 - ((2*alpha**2)/(pi*(1 + alpha**2)))))
  mu <-  4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha**2)))
  weight <-  (2/sigma)*dnorm((log10(viralload) - mu)/sigma)*pnorm(alpha*(log10(viralload) - mu)/sigma)
  return(weight)
} 

# Simulates donor population of length, within specified boundaries according to probability distribution
# above

SimDonor <- function(n, transmitter_min, transmitter_max){
  
  spvl_vec <- seq(transmitter_min, transmitter_max, length.out = 10000)
  
  sim_probs <- sapply(spvl_vec, function(x) WeightPDF(x)/sum(WeightPDF(spvl_vec)))
  
  sim_logspvl <- sample(spvl_vec, size = n, prob = sim_probs, replace = F) %>%
    cbind.data.frame('transmitter_log10SPVL' = .)
  
  return(sim_logspvl)
}





