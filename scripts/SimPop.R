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





