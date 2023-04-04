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
# Duration of Infection
D <- function(viralload){
  Dmax = 25.4
  Dfifty = 3058
  Dk = 0.41
  
  out <- Dmax*(Dfifty**Dk)/(viralload**Dk + Dfifty**Dk)
  
  
  return(out)
  
}


# The fraction of individuals at seroconversion (the point at which HIV-1 antibodies develop and become
# detectable) with each log(SPVL), vc, in the population at any given time is described by 'g'
WeightSerocons <- function(viralload){
  alpha = -3.55
  sigma <- 0.78/(sqrt(1 - ((2*alpha**2)/(pi*(1 + alpha**2)))))
  mu <-  4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha**2)))
  weight <-  (2/sigma)*dnorm((log10(viralload) - mu)/sigma)*pnorm(alpha*(log10(viralload) - mu)/sigma)
  return(weight)
} 


# Infectiousness
B <- function(viralload){
  Bmax = 0.317
  Bfifty = 13938
  Bk = 1.02
  
  out <- (Bmax*viralload**Bk)/(viralload**Bk + Bfifty**Bk)
  
  return(out)
}


# Probability of observing at least one infection for a given viral load
Obs <- function(viralload){
  d <-  D(viralload)
  
  tp <- B(viralload) * d
  
  p_0 <- ((tp*d)**0*exp(-tp*d))/factorial(tp)
  
  return(1-p_0)
}


## 4.74/0.61 (Zambia - HSX)
## 4.35/0.47 (Amsterdam - MSM)


# Simulates donor population of length n according to probability distributions above
SimDonor <- function(n){
  
  spvl_vec <- spvl_vec <- seq(2, 7, length.out = 10000) %>% raise_to_power(10,.)
  
  p_spvl <- sapply(spvl_vec, function(x) WeightSerocons(x)/sum(WeightSerocons(spvl_vec)))
  
  p_obs <- Obs(spvl_vec)
  
  sim_logspvl <- sample(spvl_vec, size = n, prob = p_spvl*p_obs, replace = F) 
  
  return(sim_logspvl)
}


##TO DO: Estimate Covariance matrix for these features from SHCS so we can simulate
RecipChars <- function(n){
  
  rg <- sample(c('M-M', 'F-M', 'M-F'), n, replace = T, prob = c(0.5, 0.15, 0.35))
  age <- sample(16:60, n, replace = T)
  
  df <- cbind.data.frame(riskgroup = rg, age_r = age) %>%

    # Infer recipient sex from risk group
    separate(riskgroup, c('sex_t', 'sex_r'),  remove = FALSE) 
  
  return(df)
    
}



