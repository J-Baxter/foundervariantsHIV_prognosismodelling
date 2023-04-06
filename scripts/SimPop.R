InitPop <- function(N, H2, transmitter_SpVL, recipient_SpVL){
  require(tidyverse)
  
  matrix <- matrix(c(1, 1-(1/H2),
                     1-(1/H2), 1), ncol = 2)
  
  d <-diag(c(transmitter_SpVL[["sd"]], recipient_SpVL[["sd"]]))
  
  S <- d * matrix * d
  
  paired_SpVL <- MASS::mvrnorm(n = N, 
                               mu = c('transmitter_log10SpVL' = transmitter_SpVL[["mean"]],
                                      'recipient_log10SpVL' = recipient_SpVL[["mean"]]), 
                               Sigma = S) %>% 
    data.frame() %>%
    mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SpVL')}"))
  
  return(paired_SpVL)
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
# detectable) with each log(SpVL), vc, in the population at any given time is described by 'g'
WeightSerocons <- function(viralload){
  alpha = -3.55
  sigma <- 0.78/(sqrt(1 - ((2*alpha**2)/(pi*(1 + alpha**2)))))
  mu <-  4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha**2)))
  weight <-  (2/sigma)*dnorm((log10(viralload) - mu)/sigma)*pnorm(alpha*(log10(viralload) - mu)/sigma)
  return(weight)
} 

logSpVL <- seq(1, 7, by = 0.01)
Prob <- WeightPDF(logSpVL)/sum(WeightPDF(logSpVL))
mydonorspVL <- sample(logSpVL, prob = Prob) 

fitted <- brms::brm(formula = recipient_log10SPVL ~ transmitter_log10SPVL, data = pop)
preds <- posterior_predict(fitted , cbind.data.frame(transmitter_log10SPVL = mydonorspVL), ndraws = 10) %>% apply(., 2, sample, 1)



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
  
  SpVL_vec <- SpVL_vec <- seq(2, 7, length.out = 10000) %>% raise_to_power(10,.)
  
  p_SpVL <- sapply(SpVL_vec, function(x) WeightSerocons(x)/sum(WeightSerocons(SpVL_vec)))
  
  p_obs <- Obs(SpVL_vec)
  
  sim_SpVL <- sample(SpVL_vec, size = n, prob = p_SpVL*p_obs, replace = F) 
  
  return(sim_SpVL)
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



