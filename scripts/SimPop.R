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


SimPop <- function(initpop, model, n){
  
  transmitter_prob <- sapply(pop$transmitter_log10SPVL, function(x) WeightPDF(x)/sum(WeightPDF(pop$transmitter_log10SPVL)))
  
  hist(transmitter_prob) #visual check - should look normal(ish) as already log
  
  sim_transmitter_log10SPVL <- InitSimTransmitter(pop_size = NPAIRS,
                                                  transmitter_min = 1, 
                                                  transmitter_max = 7, 
                                                  sample_prob = transmitter_prob)
  
  hist(sim_transmitter_log10SPVL) #visual check - should look uniform
  
  sim_spvl <- predict(h2_model, newdata= data.frame(transmitter_log10SPVL = sim_transmitter_log10SPVL)) %>% 
    cbind.data.frame(sim_recipient_log10SPVL = ., sim_transmitter_log10SPVL = sim_transmitter_log10SPVL) %>%
    mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SPVL')}")) %>%
    `colnames<-` (str_remove(colnames(.), 'sim_')) %>% mutate(sim = 'linear')
}

