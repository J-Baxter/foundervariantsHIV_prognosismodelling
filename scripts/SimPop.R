WeightPDF <- function(viralload){
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



InitSimTransmitter <- function(pop_size, transmitter_min, transmitter_max, sample_prob){
  sim_range <- seq(transmitter_min, transmitter_max, length.out = pop_size)
  sim_logspvl <- sample(sim_range, size = pop_size, prob = sample_prob, replace = T)
  
  return(sim_logspvl)
}




SimPop <- function(initpop, model, n, modeltype = 'linear'){
  
  stopifnot('transmitter_log10SPVL' %in% colnames(initpop))
  t <- pop[['transmitter_log10SPVL']]
  
  transmitter_prob <- sapply(t, function(x) WeightPDF(x)/sum(WeightPDF(t)))
  
  hist(transmitter_prob) #visual check - should look normal(ish) as already log
  
  sim_transmitter_log10SPVL <- InitSimTransmitter(pop_size = NPAIRS,
                                                  transmitter_min = 1, 
                                                  transmitter_max = 7, 
                                                  sample_prob = transmitter_prob)
  
  hist(sim_transmitter_log10SPVL) #visual check - should look uniform
  
  sim_spvl <- predict(model, newdata= data.frame(transmitter_log10SPVL = sim_transmitter_log10SPVL)) %>% 
    cbind.data.frame(sim_recipient_log10SPVL = ., sim_transmitter_log10SPVL = sim_transmitter_log10SPVL) %>%
    mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SPVL')}")) %>%
    `colnames<-` (str_remove(colnames(.), 'sim_')) %>% mutate(model = modeltype)
  
  return(sim_spvl)
}

