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

