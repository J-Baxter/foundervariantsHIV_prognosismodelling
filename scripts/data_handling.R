t <- cbind.data.frame(transmitter = pop$transmitter, prob_recip_multiple) %>% 
  rename(probTransmissionPerSexAct_average = probTransmissionPerSexAct) %>%
  rename_with(function(x) gsub('\\.(?<=\\V)', '\\1_average.\\2',x, perl= TRUE), 
              .cols = !contains('chronic') & contains('V')) %>%
  pivot_longer(cols = c(contains('chronic'),contains('average')), 
               names_to = c('type', 'stage', 'variants'), 
               names_pattern = "(^[^_]+)(_[^.]*)(.*)", 
               values_to = 'p') %>% 
  #regex wrong above, but tidied up below
  mutate(stage = gsub('^.*_', '', stage)) %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric()) 

str(t)