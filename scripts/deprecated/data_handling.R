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

ggplot(filter(t, variants == 1), 
       aes(x = transmitter, 
           y = 1 - p, 
           colour = stage))+
  geom_point(shape = 4, size = 3)+
  scale_x_log10(name = expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0,0),
                     limits = c(0,0.6),
                     breaks = seq(0, 0.6, by = 0.1)) +
  annotation_logticks(sides = 'b') +
  my_theme