# Sensitivity Analysis
# Robustness of model system to proportion of standard error incorporated in Heritability
# Model (linear, unweighted only)

linear_uw_pop_SA <- sim_donor %>% 
  cbind.data.frame(transmitter = sim_donor) %>% 
  mutate(transmitter_log10SpVL= log10(transmitter)) %>%
  
  # Bind simulated recipient characteristics
  cbind.data.frame(sim_recip_chars) %>%
  
  # Predict recipient SpVL according to heritability model
  predict(linear_model_uw, .) %>%
  cbind.data.frame(recipient_log10SpVL = .) %>%
  
  # Because A) we are predicting out-of-sample and B) our focus is the population level
  # we incorporate residual standard error into our predictions
  mutate(recipient_log10SpVL_exact = recipient_log10SpVL + rnorm(NSIM, sd = (summary(linear_model_uw)$sigma)/1)) %>%
  mutate(recipient_log10SpVL_2 = recipient_log10SpVL + rnorm(NSIM, sd = (summary(linear_model_uw)$sigma)/2)) %>%
  mutate(recipient_log10SpVL_4 = recipient_log10SpVL + rnorm(NSIM, sd = (summary(linear_model_uw)$sigma)/4)) %>%
  mutate(recipient_log10SpVL_6 = recipient_log10SpVL + rnorm(NSIM, sd = (summary(linear_model_uw)$sigma)/6)) %>%
  mutate(recipient_log10SpVL_8 = recipient_log10SpVL + rnorm(NSIM, sd = (summary(linear_model_uw)$sigma)/8)) %>%
  mutate(recipient_log10SpVL_10 = recipient_log10SpVL + rnorm(NSIM, sd = (summary(linear_model_uw)$sigma)/10)) %>%
  
  # Bind predicted recipient SpVl with transmission pair characteristics
  cbind.data.frame(transmitter_log10SpVL= log10(sim_donor)) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SpVL')}")) %>%
  `colnames<-` (str_remove(colnames(.), 'sim_')) %>% 
  cbind.data.frame(sim_recip_chars) %>%
  mutate(model = 'linear_uw') %>%
  mutate(transmitter = round(transmitter)) %>%
  
  # Join Transmission Model predictions for P(MV)
  left_join(linear_uw_variants %>% filter(variants == 1) %>% select(transmitter, p)) %>%
  select(matches('recipient_[[:digit:]]') | matches('recipient_exact'), p) %>%
  pivot_longer(cols = contains('recipient'), values_to = 'recipient', names_to = 'SE') %>% 
  mutate(SE = factor(SE, levels = c('recipient_exact', 'recipient_2', 'recipient_4', 'recipient_6', 'recipient_8', 'recipient_10'), 
                     labels =  c('100%', '50%', '25%', '16.7%', '12.5%', '10%') ))


# Plot
plt_s4 <- ggplot(linear_uw_pop_SA, 
                 aes(x = 1 - p, 
                     y = recipient #,
                     #colour = as.factor(w)
                 ))+
  geom_point(shape= 4, size = 4, colour = '#ef654a') +
  #scale_colour_brewer(palette = 'OrRd') +
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**3, 10**6),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  annotation_logticks(sides = 'l',axis.title.y = element_text(family = 'sans')) +
  facet_wrap(.~SE, labeller = )+
  my_theme + theme(legend.position = 'none')


## END ##

