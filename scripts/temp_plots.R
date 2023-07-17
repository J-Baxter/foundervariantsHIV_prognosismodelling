combined_data_CD4_PMV_wide_long <- combined_data_CD4_PMV %>%
  filter(w == 1) %>%
  pivot_wider(names_from = partner, values_from = c(sex, age.inf, riskgroup, SpVL, log10_SpVL, age.inf_category, delta_CD4, starts_with('p_')), names_sep = '_') %>%
  pivot_longer(cols = starts_with('p_'), names_to = 'x', values_to = 'probability') %>%
  separate(x, sep = '_', c(NA, 'type', 'x', 'partner'))

ggplot(combined_data_CD4_PMV_wide %>% filter(dataset == 'shcs_predicted'), 
      aes(x = 1 - p_variants_1_1, 
          y = SpVL_2 #,
          #colour = as.factor(w)
      ))+
  geom_point(shape= 4, size = 4, colour = '#ef654a') +
  #scale_colour_brewer(palette = 'OrRd') +
  scale_y_log10(name = expression(paste("Partner 2 SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**6),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_x_continuous(name = 'P(Multiple Variant Partner 1)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  annotation_logticks(sides = 'l') +
  my_theme

ggplot(combined_data_CD4_PMV_wide %>% filter(dataset == 'shcs_predicted'), 
       aes(x = 1 - p_variants_1_2, 
           y = delta_CD4_1 #,
           #colour = as.factor(w)
       ))+
  geom_point(shape= 4, size = 4, colour = '#ef654a') +
  #scale_colour_brewer(palette = 'OrRd') +
  scale_y_continuous(name = expression(paste('Partner 1', Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)),
                     expand = c(0,0),
                     limits = c(-0.6,0)) +
  scale_x_continuous(name = 'P(Multiple Variants Partner 2)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  geom_smooth(method = 'lm')+
  my_theme


plt_2e <-  ggplot(combined_data_CD4_PMV_wide_long %>% filter(type == 'variants'), aes(y = SpVL_1, x = x, fill = probability), data =  linear_uw_variants ) + 
  geom_tile()+
  scale_fill_distiller(palette = 8, 
                       direction = 1, 
                       values = c(0.0005, 0.001, 0.01, 0.1, 0.5, 0.6, 0.7, 0.8,  0.85),
                       limits =c( 0.0005,0.88), 
                       labels = c(0.001, 0.01, 0.1, 0.5, 0.8),
                       na.value = "white",
                       'P(X=X)') + 
  scale_y_log10(limits = c(10**4, 10**5.5),
                expand = c(0,0),
                name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_x_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0.1), name = 'Variants') + 
  my_theme + 
  annotation_logticks(sides = 'l') +
  theme(legend.position = 'none')


plt_2f <- ggplot(aes(y = recipient_rounded, x = nparticles, fill = mean_p_virion), data = linear_uw_virions ) + 
  geom_tile()+
  scale_fill_distiller(palette = 8, 
                       direction = 1, 
                       values = c(0.0005, 0.001, 0.01, 0.1, 0.5, 0.6, 0.7, 0.8,  0.85) , 
                       labels = c(0.001, 0.01, 0.1, 0.5, 0.8),
                       limits =c( 0.0005,0.88),
                       na.value = "white",
                       'P(X=X)') + 
  scale_y_log10(limits = c(10**4, 10**5.5),
                expand = c(0,0),
                name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_x_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0.1), name = 'Particles') + 
  my_theme + 
  annotation_logticks(sides = 'l') +
  theme(legend.position = 'none')