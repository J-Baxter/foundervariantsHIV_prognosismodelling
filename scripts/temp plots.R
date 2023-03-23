linear_pred_mp 

plt_2a <- ggplot(pop_tm %>% 
         filter(variants == 1) %>%
         filter(w == 1 ), 
       aes(x = transmitter, 
           y = 1 - p#,
           #colour = as.factor(w)
       ))+
  geom_point(shape= 4, size = 4, colour = '#ef654a') +
  #scale_colour_brewer(palette = 'OrRd') +
  scale_x_log10(name = expression(paste("Transmitter SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  annotation_logticks(sides = 'b') +
  my_theme + theme(legend.position = 'none')

plt_2b <- ggplot(linear_pred_variants %>% 
                   filter(variants == 1) %>%
                   filter(w == 1 ), 
                 aes(x = recipient, 
                     y = 1 - p#,
                     #colour = as.factor(w)
                 ))+
  geom_point(shape= 4, size = 4, colour = '#ef654a') +
  #scale_colour_brewer(palette = 'OrRd') +
  scale_x_log10(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  annotation_logticks(sides = 'b') +
  #coord_flip() + 
  my_theme + theme(legend.position = 'none')


panel_2 <- plot_grid(plt_2a, plt_2b,  align = 'hv', ncol =2)
ggsave('panel_2.eps', device = cairo_ps , plot = panel_2, width = 16, height = 10)


###########################################################################
plt_3a <- ggplot() + 
  geom_point(data = pop_pred_virions, aes(x = transmitter, y = exp.virions, colour = 'virions'), shape= 4,size = 4, colour = '#ef654a') +
  scale_x_log10(limits = c(1, 10**10),
                expand = c(0,0),
                name = expression(paste("Transmitter SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0), name = 'Virions') + 
  my_theme + 
  geom_smooth(method='lm', data = pop_pred_virions, aes(x = transmitter, y = exp.virions, colour = '#ef654a')) + 
  theme(legend.position = 'none')

plt_3b <- ggplot() + 
  geom_point(data = pop_pred_variants, aes(x = transmitter, y = exp.var, colour = 'Variants'), shape = 4,size = 4, colour = '#ef654a') + 
  scale_x_log10(limits = c(1, 10**10),
                expand = c(0,0),
                name = expression(paste("Transmitter SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0), name = 'Variants') + 
  my_theme + 
  geom_smooth(method='lm', data = pop_pred_variants, aes(x = transmitter, y = exp.var, colour = '#ef654a')) + 
  theme(legend.position = 'none')

plt_3c <- ggplot() + 
  geom_point(data = linear_pred_virions, aes(x = recipient, y = exp.virions, colour = 'Virions'), shape= 4,size = 4, colour = '#ef654a') +
  scale_x_log10(limits = c(1, 10**10),
                expand = c(0,0),
                name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0), name = 'Virions') + 
  my_theme + 
  geom_smooth(method='lm', data = linear_pred_virions, aes(x = recipient, y = exp.virions, colour = '#ef654a')) + 
  theme(legend.position = 'none')

plt_3d <- ggplot() + 
  geom_point(data = linear_pred_variants, aes(x = recipient, y = exp.var, colour = 'Variants'), shape = 4,size = 4, colour = '#ef654a') + 
  scale_x_log10(limits = c(1, 10**10),
                expand = c(0,0),
                name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0), name = 'Variants') + 
  my_theme + 
  geom_smooth(method='lm', data = linear_pred_variants, aes(x = recipient, y = exp.var, colour = '#ef654a')) + 
  theme(legend.position = 'none')



variant_decline <- linear_pred_variants %>% 
  summarise(delta_cd4_median = median(cd4_decline), .by = exp.var) %>% 
  rowwise() %>% 
  mutate(t_200 = SimCD4Decline(delta_cd4_median ,init_cd4=1000))%>%
  mutate(exp.var = factor(exp.var, levels = 1:33 ))

virion_decline <- linear_pred_virions %>% summarise(delta_cd4_median = median(cd4_decline), .by = exp.virions) %>% 
  rowwise() %>% 
  filter(exp.virions < 12)%>%
  mutate(t_200 = SimCD4Decline(delta_cd4_median ,init_cd4=1000)) %>%
  mutate(exp.virions = factor(exp.virions, levels = 1:33 )) 

plt_3f <- ggplot() + geom_blank() +
  geom_abline(data = variant_decline, aes(slope = delta_cd4_median*365.25, intercept = 1000, colour = exp.var)) + 
  scale_x_continuous(limits = c(0,14), expand = c(0,0), name = 'Time (Years)', breaks = seq(0,14, by = 2)) + 
  scale_y_continuous(limits = c(0,1200), expand = c(0,0), name = expression(paste("CD4", ' (' ," cells ", mm**-3, ')')), breaks = seq(0,1500, by =200))+ 
  scale_color_brewer(palette = 'OrRd', name = 'Variants') +
  my_theme + 
  theme(legend.position = c(0.8,0.8))

plt_3e <- ggplot() + geom_blank() +
  geom_abline(data = virion_decline, aes(slope = delta_cd4_median*365.25, intercept = 1000, colour = exp.virions)) + 
  scale_x_continuous(limits = c(0,14), expand = c(0,0), name = 'Time (Years)', breaks = seq(0,14, by = 2)) + 
  scale_y_continuous(limits = c(0,1200), expand = c(0,0), name = expression(paste("CD4", ' (' ," cells ", mm**-3, ')')), breaks = seq(0,1500, by =200))+ 
  scale_color_brewer(palette = 'OrRd', name = 'Virions') +
  my_theme + 
  theme(legend.position = c(0.8,0.8))


panel_3 <- plot_grid(plt_3a, plt_3b, plt_3c, plt_3d,plt_3e, plt_3f, align = 'hv', ncol =3, byrow = F)
ggsave('panel_3.eps', device = cairo_ps , plot = panel_3 , width = 24, height = 12)




panel_2 <- plot_grid(plt_2a, plt_2b,  align = 'hv', ncol =2)
ggsave('panel_2.eps', device = cairo_ps , plot = panel_2, width = 16, height = 10)

# delta cd4 ~ against 



###################################################################################################
# Panel 4,5,6: timing of transmission
plt_4a <- ggplot(timing_tm %>% 
                   filter(variants == 1), 
                 aes(x = transmitter, 
                     y = 1 - p,
                     colour = as.factor(w)
                 ))+
  geom_point(shape= 4, size = 4) +
  scale_colour_brewer(palette = 'OrRd', name = 'Acute Infection Weighting') +
  scale_x_log10(name = expression(paste("Transmitter SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  annotation_logticks(sides = 'b') +
  my_theme + theme(legend.position = 'none')

plt_4b <- ggplot(timing_pred_variants %>% 
                   filter(variants == 1) , 
                 aes(x = recipient, 
                     y = 1 - p,
                     colour = as.factor(w)
                 ))+
  geom_point(shape= 4, size = 4) +
  scale_colour_brewer(palette = 'OrRd', name = 'Acute Infection Weighting') +
  scale_x_log10(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  annotation_logticks(sides = 'b') +
  #coord_flip() + 
  my_theme + theme(legend.position = 'none')

legend_4 <- get_legend(
  plt_4b  + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

prow_4 <- plot_grid(plt_4a, plt_4b,  align = 'hv', ncol =2)
panel_4 <- plot_grid(prow_4, legend_4 ,  align = 'hv', ncol = 1, rel_heights = c(1, .1))
ggsave('panel_4.eps', device = cairo_ps , plot = panel_4, width = 18, height = 10)


###################################################################################################


plt_5b <- ggplot() + 
  geom_point(data = transmitter_timing_variants, aes(x = factor(w), y = variants, size =p), shape= 16,colour = '#ef654a') +
  scale_x_discrete(name = 'Weighting') + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0), name = 'Variants') + 
  scale_size(range = c(0,15), name = 'P(X=x)') + 
  my_theme + 
  theme(legend.position = 'none')

plt_5a <- ggplot() + 
  geom_point(data = transmitter_timing_virions, aes(x = factor(w), y = nparticles, size =p_virion), shape= 16,colour = '#ef654a') +
  scale_x_discrete(name = 'Weighting') + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0), name = 'Virions') + 
  scale_size(range = c(0,15), name = 'P(X=x)') + 
  my_theme + 
  theme(legend.position = 'none')

plt_5d <- ggplot() + 
  geom_point(data = timing_pred_variants, aes(x = factor(w), y = variants, size =p), shape= 16,colour = '#ef654a') +
  scale_x_discrete( name = 'Weighting') + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0), name = 'Variants') + 
  scale_size(range = c(0,15), name = 'P(X=x)') + 
  my_theme + 
  theme(legend.position = 'none')

plt_5c <- ggplot() + 
  geom_point(data = timing_pred_virions, aes(x = factor(w), y = nparticles, size =p_virion), shape= 16,colour = '#ef654a') +
  scale_x_discrete(name = 'Weighting') + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0), name = 'Virions') + 
  scale_size(range = c(0,15), name = 'P(X=x)') + 
  my_theme + 
  theme(legend.position = 'none')


legend_5 <- get_legend(
  plt_5b  + 
    guides(color = guide_legend(nrow = 1)) 
)

prow_5 <- plot_grid(plt_5a, plt_5b,plt_5c, plt_5d,  ncol =2, byrow = F)
panel_5 <- plot_grid(prow_5, legend_5 ,  align = 'hv', ncol = 2, rel_widths = c(1, .1))
ggsave('panel_5.eps', device = cairo_ps , plot = panel_5, width = 18, height = 10)

###################################################################################################
variant_decline_timing <- timing_pred_variants %>% 
  summarise(delta_cd4_median = median(cd4_decline), .by = c(exp.var, w)) %>% 
  rowwise() %>% 
  mutate(t_200 = SimCD4Decline(delta_cd4_median ,init_cd4=1000))%>%
  mutate(exp.var = factor(exp.var, levels = 1:33 ))

virion_decline_timing <- timing_pred_virions%>% summarise(delta_cd4_median = median(cd4_decline), .by = c(exp.virions, w))%>% 
  rowwise() %>% 
  filter(exp.virions < 12)%>%
  mutate(t_200 = SimCD4Decline(delta_cd4_median ,init_cd4=1000)) %>%
  mutate(exp.virions = factor(exp.virions, levels = 1:33 )) 


plt_6a <- ggplot() + geom_blank() +
  geom_abline(data = variant_decline_timing, aes(slope = delta_cd4_median*365.25, intercept = 1000, colour = exp.var)) + 
  scale_x_continuous(limits = c(0,14), expand = c(0,0), name = 'Time (Years)', breaks = seq(0,14, by = 2)) + 
  scale_y_continuous(limits = c(0,1200), expand = c(0,0), name = expression(paste("CD4", ' (' ," cells ", mm**-3, ')')), breaks = seq(0,1500, by =200))+ 
  scale_color_brewer(palette = 'OrRd', name = 'Variants') +
  facet_wrap(~w, nrow= 1) + 
  my_theme + 
  theme(panel.spacing = unit(2, "lines"))

plt_6b <- ggplot() + geom_blank() +
  geom_abline(data = virion_decline_timing, aes(slope = delta_cd4_median*365.25, intercept = 1000, colour = exp.virions)) + 
  scale_x_continuous(limits = c(0,14), expand = c(0,0), name = 'Time (Years)', breaks = seq(0,14, by = 2)) + 
  scale_y_continuous(limits = c(0,1200), expand = c(0,0), name = expression(paste("CD4", ' (' ," cells ", mm**-3, ')')), breaks = seq(0,1500, by =200))+ 
  scale_color_brewer(palette = 'OrRd', name = 'Virions') +
  facet_wrap(~w, nrow= 1) + 
  my_theme + 
  theme(panel.spacing = unit(2, "lines"))

panel_6 <- plot_grid( plt_6b + theme(strip.background = element_blank(), strip.text = element_text(size = 26)), plt_6a + theme(strip.text.x = element_blank()), align = 'hv', nrow = 2)
ggsave('panel_6.eps', device = cairo_ps , plot = panel_6, width = 18, height = 10)

