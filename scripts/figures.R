# Figures

# Dependencies
source('./scripts/dependencies.R')
font_add_google("Questrial", "Questrial")
set_null_device(cairo_pdf)
showtext_auto()


############################################## Panel 1 ##############################################
#3 model outputs (tolerance, herit, jp transmission)

# For Plot 1a, run ./scripts/model_dag.R
source('./scripts/model_dag.R')

plt_1b <- ggplot(pop, aes(x = transmitter, recipient)) +
  geom_point( #'#CB6015' #'#66c2a4','#2ca25f','#006d2c'
    colour = '#ef654a',
    shape = 4, size = 3) +
  scale_x_log10(limits = c(10**2, 10**7),
                expand = c(0.05,0),
                name = expression(paste("Transmitter SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  geom_smooth(method = 'lm', colour = '#ef654a' ) + 
  #geom_function() + Incorporate lm with fixed effects (with #Add 95% CIs to lines (geom_ribbon))
  annotation_logticks() +
  my_theme + 
  theme(legend.position = 'none')

plt_1c <- ggplot(data, aes(x = SpVL, y=CD4.decline))+
  geom_point(colour = '#ef654a', shape= 4, alpha = 0.3) +
  scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)),
                     expand = c(0,0),
                     limits = c(-2,2))+
  scale_x_continuous(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                     expand = c(0,0))+
  geom_function(fun = function(x) ((-5.6e-3) + (-1.6e-4)*20) * x**2, colour = '#fdbb84', linewidth = 1.2, xlim = c(0,6.3)) +
  geom_text(aes(x = 6.65, y = -0.35), label = "20 yrs", colour = "#fdbb84", size = 6) +
  geom_function(fun = function(x) ((-5.6e-3) + (-1.6e-4)*40) * x**2, colour = '#ef6548', linewidth = 1.2, xlim = c(0,6.3))+
  geom_text(aes(x = 6.65, y = -0.5, label = "40 yrs"), colour = "#ef6548", size = 6) +
  geom_function(fun = function(x) ((-5.6e-3) + (-1.6e-4)*60) * x**2, colour = '#b30000', linewidth = 1.2, xlim = c(0,6.3)) +
  geom_text(aes(x = 6.65, y = -0.62, label = "60 yrs"), colour = "#b30000", size = 6) + #Add 95% CIs to lines (geom_ribbon)
  my_theme + theme(axis.title.y = element_text(family = 'sans')) + 
  annotation_logticks(sides = 'b') 


plt_1d <- ggplot(joint_probs, aes(y = variants, x = nparticles))+
  geom_point(aes(size = p), colour = '#ef654a')+
  scale_size(range = c(0,15), name = 'P(Variants \u2229 Virions)')+
  scale_y_continuous(name = 'Variants', expand = c(0,0), limits = c(0.5,10.5), breaks = 1:10)+ 
  scale_x_continuous(name = 'Virions', expand = c(0,0), limits = c(0.5,10.5), breaks = 1:10)+
  my_theme+
  theme(legend.position =  c(0.25,0.77),
        legend.title = element_text(family = 'sans'))


panel_1 <- plot_grid(plt_1a, plt_1b, plt_1c, plt_1d, ncol = 2, labels = 'AUTO')
cat('Panel 1 complete. \n')


############################################## Panel 2 ##############################################
#Confounder
plt_2a <- ggplot(shcs_transmitters  %>% 
                   filter(variants == 1) %>%
                   filter(w == 1 ), 
                 aes(x = transmitter, 
                     y = 1 - p#,
                     #colour = as.factor(w)
                 ))+
  geom_point(shape= 4, size = 4, colour = '#ef654a') +
  #scale_colour_brewer(palette = 'OrRd') +
  scale_x_log10(name = expression(paste("Transmitter SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  annotation_logticks(sides = 'b') +
  my_theme + theme(legend.position = 'none')

plt_2b <- ggplot(linear_uw_variants %>% 
                   filter(variants == 1) %>%
                   filter(w == 1 ), 
                 aes(x = recipient, 
                     y = 1 - p#,
                     #colour = as.factor(w)
                 ))+
  geom_point(shape= 4, size = 4, colour = '#ef654a') +
  #scale_colour_brewer(palette = 'OrRd') +
  scale_x_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
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



plt_2c <- ggplot(aes(x = recipient_rounded, y = nparticles, size = mean_p_virion, colour = mean_p_virion), data = linear_uw_virions ) + 
  geom_point() + 
  scale_color_distiller(palette = 'OrRd') + 
  scale_size(range = c(0.00001,10), name = 'P(X=x)') + 
  scale_x_log10(limits = c(10**2, 10**7),
                expand = c(0,0),
                name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0.1), name = 'Virions') + 
  my_theme + 
  annotation_logticks(sides = 'b') +
  #geom_smooth(method='lm', data = linear_pred_virions, aes(x = recipient, y = exp.virions, colour = '#ef654a')) + 
  theme(legend.position = 'none')

plt_2d <-  ggplot(aes(x = recipient_rounded, y = variants, size = mean_p, colour = mean_p), data = linear_uw_variants ) + 
  geom_point() + 
  scale_color_distiller(palette = 'OrRd') + 
  scale_size(range = c(0,10), name = 'P(X=x)') +
  scale_x_log10(limits = c(10**2, 10**7),
                expand = c(0,0),
                name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0.1), name = 'Variants') + 
  my_theme + 
  annotation_logticks(sides = 'b') +
  #geom_smooth(method= VGAM::vglm , method.args = list(family = "poisson"), data = linear_pred_variants, aes(x = recipient, y = exp.var, colour = '#ef654a')) + # should be truncated poisson
  theme(legend.position = 'none')


panel_2 <- plot_grid(plt_2a, plt_2b, plt_2c, plt_2d, ncol = 2, labels = 'AUTO')
cat('Panel 2 complete. \n')
############################################## Panel 4 ##############################################
#Non-Linear

models = c(concave_uw ="e^{X}",
           convex_uw = 'ln(X)',
           linear_uw = 'X')

all_pops <- rbind(linear_uw_pop, concave_uw_pop, convex_uw_pop) %>%
  mutate(model = str_replace_all(model, models ))


all_vars <- rbind(linear_uw_variants , concave_uw_variants , convex_uw_variants) %>%
  filter(w == 1) %>%
  mutate(model = str_replace_all(model, models ))

all_virions <- rbind(linear_uw_virions , concave_uw_virions, convex_uw_virions) %>%
  filter(w == 1) %>%
  mutate(model = str_replace_all(model, models ))

plt_3a <- ggplot(all_pops, aes(x = transmitter, recipient)) +
  geom_point(
    shape = 4, size = 3, colour = '#ef654a') +
  scale_x_log10(limits = c(10**2, 10**7),
                expand = c(0,0),
                name = expression(paste("Transmitter SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  facet_wrap(.~model, labeller = label_parsed) +
  annotation_logticks() +
  my_theme + 
  theme(legend.position = 'none', 
        panel.spacing = unit(2, "lines"), 
        strip.background = element_blank()) 

plt_3b <- ggplot(all_vars %>% 
                   filter(variants == 1) , 
                 aes(x = recipient, 
                     y = 1 - p
                 ))+
  geom_point(shape= 4, size = 4, colour = '#ef654a' ) +
  #scale_colour_brewer(palette = 'OrRd') +
  scale_x_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  annotation_logticks(sides = 'b') +  facet_wrap(.~model, labeller = label_parsed) +
  #coord_flip() + 
  my_theme + theme(legend.position = 'none', 
                   panel.spacing = unit(2, "lines"), 
                   strip.background = element_blank())

ppt_panel_3 <- cowplot::plot_grid(plt_3a, plt_3b , nrow = 2, labels = 'AUTO', align = 'HV')
ggsave(plot = ppt_panel_3, filename = paste(figs_dir,sep = '/', "ppt_panel_3.jpeg"), device = jpeg, width = 14, height = 14) # Non - Linear


plt_3c <-  ggplot(aes(x = recipient_rounded, y = variants, size = mean_p, colour = mean_p), data = all_vars ) + 
  geom_point() + 
  scale_color_distiller(palette = 'OrRd') + 
  scale_size(range = c(0,10), name = 'P(X=x)') +
  scale_x_log10(limits = c(10**2, 10**7),
                expand = c(0,0),
                name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0.1), name = 'Variants') + 
  my_theme + 
  annotation_logticks(sides = 'b') +  
  facet_wrap(.~model,labeller = label_parsed) +
  theme(legend.position = 'none')

plt_3d <-  ggplot(aes(x = recipient_rounded, y = nparticles, size = mean_p_virion, colour = mean_p_virion), data = all_virions ) + 
  geom_point() + 
  scale_color_distiller(palette = 'OrRd') + 
  scale_size(range = c(0,10), name = 'P(X=x)') +
  scale_x_log10(limits = c(10**2, 10**7),
                expand = c(0,0),
                name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0.1), name = 'Virions') + 
  my_theme + 
  annotation_logticks(sides = 'b') +  
  facet_wrap(.~model,labeller = label_parsed) +
  theme(legend.position = 'none')

ppt_panel_4 <- cowplot::plot_grid(plt_3c, plt_3d , nrow = 2, labels = 'AUTO', align = 'HV')
ggsave(plot = ppt_panel_4, filename = paste(figs_dir,sep = '/', "ppt_panel_4.jpeg"), device = jpeg, width = 14, height = 14) # Non - Linear

plt_3e <-  ggplot(aes(x = model, y =cd4_decline), data = all_vars) + 
  geom_boxplot(colour = '#ef654a') + 
  scale_size(range = c(0,10), name = 'P(X=x)') +
  scale_x_discrete(name = 'Model') + 
  scale_y_continuous(limits = c(-1,0), expand = c(0,0.1), name =expression(paste(Delta, ' CD4'))) + 
  my_theme + 
  theme(axis.title.y = element_text(family = 'sans'))






cat('Panel 3 complete. \n')
############################################## Panel 5 ##############################################
# Timing
timing_all_vars <- rbind(linear_uw_variants_timing , concave_uw_variants_timing , convex_uw_variants_timing) %>%
  mutate(model = str_replace_all(model, models ))

my_palette <- brewer.pal(name="OrRd",n=9)[4:9]

plt_4a <- ggplot(timing_all_vars   %>% 
                   filter(variants == 1) , 
                 aes(x = recipient, 
                     y = 1 - p,
                     colour = as.factor(w)
                 ))+
  geom_point(shape= 4, size = 4) +
  scale_colour_manual(values = my_palette ) +
  scale_x_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  annotation_logticks(sides = 'b') +
  facet_wrap(.~model)+
  #coord_flip() + 
  my_theme + theme(legend.position = 'none', 
                   panel.spacing = unit(2, "lines"), 
                   strip.background = element_blank())

plt_4b <- ggplot(timing_linearonly_vars  %>% 
                   filter(variants == 1) , 
                 aes(x = recipient, 
                     y = 1 - p,
                     colour = w
                 ))+
  geom_point(shape= 4, size = 4) +
  scale_colour_brewer(palette = 'OrRd') +
  scale_x_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  annotation_logticks(sides = 'b') +  facet_wrap(.~model) +
  #coord_flip() + 
  my_theme + theme(legend.position = 'none', 
                   panel.spacing = unit(2, "lines"), 
                   strip.background = element_blank())


plt_4c <-  ggplot(aes(x = recipient_rounded, y = variants, size = mean_p), data = timing_all_vars  ) + 
  geom_point(colour = '#ef654a') + 
  scale_size(range = c(0,10), name = 'P(X=x)') +
  scale_x_log10(limits = c(1, 10**8),
                expand = c(0,0),
                name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_y_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0.1), name = 'Variants') + 
  my_theme + 
  annotation_logticks(sides = 'b') +  
  facet_wrap(.~w) +
  theme(legend.position = 'none')

# CD4 decline


cat('Panel 4 complete. \n')
cat('All plots complete. \n')