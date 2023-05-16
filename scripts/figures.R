# Figures

# Dependencies
source('./scripts/dependencies.R')
font_add_google("Questrial", "Questrial")
set_null_device(cairo_pdf)
showtext_auto()


############################################## Panel 1 ##############################################
#3 model framework and individual model outputs (tolerance, herit, jp transmission)
plt_1a <- ggdraw() +
  draw_image("dag.svg")

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
  scale_size(range = c(0,15), name = 'P(Variants \u2229 Particles)')+
  scale_y_continuous(name = 'Variants', expand = c(0,0), limits = c(0.5,10.5), breaks = 1:10)+ 
  scale_x_continuous(name = 'Particles', expand = c(0,0), limits = c(0.5,10.5), breaks = 1:10)+
  my_theme+
  theme(legend.position =  c(0.25,0.77),
        legend.title = element_text(family = 'sans'))


panel_1 <- plot_grid(plt_1a, plt_1b, plt_1c, plt_1d, ncol = 2, align = 'hv', labels = 'AUTO')
cat('Panel 1 complete. \n')


############################################## Panel 2 ##############################################
# What do we expect to observe?

plt_2a <- ggplot(linear_uw_variants %>% 
                   filter(variants == 1) %>%
                   filter(w == 1 ), 
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
  #coord_flip() + 
  my_theme + theme(legend.position = 'none')


plt_2b <- ggplot(linear_uw_variants %>% 
                   filter(variants == 1) %>%
                   filter(w == 1 ), 
                 aes(x = 1 - p, 
                     y = cd4_decline #,
                     #colour = as.factor(w)
                 ))+
  geom_point(shape= 4, size = 4, colour = '#ef654a') +
  #scale_colour_brewer(palette = 'OrRd') +
  scale_y_continuous(limits = c(-0.5,0), 
                     expand = c(0,0.01), 
                     breaks = seq(-0.5, 0, by = 0.1),
                     name =expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1))) + 
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  #coord_flip() + 
  my_theme + theme(legend.position = 'none', axis.title.y = element_text(family = 'sans'))


plt_2c <- ggplot(linear_uw_virions %>% 
                   filter(nparticles == 1) %>%
                   filter(w == 1 ), 
                 aes(x = 1 - p_virion, 
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
  scale_x_continuous(name = 'P(Multiple Particle Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,1),
                     breaks = seq(0, 1, by = 0.1)) +
  annotation_logticks(sides = 'l') +
  #coord_flip() + 
  my_theme + theme(legend.position = 'none')


plt_2d <- ggplot(linear_uw_virions %>% 
                   filter(nparticles == 1) %>%
                   filter(w == 1 ), 
                 aes(x = 1 - p_virion, 
                     y = cd4_decline #,
                     #colour = as.factor(w)
                 ))+
  geom_point(shape= 4, size = 4, colour = '#ef654a') +
  #scale_colour_brewer(palette = 'OrRd') +
  scale_y_continuous(limits = c(-0.5,0), 
                     expand = c(0,0.01), 
                     breaks = seq(-0.5, 0, by = 0.1),
                     name =expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1))) +
  scale_x_continuous(name = 'P(Multiple Particle Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,1),
                     breaks = seq(0, 1, by = 0.1)) +
  annotation_logticks(sides = 'l') +
  #coord_flip() + 
  my_theme + theme(legend.position = 'none', axis.title.y = element_text(family = 'sans'))


plt_2e <-  ggplot(aes(y = recipient_rounded, x = variants, fill = mean_p), data =  linear_uw_variants ) + 
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


panel2_tile_legend <-  get_legend(plt_2e+ theme(legend.position = 'bottom', legend.key.width=unit(2,"cm"))) 
panel_2_upper <- plot_grid(plt_2a, plt_2b, plt_2c, plt_2d, plt_2e, plt_2f,  ncol = 2, align = 'hv', labels = 'AUTO')

panel_2 <- plot_grid(panel_2_upper, panel2_tile_legend, rel_heights = c(1,0.1),nrow = 2)


cat('Panel 2 complete. \n')
############################################## Panel 3 ##############################################
#Non-Linear

models = c(concave_uw ="e^{X}",
           convex_uw = 'ln(X)',
           linear_w = 'X+P(MV)')

all_pops <- rbind(linear_w_pop, concave_uw_pop, convex_uw_pop) %>%
  mutate(model = str_replace_all(model, models ))


all_vars <- rbind(linear_w_variants , concave_uw_variants , convex_uw_variants) %>%
  filter(w == 1) %>%
  mutate(model = str_replace_all(model, models ))

all_virions <- rbind(linear_w_virions , concave_uw_virions, convex_uw_virions) %>%
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
                 aes(y = recipient, 
                     x = 1 - p
                 ))+
  geom_point(shape= 4, size = 4, colour = '#ef654a' ) +
  #scale_colour_brewer(palette = 'OrRd') +
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  annotation_logticks(sides = 'l') +
  facet_wrap(.~model, labeller = label_parsed) +
  my_theme + theme(legend.position = 'none', 
                   panel.spacing = unit(2, "lines"), 
                   strip.background = element_blank())

plt_3c <-  ggplot(all_vars %>% 
                    filter(variants == 1) , 
                  aes(x = 1 - p, 
                      y = cd4_decline
                  )) + 
  geom_point(colour = '#ef654a', shape= 4, size = 3) + 
  scale_x_continuous(limits = c(0.1, 0.5),
                     expand = c(0,0),
                     name = 'P(Multiple Variant Recipient)')+
  
  scale_y_continuous(limits = c(-0.5,0), 
                     expand = c(0,0.01), 
                     breaks = seq(-0.5, 0, by = 0.1),
                     name =expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1))) +
  facet_wrap(.~model,labeller = label_parsed) + 
  my_theme + 
  theme(axis.title.y = element_text(family = 'sans'),panel.spacing = unit(3, "lines"))

panel_3 <- cowplot::plot_grid(plt_3a, plt_3b , plt_3c,  nrow = 3, labels = 'AUTO', align = 'HV')
cat('Panel 3 complete. \n')


############################################## Panel 4 ##############################################
# Timing
my_palette <- brewer.pal(name="OrRd",n=9)[4:9]

plt_4a <- ggplot(linear_uw_variants_timing %>% 
                   filter(variants == 1) , 
                 aes(x = 1 - p, 
                     y = recipient,
                     colour = as.factor(w)
                 ))+
  geom_point(size = 3, shape = 4 ) +
  #scale_shape_manual(values = c(0,1,2,3)) + 
  scale_colour_manual(values = my_palette, 'Weight to Early Infection') +
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  my_theme + theme(legend.position = 'none', 
                   panel.spacing = unit(2, "lines"), 
                   strip.background = element_blank())

plt_4b <- ggplot(linear_uw_variants_timing %>% filter(variants == 1), 
                 aes(y = cd4_decline, 
                     x = 1 - p,
                     colour =  as.factor(w)
                 ))+
  geom_point(size = 3, shape = 4  ) + #colour = '#ef654a'
  #scale_shape_manual(values = c(0,1,2,3)) + 
  scale_colour_manual(values = my_palette, 'Weight to Early Infection') +
  scale_y_continuous(limits = c(-0.5,0), 
                     expand = c(0,0.01), 
                     breaks = seq(-0.5, 0, by = 0.1),
                     name =expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)))  + 
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.5, by = 0.1)) +
  my_theme + theme(legend.position = 'none',
                   panel.spacing = unit(2, "lines"), 
                   strip.background = element_blank(),
                   axis.title.y = element_text(family = 'sans'))



# Mean p calculated over all w - need to adjust to plot these
plt_4c <-  ggplot(aes(y = recipient_rounded, x = variants, fill = mean_p), data =  linear_uw_variants_timing) + 
  geom_tile()+
  scale_fill_distiller(palette = 8, 
                       direction = 1, 
                       values = c(0.0005, 0.001, 0.01, 0.1, 0.5, 0.6, 0.7, 0.8,  0.85, 0.9, 1) , 
                       labels = c(0.001, 0.01, 0.1, 0.5, 0.99),
                       limits =c( 0.0005,1),
                       na.value = "white",
                       'P(X=X)') + 
  scale_y_log10(limits = c(10**4, 10**5.5),
                expand = c(0,0),
                name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_x_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0.1), name = 'Variants') + 
  facet_wrap(~w, nrow = 1, scales = "fixed")+
  my_theme + 
  annotation_logticks(sides = 'l') +
  theme(legend.position = 'none')


plt_4d <- ggplot(aes(y = recipient_rounded, x = nparticles, fill = mean_p_virion), data = linear_uw_virions_timing ) + 
  geom_tile()+
  scale_fill_distiller(palette = 8, 
                       direction = 1, 
                       values = c(0.0005, 0.001, 0.01, 0.1, 0.5, 0.6, 0.7, 0.8,  0.85, 0.9, 1) , 
                       labels = c(0.001, 0.01, 0.1, 0.5, 0.99),
                       limits =c( 0.0005,1),
                       na.value = "white",
                       'P(X=X)')  + 
  scale_y_log10(limits = c(10**4, 10**5.5),
                expand = c(0,0),
                name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_x_continuous(limits = c(0,10), breaks = 1:10, expand = c(0,0.1), name = 'Particles') + 
  facet_wrap(~w, nrow = 1, scales = "fixed")+
  my_theme + 
  annotation_logticks(sides = 'l') +
  theme(legend.position = 'none')

panel4_bottom_legend <-  get_legend(plt_4d + theme(legend.position = 'bottom', legend.key.width=unit(2,"cm"))) 
panel4_upper_legend <-  get_legend(plt_4a+ theme(legend.position = 'bottom')) 

panel_4_upper <- plot_grid(plt_4a, plt_4b, align = 'hv', labels = "AUTO", ncol = 2)
panel_4 <- plot_grid(panel_4_upper,panel4_upper_legend , plt_4c, plt_4d, panel4_bottom_legend,
                     align = 'hv', labels = c(NA, NA, 'C', 'D', NA), nrow = 5,
                     rel_heights = c(1,0.1,1,1,0.1))
panel_4


cat('Panel 4 complete. \n')
cat('All plots complete. \n')