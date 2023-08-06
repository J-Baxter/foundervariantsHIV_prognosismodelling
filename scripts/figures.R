# Figures

# Dependencies
source('./scripts/dependencies.R')
source('./scripts/global_theme.R')
#set_null_device(cairo_pdf)

############################################## Panel 1 ##############################################
# From observe sig effect

# Normal distributions of recipient phenotype, multiple variant and single variant
facet.labs <- c("Recipient Population", "Multiple Variant Recipients", "Single Variant Recipients") %>%
  setNames(c('recipient_mean', 'recipient_mv', 'recipient_sv'))

plt1a <- ggplot(vls ) + 
  geom_density(aes(x = vl ,fill = stage), alpha = 0.5, colour = 'white')+
  scale_x_continuous(expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')), expand= c(0,0))+
  scale_y_continuous('Density', expand= c(0,0))+
  scale_fill_brewer(palette = 'OrRd')+
  geom_vline(xintercept = decomp_vl$mv['mean'],  linetype = 2, colour = '#fdbb84')+
  geom_vline(xintercept = 4.74 ,  linetype = 2, colour = '#fee8c8')+
  geom_vline(xintercept = decomp_vl$sv['mean'], linetype = 2, colour = '#e34a33')+
  facet_wrap(.~stage, nrow = 3,labeller = labeller(stage = facet.labs))+
  my_theme


# Plotting the proportion of tests where a significant result is returned under different 
# conditions of effect size, sample size and p(mv)
plt1b <- ggplot(simsignificances ) +
  geom_raster(aes(x = p_mv, y = effect_size, fill = value))+    
  geom_point(x = 0.32, y = 0.293, shape = 2, colour = 'white') + 
  
  geom_point(x = 0.25, y = 0.372, shape = 5, colour = 'white') + 
  
  geom_vline(xintercept = 0.21, colour = "white", linetype = 2) + 
  annotate("text", label = "MF",
           x = 0.23, y = 0.95, size = 4, colour = "white"
  )+
  geom_vline(xintercept = 0.13, colour = "white", linetype = 2)+
  annotate("text", label = "FM",
           x = 0.15, y = 0.95, size = 4, colour = "white"
  )+
  geom_vline(xintercept = 0.3, colour = "white", linetype = 2) + 
  annotate("text", label = "MM",
           x = 0.32, y = 0.95, size = 4, colour = "white"
  )+
  facet_wrap(.~sample_size, scales = "free_x") + 
  scale_y_continuous('Increase in SpVL due to Multiple Variants', expand= c(0,0)) +
  scale_x_continuous('Cohort P(Multiple Variants)', expand= c(0,0))+
  scale_fill_distiller(palette = 'OrRd', 'P(SpVL ~ P(MV) is Significant)', direction = 1) +
  my_theme 




cowplot::plot_grid(plt1a, plt1b, align = 'hv', nrow = 1, labels = 'AUTO', rel_widths = c(0.4,0.6))

############################################## Panel 2 ##############################################
# Model overview, component models

############################################## Panel 3 ##############################################
# Swiss HIV data, TM applied to SHCS data
shcs_plt1a <- ggplot(shcs_data , aes(x = SpVL_1, SpVL_2)) +
  geom_point( #'#CB6015' #'#66c2a4','#2ca25f','#006d2c'
    colour = '#ef654a',
    shape = 4, size = 3) +
  scale_x_log10(limits = c(10**2, 10**7),
                expand = c(0.05,0),
                name = expression(paste("SpVL Partner 1", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_log10(name = expression(paste("SpVL Partner 2", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  annotation_logticks() +
  my_theme

# Sex
shcs_plt1b <- ggplot(data_model , aes(x = sex, fill = partner)) +
  geom_bar() + 
  scale_y_continuous(expand = c(0,0), name = 'Count') + 
  scale_x_discrete(expand = c(0,0), name = 'Sex') +
  scale_fill_manual(values = c('#ef654a','#fdd49e'), name = 'Partner')+
  my_theme


#Risk Group
shcs_plt1b <- ggplot(data_model , aes(x =  riskgroup)) +
  geom_bar( fill = '#ef654a') + 
  scale_y_continuous(expand = c(0,0), name = 'Count') + 
  scale_x_discrete(expand = c(0,0), name = 'Risk Group') +
  #scale_fill_brewer(palette = 'OrRd', name = 'Sex')+
  my_theme


# Age
shcs_plt1c <-   #scale_fill_manual(values = c('#ef654a','#fdd49e'), name = 'Partner')+
  #scale_fill_manual(values = c('#ef654a','#fdd49e'), name = 'Partner')+
  my_theme

shcs_plt1c <- ggplot(data_model) + 
  geom_histogram(aes(x = age.inf, fill = partner)) + 
  scale_y_continuous(expand = c(0.01,0.01), name = 'Count') + 
  scale_x_continuous(expand = c(0.01,0.01), name = 'Age at Infection',
                     breaks = seq(0,12, by = 1)) +
  scale_fill_manual(values = c('#ef654a','#fdd49e'), name = 'Partner')+
  my_theme



############################################## Panel 4 ##############################################
# Adjsuted for heritability model

############################################## Panel 5 ##############################################
# Results from stratified cohorts

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

combined_data_CD4_PMV_wide <- combined_data_CD4_PMV %>%
  filter(w == 1) %>%
  pivot_wider(names_from = partner, values_from = c(sex, age.inf, riskgroup, SpVL, log10_SpVL, age.inf_category, delta_CD4, starts_with('p_')), names_sep = '_') 

combined_data_CD4_PMV_wide_long <- combined_data_CD4_PMV %>%
  filter(w == 1) %>%
  pivot_wider(names_from = partner, values_from = c(sex, age.inf, riskgroup, SpVL, log10_SpVL, age.inf_category, delta_CD4, starts_with('p_')), names_sep = '_') %>%
  pivot_longer(cols = starts_with('p_'), names_to = 'x', values_to = 'probability') %>%
  separate(x, sep = '_', c(NA, 'type', 'x', 'partner')) #Incorrect output - need partner status separated for transmitter/recipient roles



plt_2a <- ggplot(combined_data_CD4_PMV_wide %>% filter(dataset == 'shcs_empirical'), 
                 aes(x = 1 - p_variants_1_1, 
                     y = SpVL_2 #,
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


plt_2b <- ggplot(combined_data_CD4_PMV_wide %>% filter(dataset %in% c('shcs_empirical',  'shcs_predicted')), #c('stratified_HET',  'stratified_MSM', 'stratified_PWID')
                 aes(x = 1 - p_variants_1_1, 
                     y = delta_CD4_2
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
  my_theme + theme(legend.position = 'none', axis.title.y = element_text(family = 'sans')) +
  facet_wrap(.~dataset)
  #stat_smooth(method = 'lm')


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