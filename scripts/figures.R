# Figures

# Dependencies
source('./scripts/dependencies.R')
source('./scripts/global_theme.R')
#set_null_device(cairo_pdf)

############################################## Import Results ############################################## 

cd4_data <- read_table('./data/pbio.1001951.s006.tsv') %>%
  rename(SpVL = spVL)

resultsfiles <- list.files('./results/07Aug23/', full.names = T)

results <- lapply(resultsfiles, read_csv)


############################################## Panel 1 ##############################################
# From observe sig effect

# Normal distributions of recipient phenotype, multiple variant and single variant
facet.labs <- c("Recipient Population", "Multiple Variant Recipients", "Single Variant Recipients") %>%
  setNames(c('recipient_mean', 'recipient_mv', 'recipient_sv'))

plt1a <- ggplot(vls) + 
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
plt1b <- ggplot(simsignificances) +
  geom_raster(aes(x = p_mv, y = effect_size, fill = value))+    
  geom_point(x = 0.32, y = 0.293, shape = 2, colour = 'white') + 
  
  geom_point(x = 0.25, y = 0.372, shape = 5, colour = 'white') + 
  
  geom_vline(xintercept = 0.21, colour = "white", linetype = 2) + 
  annotate("text", label = "MF",
           x = 0.23, y = 0.95, size = 4, colour = "white")+
  
  geom_vline(xintercept = 0.13, colour = "white", linetype = 2)+
  annotate("text", label = "FM",
           x = 0.15, y = 0.95, size = 4, colour = "white")+
  
  geom_vline(xintercept = 0.3, colour = "white", linetype = 2) + 
  annotate("text", label = "MM",
           x = 0.32, y = 0.95, size = 4, colour = "white")+
  
  facet_wrap(.~sample_size, scales = "free_x") + 
  scale_y_continuous('Increase in SpVL due to Multiple Variants', expand= c(0,0)) +
  scale_x_continuous('Cohort P(Multiple Variants)', expand= c(0,0))+
  scale_fill_distiller(palette = 'OrRd', 'P(SpVL ~ P(MV) is Significant)', direction = 1) +
  my_theme 

plt_1 <- cowplot::plot_grid(plt1a, plt1b, align = 'hv', nrow = 1, labels = 'AUTO', rel_widths = c(0.4,0.6))

ggsave(plot = plt_1, filename = paste(figs_dir,sep = '/', "plt_s1.jpeg"), 
       device = jpeg, width = 170, height = 170, units = 'mm')


############################################## Panel 2 ##############################################
# Model overview, component models and SHCS data

# Transmission Model 
plt_2b_probabilities <- TransmissionModel2(sp_ViralLoad = 10**6, 
                                  PerVirionProbability = 8.779E-07, 
                                  PropExposuresInfective = 0.14337)  %>%
  setNames(nm = c('variant_distribution','probTransmissionPerSexAct','transmitter')) %>%
  .[['variant_distribution']] %>%
  cbind.data.frame() %>% 
  pivot_longer(cols = starts_with('V'),
               names_to = 'variants',
               values_to = 'p') %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% 
           as.numeric()) %>%
  filter(p > 0)

plt_2b <- ggplot(plt_2b_probabilities, aes(y = variants, x = nparticles))+
  geom_point(aes(size = p), colour = '#ef654a')+
  scale_size(range = c(0,5), name = 'P(Variants \u2229 Virions)')+
  scale_y_continuous(name = 'Variants', expand = c(0,0), limits = c(0.5,10.5), breaks = 1:10)+ 
  scale_x_continuous(name = 'Virions', expand = c(0.01,0), limits = c(0.5,10.5), breaks = 1:10)+
  my_theme+
  theme(legend.position =  c(0.5,0.8),
        legend.background = element_rect(fill = NA))


# Heritability Model (Frequentist vis)
freq_model <- lme4::lmer(log10_SpVL ~  + partner + sex + age.inf_category + riskgroup + (1|log10_SpVL_couplemean), 
                         data = shcs_data_long_transmitterML)

plt_2c <- mutate(shcs_data_long_transmitterML,
                 .pred = predict(freq_model)) %>%
  data.table::data.table() %>%
  mltools::one_hot(cols = c('age.inf_category', 'partner', 'sex', 'riskgroup')) %>%
  as_tibble() %>%
  select(contains(c('partner_recipient', 'sex_F', 'riskgroup', 'log10_SpVL_couplemean', 'age.inf_category' , "ID_pair", '.pred'))) %>%
  select(-c('age.inf_category_15-24', 'riskgroup_HET')) %>%
  mutate(partner_recipient = case_when(partner_recipient == 1 ~ -0.15643,
                                       .default = 0)) %>%
  mutate(sex_F = case_when(sex_F == 1 ~ -0.0351,
                           .default = 0)) %>%
  mutate(`age.inf_category_25-29` = case_when(`age.inf_category_25-29` == 1 ~  0.08640,
                                              .default = 0)) %>%
  mutate(`age.inf_category_30-39` = case_when(`age.inf_category_30-39` == 1 ~ 0.12095,
                                              .default = 0)) %>%
  mutate(`age.inf_category_40-80` = case_when(`age.inf_category_40-80` == 1 ~  0.19728,
                                              .default = 0)) %>%
  mutate(riskgroup_MSM = case_when(riskgroup_MSM == 1 ~ 0.01874,
                                   .default = 0)) %>%
  mutate(riskgroup_PWID = case_when(riskgroup_PWID == 1 ~ 0.01923,
                                    .default = 0)) %>%
  mutate(riskgroup_OTHER = case_when(riskgroup_OTHER == 1 ~ -0.03739,
                                     .default = 0)) %>%
  mutate(riskgroup_UNKNOWN = case_when(riskgroup_UNKNOWN == 1 ~ -0.23424,
                                       .default = 0)) %>%
  mutate(slope = rowSums(select(., -c(ID_pair, log10_SpVL_couplemean, .pred)))) %>%
  arrange(log10_SpVL_couplemean) %>%
  mutate(intercept = rep(coef(freq_model)$log10_SpVL_couplemean[,1], each = 2)) %>%
  ggplot() +
  geom_abline(aes(intercept =intercept, slope = slope), alpha = 0.2, colour = '#fdbb84') +
  geom_point(aes(x = log10_SpVL_couplemean, y = log10_SpVL), 
             data = shcs_data_long_transmitterML, 
             colour = '#ef654a', shape = 4, alpha = 0.4, size = 1) + 
  scale_x_continuous(name = expression(paste(SpVL[ij], ' (', Log[10], " copies ", ml**-1, ')')),
                     limits = c(1, 7),
                     expand = c(0.05,0)) + 
  scale_y_continuous(name = expression(paste("Mean Pair SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                     limits = c(1, 7),
                     expand = c(0.05,0))+ 
  my_theme


# Tolerance Model
plt_2d <- ggplot(cd4_data , aes(x = SpVL, y=CD4.decline))+
  geom_point(colour = '#ef654a', shape= 4, alpha = 0.4, size = 1) +
  scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)),  #
                     expand = c(0,0),
                     limits = c(-2,2))+
  scale_x_continuous(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                     expand = c(0,0),
                     limits = c(0,7.5),
                     breaks = seq(0,6, by = 2))+
  geom_function(fun = function(x) ((-5.6e-3) + (-1.6e-4)*20) * x**2, colour = '#fdbb84', linewidth = 1.2, xlim = c(0,6.3)) +
  geom_text(aes(x = 6.9, y = -0.35), label = "20 yrs", colour = "#fdbb84", size = 2) +
  geom_function(fun = function(x) ((-5.6e-3) + (-1.6e-4)*40) * x**2, colour = '#ef6548', linewidth = 1.2, xlim = c(0,6.3))+
  geom_text(aes(x = 6.9, y = -0.5, label = "40 yrs"), colour = "#ef6548", size = 2) +
  geom_function(fun = function(x) ((-5.6e-3) + (-1.6e-4)*60) * x**2, colour = '#b30000', linewidth = 1.2, xlim = c(0,6.3)) +
  geom_text(aes(x = 6.9, y = -0.62, label = "60 yrs"), colour = "#b30000", size = 2) + #Add 95% CIs to lines (geom_ribbon)
  my_theme

# SHCS
plt_2e <- ggplot(shcs_data , aes(x = SpVL_1, SpVL_2)) +
  geom_point( #'#CB6015' #'#66c2a4','#2ca25f','#006d2c'
    colour = '#ef654a',
    shape = 4, alpha = 0.4, size = 1) +
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
  my_theme


#Risk Group
plt_2f <- ggplot(shcs_data_long_transmitterML , aes(x =  riskgroup)) +
  geom_bar( fill = '#ef654a') + 
  scale_y_continuous(expand = c(0,0), name = 'Count') + 
  scale_x_discrete(expand = c(0,0), name = 'Risk Group', labels = c('HET', 'MSM', 'PWID', 'OTHER', 'NA')) +
  #scale_fill_brewer(palette = 'OrRd', name = 'Sex')+
  my_theme


# Age
plt_2g <- ggplot(shcs_data_long_transmitterML) + 
  geom_histogram(aes(x = age.inf), fill = '#ef654a') + 
  scale_y_continuous(expand = c(0.01,0.01), name = 'Count') + 
  scale_x_continuous('Age at Infection', limits = c(15,80), breaks = seq(16,80, by = 8), expand = c(0,0))+
  
  my_theme

panel_2 <- plot_grid(plt_2b, plt_2c, plt_2d, plt_2e, plt_2f, plt_2g,  ncol = 3, align = 'hv', label_size = 9, labels = 'AUTO')

ggsave(plot = panel_2 , filename = paste(figs_dir,sep = '/', "panel_2.jpeg"), 
       device = jpeg,  width = 170, height = 140,  units = 'mm')


############################################## Panel 3 ##############################################
shcs_empirical_results <- results[[1]]
shcs_predicted_results <- results[[2]]

plt_3a <- ggplot(shcs_empirical_results %>% filter(transmitterallocation == 'ML'))+
  geom_point(aes(y = SpVL_recipient, x = 1-p_variants_1, colour = riskgroup_recipient),size = 1, shape = 4)+
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_color_brewer(palette = 'OrRd', 'Recipient Riskgroup') +
  my_theme+
  theme(axis.title = element_text(size = 7))

plt_3b <- ggplot(shcs_empirical_results %>% filter(transmitterallocation == 'ML'))+
  geom_point(aes(y = SpVL_recipient, x = 1-p_particles_1, colour = riskgroup_recipient),size = 1, shape = 4)+
  scale_x_continuous(name = 'P(Multiple Particle Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_color_brewer(palette = 'OrRd') +
  my_theme+
  theme(axis.title = element_text(size = 7))

plt_3c <- ggplot(shcs_empirical_results %>% filter(transmitterallocation == 'ML'))+
  geom_point(aes(y = delta_CD4_recipient, x = 1-p_variants_1, colour = riskgroup_recipient),size = 1, shape = 4)+
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)),  #
                     expand = c(0,0),
                     limits = c(-1,0)) +
  scale_color_brewer(palette = 'OrRd') +
  my_theme +
  theme(axis.title = element_text(size = 7))

plt_3d <- ggplot(shcs_empirical_results %>% filter(transmitterallocation == 'ML'))+
  geom_point(aes(y = delta_CD4_recipient, x = 1-p_particles_1, colour = riskgroup_recipient),size = 1, shape = 4)+
  scale_x_continuous(name = 'P(Multiple Particle Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)),  #
                     expand = c(0,0),
                     limits = c(-1,0)) +
  scale_color_brewer(palette = 'OrRd') +
  my_theme +
  theme(axis.title = element_text(size = 7))


plt_3e <- ggplot(shcs_predicted_results %>% filter(transmitterallocation == 'ML'), 
                 aes(y = SpVL_recipient, x = 1-p_variants_1))+
  geom_density_2d_filled()+
  #scale_fill_brewer(palette = 'OrRd', direction = 1)+
  scale_fill_manual(values = my_palette) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**1, 10**7.5),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  my_theme+
  theme(axis.title = element_text(size = 7))

plt_3f <- ggplot(shcs_predicted_results %>% filter(transmitterallocation == 'ML'),
                 aes(y = SpVL_recipient, x = 1-p_particles_1))+
  geom_density_2d_filled()+
  #scale_fill_brewer(palette = 'OrRd', direction = 1)+
  scale_fill_manual(values = my_palette) +
  scale_x_continuous(name = 'P(Multiple Particle Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**1, 10**7.5),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  my_theme +
  theme(axis.title = element_text(size = 7))


plt_3g <- ggplot(shcs_predicted_results %>% filter(transmitterallocation == 'ML'),
                 aes(y = delta_CD4_recipient, x = 1-p_variants_1))+
  geom_density_2d_filled()+
  #scale_fill_brewer(palette = 'OrRd', direction = 1)+
  scale_fill_manual(values = my_palette) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_continuous(limits = c(-0.5,0), 
                     expand = c(0,0.01), 
                     breaks = seq(-0.5, 0, by = 0.1),
                     name =expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1))) +
  my_theme+
  theme(axis.title = element_text(size = 7))


plt_3h <- ggplot(shcs_predicted_results %>% filter(transmitterallocation == 'ML'),
                 aes(y = delta_CD4_recipient, x = 1-p_particles_1))+
  geom_density_2d_filled()+
  #scale_fill_brewer(palette = 'OrRd', direction = 1)+
  scale_fill_manual(values = my_palette) +
  scale_x_continuous(name = 'P(Multiple Particle Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_continuous(limits = c(-0.5,0), 
                     expand = c(0,0.01), 
                     breaks = seq(-0.5, 0, by = 0.1),
                     name =expression(paste(Delta, CD4,'+ ' , mu, l**-1, ' ', day**-1))) +
  my_theme+
  theme(axis.title = element_text(size = 7))

panel3_top <- plot_grid(plt_3a, plt_3b, plt_3c,plt_3d, nrow= 1, labels= 'AUTO', label_size = 9, align = 'hv')
panel3_legend <- get_legend(plt_3a + theme(legend.position = 'bottom')) 
panel3_bottom <- plot_grid(plt_3e, plt_3f, plt_3g, plt_3h, nrow= 1, labels= c('E', 'F', 'G', 'H'), label_size = 9, align = 'hv')

panel_3 <- plot_grid(panel3_legend, panel3_top,  panel3_bottom, nrow= 3, labels= c('', '', '','E', 'F', 'G', 'H'), label_size = 9, rel_heights =  c(0.1,1,1))

ggsave(plot = panel_3 , filename = paste(figs_dir,sep = '/', "panel_3.jpeg"), 
       device = jpeg,  width = 200, height = 130,  units = 'mm')


############################################## Panel 4 ##############################################

sim_results <-  bind_rows(results[3:5]) %>% filter(transmitterallocation == 'ML')
my_palette <- colorRampPalette(brewer.pal(9, "OrRd"))(13)

plt_4a <- ggplot(sim_results,aes(y = SpVL_recipient, x = 1-p_variants_1))+
  geom_density_2d_filled()+
  #scale_fill_brewer(palette = 'OrRd', direction = 1)+
  scale_fill_manual(values = my_palette) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**1, 10**7.5),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  facet_grid(cols= vars(riskgroup_recipient), switch = 'y')+
  my_theme


plt_4b <- ggplot(sim_results,aes(y = SpVL_recipient, x = 1-p_particles_1))+
  geom_density_2d_filled()+
  #scale_fill_brewer(palette = 'OrRd', direction = 1)+
  scale_fill_manual(values = my_palette) +
  scale_x_continuous(name = 'P(Multiple Particle Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**1, 10**7.5),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  facet_grid(cols = vars(riskgroup_recipient), switch = 'y')+
  my_theme


plt_4c <- ggplot(sim_results,aes(y = delta_CD4_recipient, x = 1-p_variants_1))+
  geom_density_2d_filled()+
  #scale_fill_brewer(palette = 'OrRd', direction = 1)+
  scale_fill_manual(values = my_palette) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_continuous(limits = c(-0.5,0), 
                     expand = c(0,0.01), 
                     breaks = seq(-0.5, 0, by = 0.1),
                     name =expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1))) +
  facet_grid(cols = vars(riskgroup_recipient), switch = 'y')+
  my_theme


plt_4d <- ggplot(sim_results,aes(y = delta_CD4_recipient, x = 1-p_particles_1))+
  geom_density_2d_filled()+
  #scale_fill_brewer(palette = 'OrRd', direction = 1)+
  scale_fill_manual(values = my_palette) +
  scale_x_continuous(name = 'P(Multiple Particle Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_continuous(limits = c(-0.5,0), 
                     expand = c(0,0.01), 
                     breaks = seq(-0.5, 0, by = 0.1),
                     name =expression(paste(Delta, CD4,'+ ' , mu, l**-1, ' ', day**-1))) +
  facet_grid(cols = vars(riskgroup_recipient), switch = 'y')+
  my_theme

panel_4 <- cowplot::plot_grid(plt_4a, plt_4b, plt_4c, plt_4d, align = 'hv', nrow = 4, labels = 'AUTO', label_size = 9)

ggsave(plot = panel_4 , filename = paste(figs_dir,sep = '/', "panel_4.jpeg"), 
       device = jpeg,  width = 170, height = 250,  units = 'mm')
cat('All plots complete. \n')