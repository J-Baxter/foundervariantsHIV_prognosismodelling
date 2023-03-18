# Run baseline models, fitting to data from Swiss HIV Cohort
# Exports plots to file as produced for publication

# Dependencies
require(tidyverse)
require(scales)
require(cowplot)
source('./scripts/populationdata_acrossVL_models.R')

# Import Theme
source('./scripts/global_theme.R')

# Sanity Check
stopifnot('data' %in% ls(envir = .GlobalEnv))



index <- c('mean'= 4.61, 'sd' = 0.63) 
secondary <- c('mean' = 4.60 , 'sd' = 0.85)

pop <- InitPop(N = NPAIRS, 
               H2 = 0.33,
               donor_vl = index, 
               recipient_vl = secondary) %>%
  cbind.data.frame(sex = sample(c('M', 'F'), nrow(.), replace = T)) %>%
  cbind.data.frame(age = sample(18:50, nrow(.), replace = T))



################################### Fit Heritability Model ###################################
h2_model <- lm(recipient_log10SPVL ~ transmitter_log10SPVL + age + transmission + sex, data = data) 

# Sanity check model (convergence, linearity)

plt_1a <- ggplot(data, aes(x = transmitter, recipient)) +
  geom_point( #'#CB6015' #'#66c2a4','#2ca25f','#006d2c'
    colour = '#ef654a',
    shape = 4, size = 3) +
  scale_x_log10(limits = c(1, 10**10),
                expand = c(0,0),
                name = expression(paste("Transmitter SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_log10(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**10),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  #geom_function() + Incorporate lm with fixed effects (with #Add 95% CIs to lines (geom_ribbon))
  annotation_logticks() +
  my_theme + 
  theme(legend.position = 'none')


################################### Fit Tolerance Model ###################################
cd4_model <- lm(CD4.decline ~ spVL + I(spVL**2), data = data)

# Sanity check model (convergence, linearity)

##################################################
# parameters from Regoes et al. PLoS Biology 2014 
##################################################
# a_f = -0.010941 
# n_m = - 0.000243 

plt_1b <- ggplot(cd4_data, aes(x = spVL, y=CD4.decline))+
  geom_point(colour = '#ef654a', shape= 4, alpha = 0.3) +
  scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)),  #
                     expand = c(0,0),
                     limits = c(-2,2))+
  scale_x_continuous(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                     expand = c(0,0))+
  geom_function(fun = function(x) ((-5.6e-3) + (-1.6e-4)*20) * x**2, colour = '#fdbb84', linewidth = 1.2, xlim = c(0,6.3)) +
  geom_text(aes(x = 6.55, y = -0.35), label = "20 yrs", colour = "#fdbb84", size = 6) +
  geom_function(fun = function(x) ((-5.6e-3) + (-1.6e-4)*40) * x**2, colour = '#ef6548', linewidth = 1.2, xlim = c(0,6.3))+
  geom_text(aes(x = 6.55, y = -0.5, label = "40 yrs"), colour = "#ef6548", size = 6) +
  geom_function(fun = function(x) ((-5.6e-3) + (-1.6e-4)*60) * x**2, colour = '#b30000', linewidth = 1.2, xlim = c(0,6.3)) +
  geom_text(aes(x = 6.55, y = -0.62, label = "60 yrs"), colour = "#b30000", size = 6) + #Add 95% CIs to lines (geom_ribbon)
  my_theme + theme(axis.title.y = element_text(family = 'sans')) + 
  annotation_logticks(sides = 'b') 


################################### Run Transmission ###################################t



base_tm_sims <- c(RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 1),
          RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 5),
          RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 10),
          RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 15),
          RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 20)) %>%
  
  # Label
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w'))

variant_baseline <- base_tm_sims %>%
  # Post-process into dataframe for P(V|Transmitter), segregated by weighting
  lapply(., cbind.data.frame) %>%
  do.call(rbind.data.frame,.) %>%
  rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
  group_by(transmitter, probTransmissionPerSexAct, w) %>%
  select(-contains('nparticles')) %>%
  summarise(across(starts_with('V'), sum)) %>%
  
  # Long format
  pivot_longer(cols = starts_with('V'),
               names_to = 'variants',
               values_to = 'p') %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric()) 



plt_1c <- ggplot(variant_baseline %>% 
                   filter(variants == 1) , 
                 aes(x = transmitter, 
                     y = 1 - p,
                     colour = as.factor(w)))+
  geom_point(shape = 3, size = 4) +
  scale_colour_brewer(palette = 'OrRd') +
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

################################### Export Plots ###################################

panel_1 <- plot_grid(plt_1a, plt_1b, plt_1c, align = 'hv', nrow = 1, labels = 'AUTO')

ggsave(plot = panel_1, filename = paste(figs_dir,sep = '/', "panel_1.eps"), device=cairo_ps, width = 10, height = 10, units= 'in')

## END ##
