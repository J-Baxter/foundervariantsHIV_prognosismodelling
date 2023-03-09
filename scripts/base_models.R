# Run baseline models, fitting to data from Swiss HIV Cohort
# Exports plots to file as produced for publication

# Dependencies
require(tidyverse)
require(scales)
require(cowplot)
source('./scripts/populationmodel_acrossVL_Environment.R')

# Import Theme
source('./scripts/plot_theme.R')

# Sanity Check
stopifnot('data' %in% ls(envir = .GlobalEnv))

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


################################### Run Transmission ###################################
transmitterspvl_variantdist <- RunParallel(populationmodel_acrossVL_Environment, pop$transmitter) %>%
  
  # Post-Processing 
  do.call(cbind.data.frame, .) %>% 
  t() %>%
  cbind.data.frame(transmitter = pop$transmitter, .) %>%
  rename(probTransmissionPerSexAct_average = probTransmissionPerSexAct) %>%
  rename_with(function(x) gsub('\\.(?<=\\V)', '\\1_average.\\2',x, perl= TRUE), 
              .cols = !contains('chronic') & contains('V')) %>%
  
  # Long Format
  pivot_longer(cols = c(contains('chronic'),contains('average')), 
               names_to = c('type', 'stage', 'variants'), 
               names_pattern = "(^[^_]+)(_[^.]*)(.*)", 
               values_to = 'p') %>% 
  mutate(stage = gsub('^.*_', '', stage)) %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric()) 

plt_1c <- ggplot(transmitterspvl_variantdist %>% 
                   filter(variants == 1) , 
                 aes(x = transmitter, 
                     y = 1 - p,
                     colour = stage))+
  geom_point(shape = 3, size = 4) +
  scale_colour_brewer(palette = 'RdOr') +
  scale_x_log10(name = expression(paste("Transmitter SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.6),
                     breaks = seq(0, 0.6, by = 0.1)) +
  annotation_logticks(sides = 'b') +
  my_theme + theme(legend.position = 'none')

################################### Export Plots ###################################

panel_1 <- plot_grid(plt_1a, plt_1b, plt_1c, align = 'hv', nrow = 1, labels = 'AUTO')

ggsave(plot = panel_1, filename = paste(figs_dir,sep = '/', "panel_1.eps"), device=cairo_ps, width = 10, height = 10, units= 'in')

## END ##
