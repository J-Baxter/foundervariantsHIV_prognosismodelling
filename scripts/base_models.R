# Run baseline models, 
# Exports plots to file as produced for publication

# Dependencies
require(tidyverse)
require(scales)
require(cowplot)


# Import Theme
source('./scripts/global_theme.R')

# Sanity Check
stopifnot('data' %in% ls(envir = .GlobalEnv))


# Initialise demo population from fitted distributions of Rakkai CC
index <- c('mean'= 4.61, 'sd' = 0.63) 
secondary <- c('mean' = 4.60 , 'sd' = 0.85)
NPAIRS = 100
pop <- InitPop(N = NPAIRS, 
               H2 = 0.33,
               transmitter_SpVL = index, 
               recipient_SpVL = secondary) 


################################### Fit Heritability Model ###################################
cat('Fitting heritablity model...\n')
heritability_model <- lm(recipient_log10SpVL ~  transmitter_log10SpVL, data = pop)

# Scatter plot of transmitter ~ recipient viral loads
plt_1a <- ggplot(pop, aes(x = transmitter, recipient)) +
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


################################### Fit Tolerance Model ###################################
# Fits model as described by Regoes et al. 2014
cat('Fitting tolerance model...\n')
tolerance_model <- lm(CD4.decline ~ SpVL + I(SpVL**2), data = data)

##################################################
# parameters from Regoes et al. PLoS Biology 2014 
##################################################
# a_f = -0.010941 
# n_m = - 0.000243 

plt_1b <- ggplot(data, aes(x = SpVL, y=CD4.decline))+
  geom_point(colour = '#ef654a', shape= 4, alpha = 0.3) +
  scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)),  #
                     expand = c(0,0),
                     limits = c(-2,2))+
  scale_x_continuous(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
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
cat('Simulating transmission model...\n')
baseline_tm <- populationmodel_acrossVL_Environment(sp_ViralLoad = 10**6, w = 1)  %>%
  setNames(nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w'))
 

# Format for marginal probability of the number of variants 
marginal_probs <- baseline_tm %>%
  cbind.data.frame() %>% 
  rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
  group_by(transmitter, probTransmissionPerSexAct, w) %>%
  select(-contains('nparticles')) %>%
  summarise(across(starts_with('V'), sum)) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with('V'),
               names_to = 'variants',
               values_to = 'p') %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% 
           as.numeric())  %>%
  filter(variants <= 10)


# Format for joint probability of the number of variants and virions
joint_probs <- baseline_tm %>%
  cbind.data.frame() %>% 
  rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
  pivot_longer(cols = starts_with('V'),
               names_to = 'variants',
               values_to = 'p') %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% 
           as.numeric()) %>%
  filter(p > 0)

# Plots
plt_1c <-ggplot(marginal_probs)+
  geom_bar(aes(x = variants, y = p), stat = 'identity', fill = '#ef654a') +
  scale_y_continuous(name = 'P(MV)', expand = c(0,0), limits = c(0,0.8), breaks = seq(0,0.8, by = 0.2))+ 
  scale_x_continuous(name = 'Variants', expand = c(0,0), breaks = 1:10 )+
  my_theme

plt_1d <- ggplot(joint_probs, aes(y = variants, x = nparticles))+
  geom_point(aes(size = p), colour = '#ef654a')+
  scale_size(range = c(0,15), name = 'P(Variants \u2229 Virions)')+
  scale_y_continuous(name = 'Variants', expand = c(0,0), limits = c(0.5,10.5), breaks = 1:10)+ 
  scale_x_continuous(name = 'Virions', expand = c(0,0), limits = c(0.5,10.5), breaks = 1:10)+
  my_theme+
  theme(legend.position =  c(0.2,0.85),
        legend.title = element_text(family = 'sans'))

################################### Export Plots ###################################
cat('Preparing outputs...\n')
#ggsave(plot = plt_1a, filename = paste(figs_dir,sep = '/', "plt_1a.eps"), device=cairo_ps, width = 10, height = 10, units= 'in')
#ggsave(plot = plt_1b, filename =  paste(figs_dir,sep = '/',"plt_1b.eps"), device=cairo_ps, width = 10, height = 10, units= 'in')
#ggsave(plot = plt_1c, filename =  paste(figs_dir,sep = '/',"plt_1c.eps"), device=cairo_ps, width = 10, height = 10, units= 'in')
#ggsave(plot = plt_1d, filename =  paste(figs_dir,sep = '/',"plt_1d.eps"), device=cairo_ps, width = 10, height = 10, units= 'in')
 
#ggsave(plot = plt_1a, filename = paste(figs_dir,sep = '/', "plt_1a.jpeg"), device=cairo_ps, width = 10, height = 10, units= 'in')
#ggsave(plot = plt_1b, filename =  paste(figs_dir,sep = '/',"plt_1b.jpeg"), device=cairo_ps, width = 10, height = 10, units= 'in')
#ggsave(plot = plt_1c, filename =  paste(figs_dir,sep = '/',"plt_1c.jpeg"), device=cairo_ps, width = 10, height = 10, units= 'in')
#ggsave(plot = plt_1d, filename =  paste(figs_dir,sep = '/',"plt_1d.jpeg"), device=cairo_ps, width = 10, height = 10, units= 'in')

cat('base_models.R complete.\n')
## END ##
