# Figures

################################### Dependencies ###################################
source('./scripts/dependencies.R')
source('./scripts/functions/BonhoefferEqns.R')
source('./scripts/functions/AllocateTransmitter.R')
source('./scripts/functions/simulate_cohorts_funcs.R')
source('./scripts/functions/PostProcessing.R')
source('./scripts/interpret_outputs.R')
source('./scripts/functions/populationdata_acrossVL_models_R.r')
#source('./scripts/global_theme.R')
library(tidyverse)
library(scales)
library(extrafont)




#https://www.fontsquirrel.com/fonts/latin-modern-sans
font_import(pattern = "lmsans10*") 
loadfonts()

my_theme <- theme_classic(base_family = "LM Sans 10")+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 8),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    strip.text  = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank()
  )

################################### Import Data ###################################
# Paired set point viral load (SpVL), age, riskgroup and sex data provided on request
# by the Swiss HIV Cohort Study. These data were used in the analysis of Bertels et. al
# 2018 https://academic.oup.com/mbe/article/35/1/27/4210012

# init distribtions for transmitter allocation
# mc and vc are the mean and variance of the Amsterdam seroconverter study, respectively.
mc <- 4.35
vc <- 0.47**2

t <- CalcTransmitter(mc,vc)
r <- CalcRecipient(mc,vc)

# Import data and pre-process
source('./scripts/import_data.R')


############################################## Import Results ############################################## 

cd4_data <- read_table('./data/pbio.1001951.s006.tsv') %>%
  rename(SpVL = spVL)

resultsfiles <- list.files('./results/21Feb24', full.names = T)

modelresults <- lapply(resultsfiles[which(grepl('rawresults',resultsfiles) & !grepl('UNKNOWN|OTHER',resultsfiles))], read_csv) %>%
  do.call(rbind.data.frame,.)

cd4results <- lapply(resultsfiles[which(grepl('ML_CD4resample',resultsfiles))], read_csv) %>% 
  setNames(., c("FM", "MF", "MMI", "MMR", "PWID")) %>%
  bind_rows(. , .id = 'riskgroup')

SpVLresults <- lapply(resultsfiles[which(grepl('ML_SpVLresample',resultsfiles))], read_csv) %>%
  setNames(., c("FM", "MF", "MMI", "MMR", "PWID")) %>%
  bind_rows(. , .id = 'riskgroup')

############################################## Panel 1 ##############################################
# From observe sig effect

# Normal distributions of recipient phenotype, multiple variant and single variant
facet.labs <- c("Recipient Population", "Multiple Variant Recipients", "Single Variant Recipients") %>%
  setNames(c('recipient_mean', 'recipient_mv', 'recipient_sv'))

plt1a <- ggplot(vls) + 
  geom_density(aes(x = vl ,fill = stage), alpha = 0.5, colour = 'white')+
  scale_x_continuous(expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')), expand= c(0,0))+
  scale_y_continuous('Density', expand= c(0,0))+
  scale_fill_brewer(palette = 'GnBu')+
  geom_vline(xintercept = decomp_vl$mv['mean'],  linetype = 2, colour = '#ccebc5')+
  geom_vline(xintercept = 4.74 ,  linetype = 2, colour = '#7bccc4')+
  geom_vline(xintercept = decomp_vl$sv['mean'], linetype = 2, colour = '#0868ac')+
  facet_wrap(.~stage, nrow = 3,labeller = labeller(stage = facet.labs))+
  my_theme


# Plotting the proportion of tests where a significant result is returned under different 
# conditions of effect size, sample size and p(mv)
studysizes <-c('Sagar et al.', 'Janes et al. RV144', 'Janes et al. STEP') 
names(studysizes) <- c( '156', '100', '63')
studysizes <- sort(studysizes, decreasing = T)
study_effects <- study_effects %>% rename(sample_size = size)

plt1b <- ggplot(simsignificances %>% filter(endpoint =='SpVL')) +
  geom_raster(aes(x = p_mv, 
                  y = effect_size, 
                  fill = value),
              interpolate = T)+    
  geom_point(aes(x = p_mv, y = es_spvl), size = 1.5,colour = 'white', data = study_effects) + 
  
  #geom_point(x = 0.25, y = 0.372, shape = 2, size = 1.5,colour = 'white') + 
  
  #geom_point(x = 0.57, y = 0.27, shape = 5, size = 1.5,colour = 'white') + 
  
  #geom_vline(xintercept = 0.21,
  #    colour = "white",
  #   linetype = 2) + 
  #annotate("text", label = "MF",
  # x = 0.24, 
  # y = 0.95, 
# size = 2,
#colour = "white")+

#geom_vline(xintercept = 0.13, 
#colour = "white", 
#linetype = 2)+

#annotate("text", label = "FM",
#x = 0.16, 
# y = 0.95, 
# size = 2, 
#colour = "white")+

#geom_vline(xintercept = 0.3, 
# colour = "white", 
# linetype = 2) + 

# annotate("text", label = "MM",
# x = 0.34,
# y = 0.95,
#size = 2, 
#colour = "white")+

facet_grid(#rows = vars(endpoint),
  rows = vars(sample_size), 
  scales = "free_x",
  labeller = labeller(sample_size = studysizes)) + 
  
  #scale_shape(solid = FALSE, 
  #name = 'Cohort',
  # labels = c('Janes-RV144',
  #'Sagar-Mombasa',
  # 'Janes-STEP'),
  #guide = guide_legend(override.aes = list(colour = 'darkgrey'),
  # direction = "horizontal"))+ 
  
  scale_y_continuous(expression(paste("Increase in ", Log[10], " SpVL due to Multiple Variants")), 
                     expand= c(0,0)) +
  
  scale_x_continuous('Study P(Multiple Variants)', 
                     expand= c(0,0))+
  
  scale_fill_distiller(palette = 'GnBu',
                       'Probability of Observing a True Effect',
                       direction = 1,
                       guide = guide_colourbar(
                         title.position="top",
                         barwidth = unit(40, 'mm'),
                         title.hjust = 0.5,
                         direction = "horizontal"
                       )) +
  my_theme +
  theme(legend.position = 'bottom',
        legend.box="vertical")

#plt_1 <- cowplot::plot_grid(plt1a, plt1b, align = 'hv', nrow = 1, labels = 'AUTO', rel_widths = c(0.4,0.6))


ggsave("figure1.eps", device=cairo_ps,  height = 180, width = 70, units = 'mm')
Sys.sleep(0.5)
plt1b
dev.off()

############################################## Panel 2 ##############################################
# Model overview, component models and SHCS data

require(ggdag)
require(dagitty)


coords <- list(
  x = c(SpVLr = 8, SpVLt = 3, deltaCD4 = 8, PMV = 3),
  y = c(SpVLr = 8, SpVLt = 8, deltaCD4 = 2, PMV = 2)) %>% 
  coords2df() %>%
  coords2list()

test_dag <- dagify(
  deltaCD4 ~ SpVLr,
  SpVLr ~ SpVLt,
  PMV ~ SpVLt )

dagitty::coordinates(test_dag) <- coords


plt_2a <- test_dag %>% 
  ggplot(., aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(size = 14, col = "white", fill = "white") +
  geom_dag_edges() +
  geom_dag_text(col = "black",  size = 3, 
                label = c('Probability of \n Multiple Variant Infection',
                          'Recipient \nSpVL',
                          'Transmitter \nSpVL',
                          'CD4 T Cell\n decline')) +
  annotate("text", x = 5.5, y = 8.3, label = "Heritability", family ="LM Sans 10", size = 3) + 
  annotate("text", x = 2.5, y = 5, label = "Transmission", angle = 90, family ="LM Sans 10", size = 3) + 
  annotate("text", x = 8.5, y = 5, label = "Tolerance", angle = 270, family ="LM Sans 10", size = 3) + 
  scale_x_continuous(limits = c(1.5,9.5), expand = c(0.01,0.01)) + 
  scale_y_continuous(limits = c(1,9), expand = c(0.01,0.01)) + 
  theme_classic(base_family = "LM Sans 10") +
  theme(
    text = element_text(size=10),
    #axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    #axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 8),
    #axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    strip.text  = element_text(size = 8),
    #legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank(),
    axis.line = element_blank(), 
    axis.ticks = element_blank(), 
    axis.text = element_blank(),
    axis.title = element_blank(),
    #axis.ticks.length = unit(0, "lines"), # Error 
    axis.ticks.margin = unit(c(0,0,0,0), "lines"), 
    legend.position = "none", 
    #panel.background = element_rect(fill = "gray"), 
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.margin = unit(c(0,0,0,0), "lines"), 
    #plot.background = element_rect(fill = "blue"),
    plot.margin = unit(c(0,0,0,0), "lines")
  )

ggsave('figure2.eps', device=cairo_ps,  height = 110, width = 140, units = 'mm')
Sys.sleep(0.5)
plt_2a 
dev.off()

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
  geom_abline(aes(intercept =intercept, slope = slope), alpha = 0.2, colour = '#7bccc4') +
  geom_point(aes(x = log10_SpVL_couplemean, y = log10_SpVL), 
             data = shcs_data_long_transmitterML, 
             colour = '#084081', shape = 4, alpha = 0.4, size = 1) + 
  scale_y_continuous(name = expression(paste(SpVL[ij], ' (', Log[10], " copies ", ml**-1, ')')),
                     limits = c(1, 7),
                     expand = c(0.05,0)) + 
  scale_x_continuous(name = expression(paste("Mean Pair SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                     limits = c(1, 7),
                     expand = c(0.05,0))+ 
  my_theme


# Tolerance Model
plt_2d <- ggplot(cd4_data , aes(x = SpVL, y=CD4.decline))+
  geom_point(colour = '#084081', shape= 4, alpha = 0.4, size = 1) +
  scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)),  #
                     expand = c(0,0),
                     limits = c(-2,2))+
  scale_x_continuous(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                     expand = c(0,0),
                     limits = c(0,7.5),
                     breaks = seq(0,6, by = 2))+
  geom_function(fun = function(x) ((-5.6e-3) + (-1.6e-4)*20) * x**2, colour = '#ccebc5', linewidth = 1.2, xlim = c(0,6.3)) +
  geom_text(aes(x = 6.9, y = -0.35), label = "20 yrs", colour = "#ccebc5", size = 2) +
  geom_function(fun = function(x) ((-5.6e-3) + (-1.6e-4)*40) * x**2, colour = '#7bccc4', linewidth = 1.2, xlim = c(0,6.3))+
  geom_text(aes(x = 6.9, y = -0.5, label = "40 yrs"), colour = "#7bccc4", size = 2) +
  geom_function(fun = function(x) ((-5.6e-3) + (-1.6e-4)*60) * x**2, colour = '#2b8cbe', linewidth = 1.2, xlim = c(0,6.3)) +
  geom_text(aes(x = 6.9, y = -0.62, label = "60 yrs"), colour = "#2b8cbe", size = 2) + #Add 95% CIs to lines (geom_ribbon)
  my_theme


# Transmission Model ~ Risk Group
mf <- sapply(raise_to_power(10, 2:7),
             populationmodel_acrossVL_Env, 
             PerVirionProbability = 1.765e-06,
             PropExposuresInfective = 0.013296)

fm <- sapply(raise_to_power(10, 2:7), 
             populationmodel_acrossVL_Env, 
             PerVirionProbability = 8.779e-07, 
             PropExposuresInfective = 0.14337)

mm_ra <- sapply(raise_to_power(10, 2:7),
                populationmodel_acrossVL_Env, 
                PerVirionProbability = 3.190e-06,
                PropExposuresInfective = 0.08923)

mm_ia <- sapply(raise_to_power(10, 2:7),
                populationmodel_acrossVL_Env,
                PerVirionProbability = 3.114e-06, 
                PropExposuresInfective = 0.008839)

pwid <-  sapply(raise_to_power(10, 2:7), 
                populationmodel_acrossVL_Env,
                PerVirionProbability = 5.566e-06, 
                PropExposuresInfective = 0.02496)

av_pmv <- rbind.data.frame(populationmodel_acrossVL_Env(PerVirionProbability = 1.765e-06, PropExposuresInfective = 0.013296), 
                           populationmodel_acrossVL_Env(PerVirionProbability = 8.779e-07, PropExposuresInfective = 0.14337),
                           populationmodel_acrossVL_Env(PerVirionProbability = 3.190e-06, PropExposuresInfective = 0.08923),
                           populationmodel_acrossVL_Env(PerVirionProbability = 3.114e-06, PropExposuresInfective = 0.008839),
                           populationmodel_acrossVL_Env(PerVirionProbability = 5.566e-06, PropExposuresInfective = 0.02496))

av_pmv$Riskgroup <- c('MF', 'FM', 'MM:RA', 'MM:IA', 'PWID')
colnames(av_pmv)[1] <- c('multiple_founder_proportion')

df <- bind_rows(list('MF' = as_tibble(cbind(mf,  rownames(mf))),
                     'FM' =  as_tibble(cbind(fm, rownames(mf))),
                     'MM:RA' =  as_tibble(cbind(mm_ra, rownames(mf))),
                     'MM:IA' =  as_tibble(cbind(mm_ia, rownames(mf))),
                     'PWID' =  as_tibble(cbind(pwid,  rownames(mf)))), .id = 'Riskgroup') %>%
  setNames(c('Riskgroup',raise_to_power(10, 2:7) %>% as.character(), 'time')) %>%
  mutate(time = as.character(time)) %>%
  
  left_join(., pivot_longer(av_pmv, cols = contains('multiple'), values_to = 'average', names_to = 'time'), by = join_by(Riskgroup, time)) %>%
  pivot_longer(cols = contains(c('1')), values_to = 'P_MV', names_to = 'SpVL') %>%
  mutate(time = factor(time, levels = c('multiple_founder_proportion',
                                        'multiple_founder_proportion_primary',
                                        'multiple_founder_proportion_chronic',
                                        'multiple_founder_proportion_preaids')))

plt_2e <- ggplot(df %>% filter(time == 'multiple_founder_proportion')) +
  geom_bar(aes(x = as.numeric(SpVL), y = as.numeric(P_MV), fill = as.factor(SpVL)), stat = 'identity') +
  geom_hline(aes(yintercept = average), linetype = 'dashed')+
  my_theme + 
  scale_x_log10(expand = c(0,0), expression(paste("Transmitter SpVL", ' (', Log[10], " copies ", ml**-1, ')')),  
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", label_math(.x))) +
  scale_y_continuous(expand = c(0,0), 'P(Multiple Variants)', limits = c(0,1), breaks = seq(0,1 ,by = 0.2))+
  scale_fill_brewer(palette = 'GnBu')+
  facet_grid(cols = vars(Riskgroup), switch = 'y',
             labeller = as_labeller( c('MF' = 'Male to Female',
                                       'PWID' = 'PWID',
                                       'FM' = 'Female to Male',
                                       'MM:IA' = 'MSM Insertive',
                                       'MM:RA' = 'MSM Receptive'))) + 
  theme(strip.placement = 'outside')


panel_2top <- plot_grid(plt_2c, plt_2d, align = 'hv', label_size = 9, labels = 'AUTO')

panel_2 <- plot_grid(panel_2top, plt_2e,   nrow = 2, label_size = 9, labels = c('', 'C'))

ggsave('figure3.eps', device=cairo_ps,  height = 140, width = 160, units = 'mm')
Sys.sleep(0.5)
panel_2
dev.off()


############################################## Panel 3 ##############################################

#my_palette <- colorRampPalette(brewer.pal(9, "OrRd"))(20)


plt_3a <- ggplot(modelresults %>% filter(transmitterallocation == 'ML'),aes(y = SpVL_recipient, x = 1-p_variants_1))+
  #stat_density_2d(aes(fill = (..density..)**(1/3)), geom = "raster", contour = FALSE) +
  geom_hex(aes(fill = (..density..)), binwidth = c(0.05, 0.2)) + 
  scale_fill_distiller(palette = 'GnBu',
                       trans = 'sqrt',
                       direction = -1) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,1),
                     breaks = seq(0, 1, by = 0.25))+
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**1, 10**7),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  facet_grid(cols= vars(dataset_id), switch = 'y',
             labeller = as_labeller( c('MF' = 'Male to Female',
                                       'PWID' = 'PWID',
                                       'FM' = 'Female to Male',
                                       'MMI' = 'MSM Insertive',
                                       'MMR' = 'MSM Receptive')))+
  my_theme 


plt_3b <- ggplot(modelresults %>% filter(transmitterallocation == 'ML'),aes(y = delta_CD4_recipient, x = 1-p_variants_1))+
  geom_hex(aes(fill = (..density..))) + 
  #stat_density_2d(aes(fill = (..density..)), geom = "raster", contour = FALSE) +
  #geom_density_2d_filled(bins = 20)+
  scale_fill_distiller(palette = 'GnBu',
                       direction = -1,
                       trans = 'sqrt',
                       'Kernel Density \nEstimate') +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0,0),
                     limits = c(0,1),
                     breaks = seq(0, 1, by = 0.25))+
  scale_y_continuous(limits = c(-0.65,0), 
                     expand = c(0,0), 
                     breaks = seq(-0.6, 0, by = 0.1),
                     name =expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1))) +
  facet_grid(cols = vars(dataset_id), switch = 'y',
             labeller = as_labeller( c('MF' = 'Male to Female',
                                       'PWID' = 'PWID',
                                       'FM' = 'Female to Male',
                                       'MMI' = 'MSM Insertive',
                                       'MMR' = 'MSM Receptive')))+
  my_theme 

legend <- get_legend(
  # create some space to the left of the legend
  plt_3b + theme(legend.position = 'right', legend.box.margin = margin(10, 0, 15, 10))
)
  
panel_3_main <- cowplot::plot_grid(plt_3a,NA,
                              plt_3b + theme(legend.key.width=unit(1,"cm")), 
                              align = 'hv',
                              nrow = 3, 
                              rel_heights = c(1,0.1, 1), 
                              labels = c('A','', 'B'), label_size = 9)

panel_3 <- cowplot::plot_grid(panel_3_main, legend,  rel_widths = c(0.85, .15))

ggsave('figure4.eps', device=cairo_ps,   width = 260, height = 140,  units = 'mm')
Sys.sleep(0.5)
panel_3
dev.off()
############################################## Panel 4 ##############################################

plt_4a <- ggplot(SpVLresults) +
  geom_boxplot(aes(x = factor(multiplicity, levels =c('log10_SpVL_single', 'log10_SpVL_multiple')), 
                   y = log10_SpVL,
                   fill = multiplicity),  
               outlier.size = 0.5,
               notchwidth = 0.5)+
  scale_x_discrete('Multiplicity',
                   labels = c('log10_SpVL_multiple' = 'Multiple',  'log10_SpVL_single' = 'Single') ) +
  scale_y_continuous('Mean Log10 SpVL') +
  scale_fill_brewer(palette = 'BuGn') +
  #scale_colour_brewer(palette = 'OrRd') +
  #coord_cartesian(xlim = c(0,365*10))+ #cut at 10 years
  facet_grid(cols = vars(riskgroup),
             labeller = as_labeller( c('MF' = 'Male to Female',
                                       'PWID' = 'PWID',
                                       'FM' = 'Female to Male',
                                       'MMI' = 'MSM Insertive',
                                       'MMR' = 'MSM Receptive')),
             
             switch = 'y')+
  my_theme

splv_tb <- SpVLresults %>% summarise(value = median(log10_SpVL), n = n(), sd = sd(log10_SpVL), .by = c('riskgroup', 'multiplicity')) %>%
  mutate(multiplicity = gsub('log10_SpVL_', '', multiplicity)) %>%
  rename(log10_SpVL = value) %>%
  pivot_wider(names_from = multiplicity, values_from = c('log10_SpVL','sd')) %>%
  group_by(riskgroup) %>%
  mutate(diff = log10_SpVL_multiple - log10_SpVL_single)

t.test(splv_tb$sd_multiple, splv_tb$sd_single, alternative = 'g', paired = T)

plt_4b <- ggplot(cd4results) +
  geom_ribbon(aes(x = time,
                  ymin = `0.01`,
                  ymax = `0.99`,
                  fill = multiplicity), 
              alpha = 0.7)+
  # geom_line(aes(x = time, y= `0.5`, 
  #      colour = multiplicity),
  #   linetype= 'dashed', 
  #  linewidth = 0.6) +
  scale_x_continuous('Days Post Infection', 
                     expand = c(0,0)) + 
  scale_y_continuous(str_wrap('Proportion of Cohort with < 350 CD4 mm3', 25),
                     expand = c(0,0)) +
  #scale_fill_brewer(palette = 'BuGn') + 
  scale_fill_manual(values = c('#2b8cbe', '#ccebc5')) + 
  scale_colour_manual(values = c('#084081', '#238b45')) + 
  coord_cartesian(xlim = c(365*2,365*7),
                  ylim = c(0.7, 1.05))+ #cut at 10 years
  facet_grid(cols = vars(riskgroup),
             switch = 'y',
             labeller = as_labeller( c('MF' = 'Male to Female',
                                       'PWID' = 'PWID',
                                       'FM' = 'Female to Male',
                                       'MMI' = 'MSM Insertive',
                                       'MMR' = 'MSM Receptive')))+
  my_theme

legend <- get_legend(
  # create some space to the left of the legend
  plt_4b + theme(legend.position = 'right', legend.box.margin = margin(10, 0, 15, 10))
)


panel_4_main <- cowplot::plot_grid(plt_4a, NA,plt_4b, 
                              align = 'hv',
                              nrow = 3, 
                              rel_heights = c(1, 0.1, 1),
                              labels = c('A', '', 'B'),
                              label_size = 9)

panel_4 <- cowplot::plot_grid(panel_4_main, legend,  rel_widths = c(0.85, .15))

ggsave('figure5.eps', device=cairo_ps,   width = 260, height = 140, units = 'mm')
Sys.sleep(0.5)
panel_4
dev.off()

###################################################################################################
# Supplementary figures

################################### Fig S1 #################################
# Change in distribution of SpVL over transmission cycle
MC <- 4
VC <- 0.5

t <- CalcTransmitter(MC, VC)
r <- CalcRecipient(MC, VC)

s1_data <- tibble(Carrier = rnorm(100000, MC, sqrt(VC)), 
                  Transmitter = rnorm(100000, t[['mean']], sqrt(t[['var']])),
                  Recipient = rnorm(100000, r[['mean']], sqrt(r[['var']]))) %>%
  pivot_longer(cols = everything(), names_to = 'stage', values_to = 'vl') %>%
  mutate(stage = forcats::fct_relevel(stage, c('Carrier', 'Transmitter', 'Recipient')))

plt_s1 <- ggplot(s1_data) + 
  geom_density(aes(x = vl,
                   fill = stage),
               alpha = 0.5,
               colour = 'white')+
  scale_x_continuous(expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')), 
                     expand= c(0,0), 
                     limits = c(0,7))+
  scale_y_continuous('Density', 
                     expand= c(0,0))+
  scale_fill_brewer(palette = 'GnBu')+
  geom_vline(xintercept = t[['mean']],
             linetype = 2, 
             colour = '#7bccc4')+
  geom_vline(xintercept = MC , 
             linetype = 2, 
             colour = '#e0f3db')+
  geom_vline(xintercept = r[['mean']],
             linetype = 2, 
             colour = '#084081')+
  facet_wrap(.~stage, ncol = 3)+
  my_theme

# A4: 210 x 297 mm, 
# 170 x 100
ggsave(plot = plt_s1,
       filename = paste(figs_dir,sep = '/', "panel_s1.jpeg"), 
       device = jpeg, 
       width = 170, 
       height = 60,
       units = 'mm')



################################### Fig S1 #################################
# source(./scripts/observing_sigeffects.R)
ggsave('figureS2.eps', device=cairo_ps,   width = 140, height = 110, units = 'mm')
Sys.sleep(0.5)
ggpubr::ggqqplot(vls$vl) + my_theme
dev.off()


################################### Fig S2 #################################
# Transmitter allocation strategies

vlpairs <- shcs_data %>%
  select(-contains('allocate')) %>%
  rowwise() %>%
  mutate(transmitter_random = AllocateTransmitter(log10_SpVL_1, log10_SpVL_2, transmitter = t, recipient = r, method = 'random')) %>%
  mutate(transmitter_max = AllocateTransmitter(log10_SpVL_1, log10_SpVL_2, transmitter = t, recipient = r, method = 'max')) %>%
  mutate(transmitter_ml = AllocateTransmitter(log10_SpVL_1, log10_SpVL_2, transmitter = t, recipient = r, method = 'ml')) %>%
  mutate(transmitter_bayes = AllocateTransmitter(log10_SpVL_1, log10_SpVL_2, transmitter = t, recipient = r, method = 'bayes')) %>%
  select(c(log10_SpVL_1, log10_SpVL_2, starts_with('transmitter'))) %>%
  pivot_longer(cols = starts_with('log10_SpVL'), values_to = 'log10SpVL') %>%
  mutate(across(starts_with('transmitter'), .fns = ~ case_when(name == 'log10_SpVL_1' &  .x == 1 ~ 't',
                                                               name == 'log10_SpVL_2' &  .x == 2 ~ 't',
                                                               .default = 'r'))) %>%
  pivot_longer(cols = starts_with('transmitter'), values_to = 'role', names_to = 'allocation') %>%
  mutate(allocation = str_remove(allocation, 'transmitter_')) %>%
  mutate(allocation = str_to_title(allocation)) %>%
  mutate(allocation = case_when(allocation == 'Ml' ~ 'ML', 
                                .default = allocation))

panel_s2 <- ggplot(vlpairs)+
  geom_histogram(aes(x = log10SpVL,
                     fill = role),
                 binwidth = 0.25) +
  scale_fill_brewer(palette = 'GnBu', 
                    'Role Allocated',
                    labels = c('Recipient', 'Transmitter'))+
  scale_y_continuous('Count',
                     expand = c(0,0),
                     limits = c(0,70), 
                     breaks = seq(0,70,by=10))+
  scale_x_continuous(expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')), 
                     expand =  c(0,0),
                     limits = c(1,7))+
  facet_wrap(.~allocation)+
  my_theme + 
  theme(legend.position = 'right')




ggsave('figureS3.eps', device=cairo_ps,   width = 140, height = 110, units = 'mm')
Sys.sleep(0.5)
panel_s2 
dev.off()




################################### Fig S3 #################################
# Validate Simulation Plots
source('/scripts/validate_simcohorts.R')

bias_covmats <- bias_covmats %>%
  mutate(across(starts_with('Var'), .fns = ~ case_when(.x == 'age.inf_transmitter' ~ 'Transmitter Age',
                                                       .x == 'age.inf_recipient' ~ 'Recipient Age',
                                                       .x == 'log10_SpVL_couplemean' ~ 'Log10SpVL Pair Mean',
                                                       .x == 'sex' ~ 'Sex'))) 

bias_means <- bias_means %>%
  mutate(variable = case_when(variable == 'age.inf_transmitter' ~ 'Transmitter Age',
                              variable == 'age.inf_recipient' ~ 'Recipient Age',
                              variable == 'log10_SpVL_couplemean' ~ 'Log10SpVL Pair Mean',
                              variable == 'sex' ~ 'Sex'))


plt_s3a <- ggplot(bias_means) +
  geom_col(aes(x = variable, y = bias, fill = abs(bias)))+ 
  facet_wrap(.~riskgroup) +
  my_theme + 
  scale_fill_distiller(palette = 'GnBu', direction = 1) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_blank())+
  scale_y_continuous(expression(paste(bar(x)["simulated"] - bar(x)['empirical'])), 
                     limits = c(-0.8, 0.8), 
                     breaks = seq(-0.8, 0.8, by = 0.2) %>% round(digits  = 2)
  )


plt_s3b <- ggplot(bias_covmats,
                  aes(x = Var1, 
                      y = Var2,
                      fill = abs(value))) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  geom_text(aes(label = signif(value, 3)), color = "black", size = 2) +
  scale_fill_distiller(palette = 'GnBu', direction = 1) + 
  facet_wrap(.~riskgroup)+ 
  my_theme+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank())

validate_df <- list(Simulated = stratified_data, Empirical = shcs_data_long_transmitterrandom) %>%
  bind_rows(.id = 'id') %>%
  select(sex, 
         riskgroup, 
         age.inf, 
         id,
         SpVL_couplemean) %>% 
  filter(riskgroup %in% c('HET', "MSM", 'PWID'))

#ggplot(validate_df)+
#geom_histogram(aes(x = SpVL_couplemean))+
#scale_x_log10()+
#facet_grid(id~riskgroup, switch = 'y')

plt_s3c <- ggplot(validate_df ) +
  geom_histogram(aes(x= age.inf, 
                     y = after_stat(density)), 
                 binwidth = 2, 
                 fill = '#7bccc4')+
  scale_x_continuous('Age at Infection', 
                     limits = c(15,80), 
                     breaks = seq(16,80, by = 8),
                     expand = c(0,0))+
  scale_y_continuous('Density', 
                     limits = c(0,0.08),
                     expand = c(0,0))+
  my_theme +
  facet_grid(id~riskgroup, 
             switch = 'y')+
  theme(strip.placement = 'outside')


panel_s3 <- cowplot::plot_grid(plt_s3a,
                               plt_s3b,
                               plt_s3c,
                               nrow = 3,
                               ncol = 1,
                               align = 'hv',
                               labels = 'AUTO',
                               label_size = 10)

ggsave('figureS4.eps', device=cairo_ps,   width = 150, height = 200, units = 'mm')
Sys.sleep(0.5)
panel_s3
dev.off()
################################### Fig S4 #################################
# Heritability Model Plots

# Prior-predictive check
color_scheme_set("brewer-GnBu") 
s4a_plt <- pp_check(priorpredictivesim, ndraws = 1000) + my_theme 
s4a_plt 

# Posterior-predictive check
s4b_plt <- pp_check(heritability_model_transmitterrandom, ndraws = 1000) + my_theme 

s4c_plt <- pp_check(heritability_model_transmitterML, ndraws = 1000) + my_theme 


# Check Residuals
# Plot residuals with associated uncertainty
s4d_plt <- shcs_data_long_transmitterrandom %>%
  add_residual_draws(heritability_model_transmitterrandom) %>%
  ggplot(aes(x = .row, y = .residual)) + 
  scale_y_continuous('Residual') +
  scale_x_continuous(expression('Individual SpVL' ['ij'])) +
  stat_pointinterval(colour = '#7bccc4', 
                     shape = 4,
                     size = 0.5, 
                     interval_alpha = 0.7) + 
  my_theme


s4e_plt <- shcs_data_long_transmitterML %>%
  add_residual_draws(heritability_model_transmitterML) %>%
  ggplot(aes(x = .row, y = .residual)) + 
  scale_y_continuous('Residual') +
  scale_x_continuous(expression('Individual SpVL' ['ij'])) +
  stat_pointinterval(colour= '#7bccc4', 
                     shape = 4, 
                     size = 0.5, 
                     interval_alpha = 0.7) + 
  my_theme


# QQ plot
s4f_plt <- shcs_data_long_transmitterrandom %>%
  add_residual_draws(heritability_model_transmitterrandom) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) + 
  geom_qq(colour = '#7bccc4') + 
  geom_qq_line()+ 
  scale_y_continuous('Sample Quantiles') +
  scale_x_continuous('Theoretical Quantiles') +
  my_theme

s4g_plt <- shcs_data_long_transmitterML %>%
  add_residual_draws(heritability_model_transmitterML) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) + 
  geom_qq(colour = '#7bccc4') + 
  geom_qq_line()+ 
  scale_y_continuous('Sample Quantiles') +
  scale_x_continuous('Theoretical Quantiles') +
  my_theme

panel_s4_top <- plot_grid(s4a_plt, 
                          s4b_plt,
                          s4c_plt, 
                          ncol = 3, 
                          align = 'hv', 
                          labels = 'AUTO', 
                          label_size = 8)

panel_s4_middle <- plot_grid(s4d_plt, 
                             s4e_plt, 
                             nrow = 2, 
                             align = 'hv', 
                             labels = c('D', 'E'),
                             label_size = 8)

panel_s4_bottom <- plot_grid(s4f_plt,
                             s4g_plt, 
                             ncol = 2, 
                             align = 'hv', 
                             labels = c('F', 'G'), 
                             label_size = 8)

panel_s4_top_legend <- get_legend(s4a_plt + theme(legend.position = 'bottom')) 

panel_s4 <- plot_grid(panel_s4_top,
                      panel_s4_top_legend, 
                      panel_s4_middle, 
                      panel_s4_bottom, 
                      nrow = 4, 
                      rel_heights = c(1,0.1,1.5,1))

ggsave('figureS5.eps', device=cairo_ps, height = 220, width = 150, units = 'mm')
Sys.sleep(0.5)
panel_s4
dev.off()

################################### Fig S5 #################################
# Plot posterior chains
h2_chains <- lapply(model_list, ggs) %>% # Warning message In custom.sort(D$Parameter) : NAs introduced by coercion
  do.call(rbind.data.frame, .) %>%
  mutate(model = sapply(str_split(rownames(.), '\\.'), '[', 1)) %>%
  filter(Parameter %in% c('b_Intercept',  'b_sexF', 'b_partnerrecipient',
                          'b_age.inf_category25M29', 'b_age.inf_category30M39',
                          'b_age.inf_category40M80', 'b_riskgroupMSM', 
                          'b_riskgroupPWID', 'b_riskgroupOTHER', 'b_riskgroupUNKNOWN',
                          'sd_log10_SpVL_couplemean__Intercept', 'sigma')) %>%
  filter(Iteration > 1000) %>%
  mutate(Parameter = factor(Parameter,
                            levels = c('b_Intercept',  'b_sexF', 'b_partnerrecipient',
                                       'b_age.inf_category25M29', 'b_age.inf_category30M39',
                                       'b_age.inf_category40M80', 'b_riskgroupMSM', 
                                       'b_riskgroupPWID', 'b_riskgroupOTHER', 'b_riskgroupUNKNOWN',
                                       'sd_log10_SpVL_couplemean__Intercept', 'sigma'),
                            ordered = TRUE,
                            labels = c(expression(paste(beta[0])),
                                       expression(paste(beta['sexF'])),
                                       expression(paste(beta['recipient'])),
                                       expression(paste(beta['25-29'])),
                                       expression(paste(beta['30-39'])),
                                       expression(paste(beta['40-80'])),
                                       expression(paste(beta['MSM'])),
                                       expression(paste(beta['PWID'])),
                                       expression(paste(beta['OTHER'])),
                                       expression(paste(beta['UNKNOWN'])),
                                       expression(paste(sigma[mu])),
                                       expression(paste(sigma[epsilon]))
                            )
  ))

plt_s5 <- ggplot(h2_chains,
                 aes(x = Iteration,
                     y = value, 
                     col = as.factor(Chain))
                 )+
  geom_line(alpha = 0.7)+
  facet_grid(Parameter ~ model,
             scale  = 'free_y',
             switch = 'y',
             labeller = label_parsed)+
  scale_colour_brewer(palette = 'GnBu', 'Chains') +
  my_theme + 
  theme(legend.position = 'bottom', 
        strip.text.y = element_text(family = 'sans'))

ggsave('figureS6.eps', device=cairo_ps, height = 250, width = 150, units = 'mm')
Sys.sleep(0.5)
plt_s5
dev.off()

################################### Fig S6 #################################
# Plot Prior ~ Posterior distributions
model_list <- list(random = heritability_model_transmitterrandom, 
                   ml = heritability_model_transmitterML)


r2_bayes <- lapply(model_list, r2_posterior) %>%
  lapply(., function(x) x$R2_Bayes)%>%
  do.call(cbind.data.frame, .) %>%
  pivot_longer(cols = everything(), names_to = 'model', values_to = 'R2')

plt_s6 <- ggplot(r2_bayes,
                 aes(x = R2,
                     fill = as.factor(model)))+
  geom_density(aes(y = after_stat(density)),
               colour = NA, 
               alpha = 0.7)+
  scale_fill_brewer(palette = 'GnBu', 
                    'Allocation', 
                    labels= c('Maximum Likelihood', 'Random') ) +
  scale_y_continuous(expand = c(0,0),
                     'Frequency')+
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,0.5), 
                     expression(R^2))+
  facet_wrap(.~model, 
             labeller = as_labeller(c('ml' = 'ML', 'random' = 'Random')))+
  my_theme #+ 
# theme(legend.position = c(0.75,0.75))

ggsave('figureS7.eps', device=cairo_ps, height = 70, width = 140, units = 'mm')
Sys.sleep(0.5)
plt_s6
dev.off()

################################### Fig S7 #################################
# Transmission Model Plots
# This script plots extracted estimates/data from the Thompson transmission model to produce 
# plots describing the within-host processes encoded within the model
my_palette <- RColorBrewer::brewer.pal(name="GnBu",n=9)[4:9]

# For a given SpVL, plot the probabilities of transmission over time since infection (Thompson)
spvl_dur <- read_csv('./data/spvl_duration.csv') %>%
  mutate(stage = forcats::fct_relevel(stage, 'primary', 'not_primary'))


plt_s7a <- spvl_dur %>%
  ggplot() +
  geom_line(aes(x = time, 
                y = vl, 
                colour = as.factor(spvl)), 
            size = 1.5) +
  scale_colour_brewer(palette = 'GnBu', name = 'SpVL')+
  my_theme+
  scale_y_log10('Viral Load',
                limits = c(10**0, 10**8),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  facet_grid(rows= vars(model), cols = vars(stage), scales = 'free_x', switch = 'y', 
             labeller = as_labeller(c('new' = 'Updated', 
                                      'thompson' = 'Thompson et al.', 
                                      'primary' = 'Primary Infection', 
                                      'not_primary' = 'Chronic & Pre-AIDS'))) + 
  scale_x_continuous(expand = c(0.05,0), 'Time (Years)', 
                     limits = ~ c(min(.x), max(.x))) +
  theme(strip.placement = 'outside',
        panel.spacing = unit(0, "lines"))

# Viral Loads in pimary and acute stages
feibig_vls <-matrix(c( 97,     42084,   64344,    3878,     1714,     # 0.1 (w = 0.1)
                       421,     92159,   136380,   27524,     7214,    # 0.25 (w = 0.15)
                       2742,    262131,   491492,   88507,    25176,   # 0.5 (w = 0.5)
                       22626,   967867,  1216484,   427765,   105021,  # 0.75 (w = 0.15)
                       42084,  3574304,  4492440,   1341713,   471801), ncol = 5) %>%
  as_tibble() %>%
  rename(.,all_of(c(`0.1` = 'V1', `0.25` = 'V2', `0.5` = 'V3', `0.75` = 'V4', `0.9` = 'V5'))) %>%
  mutate(stage = c('I', 'II', 'III', "IV", 'V')) 

plt_s7b <- ggplot(feibig_vls , 
                  aes(x = as.factor(stage),
                      fill =  as.factor(stage))) +
  geom_boxplot(aes(
    lower = `0.25`, 
    upper = `0.75`, 
    middle = `0.5`, 
    ymin = `0.1`, 
    ymax = `0.9`),
    stat = "identity")+
  scale_y_log10(expression(paste("Viral Load", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**0, 10**8),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  scale_x_discrete('Fiebig Stage')+
  scale_fill_brewer(palette = 'GnBu')+
  my_theme


# The different VLs during preAIDS (Mellors et al 1995)
fitted_dist <- fitdistrplus::fitdist(log(c(20190, 62530, 103900, 188000, 1235000),10), distr = "norm", method = "mle") 
preaids_vals <- rnorm(n = 10000, mean = fitted_dist$estimate[1], sd = fitted_dist$estimate[2]) %>%
  cbind.data.frame() %>%
  mutate(vl = raise_to_power(10, `.`))

plt_s7c <- ggplot(preaids_vals)+
  geom_histogram(aes(x = vl, 
                     y =  after_stat(ndensity)), 
                 binwidth = 0.2, 
                 fill = '#7bccc4')+
  scale_x_log10(expression(paste("Viral Load", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**8),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  scale_y_continuous('Density', 
                     expand = c(0,0),
                     breaks = seq(0, 1, by = 0.2),
                     limits = c(0,1.1))+
  my_theme


p_aq <- data.frame(route = c('FM', 'MF', 'MSM (IA)', 'MSM (RA)', 'PWID', 'MTC'),
                   p_aq.ci.lower = c(1/10000, 6/10000, 4/10000, 102/10000, 41/10000, 1700/10000),
                   p_aq = c(4/10000, 8/10000, 11/10000, 138/10000, 63/10000, 2260/10000),
                   p_aq.ci.upper = c(14/10000, 11/10000, 28/10000, 186/10000, 92/10000, 2900/10000)) %>% 
  pivot_longer(!route, names_to = 'type', values_to ='probability') %>%
  mutate(route = factor(route, levels = c('FM', 'MF', 'MSM (IA)', 'PWID', 'MSM (RA)','MTC')))



p_mv_sysreview <- data.frame('exposure' = factor(c('HSX:mf', 'HSX:fm', 'HSX:unknown', 'MSM', 'MTC:pre', 'MTC:intra', 'MTC:post', 'MTC:unknown', 'PWID'), 
                                                 levels = c('HSX:mf', 'HSX:fm', 'HSX:unknown', 'MSM', 'MTC:pre', 'MTC:intra', 'MTC:post', 'MTC:unknown', 'PWID')),
                             p_mv.ci.lower = c(0.14, 0.08, 0.12, 0.22, 0.08, 0.14, 0.03, 0.08, 0.24),
                             p_mv = c(0.21, 0.13, 0.32, 0.30, 0.17, 0.27, 0.18, 0.18, 0.37),
                             p_mv.ci.upper = c(0.31, 0.21, 0.61, 0.4, 0.33, 0.45, 0.57, 0.35, 0.53))

#write_csv(p_mv_sysreview , 'p_mv_sysreview.csv')

plt_s7d <- p_aq %>%
  pivot_wider(names_from = type, 
              values_from = probability) %>%
  mutate(exposure = c('HSX:fm', 'HSX:mf', 'MSM', 'MSM', 'PWID', 'MTC:unknown')) %>%
  right_join(p_mv_sysreview, 
             by = 'exposure') %>%
  filter(!is.na(route)) %>%
  ggplot() +
  geom_point(aes(x = p_aq, 
                 y = p_mv,
                 color = exposure),
             key_glyph = draw_key_rect)+
  geom_linerange(aes(xmin = p_aq.ci.lower,
                     xmax = p_aq.ci.upper, 
                     y = p_mv,
                     color = exposure ), 
                 key_glyph = draw_key_rect)+
  geom_linerange(aes(ymin = p_mv.ci.lower, 
                     ymax = p_mv.ci.upper,
                     x = p_aq, 
                     color = exposure),
                 key_glyph = draw_key_rect) +
  scale_x_log10(name = 'Probability of Acquisition', 
                expand = c(0,0), 
                limits = c(0.0001, 1), 
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", label_math())) + 
  scale_y_continuous(name = 'Probability of Multiple Founders', 
                     expand = c(0,0), 
                     limits = c(0,0.8)) +
  scale_colour_manual(values = c('#ccebc5','#a8ddb5','#7bccc4','#43a2ca','#0868ac'), 'Exposure', 
                      labels = as_labeller(c('HSX:fm' = 'HSX:FM',
                                             'HSX:mf' = 'HSX:MF',
                                             'MSM' = 'MSM',
                                             'MTC:unknown' = 'MTC',
                                             'PWID' = 'PWID'))) +
  #annotation_logticks()+
  my_theme +
  guides(color = guide_legend(override.aes = list(size = 2)))+
  theme(legend.position = c(0.85,0.7))+ coord_flip()

# Plot the proportion of the xth most common variants over time since infection
var_t <- read_csv('./data/variant_distribution_time.csv') %>% 
  mutate(x_var = as.factor(x_var)) %>% 
  rename(time = T)

plt_s7e <- ggplot(var_t, aes(x = time, y = P, fill = x_var)) + 
  geom_area() + 
  scale_x_continuous(name = 'Time', expand = c(0,0))+
  scale_y_continuous(name = 'Proportion of Xth Most Common Variant', expand = c(0,0))+
  scale_fill_brewer(palette = 'GnBu', na.value = 'black')+
  my_theme+
  theme(legend.position = 'none')

panel_s7_bottom <- cowplot::plot_grid(plt_s7b, 
                                      plt_s7c,
                                      plt_s7d, 
                                      plt_s7e,
                                      nrow = 2, 
                                      align = 'hv', 
                                      labels = c('B', 'C', 'D', 'E'), 
                                      label_size = 8)

panel_s7 <- cowplot::plot_grid(plt_s7a, 
                               panel_s7_bottom, 
                               labels= c('A', ''), 
                               rel_heights = c(0.75, 1), 
                               nrow = 2,
                               label_size = 8)

ggsave('figureS8.eps', device=cairo_ps, height = 250, width = 150, units = 'mm')
Sys.sleep(0.5)
panel_s7
dev.off()

################################### Fig S8 #################################
# Panels S8 and S9 are imported from Julián's model fitting script and are not described here


################################### Fig S9 #################################
# Panels S8 and S9 are imported from Julián's model fitting script and are not described here
 load('./data/posteriors_by_exposure.RData')

posterior_fp <- list(fm = ftm, 
                     mf = mtf,
                     mmi = msmi, 
                     mmr = msmr, 
                     pwid = pwid) %>%
  bind_rows(., .id = 'riskgroup') %>%
  as_tibble()%>%
  pivot_longer(cols = c(ppp, fEnv), values_to = 'value', names_to = 'parameter')

plt_s9 <- ggplot(posterior_fp) + 
  geom_histogram(aes(x = value,  y = after_stat(count / sum(count))), fill = "#4EB3D3") +
  scale_x_continuous('Parameter value', expand = c(0,0))+
  scale_y_continuous('Frequency', expand = c(0,0))+
  facet_grid2(rows = vars(parameter), 
             cols = vars(riskgroup),
             scales = 'free', independent = 'x',
             labeller =  as_labeller(c('fm' = 'Female-to-male',
                                       'mf' = 'Male-to-female',
                                       'mmi' = 'MSM Insertive',
                                       'mmr' = 'MSM Receptive',
                                       'pwid' = 'PWID',
                                       'ppp' = 'p',
                                       'fEnv' = 'f'))) +  
  my_theme

ggsave('figureS9.eps', device=cairo_ps, width = 250, height  = 110, units = 'mm')
Sys.sleep(0.5)
plt_s9
dev.off()



###########################################################################################
