# Figures

################################### Dependencies ###################################
source('./scripts/dependencies.R')
source('./scripts/BonhoefferEqns.R')
source('./scripts/AllocateTransmitter.R')
source('./scripts/simulate_cohorts_funcs.R')
source('./scripts/PostProcessing.R')
source('./scripts/interpret_outputs.R')
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
source('./scripts/deprecated/populationdata_models.R')
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


panel_2top <- plot_grid(plt_2a, plt_2c, plt_2d, ncol = 3, align = 'hv', label_size = 9, labels = 'AUTO')

panel_2 <- plot_grid(panel_2top, plt_2e,   nrow = 2, label_size = 9, labels = c('', 'D'))

ggsave(plot = panel_2 , filename = paste(figs_dir,sep = '/', "panel_2.jpeg"), 
       device = jpeg,  width = 170, height = 140,  units = 'mm')


############################################## Panel 3 ##############################################

#my_palette <- colorRampPalette(brewer.pal(9, "OrRd"))(20)


plt_3a <- ggplot(modelresults %>% filter(transmitterallocation == 'ML'),aes(y = SpVL_recipient, x = 1-p_variants_1),aes(y = y, x = 1-p_variants_1))+
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
  facet_grid(cols= vars(dataset_id), switch = 'y')+
  my_theme 


plt_3b <- ggplot(modelresults %>% filter(transmitterallocation == 'ML'),aes(y = delta_CD4_recipient, x = 1-p_variants_1))+
  geom_hex(aes(fill = (..density..))) + 
  #stat_density_2d(aes(fill = (..density..)), geom = "raster", contour = FALSE) +
  #geom_density_2d_filled(bins = 20)+
  scale_fill_distiller(palette = 'GnBu',
                       direction = -1,
                       trans = 'sqrt',
                       'Kernel Density Estimate') +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0,0),
                     limits = c(0,1),
                     breaks = seq(0, 1, by = 0.25))+
  scale_y_continuous(limits = c(-0.65,0), 
                     expand = c(0,0), 
                     breaks = seq(-0.6, 0, by = 0.1),
                     name =expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1))) +
  facet_grid(cols = vars(dataset_id), switch = 'y')+
  my_theme + 
  theme(legend.position = 'bottom')

panel_3 <- cowplot::plot_grid(plt_3a,
                              plt_3b + theme(legend.key.width=unit(1,"cm")), 
                              align = 'hv',
                              nrow = 2, 
                              rel_heights = c(0.8, 1), 
                              labels = c('A', 'B'), label_size = 9)

ggsave(plot = panel_3 , filename = paste(figs_dir,sep = '/', "panel_3.jpeg"), 
       device = jpeg,  width = 170, height = 140,  units = 'mm')


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
             switch = 'y')+
  my_theme+ 
  theme(legend.position = 'bottom')

panel_4 <- cowplot::plot_grid(plt_4a, plt_4b, 
                              align = 'hv',
                              nrow = 2, 
                              rel_heights = c(0.8,1),
                              labels = 'AUTO',
                              label_size = 9)

ggsave(plot = panel_4 , filename = paste(figs_dir,sep = '/', "panel_4.jpeg"), 
       device = jpeg,  width = 180, height = 120,  units = 'mm')



cat('All plots complete. \n')