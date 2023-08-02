###################################################################################################
# Supplementary Figures


################################### Dependencies ###################################
require(tidyverse)
require(showtext)
require(showtextdb)
require(cowplot)
require(RColorBrewer)

source('./scripts/BonhoefferEqns.R')


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
  geom_density(aes(x = vl ,fill = stage), alpha = 0.5, colour = 'white')+
  scale_x_continuous(expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')), expand= c(0,0), limits = c(0,7))+
  scale_y_continuous('Density', expand= c(0,0))+
  scale_fill_brewer(palette = 'OrRd')+
  geom_vline(xintercept = MD,  linetype = 2, colour = '#fdbb84')+
  geom_vline(xintercept = MC ,  linetype = 2, colour = '#fee8c8')+
  geom_vline(xintercept = MR, linetype = 2, colour = '#e34a33')+
  facet_wrap(.~stage, ncol = 3)+
  my_theme

# A4: 210 x 297 mm, 
# 170 x 100
ggsave(plot = plt_s1, filename = paste(figs_dir,sep = '/', "plt_s1.jpeg"), 
       device = jpeg, width = 170, height = 60, units = 'mm')


################################### Fig S2 #################################
# Densities and correlations of simulated variables.


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


plt_s2a <- ggplot(bias_means) +
  geom_col(aes(x = variable, y = bias, fill = abs(bias)))+ 
  facet_wrap(.~riskgroup) +
  my_theme + 
  scale_fill_distiller(palette = 'OrRd', direction = 1) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(expression(paste(bar(x)["simulated"] - bar(x)['empirical'])), 
                     limits = c(-0.8, 0.8), 
                     breaks = seq(-0.8, 0.8, by = 0.2) %>% round(digits  = 2)
                     )


plt_s2b <- ggplot(bias_covmats, aes(x = Var1, y = Var2, fill = abs(value))) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  geom_text(aes(label = signif(value, 3)), color = "black", size = 3) +
  scale_fill_distiller(palette = 'OrRd', direction = 1) + 
  facet_wrap(.~riskgroup)+ 
  my_theme+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank())

validate_df <- list(Simulated = stratified_data, Empirical = shcs_data_long) %>%
  bind_rows(.id = 'id') %>%
  select(sex, riskgroup, age.inf,  id) %>% filter(riskgroup %in% c('HET', "MSM", 'PWID'))

plt_2c <- ggplot(validate_df ) +
  geom_histogram(aes(x= age.inf, y = after_stat(density)), binwidth = 2, fill = '#e34a33')+
  facet_wrap(id~riskgroup, ncol = 3)+
  scale_x_continuous('Age at Infection', limits = c(15,80), breaks = seq(16,80, by = 8), expand = c(0,0))+
  scale_y_continuous('Density', limits = c(0,0.08), expand = c(0,0))+
  my_theme


panel_s2 <- cowplot::plot_grid(plt_s2a, plt_s2b, plt_2c,nrow = 3, ncol = 1,align = 'hv', labels = 'AUTO')

ggsave(plot = panel_s2  , filename = paste(figs_dir,sep = '/', "panel_s2_new.jpeg"), device = jpeg, width = 15, height = 20)


################################### Fig S4 #################################
# Heritability Model Plots

# Check Residuals
# Plot residuals with associated uncertainty
plt_residuals <- shcs_data_long %>%
  add_residual_draws(heritability_model) %>%
  ggplot(aes(x = .row, y = .residual)) + 
  stat_pointinterval() + 
  my_theme

# QQ plot
plt_qq <- shcs_data_long %>%
  add_residual_draws(heritability_model) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) + 
  geom_qq() + 
  geom_qq_line()+ 
  my_theme


# Check Residuals
# Plot residuals with associated uncertainty
plt_residuals <- shcs_data_long_transmitterallocated %>%
  add_residual_draws(heritability_model_transmittermax) %>%
  ggplot(aes(x = .row, y = .residual)) + 
  stat_pointinterval() + 
  my_theme

# QQ plot
plt_qq <- shcs_data_long_transmitterallocated %>%
  add_residual_draws(heritability_model_transmittermax) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) + 
  geom_qq() + 
  geom_qq_line()+ 
  my_theme


################################### Fig S5 #################################

################################### Fig S6 #################################
# Plot Prior ~ Posterior distributions
model_list <- list(random = heritability_model_randomallocation, 
                   max.spvl = heritability_model_transmittermax)


r2_bayes <- lapply(model_list, r2_posterior) %>%
  lapply(., function(x) x$R2_Bayes)%>%
  do.call(cbind.data.frame, .) %>%
  pivot_longer(cols = everything(), names_to = 'model', values_to = 'R2')
  
plot_r2_density <- ggplot(r2_bayes,
                          aes(x = R2,
                              fill = as.factor(model)))+
  geom_histogram(aes(y = after_stat(count / max(count))), colour = NA)+
  scale_fill_brewer(palette = 'OrRd', 'Chains') +
  scale_y_continuous(expand = c(0,0), 'Frequency')+
  scale_x_continuous(expand = c(0,0), limits = c(0,1), expression(R^2))+
  my_theme + 
  theme(legend.position = 'right', 
        strip.text.y = element_text(family = 'sans'))

ggsave(plot = plot_r2_density , filename = paste(figs_dir,sep = '/', "plot_r2_density.jpeg"), device = jpeg, width = 14, height = 14)


################################### Fig S7 #################################
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

plot_chains <- ggplot(h2_chains,
                      aes(x   = Iteration,
                          y   = value, 
                          col = as.factor(Chain)))+
  geom_line()+
  facet_grid(Parameter ~ model,
             scale  = 'free_y',
             switch = 'y',
             labeller = label_parsed)+
  scale_colour_brewer(palette = 'OrRd', 'Chains') +
  my_theme + 
  theme(legend.position = 'bottom', 
        strip.text.y = element_text(family = 'sans'))

ggsave(plot = plot_chains, filename = paste(figs_dir,sep = '/', "panel_chains.jpeg"), device = jpeg, width = 12, height = 19)


################################### Fig S8 #################################
# Transmission Model Plots
# This script plots extracted estimates/data from the Thompson transmission model to produce 
# plots describing the within-host processes encoded within the model


my_palette <- RColorBrewer::brewer.pal(name="OrRd",n=9)[4:9]

# For a given SpVL, plot the probabilities of transmission over time since infection (Thompson)
spvl_infdur <- read_csv('./data/spvl_infectionduration.csv')
plt_tda <- 
  ggplot(spvl_infdur)+
  geom_line(aes(x = t, y = p, colour = as.factor(spvl)), size = 1.5)+
  scale_color_manual(values = my_palette, name = 'SpVL')+
  my_theme+
  scale_x_continuous(expand = c(0,0), limits = c(0,17), 'Time') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.04), 'P(Acquisition)')+
  theme(legend.position = c(0.92,0.88))



# Primary Infection VLs (Feibig et al. 2003)
fiebig_q <- c(0.1,   0.15,      0.5,    0.15,    0.1)
# I 0-5 days\\II 5-10.3 days\\III 10.3-13.5 days\\IV 13.5-19.1 days\\V 19.1-88.6 days

#            I\\      II\\      III\\    IV\\       V 
npVals <- c( 97,     42084,   64344,    3878,     1714,     # 0.1 (w = 0.1)
            421,     92159,   136380,   27524,     7214,    # 0.25 (w = 0.15)
            2742,    262131,   491492,   88507,    25176,   # 0.5 (w = 0.5)
            22626,   967867,  1216484,   427765,   105021,  # 0.75 (w = 0.15)
            42084,  3574304,  4492440,   1341713,   471801) # 0.9 (w = 0.1)

fiebig_t <- c(5,5.3,3.2,5.6,69.5)/365 #Stage duration (Fiebig et al 2003) (Total time in acute infection = 0.243)
fiebig_q <- c(0.1,0.15,0.5,0.15,0.1) #Quantiles 
gp <- as.vector(outer(fiebig_t/sum(fiebig_t), fiebig_q))

# preAIDS  Infection VLs (Mellors et al 1995)
naVals <- round(c(10^4, 10^4.5, 10^5, 10^5.5, 10^6))
ga <- c(0.06053491, 0.23398873, 0.34935005, 0.25220102, 0.10392529) # Probability density function of Mellors VLs (ie. estimate is averaged over these vls for fixed period 0.75)




tauc = array(0, length(sp_ViralLoad))
taup = sum(fiebig_t)
taua = 0.75


plt_tda <- 
  ggplot(spvl_infdur)+
  geom_line(aes(x = t, y = p, colour = as.factor(spvl)), size = 1.5)+
  scale_color_manual(values = my_palette, name = 'SpVL')+
  my_theme+
  scale_x_continuous(expand = c(0,0), limits = c(0,17), 'Time') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.04), 'P(Acquisition)')+
  theme(legend.position = c(0.92,0.88))


# Plot the proportion of the xth most common variants over time since infection
var_t <- read_csv('./data/variant_distribution_time.csv') %>% 
  mutate(x_var = as.factor(x_var)) %>% 
  rename(time = T)

plot_var <- ggplot(var_t, aes(x = time, y = P, fill = x_var)) + 
  geom_area() + 
  scale_x_continuous(name = 'Time', expand = c(0,0))+
  scale_y_continuous(name = 'Proportion of Xth Most Common Variant', expand = c(0,0))+
  scale_fill_brewer(palette = 'OrRd', na.value = 'black')+
  my_theme+
  theme(legend.position = 'none')

# Values of F and P


# P(Aq) ~ P(MV)

panel_s1 <- cowplot::plot_grid( plot_var,NA, plt_tda, nrow = 1, align = 'hv', labels = c('A', '', 'B'), rel_widths = c(1,0.2,1))

## END ##
