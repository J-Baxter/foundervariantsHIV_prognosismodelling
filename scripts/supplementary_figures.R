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



ggsave(plot = panel_s2, 
       filename = paste(figs_dir,sep = '/', "panel_s2_new.jpeg"), 
       device = jpeg, 
       width = 170, 
       height = 140, 
       units = 'mm')


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

ggsave(plot = panel_s3, 
       filename = paste(figs_dir,sep = '/', "panel_s3_new.jpeg"), 
       device = jpeg,  
       width = 170, 
       height = 230,  
       units = 'mm')

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

ggsave(plot = panel_s4,
       filename = paste(figs_dir,sep = '/', "panel_s4.jpeg"), 
       device = jpeg,  
       width = 170, 
       height = 200, 
       units = 'mm')


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
                     col = as.factor(Chain)))+
  geom_line()+
  facet_grid(Parameter ~ model,
             scale  = 'free_y',
             switch = 'y',
             labeller = label_parsed)+
  scale_colour_brewer(palette = 'GnBu', 'Chains') +
  my_theme + 
  theme(legend.position = 'bottom', 
        strip.text.y = element_text(family = 'sans'))

ggsave(plot = plt_s5, 
       filename = paste(figs_dir,sep = '/', "panel_s5.jpeg"),
       device = jpeg,
       width = 170, 
       height = 250,
       units = 'mm')

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

ggsave(plot = plt_s6,
       filename = paste(figs_dir,sep = '/', "panel_s6.jpeg"), 
       device = jpeg,
       width = 170, 
       height = 80,  
       units = 'mm')


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
  theme(strip.placement = 'outside')

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

ggsave(plot = panel_s7, 
       filename = paste(figs_dir,sep = '/', "panel_s7.jpeg"), 
       device = jpeg,   
       width = 170, 
       height = 250, 
       units = 'mm')

################################### Fig S8 #################################
# Panels S8 and S9 are imported from Julián's model fitting script and are not described here


################################### Fig S9 #################################
# Panels S8 and S9 are imported from Julián's model fitting script and are not described here


################################### Fig S10 #################################
# Panel S10 - Comparing transmission model estimates of P(MV) (and across VLs)
source('./scripts/populationdata_acrossVL_models_R.r')

mf <- sapply(raise_to_power(10, 2:7), 
             populationmodel_acrossVL_Env, 
             PerVirionProbability = 1.765e-06, 
             PropExposuresInfective = 0.013296)

av_pmv <- populationmodel_acrossVL_Env(PerVirionProbability = 1.765e-06, 
                                        PropExposuresInfective = 0.013296) %>%
  cbind.data.frame() %>%
  pivot_longer(cols = contains('multiple'), 
               values_to = 'average',
               names_to = 'time') %>%
  mutate(Riskgroup = 'MF') 


df <- as_tibble(cbind(mf,  rownames(mf))) %>%
  mutate(Riskgroup = 'MF') %>%
  setNames(c(raise_to_power(10, 2:7) %>% as.character(), 'time', 'Riskgroup')) %>%
  mutate(time = as.character(time)) %>%
  left_join(., av_pmv, by = join_by(Riskgroup, time)) %>%
  pivot_longer(cols = contains(c('1')), values_to = 'P_MV', names_to = 'SpVL') %>%
  mutate(time = factor(time, levels = c('multiple_founder_proportion',
                                        'multiple_founder_proportion_primary',
                                        'multiple_founder_proportion_chronic',
                                        'multiple_founder_proportion_preaids')))
  

panel_s10 <- ggplot(df) +
  geom_bar(aes(x = as.numeric(SpVL), y = as.numeric(P_MV), fill = as.factor(SpVL)), stat = 'identity') +
  geom_hline(aes(yintercept = average), linetype = 'dashed')+
  my_theme + 
  scale_x_log10(expand = c(0,0), expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')),  
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", label_math(.x))) +
  scale_y_continuous(expand = c(0.01,0), 'P(Multiple Variants)', limits = c(0,1), breaks = seq(0,1 ,by = 0.2))+
  scale_fill_brewer(palette = 'GnBu')+
  facet_grid(cols = vars(time), switch = 'y',
             labeller = as_labeller( c('multiple_founder_proportion' = 'Average',
                                       'multiple_founder_proportion_primary' = 'Primary',
                                       'multiple_founder_proportion_chronic' = 'Chronic' ,
                                       'multiple_founder_proportion_preaids' = 'Pre-AIDS'))) + 
  theme(strip.placement = 'outside')

ggsave(plot = panel_s10,
       filename = paste(figs_dir,sep = '/', "panel_s10.jpeg"),
       device = jpeg,  
       width = 170, 
       height = 60,  units = 'mm')



# Temporary
#ggplot(bind_rows(results) %>% filter(transmitterallocation == 'ML') %>% filter(riskgroup_recipient != 'UNKNOWN')) + 
  #geom_histogram(aes(x = 1-p_variants_1, y = after_stat(count / max(count))), binwidth = 0.1) + 
  #geom_vline(aes(xintercept = multiple_founder_proportion), data = av_pmv %>% rename(riskgroup_recipient = Riskgroup), linetype = 'dashed')+
  #scale_x_continuous('P(MV)', limits = c(0,1), breaks = seq(0,1, by = 0.2), expand = c(0,0))+
  #scale_y_continuous('Frequency', limits = c(0,1), breaks = seq(0,1, by = 0.2), expand = c(0,0))+
  #my_theme + 
  #facet_grid(cols  = vars(riskgroup_recipient)) 


################################### Fig S11 #################################
plt_s11a <-  modelresults %>% 
  filter(transmitterallocation == 'ML') %>%
  ggplot(aes(y = SpVL_recipient, 
             x = 1-p_particles_1))+
  stat_density_2d(aes(fill = (..density..)**(1/3)),
                  geom = "raster", 
                  contour = FALSE) +
  scale_fill_distiller(palette = 'BuGn',
                       direction = 1,
                       breaks=1e-6*seq(0,10,by=2)) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,1),
                     breaks = seq(0, 1, by = 0.25))+
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**1, 10**7.5),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  my_theme +
  facet_grid(rows = vars(dataset_id), 
             switch = 'y')


plt_s11b <- modelresults %>% 
  filter(transmitterallocation == 'ML') %>%
  ggplot(aes(y = delta_CD4_recipient,
             x = 1-p_particles_1))+
  stat_density_2d(aes(fill = (..density..)**(1/3)), 
                  geom = "raster",
                  contour = FALSE) +
  scale_fill_distiller(palette = 'BuGn',
                       direction = 1, 
                       breaks=1e-6*seq(0,10,by=2)) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0,0),
                     limits = c(0,1),
                     breaks = seq(0, 1, by = 0.25))+
  scale_y_continuous(limits = c(-0.5,0), 
                     expand = c(0,0), 
                     breaks = seq(-0.5, 0, by = 0.1),
                     name =expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1))) +
  my_theme +
  facet_grid(cols = vars(dataset_id), 
             switch = 'y') 




################################### Fig S12 #################################

plt_s12a <- ggplot(modelresults ,aes(y = SpVL_recipient, x = 1-p_variants_1))+
  stat_density_2d(aes(fill = (..density..)**(1/3)),
                  geom = "raster", 
                  contour = FALSE) +
  scale_fill_distiller(palette = 'BuGn',
                       direction = 1,
                       breaks=1e-6*seq(0,10,by=2)) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,1),
                     breaks = seq(0, 1, by = 0.25))+
  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**1, 10**7.5),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  my_theme +
  facet_grid(rows = vars(dataset_id), 
             cols = vars(transmitterallocation),
             switch = 'y')


plt_s12b <- ggplot(modelresults,
                   aes(y = delta_CD4_recipient,
                       x = 1-p_variants_1))+
  stat_density_2d(aes(fill = (..density..)**(1/3)), 
                  geom = "raster",
                  contour = FALSE) +
  scale_fill_distiller(palette = 'BuGn',
                       direction = 1, 
                       breaks=1e-6*seq(0,10,by=2)) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0,0),
                     limits = c(0,1),
                     breaks = seq(0, 1, by = 0.25))+
  scale_y_continuous(limits = c(-0.5,0), 
                     expand = c(0,0), 
                     breaks = seq(-0.5, 0, by = 0.1),
                     name =expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1))) +
  my_theme + 
  facet_grid(rows = vars(dataset_id), 
             cols = vars(transmitterallocation),
             switch = 'y')

  
# SHCS
#plt_2e <- ggplot(shcs_data , aes(x = SpVL_1, SpVL_2)) +
# geom_point( #'#CB6015' #'#66c2a4','#2ca25f','#006d2c'
# colour = '#ef654a',
#  shape = 4, alpha = 0.4, size = 1) +
#scale_x_log10(limits = c(10**2, 10**7),
#              expand = c(0.05,0),
#              name = expression(paste("SpVL Partner 1", ' (', Log[10], " copies ", ml**-1, ')')),
#              breaks = trans_breaks("log10", function(x) 10**x),
#              labels = trans_format("log10", math_format(.x))) +
#scale_y_log10(name = expression(paste("SpVL Partner 2", ' (', Log[10], " copies ", ml**-1, ')')),
#              limits = c(10**2, 10**7),
#              expand = c(0.05,0),
#               breaks = trans_breaks("log10", function(x) 10**x),
#               labels = trans_format("log10", math_format(.x))) +
#my_theme


#Risk Group
#plt_2f <- ggplot(shcs_data_long_transmitterML , aes(x =  riskgroup)) +
#geom_bar( fill = '#ef654a') + 
#scale_y_continuous(expand = c(0,0), name = 'Count') + 
#scale_x_discrete(expand = c(0,0), name = 'Risk Group', labels = c('HET', 'MSM', 'PWID', 'OTHER', 'NA')) +
#my_theme


# Age
#plt_2g <- ggplot(shcs_data_long_transmitterML) + 
#geom_histogram(aes(x = age.inf), fill = '#ef654a') + 
# scale_y_continuous(expand = c(0.01,0.01), name = 'Count') + 
#scale_x_continuous('Age at Infection', limits = c(15,80), breaks = seq(16,80, by = 8), expand = c(0,0))+
#my_theme

############################################## Panel 3 ##############################################


#plt_3a <- ggplot(modelres %>% filter(transmitterallocation == 'ML'))+
#  geom_point(aes(y = SpVL_recipient, x = 1-p_variants_1, colour = riskgroup_recipient),size = 1, shape = 4)+
#  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
#                     expand = c(0.02,0.02),
#                     limits = c(0,0.4),
#                     breaks = seq(0, 0.4, by = 0.1))+
#  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
#                limits = c(10**2, 10**7),
#                expand = c(0.05,0),
#                breaks = trans_breaks("log10", function(x) 10**x),
#                labels = trans_format("log10", math_format(.x))) +
#  scale_color_brewer(palette = 'GnBu', 'Recipient Riskgroup') +
#  my_theme+
#  theme(axis.title = element_text(size = 7))

#plt_3b <- ggplot(shcs_empirical_results %>% filter(transmitterallocation == 'ML'))+
#  geom_point(aes(y = SpVL_recipient, x = 1-p_particles_1, colour = riskgroup_recipient),size = 1, shape = 4)+
#  scale_x_continuous(name = 'P(Multiple Particle Recipient)',
#                     expand = c(0.02,0.02),
#                     limits = c(0,0.4),
#                     breaks = seq(0, 0.4, by = 0.1))+
#  scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
#                limits = c(10**2, 10**7),
#                expand = c(0.05,0),
#                breaks = trans_breaks("log10", function(x) 10**x),
 # scale_color_brewer(palette = 'GnBu') +
##                labels = trans_format("log10", math_format(.x))) +
 # theme(axis.title = element_text(size = 7))

###  my_theme+
#plt_3c <- ggplot(shcs_empirical_results %>% filter(transmitterallocation == 'ML'))+
 # geom_point(aes(y = delta_CD4_recipient, x = 1-p_variants_1, colour = riskgroup_recipient),size = 1, shape = 4)+
 # scale_x_continuous(name = 'P(Multiple Variant Recipient)',
  #                   expand = c(0.02,0.02),
 #                    limits = c(0,0.4),
  #                   breaks = seq(0, 0.4, by = 0.1))+
 # scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)),  #
  #                   expand = c(0,0),
  #                   limits = c(-1,0)) +
  #scale_color_brewer(palette = 'GnBu') +
 # my_theme +
  #theme(axis.title = element_text(size = 7))

#plt_3d <- ggplot(shcs_empirical_results %>% filter(transmitterallocation == 'ML'))+
#  geom_point(aes(y = delta_CD4_recipient, x = 1-p_particles_1, colour = riskgroup_recipient),size = 1, shape = 4)+
 # scale_x_continuous(name = 'P(Multiple Particle Recipient)',
 #                    expand = c(0.02,0.02),
 #                    limits = c(0,0.4),
 #                    breaks = seq(0, 0.4, by = 0.1))+
 # scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)),  #
 #                    expand = c(0,0),
 #                    limits = c(-1,0)) +
 # scale_color_brewer(palette = 'GnBu') +
 # my_theme +
  #theme(axis.title = element_text(size = 7))


#plt_3e <- ggplot(modelresults %>% filter(transmitterallocation == 'ML'), 
 #                aes(y = SpVL_recipient, x = 1-p_variants_1))+
 # geom_density_2d_filled()+
  #scale_fill_brewer(palette = 'GnBu', direction = 1)+
  #scale_fill_manual(values = my_palette) +
 # scale_x_continuous(name = 'P(Multiple Variant Recipient)',
  #                   expand = c(0.02,0.02),
  #                   limits = c(0,0.4),
  #                   breaks = seq(0, 0.4, by = 0.1))+
 # scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
  #              limits = c(10**1, 10**7.5),
  #              expand = c(0.05,0),
  #              breaks = trans_breaks("log10", function(x) 10**x),
   #             labels = trans_format("log10", math_format(.x))) +
  #my_theme+
  #theme(axis.title = element_text(size = 7))

#plt_3f <- ggplot(shcs_predicted_results %>% filter(transmitterallocation == 'ML'),
#                 aes(y = SpVL_recipient, x = 1-p_particles_1))+
 # geom_density_2d_filled()+
  #scale_fill_brewer(palette = 'GnBu', direction = 1)+
#  scale_fill_manual(values = my_palette) +
 # scale_x_continuous(name = 'P(Multiple Particle Recipient)',
 #                    expand = c(0.02,0.02),
 #                    limits = c(0,0.4),
 #                    breaks = seq(0, 0.4, by = 0.1))+
 # scale_y_log10(name = expression(paste("Recipient SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
 #               limits = c(10**1, 10**7.5),
 #               expand = c(0.05,0),
#                breaks = trans_breaks("log10", function(x) 10**x),
#                labels = trans_format("log10", math_format(.x))) +
#  my_theme +
#  theme(axis.title = element_text(size = 7))


#plt_3g <- ggplot(shcs_predicted_results %>% filter(transmitterallocation == 'ML'),
   #              aes(y = delta_CD4_recipient, x = 1-p_variants_1))+
 # geom_density_2d_filled()+
  #scale_fill_brewer(palette = 'GnBu', direction = 1)+
  #scale_fill_manual(values = my_palette) +
 # scale_x_continuous(name = 'P(Multiple Variant Recipient)',
     #                expand = c(0.02,0.02),
    #                 limits = c(0,0.4),
  #                   breaks = seq(0, 0.4, by = 0.1))+
  #scale_y_continuous(limits = c(-0.5,0), 
    #                 expand = c(0,0.01), 
    #                 breaks = seq(-0.5, 0, by = 0.1),
   #                  name =expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1))) +
  #my_theme+
  #theme(axis.title = element_text(size = 7))


#plt_3h <- ggplot(shcs_predicted_results %>% filter(transmitterallocation == 'ML'),
    #             aes(y = delta_CD4_recipient, x = 1-p_particles_1))+
  #geom_density_2d_filled()+
  #scale_fill_brewer(palette = 'GnBu', direction = 1)+
  #scale_fill_manual(values = my_palette) +
  #scale_x_continuous(name = 'P(Multiple Particle Recipient)',
      #               expand = c(0.02,0.02),
      #               limits = c(0,0.4),
      #               breaks = seq(0, 0.4, by = 0.1))+
  #scale_y_continuous(limits = c(-0.5,0), 
   #                  expand = c(0,0.01), 
  #                   breaks = seq(-0.5, 0, by = 0.1),
   #                  name =expression(paste(Delta, CD4,'+ ' , mu, l**-1, ' ', day**-1))) +
 # my_theme+
  #theme(axis.title = element_text(size = 7))

#panel3_top <- plot_grid(plt_3a, plt_3b, plt_3c,plt_3d, nrow= 1, labels= 'AUTO', label_size = 9, align = 'hv')
#panel3_legend <- get_legend(plt_3a + theme(legend.position = 'bottom')) 
##panel3_bottom <- plot_grid(plt_3e, plt_3f, plt_3g, plt_3h, nrow= 1, labels= c('E', 'F', 'G', 'H'), label_size = 9, align = 'hv')

#panel_3 <- plot_grid(panel3_legend, panel3_top,  panel3_bottom, nrow= 3, labels= c('', '', '','E', 'F', 'G', 'H'), label_size = 9, rel_heights =  c(0.1,1,1))

#ggsave(plot = panel_3 , filename = paste(figs_dir,sep = '/', "panel_3.jpeg"), 
 #      device = jpeg,  width = 200, height = 130,  units = 'mm')
#
## END ##

