# Change in distribution of SpVL over transmission cycle
MC <- 4
VC <- 0.5
mu_0 <- 4.46  # Gaussian approximation of Transmission Potential
v_0 <- 1  # Gaussian approximation of Transmission Potential
mu_e <- 3  # Distribution of environmental effects 
v_e <- 1  # Distribution of environmental effects 
v_t <- 0.15 # variance asssociated with sampling at transmission

# Subtract environmental effects to calculate genotypic distribution
m_c <- MC - mu_e
v_c <- VC - v_e

# Apply transmission potential to infer donor genotype
m_d <- (m_c*(v_e + v_0) + (mu_0 - mu_e)*v_c)/(v_c + v_e + v_0)
v_d <- (v_c*(v_e + v_0))/(v_c + v_e + v_0)

# Infer donor phenotype
MD <- (MC*v_0 + mu_0*VC)/(v_0 + VC)
VD <- (VC*v_0)/(v_0 + VC)

# sample from donor genotypes for recipient genotypes
m_r <- m_d 
v_r <- v_d + v_t

# Re-distribute environmental variance for recipient phenotype
MR <- m_r + mu_e
VR <- v_r + v_e
s1_data <- tibble(Carrier = rnorm(100000, MC, VC), 
                  Transmitter = rnorm(100000, MD, sqrt(VD)),
                  Recipient = rnorm(100000, MR, sqrt(VR))) %>%
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
  facet_wrap(.~stage, nrow = 3)+
  my_theme

ggsave(plot = plt_s1, filename = paste(figs_dir,sep = '/', "plt_s1.jpeg"), device = jpeg, width = 6, height =8)


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
  geom_col(aes(x = variable, y = bias, fill = bias))+ 
  facet_wrap(.~riskgroup) +
  my_theme + 
  scale_fill_distiller(palette = 'OrRd') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(expression(paste(bar(x)["simulated"] - bar(x)['empirical'])), 
                     limits = c(-0.3, 0.3), 
                     breaks = seq(-0.3, 0.3, by = 0.1) %>% round(digits  = 2))


plt_s2b <- ggplot(bias_covmats, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  geom_text(aes(label = signif(value, 3)), color = "black", size = 3) +
  scale_fill_distiller(palette = 'OrRd') + 
  facet_wrap(.~riskgroup)+ 
  my_theme+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank())

plt_s2 <- cowplot::plot_grid(plt_s2a, plt_s2b, nrow = 2, ncol = 1,align = 'hv', labels = 'AUTO')
ggsave(plot = plt_s2 , filename = paste(figs_dir,sep = '/', "panel_s2.jpeg"), device = jpeg, width = 15, height = 14)


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

# Plot Prior ~ Posterior distributions
h2_df <- ggs(heritability_model_randomallocation)


# Plot posterior chains
plot_df <- h2_df %>%
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

plot_chains <- ggplot(plot_df,
                      aes(x   = Iteration,
                          y   = value, 
                          col = as.factor(Chain)))+
  geom_line()+
  facet_grid(Parameter ~ .,
             scale  = 'free_y',
             switch = 'y',
             labeller = label_parsed)+
  scale_colour_brewer(palette = 'OrRd', 'Chains') +
  my_theme + 
  theme(legend.position = 'bottom', 
        strip.text.y = element_text(family = 'sans'))

ggsave(plot = plot_chains, filename = paste(figs_dir,sep = '/', "panel_chains.jpeg"), device = jpeg, width = 12, height = 19)

plot_density <- ggplot(plot_df,
                      aes(x = value,
                          fill = as.factor(Chain)))+
  geom_density(alpha = 0.5, colour = NA)+
  facet_wrap(Parameter ~ .,
             scale  = 'free',
             switch = 'y',
             labeller = label_parsed, ncol = 2)+
  scale_fill_brewer(palette = 'OrRd', 'Chains') +
  my_theme + 
  theme(legend.position = 'bottom', 
        strip.text.y = element_text(family = 'sans'))

ggsave(plot = plot_density , filename = paste(figs_dir,sep = '/', "panel_density.jpeg"), device = jpeg, width = 14, height = 16)
