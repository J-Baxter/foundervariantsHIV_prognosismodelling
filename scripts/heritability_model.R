################################### Heritability Model ###################################
# Set priors and fit a Bayesian linear mixed model to estimate the heritability
# of the SHCS data. 
# Uses brms package which acts as an interface to Stan


################################### Dependencies ###################################
# (Should already by loaded if ./scripts/dependencies.R has been run)
require(brms)
require(tidybayes)
require(performance)


################################### Set Model Priors ###################################
# Partner2, SexF and age categories informed from Hollingsworth et al 2010 linear model 
# with random intercept for mean pair log10 SpVL. Mean intercept informed from summary 
# statistics of the SHCS data. 

# Link to Hollingsworth et al.:
# https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1000876

model_priors <- c( set_prior(coef = 'age.inf_category25M29', 
                            class = 'b', 
                            prior = 'normal(0.32, 2)'), # from Hollingsworth et al. 2010
                  
                  set_prior(coef = 'age.inf_category30M39', 
                            class = 'b', 
                            prior = 'normal(0.35, 2)'), # from Hollingsworth et al. 2010
                  
                  set_prior(coef = 'age.inf_category40M80', 
                            class = 'b', 
                            prior = 'normal(0.73, 2)'), # from Hollingsworth et al. 2010
                  
                  set_prior(coef = 'partner2', 
                            class = 'b', 
                            prior = "normal(0.42, 2)"), # from Hollingsworth et al. 2010
                
                  set_prior(coef = 'riskgroupMSM', 
                            class = 'b',
                            prior = "normal(0, 3)"), # Relatively uninformative
                  
                  set_prior(coef = 'riskgroupPWID',
                            class = 'b', 
                            prior = "normal(0,4)"), # Relatively uninformative
                  
                  set_prior(coef = 'riskgroupOTHER', 
                            class = 'b', 
                            prior = "normal(0, 4)"), # Relatively uninformative
                  
                  set_prior(coef = 'riskgroupUNKNOWN', 
                            class = 'b', 
                            prior = "normal(0, 4)"), # Relatively uninformative
                  
                  set_prior(coef = 'sexF', 
                            class = 'b', 
                            prior = "normal(0.15, 2)"), # from Hollingsworth et al. 2010
                  
                  set_prior(coef = '', # left intentionally blank
                            class = 'Intercept', 
                            prior = "normal(4, 3)"), # From mean of log10_SpVL_couplemean
                  
                  set_prior(coef = '', # left intentionally blank
                            class = 'sd', 
                            prior = "normal(0, 4)"), # Relatively uninformative
                
                 # set_prior(group = 'log10_SpVL_couplemean', 
                           # class = "sd",
                           # prior = "normal(0, 6)"), # Relatively uninformative
                  
                  set_prior(coef = 'Intercept',
                            class = "sd",
                            group = 'log10_SpVL_couplemean', 
                            prior = "normal(0, 4)"), # Relatively uninformative
                  
                  set_prior(coef = '', # left intentionally blank
                            class = 'sigma',
                            prior = "normal(0, 4)") # Relatively uninformative
                  )


################################### Fit Model ###################################
# Fit model to SHCS data (long format). This 'should' run relatively quickly.
# NB: User should set 'cores' parameter as appropriate for their device: recommend 
# 1 < cores.used ≤ available.cores - 1 (No point of exceeding cores.used  = 4)

# ESS (the number of effectively independent draws from the posterior distribution
# that the Markov chain is equivalent to) should be at least 200. tail-ESS represents
#  effective sample sizes for 5% and 95% quantiles. (our run all bulk ESS > 28000)

# R-hat (compares the between- and within-chain estimates for model parameters) should
# be less than 1.05 (our run == 1)

heritability_model <- brm(log10_SpVL ~  + partner + sex + age.inf_category + riskgroup + (1|log10_SpVL_couplemean), 
                          prior = model_priors,
                          data = shcs_data_long,
                          chains = 4,
                          iter = 10000,
                          warmup = 1000, # 10% burn in 
                          cores = 4,
                          control = list(adapt_delta = 0.95)
                          )


################################### Misc Evaluation Code ###################################
# performance(heritability_model) # tibble output of model metrics including R2, ELPD, LOOIC, RMSE

# plot(heritability_model) # default output plot of brms showing posterior distributions of
                           # parameters and MCMC chains

# prior_summary(b2_model) obtain dataframe of priors used in model.
                
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

# Plot Prior ~ Posterior distributions

# Plot Joint posterior chain



################################### END ###################################