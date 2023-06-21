# SHCS Bayesian
library(brms)

model_priors <- c(set_prior(coef = 'age.inf', 
                            class = 'b', 
                            prior = 'normal(0, 3)'), # Relatively uninformative
                  
                  set_prior(coef = 'partner2', 
                            class = 'b', 
                            prior = "normal(0.42, 2)"), # from Hollingsworth et al. 2010
                  
                  set_prior(coef = 'riskgroupMSM', 
                            class = 'b',
                            prior = "normal(0, 3)"), # Relatively uninformative
                  
                  set_prior(coef = 'riskgroupPWID',
                            class = 'b', 
                            prior = "normal(0,3)"), # Relatively uninformative
                  
                  set_prior(coef = 'riskgroupOTHER', 
                            class = 'b', 
                            prior = "normal(0, 3)"), # Relatively uninformative
                  
                  set_prior(coef = 'riskgroupUNKNOWN', 
                            class = 'b', 
                            prior = "normal(0, 3)"), # Relatively uninformative
                  
                  set_prior(coef = 'sexF', 
                            class = 'b', 
                            prior = "normal(0.15, 2)"), # from Hollingsworth et al. 2010
                  
                  set_prior(coef = '', # left intentionally blank
                            class = 'Intercept', 
                            prior = "normal(4, 3)"), # From mean of log10_SpVL_couplemean
                  
                  set_prior(coef = '', # left intentionally blank
                            class = 'sd', 
                            prior = "normal(0, 6)"), # Relatively uninformative
                
                 # set_prior(group = 'log10_SpVL_couplemean', 
                           # class = "sd",
                           # prior = "normal(0, 6)"), # Relatively uninformative
                  
                  set_prior(coef = 'Intercept',
                            class = "sd",
                            group = 'log10_SpVL_couplemean', 
                            prior = "normal(0, 6)"), # Relatively uninformative
                  
                  set_prior(coef = '', # left intentionally blank
                            class = 'sigma',
                            prior = "normal(0, 6)") # Relatively uninformative
                  )

b2_model <- brm(log10_SpVL ~  + partner + sex + age.inf + riskgroup + (1|log10_SpVL_couplemean), 
                prior = model_priors,
                data = shcs_data_long,
                chains = 4,
                iter = 4000,
                cores = 4,
                control = list(adapt_delta = 0.95)
                )

tidy_pred <- b2_model %>%
  epred_draws( newdata = shcs_data_, allow_new_levels = TRUE, sample_new_levels = "gaussian", re_formula = NULL) %>% 
  group_by(log10_SpVL_couplemean, partner) %>%
  summarise(mean.pred = mean(.epred))
