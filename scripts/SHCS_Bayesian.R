# SHCS Bayesian
library(brms)

model_priors <- c(set_prior(coef = 'age.inf', group = 'log10_SpVL_couplemean', class = 'sd', prior = 'normal(0.005, 1)'), #fit lognormal to age dat
                  set_prior(coef = 'partner2', group = 'log10_SpVL_couplemean', class = 'sd', prior = "normal(0.037, 1)"), # binomial p = 0.5
                  set_prior(coef = 'riskgroupHET', group = 'log10_SpVL_couplemean', class = 'sd', prior = "normal(-0.011,1)"),
                  set_prior(coef = 'riskgroupMSM', group = 'log10_SpVL_couplemean', class = 'sd', prior = "normal(0.014,1)"),
                  set_prior(coef = 'riskgroupOTHER', group = 'log10_SpVL_couplemean', class = 'sd', prior = "normal(0.014, 1)"),
                  set_prior(coef = 'riskgroupUNKNOWN', group = 'log10_SpVL_couplemean', class = 'sd', prior = "normal(-0.326,1)"),
                  set_prior(coef = 'sexF', group = 'log10_SpVL_couplemean', class = 'sd', prior = "normal(-0.052, 1)"),
                  set_prior( group = 'log10_SpVL_couplemean', class = "cor", prior = "lkj(2)"),
                  set_prior(class = 'Intercept', prior = "normal(4.43, 0.537)") #mean and sd of log10_SpVL_couplemean
                  )

p <- get_prior(formula = log10_SpVL ~  + partner + sex + age.inf + riskgroup + 1|log10_SpVL_couplemean, 
          data = shcs_data_long, family = gaussian())
b_model <- brm(formula = log10_SpVL ~  + partner + sex + age.inf + riskgroup + 1|log10_SpVL_couplemean, 
               data = shcs_data_long,
               prior = NULL,
               chains = 4,
               iter = 8000,
               cores = 4
               )

fit_params <- MASS::fitdistr(shcs_data_long$age.inf,"lognormal")