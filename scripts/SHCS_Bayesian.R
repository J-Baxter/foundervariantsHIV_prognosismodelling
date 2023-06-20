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




b2_model <- brm(bf (log10_SpVL ~  + partner + sex + age.inf + riskgroup + (1|log10_SpVL_couplemean), 
                sigma ~ (1|log10_SpVL_couplemean)), vg
               data = shcs_data_long,
               chains = 4,
               iter = 8000,
               cores = 4,
               control = list(adapt_delta = 0.9)
               )

ppe <- posterior_predict(b_model)
str(ppe)
ppe_simdata <- posterior_predict(b2_model, newdata = sim_data, allow_new_levels = TRUE, sample_new_levels = "gaussian", re_formula = (1|log10_SpVL_couplemean))
str(ppe_simdata)

fit_params <- MASS::fitdistr(shcs_data_long$age.inf,"lognormal")

tidy_pred <- b2_model %>%
  epred_draws( newdata = shcs_data_, allow_new_levels = TRUE, sample_new_levels = "gaussian", re_formula = NULL) %>% 
  group_by(log10_SpVL_couplemean, partner) %>%
  summarise(mean.pred = mean(.epred))



shcs_data_x <- shcs_data_long %>% mutate(log10_SpVL_couplemean = sort(log10_SpVL_couplemean))

