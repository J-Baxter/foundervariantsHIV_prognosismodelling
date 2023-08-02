# SHCS frequentist
# Linear Models

# Null Model: one coefficient, the overall mean for all individuals
null_model <- lm(log10_SpVL ~ log10_SpVL_cohortmean, data = shcs_data_long)


# Unadjusted Model: each viral load is predicted by the coefficient for the couple
unadjusted_model <- lmer(log10_SpVL ~ (1|log10_SpVL_couplemean), data = shcs_data_long) 


# Adjusted Model - No couple effect
adjustednocouple_model <- lm(log10_SpVL ~ 1 + partner + sex + age.inf + riskgroup, data = shcs_data_long) 


# Adjusted Model
# Baseline covariates: partner 1, sex M, riskgroup HET
adjusted_model <- lmer(log10_SpVL ~ (1|log10_SpVL_couplemean) + partner + sex + age.inf + riskgroup, data = shcs_data_long) 
adjusted_model <- lmer(log10_SpVL_1 ~ (log10_SpVL_2|ID_pair) +  sex_1 + age.inf_1 + riskgroup_1, data = shcs_data) 


adjusted_model %>%
  tidy("fixed") %>% 
  mutate(p_value = 2*(1 - pnorm(abs(statistic))))

adjusted_model.performance <- performance(adjusted_model)

adjusted_model.fixed <- adjusted_model %>%
  tidy("fixed") %>% 
  mutate(p_value = 2*(1 - pnorm(abs(statistic))))

adjusted_model.ranef <- adjusted_model %>%
  tidy("ran_coefs")




# AOV (Log10SpVL_CoupleMean ~ Riskgroup and Sex)




# Q: not coming up with adjusted/unadjustd R2 - suspect this is lmm vs lm. Not clear 
# how couple intercept acheived without this.

# Simulated data
shcs_data_int <- shcs_data_long %>%
  select(-starts_with('ID'), -starts_with('SpVL'), -contains('cohortmean'), -log10_SpVL) %>%
  mutate(across(where(is.factor), .fns = ~ match(.x, unique(.x))))

cum_probs <- shcs_data_int %>% 
  select(where(is.integer)) %>% 
  apply(.,2, function(x) cumsum(table(x))/length(x))

sim_data <- mvtnorm::rmvnorm(n = 100, mean = as.vector(colMeans(shcs_data_int)), sigma = cov(shcs_data_int)) %>%
  as_tibble(name_repair = NULL) %>%
  setNames(colnames(shcs_data_int)) %>%
  mutate(across(c(partner,riskgroup,sex), 
                .fns = ~cut(.x,
                            right = T,
                            include.lowest = T,
                            breaks = as.vector(c(min(.x), quantile(.x,cum_probs[[cur_column()]]))),
                            labels = 1:length(cum_probs[[cur_column()]])))) %>%
  mutate(sex = case_when(sex == 1 ~ 'M',
                         sex == 2 ~ 'F')) %>%
  mutate(riskgroup = case_when(riskgroup == 1 ~ "HET",
                               riskgroup == 2 ~ "MSM",
                               riskgroup == 3 ~ "PWID",
                               riskgroup == 4 ~ "UNKNOWN",
                               riskgroup == 5 ~ "OTHER")) %>%
  mutate(across(where(is.character), .fns = ~as.factor(.x))) %>%
  mutate(sex = relevel(sex, 'M')) %>%
  mutate(riskgroup = factor(riskgroup, levels = c('HET', 'MSM', 'PWID', 'OTHER', 'UNKNOWN'))) %>%
  filter(min(shcs_data_long$log10_SpVL_couplemean)<log10_SpVL_couplemean && log10_SpVL_couplemean<max(shcs_data_long$log10_SpVL_couplemean))

preds <-predict(adjusted_model, newdata = sim_data, re.form = ~(1|log10_SpVL_couplemean)) #newlevels detected in data - current suspicion is random effects?

