freq_model <- lme4::lmer(log10_SpVL ~  + partner + sex + age.inf_category + riskgroup + (1|log10_SpVL_couplemean), 
                         data = shcs_data_long)


predicted_df_tidy <- mutate(shcs_data_long,
                            .pred = predict(freq_model)) %>%
  data.table() %>%
  mltools::one_hot(cols = c('age.inf_category', 'partner', 'sex', 'riskgroup')) %>%
  as_tibble() %>%
  select(contains(c('partner_2', 'sex_F', 'riskgroup', 'log10_SpVL_couplemean', 'age.inf_category' , "ID_pair", '.pred'))) %>%
  select(-c('age.inf_category_15-24', 'riskgroup_HET')) %>%
  mutate(partner_2 = case_when(partner_2 == 1 ~ -0.03912252,
                               .default = 0)) %>%
  mutate(sex_F = case_when(sex_F == 1 ~ -0.03394273,
                               .default = 0)) %>%
  mutate(`age.inf_category_25-29` = case_when(`age.inf_category_25-29` == 1 ~ 0.09351565,
                               .default = 0)) %>%
  mutate(`age.inf_category_30-39` = case_when(`age.inf_category_30-39` == 1 ~ 0.1278834,
                                           .default = 0)) %>%
  mutate(`age.inf_category_40-80` = case_when(`age.inf_category_40-80` == 1 ~ 0.1901484,
                                           .default = 0)) %>%
  mutate(riskgroup_MSM = case_when(riskgroup_MSM == 1 ~ 0.02510452,
                                           .default = 0)) %>%
  mutate(riskgroup_PWID = case_when(riskgroup_PWID == 1 ~ 0.02278369,
                                   .default = 0)) %>%
  mutate(riskgroup_OTHER = case_when(riskgroup_OTHER == 1 ~ 0.03121261,
                                   .default = 0)) %>%
  mutate(riskgroup_UNKNOWN = case_when(riskgroup_UNKNOWN == 1 ~ -0.2791359,
                                   .default = 0)) %>%
  mutate(slope = rowSums(select(., -c(ID_pair, log10_SpVL_couplemean, .pred)))) %>%
  arrange(log10_SpVL_couplemean) %>%
  mutate(intercept = rep(coef(freq_model)$log10_SpVL_couplemean[,1], each = 2))
  

ggplot(predicted_df_tidy)+
  geom_point(aes(x = log10_SpVL_couplemean, y = log10_SpVL), data = shcs_data_long, colour = '#d7301f', shape = 4) + 
  geom_abline(aes(intercept =intercept, slope = slope), alpha = 0.1, colour = '#fdbb84') +
  xlim(0, 7) +
  ylim(0, 7) + theme(legend.position = 'none') + my_theme