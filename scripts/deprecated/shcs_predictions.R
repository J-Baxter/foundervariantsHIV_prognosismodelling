sigma_0 <- 0.3198
b_0 <- 4.141712
SpVL = b_0 + rnorm(1,0,0.3198)*sigma_0  + -0.037020*partner2 + -0.052022*sexF + 0.004622*age.inf + 
  0.026272*riskgroupMSM + 0.011486*riskgroupPWID + 0.025593*riskgroupOTHER + -0.314611*riskgroupUNKNOWN

# Partner 1, male, 35, MSM
b_0 + rnorm(1,0,0.3198)*sigma_0  + -0.037020*0 + -0.052022*0 + 0.004622*35 + 
  0.026272*1 + 0.011486*0 + 0.025593*0 + -0.314611*0



sim_intercepts <- b_0 + rnorm(196,0,0.3198)*sigma_0
empirical_intercepts <- coef(adjusted_model)$log10_SpVL_couplemean %>% select(`(Intercept)`) %>% unlist()

intercepts <- data.frame(sim = sim_intercepts, empirical = empirical_intercepts) %>% pivot_longer(cols = everything(), names_to = 'origin', values_to = 'intercept')

ggplot(intercepts, aes(x = origin, y = intercept)) +
  geom_boxplot()+
  stat_compare_means(method = "t.test", label.x = 1.5, label.y = 6) + 
  theme_classic(base_family = "Questrial")+
  theme(
    text = element_text(size=16),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines")
  )


ModelAsFunction <- function(x){
  sigma_1 <-0.3198
  sigma_2 <- 0.6133 #? overall regression coefficient, or the slope, between the dependent variable and the Level 2 predictor
  gamma <- 4.141712
  delta <- #? verall regression coefficient, or the slope, between the dependent variable and the Level 1 predictor.
    
  v_i <- rnorm(1, 0, sigma_0)
  t_ij <- rnorm(1, 0, )

  b_0 <- 4.141712

  x['Pred_SpVL'] <- b_0 + rnorm(1, 0,sigma_0)*sigma_0  + 
    -0.037020*x['partner_2'] + 
    -0.052022*x['sex_F'] +
    0.004622*x['age.inf'] + 
    0.026272*x['riskgroup_MSM'] + 0.011486*x['riskgroup_PWID'] + 0.025593*x['riskgroup_OTHER'] + -0.314611*x['riskgroup_UNKNOWN'] +
    rnorm(1, 0,sigma_resid)
  
  return(x)
}

shcs_data_onehot <- shcs_data_long %>% fastDummies::dummy_cols(select_columns = c('partner', 'sex', 'riskgroup'),
                                                               remove_first_dummy = T, 
                                                               remove_selected_columns = T)

shcs_data_preds <- shcs_data_onehot %>%
  rowwise() %>% 
  ModelAsFunction()


SpVL = b_0 + rnorm(1,0,0.3198)*sigma_0  + -0.037020*partner2 + -0.052022*sexF + 0.004622*age.inf + 
  0.026272*riskgroupMSM + 0.011486*riskgroupPWID + 0.025593*riskgroupOTHER + -0.314611*riskgroupUNKNOWN

cov(shcs_data$log10_SpVL_1, shcs_data$log10_SpVL_2)/var(shcs_data$log10_SpVL_1)