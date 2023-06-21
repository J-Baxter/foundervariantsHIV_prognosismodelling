# Main Analysis for "RECONCILING THE FOUNDER VARIANT MULTIPLICITY OF HIV-1 INFECTION WITH THE RATE
# OF CD4+ DECLINE"
# J Baxter

# This script combines the functions for the analysis of paired SpVL distributions, occuring in
# three main stages:

# ALL VIRAL LOADS ARE REAL NUMBERS UNLESS EXPLICITLY STATED


################################### Dependencies ###################################
source('./scripts/dependencies.R')
source('./scripts/populationdata_acrossVL_models.R')
source('./scripts/SimPop.R')
source('./scripts/global_theme.R')
source('./scripts/postprocessingfuncs.R')


# Set Seed
set.seed(4472)
options(scipen = 100) #options(scipen = 100, digits = 4)


################################### Import Data ###################################
# Paired set point viral load (SpVL), age, riskgroup and sex data provided on request
# by the Swiss HIV Cohort Study. These data were used in the analysis of Bertels et. al
# 2018 https://academic.oup.com/mbe/article/35/1/27/4210012

# NB: partner number is allocated randomly and does not assign transmitter/recipient status
shcs_data <-read_csv("data/shcs_data.csv", 
                     col_types = cols(sex.1 = readr::col_factor(levels = c("M", "F")), 
                                      sex.2 = readr::col_factor(levels = c("M", "F")), 
                                      riskgroup.1 = readr::col_factor(levels = c("HET", "MSM",
                                                                                 "IDU", "OTHER",
                                                                                 "UNKNOWN")), 
                                      riskgroup.2 = readr::col_factor(levels = c("HET", "MSM", 
                                                                                 "IDU", "OTHER", 
                                                                                 "UNKNOWN"))))%>%
  select(-contains('ID')) %>%
  rowid_to_column( "ID.pair") %>%
  mutate(across(contains('riskgroup'), 
                .fns = ~ fct_recode(.x, PWID = "IDU"))) %>%
  rename_with( ~ stri_replace_last_fixed(.x, 'spVL', 'SpVL')) %>%
  mutate(across(contains('SpVL'), 
                .fns = ~ raise_to_power(10, .x))) %>%
  mutate(SpVL.couplemean = rowMeans(across(contains('SpVL')))) %>%
  mutate(across(contains('SpVL'),
                .fns = ~log10(.x), .names = "log10_{.col}")) %>%
  rename_with( ~ stri_replace_last_fixed(.x, '.', '_')) 


# Transform to long-format data table
shcs_data_long <- shcs_data %>% 
  pivot_longer(cols = -c(ID_pair, contains('couplemean')), 
               names_to = c(".value", "partner"),
               names_pattern  = "^(.*)_([0-9])$") %>% #matches('SpVL_[[:digit:]]')
  mutate(across(ends_with('SpVL'),
                .fns = ~ mean(.x), .names = "{.col}_cohortmean")) %>%
  mutate(partner = factor(partner, levels = c("1", "2"))) %>%
  relocate(ID_pair) %>%
  relocate(age.inf, .after = sex) %>%
  relocate(ends_with('SpVL'), .after = riskgroup) %>%
  relocate(ends_with('couplemean'), .after = log10_SpVL) 


################################### Fit/Import Base Models ###################################
# Import fitted model
source('./scripts/transmission_model.R')
source('./scripts/tolerance_model.R')

# Fit Bayesian linear mixed model to SHCS data
source('./scripts/heritability_model.R')
#source('./scripts/heritability_extra_model.R') # vary assumptions of heritability


################################### Simulate Stratified Cohorts ###################################
# Simulate riskgroup stratified cohorts according to covariance matrices from SHCS, calculated from
# grouped data. Covariance and bias analysis at ./scripts/SHCS_simulation_misc.R

# Filter data and enumerate remaining categorical levels
shcs_data_int_list <- shcs_data %>%
  rowwise() %>%
  filter(riskgroup_1 == riskgroup_2) %>% # Only include pairs where both individuals share the same riskgroup
  filter(!if_any(starts_with('riskgroup'), ~ . %in% c('UNKNOWN', 'OTHER'))) %>%
  EnumerateLevels(.) #NB will drop any column with only 1 level


# Infer cumulative probabilities for each categorical level
cum_probs_list <- lapply(shcs_data_int_list, function(x) x %>%
                           select(where(is.integer)) %>% 
                           apply(2, function(i) cumsum(table(i)) / length(i), simplify = FALSE))


# Run SimCohorts - see ./scripts/SHCS_simulation.R for details
sim_data <- mapply(SimCohorts,
                   data = shcs_data_int_list,
                   probs = cum_probs_list,
                   SIMPLIFY = F) %>%
  setNames(., c('HET', 'MSM', 'PWID')) %>%
  
  # Re-organising simulated data
  bind_rows(., .id = "riskgroup") %>% 
  rowid_to_column( "ID_pair") %>%
  pivot_longer(cols = - c(contains('couplemean'), riskgroup, ID_pair), 
               names_to = c(".value", "partner"), 
               names_pattern  = "^(.*)_([0-9])$") %>% 
  mutate(partner = factor(partner, levels = c("1", "2"))) %>%
  relocate(c(partner, sex), .before = riskgroup) %>%
  relocate(ends_with('couplemean'), .after = last_col()) 


################################### Predict SpVLs ###################################
#Bind all datasets and label?

shcs_pred <- epred_draws(heritability_model,
                         newdata = shcs_data_long,
                         allow_new_levels = TRUE, # allow new random effects levels
                         sample_new_levels = "old_levels", #
                         re_formula = ~ (1|log10_SpVL_couplemean),
                         value = 'predicted_log10_SpVL') %>% # NULL retains the random effects formula in the model
  select(-c(.chain, .iteration, .draw)) %>%
  group_by(ID_pair, partner) %>%
  mutate(predicted_log10_SpVL = mean(predicted_log10_SpVL)) %>%
  select(-.row) %>%
  distinct()

segregated_pred <- epred_draws(heritability_model,
                              newdata = sim_data,
                              allow_new_levels = TRUE, # allow new random effects levels
                              sample_new_levels = "old_levels", #
                              re_formula = ~ (1|log10_SpVL_couplemean),
                              ndraws = 4000,
                              value = 'predicted_log10_SpVL') %>% # NULL retains the random effects formula in the model
  select(-c(.chain, .iteration, .draw)) %>%
  group_by(ID_pair, partner) %>%
  mutate(predicted_log10_SpVL = mean(predicted_log10_SpVL)) %>%
  select(-.row) %>%
  distinct()


################################### Predict CD4 decline ###################################
shcs_pred$delta_CD4 <- ToleranceModel(shcs_pred$log10_SpVL_couplemean, shcs_pred$age.inf, shcs_pred$sex)

segregated_pred$delta_CD4 <- ToleranceModel(segregated_pred$log10_SpVL_couplemean, segregated_pred$age.inf, segregated_pred$sex)



################################### Calculate Joint Probability Dist MV/MP ###################################



################################### Estimate Heritability Under different Assumptions ###############################
# Fit 
linear_model_uw <- heritability_model # Linear model imported from base_models

concave_model_uw <- lm(recipient_log10SpVL ~  exp(transmitter_log10SpVL), data = pop)

convex_model_uw <-lm(recipient_log10SpVL ~  log(transmitter_log10SpVL), data = pop) 

pop$p_mv <- shcs_pred_variants %>% filter(w == 1) %>% filter(variants == 1) %>% select(p) %>% mutate(1-p,.keep='none')%>%unlist()
linearweighted_model_uw <- lm(recipient_log10SpVL ~  transmitter_log10SpVL + p_mv, data = pop)


################################### Simulate Model Populations ###################################
# Generate donor population of SpVLs - sampled according to the probability that a given SpVL
# is present in the population and that at least one transmission event is observed.
sim_donor <- SimDonor(NSIM)

sim_recip_chars <- RecipChars(NSIM)

linear_uw_pop <- sim_donor %>% 
  cbind.data.frame(transmitter = sim_donor) %>% 
  mutate(transmitter_log10SpVL= log10(transmitter)) %>%
  
  # Bind simulated recipient characteristics
  cbind.data.frame(sim_recip_chars) %>%
  
  # Predict recipient SpVL according to heritability model
  predict(linear_model_uw, .) %>%
  cbind.data.frame(recipient_log10SpVL = .) %>%
  
  # Because A) we are predicting out-of-sample and B) our focus is the population level
  # we incorporate residual standard error into our predictions
  mutate(recipient_log10SpVL = recipient_log10SpVL + rnorm(NSIM, sd = (summary(linear_model_uw)$sigma)/10)) %>%
  
  # Bind predicted recipient SpVl with transmission pair characteristics
  cbind.data.frame(transmitter_log10SpVL= log10(sim_donor)) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SpVL')}")) %>%
  `colnames<-` (str_remove(colnames(.), 'sim_')) %>% 
  cbind.data.frame(sim_recip_chars) %>%
  mutate(model = 'linear_uw') 


concave_uw_pop <-  sim_donor %>% 
  cbind.data.frame(transmitter = sim_donor) %>% 
  mutate(transmitter_log10SpVL= log10(transmitter)) %>%
  
  # Bind simulated recipient characteristics
  cbind.data.frame(sim_recip_chars) %>%
  
  # Predict recipient SpVL according to heritability model
  predict(concave_model_uw, .) %>%
  cbind.data.frame(recipient_log10SpVL = .) %>%
  
  # Because A) we are predicting out-of-sample and B) our focus is the population level
  # we incorporate residual standard error into our predictions
  mutate(recipient_log10SpVL = recipient_log10SpVL + rnorm(NSIM, sd = summary(concave_model_uw)$sigma)/10) %>%
  
  # Bind predicted recipient SpVl with transmission pair characteristics
  cbind.data.frame(recipient_log10SpVL = ., transmitter_log10SpVL = log10(sim_donor)) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SpVL')}")) %>%
  `colnames<-` (str_remove(colnames(.), 'sim_'))  %>% 
  cbind.data.frame(sim_recip_chars) %>%
  mutate(model = 'concave_uw') 


convex_uw_pop <-  sim_donor %>% 
  cbind.data.frame(transmitter = sim_donor) %>% 
  mutate(transmitter_log10SpVL= log10(transmitter)) %>%
  
  # Bind simulated recipient characteristics
  cbind.data.frame(sim_recip_chars) %>%
  
  # Predict recipient SpVL according to heritability model
  predict(convex_model_uw, .) %>%
  cbind.data.frame(recipient_log10SpVL = .) %>%
  
  # Because A) we are predicting out-of-sample and B) our focus is the population level
  # we incorporate residual standard error into our predictions
  mutate(recipient_log10SpVL = recipient_log10SpVL + rnorm(NSIM, sd = summary(convex_model_uw)$sigma)/10) %>%
  
  # Bind predicted recipient SpVl with transmission pair characteristics
  cbind.data.frame(recipient_log10SpVL = ., transmitter_log10SpVL = log10(sim_donor)) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SpVL')}")) %>%
  `colnames<-` (str_remove(colnames(.), 'sim_')) %>% 
  cbind.data.frame(sim_recip_chars) %>%
  mutate(model = 'convex_uw') 

linear_w_pop <- sim_donor %>% 
  cbind.data.frame(transmitter = sim_donor) %>% 
  mutate(transmitter_log10SpVL= log10(transmitter)) %>%
  
  # Bind simulated recipient characteristics
  cbind.data.frame(sim_recip_chars) %>%
  
  # Predict recipient SpVL according to heritability model
  predict(linear_model_uw, .) %>%
  cbind.data.frame(recipient_log10SpVL = .) %>%
  
  # Because A) we are predicting out-of-sample and B) our focus is the population level
  # we incorporate residual standard error into our predictions
  mutate(recipient_log10SpVL = recipient_log10SpVL + rnorm(NSIM, sd = (summary(linear_model_uw)$sigma)/10)) %>%
  
  # Bind predicted recipient SpVl with transmission pair characteristics
  cbind.data.frame(transmitter_log10SpVL= log10(sim_donor)) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SpVL')}")) %>%
  `colnames<-` (str_remove(colnames(.), 'sim_')) %>% 
  cbind.data.frame(sim_recip_chars) %>%
  mutate(model = 'linear_w') 


################################### Run Transmission Models ##################################
linear_uw_tm <- c(RunParallel(populationmodel_acrossVL_Environment, linear_uw_pop$transmitter, w= 1),
                  RunParallel(populationmodel_acrossVL_Environment, linear_uw_pop$transmitter, w= 5),
                  RunParallel(populationmodel_acrossVL_Environment, linear_uw_pop$transmitter, w= 10),
                  RunParallel(populationmodel_acrossVL_Environment, linear_uw_pop$transmitter, w= 20)) %>%
  
  # Label
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w'))


concave_uw_tm <- c(RunParallel(populationmodel_acrossVL_Environment, concave_uw_pop$transmitter, w= 1),
                  RunParallel(populationmodel_acrossVL_Environment, concave_uw_pop$transmitter, w= 5),
                  RunParallel(populationmodel_acrossVL_Environment, concave_uw_pop$transmitter, w= 10),
                  RunParallel(populationmodel_acrossVL_Environment, concave_uw_pop$transmitter, w= 20)) %>%
  
  # Label
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w'))


convex_uw_tm <- c(RunParallel(populationmodel_acrossVL_Environment, convex_uw_pop$transmitter, w= 1),
                  RunParallel(populationmodel_acrossVL_Environment, convex_uw_pop$transmitter, w= 5),
                  RunParallel(populationmodel_acrossVL_Environment, convex_uw_pop$transmitter, w= 10),
                  RunParallel(populationmodel_acrossVL_Environment, convex_uw_pop$transmitter, w= 20)) %>%
  
  # Label
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w'))


linear_w_tm <- c(RunParallel(populationmodel_acrossVL_Environment, linear_w_pop$transmitter, w= 1),
                  RunParallel(populationmodel_acrossVL_Environment, linear_w_pop$transmitter, w= 5),
                  RunParallel(populationmodel_acrossVL_Environment, linear_w_pop$transmitter, w= 10),
                  RunParallel(populationmodel_acrossVL_Environment, linear_w_pop$transmitter, w= 20)) %>%
  
  # Label
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w'))


################################### Extract Model Outputs ##################################
# Q1: All else equal, what effect would we expect to observe between P(Multiple Variants
# and Recipient SpVL?

shcs_transmitters <- shcs_tm  %>%
  lapply(., cbind.data.frame) %>%
  do.call(rbind.data.frame,.) %>%
  rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
  dplyr::select(-contains('nparticles')) %>%
  dplyr::summarise(across(starts_with('V'), .fns =sum),.by = transmitter) %>% #by w too?
  pivot_longer(cols = starts_with('V'),
               names_to = 'variants',
               values_to = 'p') %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% 
           as.numeric())  %>%
  filter(variants <= 10) %>%
  mutate(w = 1)

shcs_pred_variants <- shcs_tm %>%
  VariantPP(pop = pop) 

shcs_pred__virions <- shcs_tm %>%
  VirionPP(pop = pop) 

linear_uw_variants <- linear_uw_tm %>%
  .[sapply(., function(x) x[['w']] == 1)] %>%
  VariantPP(pop = linear_uw_pop) %>% 
  mutate(model = 'linear_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 
  
linear_uw_virions <- linear_uw_tm %>% 
  .[sapply(., function(x) x[['w']] == 1)] %>%
  VirionPP(pop = linear_uw_pop) %>% 
  mutate(model = 'linear_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient))))


# Q2: Is the SpVL in recipient partner determined by a non-linear relationship with the SpVL in the
# transmitting partner, and how is this affected by the number of variants initiating infection in 
# the recipient partner?

concave_uw_variants <- concave_uw_tm  %>%
  .[sapply(., function(x) x[['w']] == 1)] %>%
  VariantPP(pop = concave_uw_pop) %>%
  mutate(model = 'concave_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 

concave_uw_virions <- concave_uw_tm %>%
  .[sapply(., function(x) x[['w']] == 1)] %>%
  VirionPP(pop = concave_uw_pop) %>%
  mutate(model = 'concave_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 


convex_uw_variants <- convex_uw_tm  %>%
  .[sapply(., function(x) x[['w']] == 1)] %>%
  VariantPP(pop = convex_uw_pop) %>%
  mutate(model = 'convex_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 

convex_uw_virions <- convex_uw_tm %>%
  .[sapply(., function(x) x[['w']] == 1)] %>%
  VirionPP(pop = convex_uw_pop) %>%
  mutate(model = 'convex_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 



linear_w_variants <- linear_w_tm %>%
  .[sapply(., function(x) x[['w']] == 1)] %>%
  VariantPP(pop = linear_w_pop) %>% 
  mutate(model = 'linear_w') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 

linear_w_virions <- linear_w_tm %>% 
  .[sapply(., function(x) x[['w']] == 1)] %>%
  VirionPP(pop = linear_w_pop) %>% 
  mutate(model = 'linear_w') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient))))


# Q3: How does the timing of transmission impact observations of the association between the number
# of founder variants and CD4+ T cell decline?

shcs_timing <- shcs_tm  %>%
  lapply(., cbind.data.frame) %>%
  do.call(rbind.data.frame,.) %>%
  rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
  dplyr::select(-contains('nparticles')) %>%
  dplyr::summarise(across(starts_with('V'), .fns =sum),.by = transmitter) %>%
  pivot_longer(cols = starts_with('V'),
               names_to = 'variants',
               values_to = 'p') %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% 
           as.numeric())

linear_uw_variants_timing <- linear_uw_tm %>%
  VariantPP(pop = linear_uw_pop) %>% 
  mutate(model = 'linear_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 

linear_uw_virions_timing <- linear_uw_tm %>% 
  VirionPP(pop = linear_uw_pop) %>% 
  mutate(model = 'linear_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient))))


# The following are not shown in the manuscript as the trend in timing is similar across 
# assumptions about the relationship between transmitter and recipient SpVL
concave_uw_variants_timing <- concave_uw_tm  %>%
  VariantPP(pop = concave_uw_pop) %>%
  mutate(model = 'concave_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 

concave_uw_virions_timing <- concave_uw_tm %>%
  VirionPP(pop = concave_uw_pop) %>%
  mutate(model = 'concave_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 


convex_uw_variants_timing <- convex_uw_tm  %>%
  VariantPP(pop = convex_uw_pop) %>%
  mutate(model = 'convex_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 

convex_uw_virions_timing <- convex_uw_tm %>%
  VirionPP(pop = convex_uw_pop) %>%
  mutate(model = 'convex_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 


linear_w_variants_timing <- linear_w_tm %>%
  VariantPP(pop = linear_w_pop) %>% 
  mutate(model = 'linear_w') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 

linear_w_virions_timing <- linear_w_tm %>% 
  VirionPP(pop = linear_w_pop) %>% 
  mutate(model = 'linear_w') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient))))


# Q4: Do we expect epidemic characteristics to impact the observation of the assocation between 
# multiple variant infection and CD4+ T Cell decline?

# Baseline population; characteristics sampled at random. Re-simulate with characteristics extracted
# from RV144 and HVTN502 (STEP) (representing a HSX and MSM epidemic)

rv144_chars
hvtn502_chars


################################### Sensitivity Analysis ################################### 
# Robustness of model system to proportion of standard error incorporated in Heritability
# Model (linear, unweighted only)

source('./scripts/sa_heritabilityse.R')



################################### Write to file ################################### 
source('./scripts/figures.R')

ggsave(plot = panel_1, filename = paste(figs_dir,sep = '/', "panel_1.jpeg"), device = jpeg, width = 14, height = 14) # Model Components - functional
ggsave(plot = panel_2, filename = paste(figs_dir,sep = '/', "panel_2.jpeg"), device = jpeg, width = 14, height = 18) # Confounder - functional 
ggsave(plot = panel_3, filename = paste(figs_dir,sep = '/', "panel_3.jpeg"), device = jpeg, width = 14, height = 18) # Non - Linear
ggsave(plot = panel_4, filename = paste(figs_dir,sep = '/', "panel_4.jpeg"), device = jpeg, width = 14, height = 18) # Timing of transmission

# Supplementary plots
source('./scripts/tm_withinhostprocesses.R') 
source('./scripts/transmitter_simulation.R') 
source('./scripts/mediation_analysis.R') 

ggsave(plot = panel_s1, filename = paste(figs_dir,sep = '/', "panel_s1.jpeg"), device = jpeg, width = 18, height = 12) #Within-host dynamics -functional
ggsave(plot = panel_s2, filename = paste(figs_dir,sep = '/', "panel_s2.jpeg"), device = jpeg, width = 18, height = 12) #Simulating transmitter population - functional
ggsave(plot = panel_s3, filename = paste(figs_dir,sep = '/', "panel_s3.jpeg"), device = jpeg, width = 14, height = 14) #Heritability model fitting 
ggsave(plot = plt_s4, filename = paste(figs_dir,sep = '/', "panel_s4.jpeg"), device = jpeg, width = 14, height = 14) #Sensitivity Analysis - functional

###################################################################################################

