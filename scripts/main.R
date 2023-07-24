# Main Analysis for "RECONCILING THE FOUNDER VARIANT MULTIPLICITY OF HIV-1 INFECTION WITH THE RATE
# OF CD4+ DECLINE"
# J Baxter

# This script combines the functions for the analysis of paired SpVL distributions, occuring in
# three main stages:

# ALL VIRAL LOADS ARE REAL NUMBERS UNLESS EXPLICITLY STATED


################################### Dependencies ###################################
source('./scripts/dependencies.R')
source('./scripts/simulate_cohorts_funcs.R')
source('./scripts/global_theme.R')


# Set Seed
set.seed(4472)
options(scipen = 100) #options(scipen = 100, digits = 4)


################################### Import Data ###################################
# Paired set point viral load (SpVL), age, riskgroup and sex data provided on request
# by the Swiss HIV Cohort Study. These data were used in the analysis of Bertels et. al
# 2018 https://academic.oup.com/mbe/article/35/1/27/4210012

source('./scripts/import_data.R')


################################### Fit/Import Base Models ###################################
# Import fitted model
source('./scripts/transmission_model.R')
source('./scripts/tolerance_model.R')

# Fit Bayesian linear mixed model to SHCS data (<2 mins to fit)
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
  dplyr::filter(!(riskgroup_1 == 'HET' & sex_1 == sex_2)) %>% # currently makes het cov matrix not positive definitive (ie singular)
  EnumerateLevels(.) #NB will drop any column with only 1 level

# Sanity Check
stopifnot(all(lapply(shcs_data_int_list, function(x) cov(x) %>% is.positive.definite())))


# Infer cumulative probabilities for each categorical level
cumulativeprobs_list <- lapply(shcs_data_int_list, function(x) x %>%
                           select(where(is.integer)) %>% 
                           apply(2, function(i) cumsum(table(i)) / length(i), simplify = FALSE))


# Run SimCohorts - see ./scripts/SHCS_simulation.R for details
stratified_data <- mapply(SimCohorts,
                          data = shcs_data_int_list,
                          probs = cumulativeprobs_list,
                          SIMPLIFY = F) %>%
  setNames(., c('HET', 'MSM', 'PWID')) %>%
  
  # Re-organising simulated data
  bind_rows(., .id = "riskgroup") %>% 
  rowid_to_column( "ID_pair") %>%
  separate_wider_position(cols = sex, widths = c(sex_1 = 1, sex_2 = 1)) %>%
  pivot_longer(cols = - c(contains('couplemean'), riskgroup, ID_pair), 
               names_to = c(".value", "partner"), 
               names_pattern  = "^(.*)_([0-9])$") %>% 
  mutate(partner = case_when(partner == 1 ~ 'transmitter',
                             partner == 2 ~ 'recipient')) %>%
  mutate(partner = factor(partner, levels = c("transmitter", "recipient"))) %>%
  mutate(age.inf_category = cut(age.inf, 
                                breaks = c(15,24,29,39,80), 
                                labels = c('15-24', '25-29', '30-39','40-80'))) %>%
  relocate(c(partner, sex), .before = riskgroup) %>%
  relocate(ends_with('couplemean'), .after = last_col()) %>%
  mutate(SpVL_couplemean = raise_to_power(10, log10_SpVL_couplemean)) %>%
  mutate(dataset = paste('stratified', riskgroup, sep = '_'))


# Validate composition of simulated cohorts
source('./scripts/validate_simcohorts.R')


################################### Predict SpVLs ###################################

# Random allocation of transmitter/recipient
shcs_h2preds_randomallocation <- predicted_draws(heritability_model_randomallocation,
                                newdata = shcs_data_long,
                                ndraws = 100,
                                allow_new_levels = FALSE, # allow new random effects levels
                                #sample_new_levels = "old_levels", #
                         ##re_formula = ~ (1|log10_SpVL_couplemean),
                         value = 'predicted_log10_SpVL') %>% 
  #select(-c(.chain, .iteration, .draw)) %>%
  group_by(ID_pair, partner) %>%
  mutate(predicted_SpVL = raise_to_power(10, predicted_log10_SpVL)) %>%
  select(-.row) %>%
  distinct() %>%
  ungroup()


stratified_pred_randomallocation <- predicted_draws(heritability_model_randomallocation,
                              newdata = stratified_data,
                              allow_new_levels = TRUE, # allow new random effects levels
                              sample_new_levels = "old_levels", #
                              re_formula = ~ (1|log10_SpVL_couplemean),
                              ndraws = 100,
                              value = 'predicted_log10_SpVL') %>% 
  #select(-c(.chain, .iteration, .draw)) %>%
  group_by(ID_pair, partner) %>%
  mutate(predicted_SpVL = raise_to_power(10, predicted_log10_SpVL)) %>%
  select(-.row) %>%
  distinct() %>%
  ungroup()


# Transmitter = Max(SpVLij) of couple j
shcs_h2preds_transmittermax <- predicted_draws(heritability_model_transmittermax,
                                               newdata = shcs_data_long,
                                               ndraws = 100,
                                               allow_new_levels = FALSE, # allow new random effects levels
                                               #sample_new_levels = "old_levels", #
                                               ##re_formula = ~ (1|log10_SpVL_couplemean),
                                               value = 'predicted_log10_SpVL') %>% 
  #select(-c(.chain, .iteration, .draw)) %>%
  group_by(ID_pair, partner) %>%
  mutate(predicted_SpVL = raise_to_power(10, predicted_log10_SpVL)) %>%
  select(-.row) %>%
  distinct() %>%
  ungroup()


stratified_pred_transmittermax <- predicted_draws(heritability_model_transmittermax,
                                                  newdata = stratified_data,
                                                  allow_new_levels = TRUE, # allow new random effects levels
                                                  sample_new_levels = "old_levels", #
                                                  re_formula = ~ (1|log10_SpVL_couplemean),
                                                  ndraws = 100,
                                                  value = 'predicted_log10_SpVL') %>% 
  #select(-c(.chain, .iteration, .draw)) %>%
  group_by(ID_pair, partner) %>%
  mutate(predicted_SpVL = raise_to_power(10, predicted_log10_SpVL)) %>%
  select(-.row) %>%
  distinct() %>%
  ungroup()


################################### Separate DFs ###################################
combined_data <- list(
  
  # Empirical SpVL, not adjustment from H2 model
  shcs_empirical = shcs_data_long %>%
    select(-contains('cohortmean')) %>%
    mutate(dataset = 'shcs_empirical') %>%
    mutate(transmitterallocation = 'random'),
  
  # Empirical SpVL means, individual SpVL predicted using H2 model
  # Random allocation of transmitter
  shcs_predicted = shcs_h2preds_randomallocation %>%
    select(-c(contains('cohortmean'), 'SpVL', 'log10_SpVL'))%>%
    mutate(dataset = 'shcs_predicted') %>%
    mutate(ID_pair = ID_pair + 196) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted')) %>%
    mutate(transmitterallocation = 'random'),
  
  # Simulated SpVL means for Heterosexuals, individual SpVL predicted using H2 model
  # Random allocation of transmitter
  stratified_het = stratified_pred_randomallocation  %>%
    filter(dataset == 'stratified_HET') %>%
    mutate(ID_pair = ID_pair + 196*2) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'random'),
  
  # Simulated SpVL means for MSM, individual SpVL predicted using H2 model
  # Random allocation of transmitter
  stratified_msm = stratified_pred_randomallocation  %>%
    filter(dataset == 'stratified_MSM') %>%
    mutate(ID_pair = ID_pair + 196*3) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'random'),
  
  # Simulated SpVL means for PWID, individual SpVL predicted using H2 model
  # Random allocation of transmitter
  stratified_pwid = stratified_pred_randomallocation  %>%
    filter(dataset == 'stratified_PWID') %>%
    mutate(ID_pair = ID_pair + 196*4) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'random'),
  
  
  # Empirical SpVL means, individual SpVL predicted using H2 model
  # Transmitter = max(SpVL)
  shcs_predicted = shcs_h2preds_transmittermax %>%
    select(-c(contains('cohortmean'), 'SpVL', 'log10_SpVL'))%>%
    mutate(dataset = 'shcs_predicted') %>%
    mutate(ID_pair = ID_pair + 196*5) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted')) %>%
    mutate(transmitterallocation = 'max'),
  
  # Simulated SpVL means for Heterosexuals, individual SpVL predicted using H2 model
  # Transmitter = max(SpVL)
  stratified_het = stratified_pred_transmittermax   %>%
    filter(dataset == 'stratified_HET') %>%
    mutate(ID_pair = ID_pair + 196*6) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'max'),
  
  # Simulated SpVL means for MSM, individual SpVL predicted using H2 model
  # Transmitter = max(SpVL)
  stratified_msm = stratified_pred_transmittermax   %>%
    filter(dataset == 'stratified_MSM') %>%
    mutate(ID_pair = ID_pair + 196*7) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'max'),
  
  # Simulated SpVL means for PWID, individual SpVL predicted using H2 model
  # Transmitter = max(SpVL)
  stratified_pwid = stratified_pred_transmittermax   %>%
    filter(dataset == 'stratified_PWID') %>%
    mutate(ID_pair = ID_pair + 196*8) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'max')
) %>%
  bind_rows() 


################################### Predict CD4 decline and format ###################################

combined_data_CD4 <- combined_data  %>%
  mutate(delta_CD4 = ToleranceModel(log10_SpVL,
                                    age.inf,
                                    sex)) %>%
  mutate(riskgroup = case_when(riskgroup == 'HET' & sex == 'M' & partner == 'transmitter' ~ 'MF',
                               riskgroup == 'HET' & sex == 'F' & partner == 'recipient' ~ 'MF',
                               riskgroup == 'MSM' ~ 'MM',
                               riskgroup == 'PWID' ~ 'PWID',
                               riskgroup == 'UNKNOWN' ~ 'UNKNOWN',
                               riskgroup == 'OTHER' ~ 'OTHER',
                               .default = 'FM'
                               )) %>%
  dplyr::select(-c(sex, contains('age'))) %>%
  pivot_wider(names_from = partner, values_from = c(contains('SpVL'), delta_CD4, riskgroup)) %>% 
  group_split(riskgroup_recipient) %>%
  setNames(c('FM', 'MF', 'MM', 'PWID', 'UNKNOWN'))


################################### Calculate Joint Probability Dist MV/MP ###################################
# Implement transmission model (Thompson et al. and re-parameterised by Villabona-Arenas et al. (Unpublished))
# NB: THIS IS COMPUTATIONALLY EXPENSIVE AND WILL RUN FOR MANY HOURS (159396 iterations)

combined_data_PMV <- c(
  
  # Female-to-Male
  RunParallel(TransmissionModel2,
              combined_data_CD4$FM$SpVL_transmitter,
              PerVirionProbability = 8.779E-07, 
              PropExposuresInfective = 0.14337),
                       
  # Male-to-Female
  RunParallel(TransmissionModel2,
              combined_data_CD4$MF$SpVL_transmitter,
              PerVirionProbability = 1.765E-06, 
              PropExposuresInfective = 0.13762),
  
  # MSM: currently using MSM:insertive params
  RunParallel(TransmissionModel2, 
              combined_data_CD4$MM$SpVL_transmitter,
              PerVirionProbability = 8.779E-07, 
              PropExposuresInfective = 0.14337),
  
  #PWID: Currently using MSM:receptive params
  RunParallel(TransmissionModel2, 
              combined_data_CD4$PWID$SpVL_transmitter,
              PerVirionProbability = 3.19E-06, 
              PropExposuresInfective = 0.08923),
  
  # Unknown: Currently using Male-to-Female
  RunParallel(TransmissionModel2, 
              combined_data_CD4$UNKNOWN$SpVL_transmitter,
              PerVirionProbability = 8.779E-07, 
              PropExposuresInfective = 0.14337)) %>%
                                   
  
  # Label
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) #%>%
  #mapply(function(x,y) c(x, y ),
         #x= ., 
         #y = combined_data_CD4 %>% group_split(riskgroup_recipient, SpVL), 
        # SIMPLIFY = F)


################################### Format Model Outputs ###################################
# Warning: this will take some time to run

combined_data_CD4_PMV <- combined_data_PMV %>%
  # Sum joint probabilities to estimate marginal probabilities of the number of variants and 
  # virus particles initiating infection
  
  lapply(., function(x) c(x, list('p_particles' = rowSums(x[['variant_distribution']] %>%
                                                            select(-nparticles) %>% 
                                                            as.vector())))) %>%
  lapply(., function(x) c(x, list('p_variants' = colSums(x[['variant_distribution']] %>% 
                                                           select(-c(nparticles, prob_nparticles)))))) %>%
  
  # Process model outputs to form dataframe where each row contains all results for one participant
  lapply(., function(x) bind_cols(bind_cols(x[-c(1,2,3, 17,18)]), 
                                  bind_cols(particles = 1:33, x[17]) %>% 
                                    pivot_wider(names_from = particles, values_from = p_particles, names_prefix = 'p_particles_'), 
                                  bind_cols(variants = 1:33, x[18]) %>%
                                    pivot_wider(names_from = variants, values_from = p_variants, names_prefix = 'p_variants_'))) %>%
  bind_rows() %>%
  
  # Remove high number variants/virions (probability of these events is approximately equivalent to 0)
  select(-ends_with(as.character(13:33)))


# Confirm output is a dataframe of dimensions 
stopifnot(dim(combined_data_CD4_PMV) == c(7936, 37))
  
  
# Save to file (so you don't have to run again!)
write.csv(combined_data_CD4_PMV, paste0(results_dir, '/combined_results_CD4_PM.csv'))

combined_data_CD4_PMV <- read_csv('./results/24Jun23/combined_results_CD4_PM.csv',col_types = cols(...1 = col_skip()))

combined_data_CD4_PMV_wide_long <- combined_data_CD4_PMV %>%
  filter(w == 1) %>% # No weighting for early transmission
  pivot_wider(names_from = partner, values_from = c(sex, age.inf, riskgroup, SpVL, log10_SpVL, age.inf_category, delta_CD4, starts_with('p_')), names_sep = '_') %>%
  pivot_longer(cols = starts_with('p_'), names_to = 'x', values_to = 'probability') %>%
  separate(x, sep = '_', c(NA, 'type', 'x', 'partner'))


combined_data_CD4_PMV_wide <- combined_data_CD4_PMV %>%
  filter(w == 1) %>%
  pivot_wider(names_from = partner, values_from = c(sex, age.inf, riskgroup, SpVL, log10_SpVL, age.inf_category, delta_CD4, starts_with('p_')), names_sep = '_') 

################################### Write to file ################################### 
source('./scripts/figures.R')

ggsave(plot = panel_1, filename = paste(figs_dir,sep = '/', "panel_1.jpeg"), device = jpeg, width = 14, height = 14) # Model Components - functional
ggsave(plot = panel_2, filename = paste(figs_dir,sep = '/', "panel_2.jpeg"), device = jpeg, width = 14, height = 18) # Confounder - functional 
ggsave(plot = panel_3, filename = paste(figs_dir,sep = '/', "panel_3.jpeg"), device = jpeg, width = 14, height = 18) # Non - Linear
ggsave(plot = panel_4, filename = paste(figs_dir,sep = '/', "panel_4.jpeg"), device = jpeg, width = 14, height = 18) # Timing of transmission

# Supplementary plots
# SHCS summary plots
# Simulation Bias Checks
# H2 model plots
# SA ?
source('./scripts/tm_withinhostprocesses.R') 


ggsave(plot = panel_s1, filename = paste(figs_dir,sep = '/', "panel_s1.jpeg"), device = jpeg, width = 18, height = 12) #Within-host dynamics -functional
ggsave(plot = panel_s2, filename = paste(figs_dir,sep = '/', "panel_s2.jpeg"), device = jpeg, width = 18, height = 12) #Simulating transmitter population - functional
ggsave(plot = panel_s3, filename = paste(figs_dir,sep = '/', "panel_s3.jpeg"), device = jpeg, width = 14, height = 14) #Heritability model fitting 
ggsave(plot = plt_s4, filename = paste(figs_dir,sep = '/', "panel_s4.jpeg"), device = jpeg, width = 14, height = 14) #Sensitivity Analysis - functional

###################################################################################################

