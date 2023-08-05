# Main Analysis for "RECONCILING THE FOUNDER VARIANT MULTIPLICITY OF HIV-1 INFECTION WITH THE RATE
# OF CD4+ DECLINE"
# J Baxter

# This script combines the functions for the analysis of paired SpVL distributions, occuring in
# three main stages:

# ALL VIRAL LOADS ARE REAL NUMBERS UNLESS EXPLICITLY STATED


################################### Dependencies ###################################
source('./scripts/dependencies.R')
source('./scripts/BonhoefferEqns.R')
source('./scripts/AllocateTransmitter.R')
source('./scripts/simulate_cohorts_funcs.R')
source('./scripts/PostProcessing.R')
source('./scripts/global_theme.R')


# Set Seed
set.seed(4472)
options(scipen = 100) #options(scipen = 100, digits = 4)


################################### Import Data ###################################
# Paired set point viral load (SpVL), age, riskgroup and sex data provided on request
# by the Swiss HIV Cohort Study. These data were used in the analysis of Bertels et. al
# 2018 https://academic.oup.com/mbe/article/35/1/27/4210012

# init distribtions for transmitter allocation
# mc and vc are the mean and variance of the Amsterdam seroconverter study, respectively.
mc <- 4.35
vc <- 0.478**2

t <- CalcTransmitter(mc,vc)
r <- CalcRecipient(mc,vc)

# Import data and pre-process
source('./scripts/import_data.R')


################################### Fit/Import Base Models ###################################
# Import fitted model
source('./scripts/TransmissionModel.R')
source('./scripts/ToleranceModel.R')


# Fit Bayesian linear mixed model to SHCS data (<2 mins to fit)
source('./scripts/HeritabilityModel.R')


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


################################### Validate composition of simulated cohorts ###################################
source('./scripts/validate_simcohorts.R')


################################### Predict SpVLs ###################################

# Random allocation of transmitter/recipient
shcs_h2preds_transmitterrandom <- predicted_draws(heritability_model_transmitterrandom,
                                newdata = shcs_data_long_transmitterrandom,
                                ndraws = 50,
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


stratified_pred_transmitterrandom  <- predicted_draws(heritability_model_transmitterrandom ,
                              newdata = stratified_data,
                              allow_new_levels = TRUE, # allow new random effects levels
                              sample_new_levels = "old_levels", #
                              re_formula = ~ (1|log10_SpVL_couplemean),
                              ndraws = 50,
                              value = 'predicted_log10_SpVL') %>% 
  #select(-c(.chain, .iteration, .draw)) %>%
  group_by(ID_pair, partner) %>%
  mutate(predicted_SpVL = raise_to_power(10, predicted_log10_SpVL)) %>%
  select(-.row) %>%
  distinct() %>%
  ungroup()


# Transmitter = Max(SpVLij) of couple j
shcs_h2preds_transmitterML <- predicted_draws(heritability_model_transmitterML,
                                               newdata = shcs_data_long_transmitterML,
                                               ndraws = 50,
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


stratified_pred_transmitterML <- predicted_draws(heritability_model_transmitterML,
                                                  newdata = stratified_data,
                                                  allow_new_levels = TRUE, # allow new random effects levels
                                                  sample_new_levels = "old_levels", #
                                                  re_formula = ~ (1|log10_SpVL_couplemean),
                                                  ndraws = 50,
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
  shcs_empirical_transmitterrandom = shcs_data_long_transmitterrandom %>%
    select(-contains('cohortmean')) %>%
    mutate(dataset = 'shcs_empirical') %>%
    mutate(transmitterallocation = 'random'),
  
  shcs_empirical_transmitterML = shcs_data_long_transmitterML %>%
    select(-contains('cohortmean')) %>%
    mutate(dataset = 'shcs_empirical') %>%
    mutate(transmitterallocation = 'ML'),
  
  # Empirical SpVL means, individual SpVL predicted using H2 model
  # Random allocation of transmitter
  shcs_predicted_transmitterrandom = shcs_h2preds_transmitterrandom %>%
    select(-c(contains('cohortmean'), 'SpVL', 'log10_SpVL'))%>%
    mutate(dataset = 'shcs_predicted') %>%
    mutate(ID_pair = ID_pair + 196) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted')) %>%
    mutate(transmitterallocation = 'random'),
  
  # Empirical SpVL means, individual SpVL predicted using H2 model
  # Transmitter = Maximum Likelihood(SpVL)
  shcs_predicted_transmitterML = shcs_h2preds_transmitterML %>%
    select(-c(contains('cohortmean'), 'SpVL', 'log10_SpVL'))%>%
    mutate(dataset = 'shcs_predicted') %>%
    mutate(ID_pair = ID_pair + 196*5) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted')) %>%
    mutate(transmitterallocation = 'ML'),
  
  # Simulated SpVL means for Heterosexuals, individual SpVL predicted using H2 model
  # Random allocation of transmitter
  stratified_het_transmitterrandom = stratified_pred_transmitterrandom  %>%
    filter(dataset == 'stratified_HET') %>%
    mutate(ID_pair = ID_pair + 196*2) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'random'),
  
  # Simulated SpVL means for MSM, individual SpVL predicted using H2 model
  # Random allocation of transmitter
  stratified_msm_transmitterrandom = stratified_pred_transmitterrandom  %>%
    filter(dataset == 'stratified_MSM') %>%
    mutate(ID_pair = ID_pair + 196*3) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'random'),
  
  # Simulated SpVL means for PWID, individual SpVL predicted using H2 model
  # Random allocation of transmitter
  stratified_pwid_transmitterrandom = stratified_pred_transmitterrandom  %>%
    filter(dataset == 'stratified_PWID') %>%
    mutate(ID_pair = ID_pair + 196*4) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'random'),
  
  # Simulated SpVL means for Heterosexuals, individual SpVL predicted using H2 model
  # Transmitter = Maximum Likelihood(SpVL)
  stratified_het_transmitterML = stratified_pred_transmitterML  %>%
    filter(dataset == 'stratified_HET') %>%
    mutate(ID_pair = ID_pair + 196*6) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'ML'),
  
  # Simulated SpVL means for MSM, individual SpVL predicted using H2 model
  # Transmitter = Maximum Likelihood(SpVL)
  stratified_msm_transmitterML = stratified_pred_transmitterML   %>%
    filter(dataset == 'stratified_MSM') %>%
    mutate(ID_pair = ID_pair + 196*7) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'ML'),
  
  # Simulated SpVL means for PWID, individual SpVL predicted using H2 model
  # Transmitter = Maximum Likelihood(SpVL)
  stratified_pwid_transmitterML = stratified_pred_transmitterML   %>%
    filter(dataset == 'stratified_PWID') %>%
    mutate(ID_pair = ID_pair + 196*8) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'ML')
) 


################################### Predict CD4 decline and format ###################################

combined_data_CD4 <- combined_data  %>%
  lapply(., function(x) x %>% mutate(delta_CD4 = ToleranceModel(log10_SpVL,
                                                                age.inf,
                                                                sex))) %>%
  lapply(., function(x) x %>%   mutate(riskgroup = case_when(riskgroup == 'HET' & sex == 'M' & partner == 'transmitter' ~ 'MF',
                                                             riskgroup == 'HET' & sex == 'F' & partner == 'recipient' ~ 'MF',
                                                             riskgroup == 'MSM' ~ 'MM',
                                                             riskgroup == 'PWID' ~ 'PWID',
                                                             riskgroup == 'UNKNOWN' ~ 'UNKNOWN',
                                                             riskgroup == 'OTHER' ~ 'OTHER',
                                                             .default = 'FM'))) %>%
  lapply(., function(x) x %>%  dplyr::select(-c(sex, contains('age')))) %>%
  lapply(., function(x) x %>%  pivot_wider(names_from = partner, values_from = c(log10_SpVL, SpVL, delta_CD4, riskgroup))) %>%
  bind_rows() %>%
  rowid_to_column( "index") %>% 
  group_split(riskgroup_recipient) %>%
  setNames(c('FM', 'MF', 'MM', 'OTHER', 'PWID', 'UNKNOWN')) 


################################### Calculate Joint Probability Dist MV/MP ###################################
# Implement transmission model (Thompson et al. and re-parameterised by Villabona-Arenas et al. (Unpublished))
# NB: THIS IS COMPUTATIONALLY EXPENSIVE AND WILL RUN FOR ~ 25 HOURS USING 4 CORES

# Female-to-Male
FM_results <- RunParallel(TransmissionModel2,
                          combined_data_CD4$FM$SpVL_transmitter,
                          PerVirionProbability = 8.779E-07, 
                          PropExposuresInfective = 0.14337) %>%
  
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$FM  %>%
           group_split(index), 
         SIMPLIFY = F) %>%
  bind_rows()


# Male-to-Female
MF_results <- RunParallel(TransmissionModel2,
                          combined_data_CD4$MF$SpVL_transmitter,
                          PerVirionProbability = 1.765E-06, 
                          PropExposuresInfective = 0.13762) %>%
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$FM  %>%
           group_split(index), 
         SIMPLIFY = F)%>%
  bind_rows()


# MSM: currently using MSM:insertive params
MM_results <- RunParallel(TransmissionModel2, 
                          combined_data_CD4$MM$SpVL_transmitter,
                          PerVirionProbability = 8.779E-07, 
                          PropExposuresInfective = 0.14337) %>%
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$FM  %>%
           group_split(index), 
         SIMPLIFY = F)%>%
  bind_rows()



# PWID: Currently using MSM:receptive params
PWID_results <- RunParallel(TransmissionModel2, 
                            combined_data_CD4$PWID$SpVL_transmitter,
                            PerVirionProbability = 3.19E-06, 
                            PropExposuresInfective = 0.08923) %>%
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$FM  %>%
           group_split(index), 
         SIMPLIFY = F)%>%
  bind_rows()


# OTHER: Currently using Male-to-Female
OTHER_results <- RunParallel(TransmissionModel2, 
                               combined_data_CD4$OTHER$SpVL_transmitter,
                               PerVirionProbability = 8.779E-07, 
                               PropExposuresInfective = 0.14337) %>%
lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$FM  %>%
           group_split(index), 
         SIMPLIFY = F)%>%
  bind_rows()


# UNKNOWN: Currently using Male-to-Female
UNKNOWN_results <- RunParallel(TransmissionModel2, 
                               combined_data_CD4$UNKNOWN$SpVL_transmitter,
                               PerVirionProbability = 8.779E-07, 
                               PropExposuresInfective = 0.14337) %>%
lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$FM  %>%
           group_split(index), 
         SIMPLIFY = F)%>%
  bind_rows()


################################### Write to file ################################### 

# Sort into 'dataframes' (dataset:transmitterselection)
datanames <- unique(FM_results$dataset)

combinded_results <- list(FM_results, 
                          MF_results,
                          MM_results,
                          PWID_results,
                          OTHER_results,
                          UNKNOWN_results) %>%
  do.call(bind.data.frame, .) %>%
  group_split(dataset) 

# write csv to file
filenames <- paste0(results_dir, '/', datanames, '_modelresults.csv')

mapply(write_csv, combinded_results, file = filenames)


#source('./scripts/figures.R')

#ggsave(plot = panel_1, filename = paste(figs_dir,sep = '/', "panel_1.jpeg"), device = jpeg, width = 14, height = 14) # Model Components - functional
#ggsave(plot = panel_2, filename = paste(figs_dir,sep = '/', "panel_2.jpeg"), device = jpeg, width = 14, height = 18) # Confounder - functional 
#ggsave(plot = panel_3, filename = paste(figs_dir,sep = '/', "panel_3.jpeg"), device = jpeg, width = 14, height = 18) # Non - Linear
#ggsave(plot = panel_4, filename = paste(figs_dir,sep = '/', "panel_4.jpeg"), device = jpeg, width = 14, height = 18) # Timing of transmission

# Supplementary plots
# SHCS summary plots
# Simulation Bias Checks
# H2 model plots
# SA ?
#source('./scripts/tm_withinhostprocesses.R') 


#ggsave(plot = panel_s1, filename = paste(figs_dir,sep = '/', "panel_s1.jpeg"), device = jpeg, width = 18, height = 12) #Within-host dynamics -functional
#ggsave(plot = panel_s2, filename = paste(figs_dir,sep = '/', "panel_s2.jpeg"), device = jpeg, width = 18, height = 12) #Simulating transmitter population - functional
#ggsave(plot = panel_s3, filename = paste(figs_dir,sep = '/', "panel_s3.jpeg"), device = jpeg, width = 14, height = 14) #Heritability model fitting 
#ggsave(plot = plt_s4, filename = paste(figs_dir,sep = '/', "panel_s4.jpeg"), device = jpeg, width = 14, height = 14) #Sensitivity Analysis - functional

###################################################################################################

