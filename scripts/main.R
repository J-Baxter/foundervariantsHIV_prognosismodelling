# Main Analysis for "RECONCILING THE FOUNDER VARIANT MULTIPLICITY OF HIV-1 INFECTION WITH THE RATE
# OF CD4+ DECLINE"
# J Baxter

# This script combines the functions for the analysis of paired SpVL distributions, occuring in
# three main stages:

# ALL VIRAL LOADS ARE REAL NUMBERS UNLESS EXPLICITLY STATED


################################### Dependencies ###################################
source('./scripts/dependencies.R')
source('./scripts/functions/BonhoefferEqns.R')
source('./scripts/functions/AllocateTransmitter.R')
source('./scripts/functions/simulate_cohorts_funcs.R')
source('./scripts/functions/PostProcessing.R')
source('./scripts/interpret_outputs.R')
#source('./scripts/global_theme.R')


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
vc <- 0.47**2

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
                                ndraws = 30,
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
                              ndraws = 30,
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
                                               ndraws =30,
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
                                                  ndraws = 30,
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
    mutate(ID_pair = ID_pair + 196) %>% 
    mutate(transmitterallocation = 'ML'),
  
  # Empirical SpVL means, individual SpVL predicted using H2 model
  # Random allocation of transmitter
  shcs_predicted_transmitterrandom = shcs_h2preds_transmitterrandom %>%
    select(-c(contains('cohortmean'), 'SpVL', 'log10_SpVL'))%>%
    mutate(dataset = 'shcs_predicted') %>%
    mutate(ID_pair = ID_pair + 196*2) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted')) %>%
    mutate(transmitterallocation = 'random'),
  
  # Empirical SpVL means, individual SpVL predicted using H2 model
  # Transmitter = Maximum Likelihood(SpVL)
  shcs_predicted_transmitterML = shcs_h2preds_transmitterML %>%
    select(-c(contains('cohortmean'), 'SpVL', 'log10_SpVL'))%>%
    mutate(dataset = 'shcs_predicted') %>%
    mutate(ID_pair = ID_pair + 196*3) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted')) %>%
    mutate(transmitterallocation = 'ML'),
  
  # Simulated SpVL means for Heterosexuals, individual SpVL predicted using H2 model
  # Random allocation of transmitter
  stratified_het_transmitterrandom = stratified_pred_transmitterrandom  %>%
    filter(dataset == 'stratified_HET') %>%
    mutate(ID_pair = ID_pair + 196*4) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'random'),
  
  # Simulated SpVL means for MSM, individual SpVL predicted using H2 model
  # Random allocation of transmitter
  stratified_msm_transmitterrandom = stratified_pred_transmitterrandom  %>%
    filter(dataset == 'stratified_MSM') %>%
    mutate(ID_pair = ID_pair + 196*5) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'random'),
  
  # Simulated SpVL means for PWID, individual SpVL predicted using H2 model
  # Random allocation of transmitter
  stratified_pwid_transmitterrandom = stratified_pred_transmitterrandom  %>%
    filter(dataset == 'stratified_PWID') %>%
    mutate(ID_pair = ID_pair + 196*6) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'random'),
  
  # Simulated SpVL means for Heterosexuals, individual SpVL predicted using H2 model
  # Transmitter = Maximum Likelihood(SpVL)
  stratified_het_transmitterML = stratified_pred_transmitterML  %>%
    filter(dataset == 'stratified_HET') %>%
    mutate(ID_pair = ID_pair + 196*7) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'ML'),
  
  # Simulated SpVL means for MSM, individual SpVL predicted using H2 model
  # Transmitter = Maximum Likelihood(SpVL)
  stratified_msm_transmitterML = stratified_pred_transmitterML   %>%
    filter(dataset == 'stratified_MSM') %>%
    mutate(ID_pair = ID_pair + 196*8) %>% 
    rename_with(~ gsub("predicted_", "", .x), starts_with('predicted'))%>%
    mutate(transmitterallocation = 'ML'),
  
  # Simulated SpVL means for PWID, individual SpVL predicted using H2 model
  # Transmitter = Maximum Likelihood(SpVL)
  stratified_pwid_transmitterML = stratified_pred_transmitterML   %>%
    filter(dataset == 'stratified_PWID') %>%
    mutate(ID_pair = ID_pair + 196*9) %>% 
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
  setNames(c('FM', 'MF', 'MM', 'OTHER', 'PWID', 'UNKNOWN')) %>% as.list()



################################### Allocate f and p from joint posterior ###################################
load("./data/posteriors_by_exposure.RData")
fp <- list( mtf, ftm, msmi, msmr, pwid)

combined_data_CD4$FM <- combined_data_CD4$FM %>%
  cbind.data.frame(ftm[sample(1:nrow(ftm), nrow(combined_data_CD4$FM), replace = T),])

combined_data_CD4$MF <- combined_data_CD4$MF  %>%
  cbind.data.frame(mtf[sample(1:nrow(mtf), nrow(combined_data_CD4$MF), replace = T),])

combined_data_CD4$MMI<- combined_data_CD4$MM  %>%
  cbind.data.frame(msmi[sample(1:nrow(msmi), nrow(combined_data_CD4$MM), replace = T),])

combined_data_CD4$MMR<- combined_data_CD4$MM %>%
  cbind.data.frame(msmr[sample(1:nrow(msmr), nrow(combined_data_CD4$MM), replace = T),])

combined_data_CD4$PWID<- combined_data_CD4$PWID  %>%
  cbind.data.frame(pwid[sample(1:nrow(pwid), nrow(combined_data_CD4$PWID), replace = T),])

combined_data_CD4$OTHER <- combined_data_CD4$OTHER  %>%
  cbind.data.frame(mtf[sample(1:nrow(mtf), nrow(combined_data_CD4$OTHER), replace = T),])

combined_data_CD4$UNKNOWN <- combined_data_CD4$UNKNOWN  %>%
  cbind.data.frame(mtf[sample(1:nrow(mtf), nrow(combined_data_CD4$UNKNOWN), replace = T),])


################################### Calculate Joint Probability Dist MV/MP ###################################
# Implement transmission model (Thompson et al. and re-parameterised by Villabona-Arenas et al. (Unpublished))
# NB: THIS IS COMPUTATIONALLY EXPENSIVE 

# Set up cluster (fork)
cl <- detectCores() %>% `-` (2) 

# Female-to-Male
start <- Sys.time()
print(start)
FM_results <- mcmapply(TransmissionModel2,
                       sp_ViralLoad = combined_data_CD4$FM$SpVL_transmitter,
                       PerVirionProbability = combined_data_CD4$FM$ppp, 
                       PropExposuresInfective = combined_data_CD4$FM$fEnv,
                       SIMPLIFY = FALSE,
                       mc.cores = cl,
                       mc.set.seed = FALSE) %>%
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$FM  %>%
           group_split(index), 
         SIMPLIFY = F) %>%
  bind_rows()

end <- Sys.time()
elapsed <- end-start
print(end)
print(elapsed)

# Male-to-Female
start <- Sys.time()
print(start)

MF_results <- mcmapply(TransmissionModel2,
                       sp_ViralLoad = combined_data_CD4$MF$SpVL_transmitter,
                       PerVirionProbability = combined_data_CD4$MF$ppp, 
                       PropExposuresInfective = combined_data_CD4$MF$fEnv,
                       SIMPLIFY = FALSE,
                       mc.cores = cl,
                       mc.set.seed = FALSE) %>%
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$MF  %>%
           group_split(index), 
         SIMPLIFY = F)%>%
  bind_rows()

end <- Sys.time()
elapsed <- end-start
print(end)
print(elapsed)

# MSMI
start <- Sys.time()
print(start)

MMI_results <- mcmapply(TransmissionModel2, 
                        sp_ViralLoad = combined_data_CD4$MMI$SpVL_transmitter,
                        PerVirionProbability = combined_data_CD4$MMI$ppp, 
                        PropExposuresInfective = combined_data_CD4$MMI$fEnv,
                        SIMPLIFY = FALSE,
                        mc.cores = cl,
                        mc.set.seed = FALSE) %>%
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$MM  %>%
           group_split(index), 
         SIMPLIFY = F)%>%
  bind_rows()

end <- Sys.time()
elapsed <- end-start
print(end)
print(elapsed)

# MSMR
start <- Sys.time()
print(start)

MMR_results <- mcmapply(TransmissionModel2, 
                        sp_ViralLoad = combined_data_CD4$MMR$SpVL_transmitter,
                        PerVirionProbability = combined_data_CD4$MMR$ppp, 
                        PropExposuresInfective = combined_data_CD4$MMR$fEnv,
                        SIMPLIFY = FALSE,
                        mc.cores = cl,
                        mc.set.seed = FALSE) %>%
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$MM  %>%
           group_split(index), 
         SIMPLIFY = F)%>%
  bind_rows()

end <- Sys.time()
elapsed <- end-start
print(end)
print(elapsed)

# PWID
start <- Sys.time()
print(start)

PWID_results <- mcmapply(TransmissionModel2, 
                         sp_ViralLoad = combined_data_CD4$PWID$SpVL_transmitter,
                         PerVirionProbability = combined_data_CD4$PWID$ppp, 
                         PropExposuresInfective = combined_data_CD4$PWID$fEnv,
                         SIMPLIFY = FALSE,
                         mc.cores = cl,
                         mc.set.seed = FALSE) %>%
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$PWID  %>%
           group_split(index), 
         SIMPLIFY = F)%>%
  bind_rows()

end <- Sys.time()
elapsed <- end-start
print(end)
print(elapsed)

# OTHER: Currently using Male-to-Female
start <- Sys.time()
print(start)

OTHER_results <- mcmapply(TransmissionModel2, 
                          sp_ViralLoad = combined_data_CD4$OTHER$SpVL_transmitter,
                          PerVirionProbability = combined_data_CD4$MF$ppp, 
                          PropExposuresInfective = combined_data_CD4$MF$fEnv,
                          SIMPLIFY = FALSE,
                          mc.cores = cl,
                          mc.set.seed = FALSE) %>%
lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$OTHER  %>%
           group_split(index), 
         SIMPLIFY = F)%>%
  bind_rows()

end <- Sys.time()
elapsed <- end-start
print(end)
print(elapsed)


# UNKNOWN: Currently using Male-to-Female
start <- Sys.time()
print(start)

UNKNOWN_results <- mcmapply(TransmissionModel2, 
                            sp_ViralLoad = combined_data_CD4$UNKNOWN$SpVL_transmitter,
                            PerVirionProbability = combined_data_CD4$MF$ppp, 
                            PropExposuresInfective = combined_data_CD4$MF$fEnv,
                            SIMPLIFY = FALSE,
                            mc.cores = cl,
                            mc.set.seed = FALSE) %>%
lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','SpVL')) %>%
  
  mapply(PostProcess,
         tmresult = ., 
         metadata = combined_data_CD4$UNKNOWN  %>%
           group_split(index), 
         SIMPLIFY = F)%>%
  bind_rows()

end <- Sys.time()
elapsed <- end-start
print(end)
print(elapsed)


################################### Write to file ################################### 
# Sort into 'dataframes' (dataset:transmitterselection)
datanames <- c('FM', 'MF', 'MMI', 'MMR',  'PWID')

combined_results <- list(FM_results, 
                          MF_results,
                          MMI_results,
                          MMR_results,
                          PWID_results) %>%
  setNames(., datanames)%>%
  bind_rows(., .id = 'dataset_id') %>%
  dplyr::select(-ends_with(as.character(35:200)))


combined_results_list <- combined_results %>%
  group_split(dataset_id) %>% 
  setNames(., datanames)

# write csv to file
filenames <- paste0(results_dir, '/', datanames, '_rawresults.csv')
mapply(write_csv, combined_results_list, file = filenames)


################################### Resample: Compare CD4+ ################################### 
combined_results_simonly <- combined_results  %>%
  filter(grepl('stratified', dataset)) %>%
  group_split(dataset_id, transmitterallocation)

CD4_resample <- lapply(combined_results_simonly, GetCD4Survival, replicates = 1000)


# write csv to file
datanames <- lapply(combined_results_simonly, function(x) paste(unique(x$dataset_id),
                                                                 unique(x$transmitterallocation), sep ='_')) %>% 
  unlist()
filenames <- paste0(results_dir, '/', datanames, '_CD4resample.csv')
mapply(write_csv, CD4_resample, file = filenames)


################################### Resample: Compare SpVL ################################### 
SpVL_resample <- lapply(combined_results_simonly, GetSpVLDiff)

# write csv to file
filenames <- paste0(results_dir, '/', datanames, '_SpVLresample.csv')
mapply(write_csv, SpVL_resample, file = filenames)


source('./scripts/figures.R')
#source ('./scripts/supplementary_figures.R)


###################################################################################################
############################################## END ################################################ 
###################################################################################################