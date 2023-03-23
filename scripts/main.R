# Main Analysis for "RECONCILING THE FOUNDER VARIANT MULTIPLICITY OF HIV-1 INFECTION WITH THE RATE
# OF CD4+ DECLINE"
# J Baxter

# This script combines the functions for the analysis of paired SPVL distributions, occuring in
# three main stages:

# 1. Initialisation of study population
# 2. Null hypothesis
# 3. Ecological hypotheses


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
data <- read_delim('./data/pbio.1001951.s006.tsv', delim = '\t') #Post processing

NSIM <- 100
################################### Run Base Models ###################################
source('./scripts/base_models.R')

################################### Q1 ###################################
# Does SPVL in the transmitting partner confound the relationship between the number of founder 
# variants and viral load in the recipient?

transmitter_tm <- RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 1)%>%
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w')) 

pop_tm <- transmitter_tm %>%
  lapply(., cbind.data.frame) %>%
  do.call(rbind.data.frame,.) %>%
  rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
  dplyr::select(-contains('nparticles')) %>%
  dplyr::summarise(across(starts_with('V'), .fns =sum),.by = transmitter) %>%
  pivot_longer(cols = starts_with('V'),
               names_to = 'variants',
               values_to = 'p') %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% 
           as.numeric())  %>%
  filter(variants <= 10) %>%
  mutate(w = 1)




linear_sim <- SimDonor(100, 1,7) %>%
  posterior_predict(heritability_model, .) %>%
  apply(., 2, sample, 1) %>% #Sample one value from posterior predictions per transmitter
  cbind.data.frame(recipient_log10SPVL = ., transmitter_log10SPVL = pop[['transmitter_log10SPVL']]) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SPVL')}")) %>%
  `colnames<-` (str_remove(colnames(.), 'sim_')) %>% 
  mutate(model = 'linear') 

pop_pred_variants <- transmitter_tm  %>%
  VariantPP(pop = pop) 

pop_pred_virions <- transmitter_tm %>%
  VirionPP(pop = pop) 


linear_pred <- RunParallel(populationmodel_acrossVL_Environment, linear_sim$transmitter, w= 1) %>%
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w'))

linear_pred_variants <- linear_pred %>%
  VariantPP(pop = linear_sim) %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(spVL = log10(recipient)))) 
  
linear_pred_virions <- linear_pred %>% 
  VirionPP(pop = linear_sim) %>% 
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(spVL = log10(recipient))))


  
################################### Q2 ###################################
# Is the SPVL in recipient partner determined by a non-linear relationship with the SPVL in the
# transmitting partner, and how is this affect by the number of variants initiating infection in 
# the recipient partner?

quad_model_uw <- lm(recipient_log10SPVL ~ transmitter_log10SPVL+ I(transmitter_log10SPVL**2), data = pop) 

exp_model_uw <- lm(recipient_log10SPVL ~ exp(transmitter_log10SPVL), data = pop) 

sigmoid_model_uw <-  

# Check Fit

# Simulate Populations
quad_sim_uw  <- SimPop(data$transmitter, quad_model_uw , NSIM)

exp_sim_uw  <- SimPop(data$transmitter, exp_model_uw , NSIM)

sigmoid_sim_uw  <- SimPop(data$transmitter, sigmoid_model_uw , NSIM)

# Run Transmission Model on Simulations
quad_pred_uw  <- RunTM4Sim(quad_sim_uw , modeltype = 'quad')

exp_pred_uw  <- RunTM4Sim(exp_sim_uw , modeltype = 'exp')

sigmoid_pred_uw  <- RunTM4Sim(sigmoid_sim_uw , modeltype = 'sigmoid')


# Models with weights



################################### Q3 ###################################
# How does the timing of transmission impact observations of the association between the number
# of founder variants and CD4+ T cell decline?

# weighting either random or fitted from network model
transmitter_timing <- c(RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 1),
                  RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 5),
                  RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 10),
                  #RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 15),
                  RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 20)) %>%
  
  # Label
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w'))


timing_tm <- transmitter_timing  %>%
  lapply(., cbind.data.frame) %>%
  do.call(rbind.data.frame,.) %>%
  rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
  dplyr::select(-contains('nparticles')) %>%
  dplyr::summarise(across(starts_with('V'), .fns =sum),.by = c(transmitter, w)) %>%
  pivot_longer(cols = starts_with('V'),
               names_to = 'variants',
               values_to = 'p') %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% 
           as.numeric()) 

transmitter_timing_variants <- transmitter_timing  %>%
  VariantPP(pop = linear_sim) 

transmitter_timing_virions <- transmitter_timing %>%
  VirionPP(pop = linear_sim) 


timing_pred <- c(RunParallel(populationmodel_acrossVL_Environment, linear_sim$transmitter, w= 1),
                 RunParallel(populationmodel_acrossVL_Environment, linear_sim$transmitter, w= 5),
                 RunParallel(populationmodel_acrossVL_Environment, linear_sim$transmitter, w= 10),
                 #RunParallel(populationmodel_acrossVL_Environment, linear_sim$transmitter, w= 15),
                 RunParallel(populationmodel_acrossVL_Environment, linear_sim$transmitter, w= 20)) %>%
  
  # Label
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w'))

timing_pred_variants <- timing_pred %>%
  VariantPP(pop = linear_sim) %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(spVL = log10(recipient)))) 

timing_pred_virions <- timing_pred %>% 
  VirionPP(pop = linear_sim) %>% 
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(spVL = log10(recipient))))
################################### Write to file ################################### To be changed

# transmitter SPVL & variant distribution
write_csv(transmitterspvl_variantdist, file = paste(results_dir, 'transmitterspvl_variantdist.csv', sep = '/'))

# recipient SPVL & variant distribution
write_csv(sim_recipientspvl_variantdist, file = paste(results_dir, 'recipientspvl_variantdist.csv', sep = '/'))


###################################################################################################
################################### Non-Linear Models (Probability MV weights) ###################################
# Four groups of nonlinear relationships trialled: Polynomial, Concave/Convex Curves, Sigmoidal, Curves with Max/Min

# Calculate weights based on transmitter SPVL


poly_model <- 
  
concave_model <- 
  
convex_model <- 
  
sigmoid_model <- 
  
max_model <- 
  


