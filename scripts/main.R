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
data <- read_delim('./data/pbio.1001951.s006.tsv', delim = '\t') %>%
  rename(SpVL = spVL) #Post processing

NSIM <- 100
################################### Run Base Models ###################################
source('./scripts/base_models.R')


################################### Estimate Heritability Under different Assumptions ###############################
# Fit 

linear_model_uw <- heritability_model # Linear model imported from base_models


prior2 <- prior(normal(1, 2), nlpar = "b1") +
  prior(normal(0, 2), nlpar = "b2")

concave_model_uw <- brms::brm(
  bf(recipient_log10SpVL ~  b1 * exp(b2 * transmitter_log10SpVL),  
     b1 + b2 ~1, 
     nl = TRUE),
  data = pop,
  prior = prior2
)


prior3 <- prior(normal(1, 2), nlpar = "b1") +
  prior(normal(0, 2), nlpar = "b2")

convex_model_uw <- brms::brm(
  bf(recipient_log10SpVL ~  b1 + b2*log(transmitter_log10SpVL),  
     b1 + b2 ~1, 
     nl = TRUE),
  data = pop,
  prior = prior3
)


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
  posterior_predict(linear_model_uw, .) %>%
  
  #Sample one value from posterior predictions per transmitter
  apply(., 2, sample, 1) %>% 
  
  # Bind predicted recipient SpVl with transmission pair characteristics
  cbind.data.frame(recipient_log10SpVL = ., transmitter_log10SpVL= log10(sim_donor),sim_recip_chars) %>%
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
  posterior_predict(linear_model_uw, .) %>%
  
  #Sample one value from posterior predictions per transmitter
  apply(., 2, sample, 1) %>% 
  
  # Bind predicted recipient SpVl with transmission pair characteristics
  cbind.data.frame(recipient_log10SpVL = ., transmitter_log10SpVL = log10(sim_donor), sim_recip_chars) %>%
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
  posterior_predict(linear_model_uw, .) %>%
  
  #Sample one value from posterior predictions per transmitter
  apply(., 2, sample, 1) %>% 
  
  # Bind predicted recipient SpVl with transmission pair characteristics
  cbind.data.frame(recipient_log10SpVL = ., transmitter_log10SpVL = log10(sim_donor)) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SpVL')}")) %>%
  `colnames<-` (str_remove(colnames(.), 'sim_')) %>% 
  cbind.data.frame(sim_recip_chars) %>%
  mutate(model = 'convex_uw') 


################################### Run Transmission Models ##################################

shcs_tm <- RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 1)%>%
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w')) 


linear_uw_tm <- c(RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 1),
                  RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 5),
                  RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 10),
                  RunParallel(populationmodel_acrossVL_Environment, pop$transmitter, w= 20)) %>%
  
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


################################### Extract Model Outputs ##################################

# Q1: Does SpVL in the transmitting partner confound the relationship between the number of founder 
# variants and viral load in the recipient?

shcs_transmitters <- shcs_tm  %>%
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
# transmitting partner, and how is this affect by the number of variants initiating infection in 
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


# Q3: How does the timing of transmission impact observations of the association between the number
# of founder variants and CD4+ T cell decline?





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



linear_uw_variants_timing <- linear_uw_tm %>%
  VariantPP(pop = linear_uw_pop) %>% 
  mutate(model = 'linear_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient)))) 

linear_uw_virions_timing <- linear_uw_tm %>% 
  VirionPP(pop = linear_uw_pop) %>% 
  mutate(model = 'linear_uw') %>%
  mutate(cd4_decline = predict(tolerance_model, newdata = data.frame(SpVL = log10(recipient))))


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


################################### Write to file ################################### 
source('./scripts/figures.R')

ggsave(plot = panel_1, filename = paste(figs_dir,sep = '/', "panel_1.jpeg"), device = jpeg, width = 14, height = 14) # Model Components
ggsave(plot = panel_2, filename = paste(figs_dir,sep = '/', "panel_2.jpeg"), device = jpeg, width = 14, height = 14) # Confounder
ggsave(plot = panel_3, filename = paste(figs_dir,sep = '/', "panel_3.jpeg"), device = jpeg, width = 14, height = 14) # Non - Linear
ggsave(plot = panel_4, filename = paste(figs_dir,sep = '/', "panel_4.jpeg"), device = jpeg, width = 14, height = 14) # Timing of transmission

# Supplementary plots
ggsave(plot = panel_s1, filename = paste(figs_dir,sep = '/', "panel_s1.jpeg"), device = jpeg, width = 14, height = 14) #Within-host dynamics
ggsave(plot = panel_s2, filename = paste(figs_dir,sep = '/', "panel_s2.jpeg"), device = jpeg, width = 14, height = 14) #Simulating transmitter population
ggsave(plot = panel_s3, filename = paste(figs_dir,sep = '/', "panel_s3.jpeg"), device = jpeg, width = 14, height = 14) #Heritability model fitting


###################################################################################################

