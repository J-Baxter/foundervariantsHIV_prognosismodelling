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
source('./scripts/RunTM4Sim.R')
source('./scripts/global_theme.R')


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
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w')) %>%
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

linear_sim <- SimPop(pop, h2_model, 100)

linear_pred <- RunParallel(populationmodel_acrossVL_Environment, linear_sim$transmitter, w= 1) %>%
  lapply(., setNames, nm = c('variant_distribution','probTransmissionPerSexAct','transmitter',  'w'))


CleanUpOnAisle6 <- function(x){
  out <- lapply(x, function(i) any(is.na(i[['variant_distribution']]))) %>% unlist() %>% which()
  return(out)
}

linear_pred_mp <- linear_pred[-CleanUpOnAisle6(linear_pred)] %>%
  lapply(., cbind.data.frame) %>%
  do.call(rbind.data.frame,.) %>%
  rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
  dplyr::select(-contains('nparticles')) %>%
  dplyr::summarise(across(starts_with('V'), .fns =sum),.by = transmitter) %>%
  mutate(recipient = pop$recipient[-CleanUpOnAisle6(linear_pred)]) %>%
  pivot_longer(cols = starts_with('V'),
               names_to = 'variants',
               values_to = 'p') %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% 
           as.numeric())  %>%
  filter(variants <= 10) %>%
  mutate(model = 'linear', w = 1) %>%
  
  # Predict CD4 decline using fitted CD4 model
  mutate(cd4_decline = predict(cd4_model, newdata = data.frame(spVL = log10(recipient)))) 
  
  

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
  
###################################################################################################
################################### Variants vs Virions ###################################



# Generate a matrix of transmitter/founding particles and associated viral loads
FoundingVars <- function(spvl, sd, p_dists){
  
  #Round Log10 SPVL to nearest .5
  spvl_disc <- round(spvl/0.5)*0.5
  
  #Identify appropriate probability distribution
  lookup <- which(spvl_disc == seq(2,7,by=0.5))
  prob_var <- p_dists[,lookup]
  
  n_var <- sample(x = 1:33, 1, replace = T, prob = prob_var)
  
  within_host_genotypes <- rnorm(n = n_var, mean = spvl, sd = sd)
  
  m <- matrix(NA, nrow = 1, ncol = 33)
  
  m[1:length(within_host_genotypes)] <- within_host_genotypes
  
  
  return(m)
  
}

spvl <- seq(2,7,by=0.5)
prob_dists <- lapply(spvl, populationmodel_acrossVL_Environment) %>% sapply(., function(x) x[2:length(x)])

variant_matrix <- sapply(pop$transmitter_log10SPVL, FoundingVars, sd = 0.5, p_dists = prob_dists) %>% t()

# Sanity check
apply(variant_matrix, 1, function(x) sum(!is.na(x))) %>%
  all() %>%
  stopifnot()


