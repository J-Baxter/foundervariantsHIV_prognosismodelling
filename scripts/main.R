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


################################### Run Base Models ###################################
source('./scripts/base_models.R')


################################### Fit Non Linear Models (unweighted) ###################################

quad_model <- lm(recipient_log10SPVL ~ transmitter_log10SPVL+ I(transmitter_log10SPVL**2), data = pop) 

exp_model <- lm(recipient_log10SPVL ~ exp(transmitter_log10SPVL), data = pop) 


# Check Fit


################################### Run Simulations ###################################
NSIM <- 200
  
linear_sim <- SimPop(data$transmitter, h2_model, NSIM)

quad_sim <- SimPop(data$transmitter, quad_model, NSIM)

exp_sim <- SimPop(data$transmitter, exp_model, NSIM)


################################### Run Transmission Model on Simulations ###################################
linear_pred <- RunTM4Sim(linear_sim, modeltype = 'linear')

quad_pred <- RunTM4Sim(quad_sim, modeltype = 'quad')

exp_pred <- RunTM4Sim(exp_sim, modeltype = 'exp')


################################### Write to file ################################### To be changed
# population only
write_csv(pop, file = paste(results_dir, 'base_population.csv', sep = '/'))

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


