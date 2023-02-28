###################################################################################################
# Script for running multiple pathogen SIR models
# Imports model functions from scripts within this directory

# Future Plans:
# 1. to test over a range of parameter values and combinations
# 2. stochastic implementation with random introduction of 'second strain' 


###################################################################################################
##### Load dependencies #####
source('./scripts/dependencies.R')
source('./scripts/parms.R')
source('./scripts/ONS.R')
source('./scripts/init.R')
source('./scripts/event_func.R')

df <- read.csv('./variantsim_app/default_parms.csv') 
library(Rcpp)
library(deSolve)

###################################################################################################
##### Compile C++ Code #####
sourceCpp("./scripts/withinhost_fraser.cpp")


###################################################################################################
#### Set seed ####
set.seed(4472)


###################################################################################################
#### Set initial states ####
init <- c(cd4 =0.66, 
             cd8 = 0.33, 
             cd4_activated = 0, 
             cd8_activated = 0, 
             infectious_active = 0, 
             infectious_latent = 0, 
             z = 0, 
             v = 0.1, 
             k_cd4 = 0, 
             k_cd8 = 0, 
             dpk_cd4 = 0, 
             dpk_cd8 =0)



###################################################################################################
#### Set parameters ####
parms <- list(lambda_cd4 = (10**(-4)/6)*5, #daily thymic production of new cd4
              lambda_cd8 = 10**(-4)/6, #daily thymic production of new cd8
              a_0 = 1/10**-4, #average rate of T cell activation per antigenic exposure
              mu = 0.01, #daily rate of non-antigen-driven homeostaatic T cell division
              x_s = 0.05, #relative t cell pool size, below which t cell activation fails due to exhaustion
              mu_a = 1, # activated t cell division rate
              p_a = 0.55, # average probability of an activated t cell successfully dividing in an HIV negative control
              beta = 1/754, # average per virion infection rate of an activated CD4+ T cell
              alpha = 1/0.5, #death rate of productively infected cell (in the absence of CTL)
              alpha_latent = 1/0.01, #rate of reactivation of latent infected cells
              frac_latent = 10**-5, #proportion of successful infections that result in latency
              p = 1/0.5, #maximum proliferation rate of anti-HIV CTLs
              d = 1/0.05, #death rate of anti-HIV CTL
              z_0 = 10**-06, #pre-infection frequency of anti-HIV CTLs
              infectious_0 = 10**-3.5,
              bc_ratio = 265, # ratio of viral production rate in productively infected cells and viral lifetime
              sigma = 1/(10**4), # max rat of CTL killing of HIV infected cells
              # N_PB = parms["N_PB"],
              theta_cd4 = 1/0.02, # average clearance rate in antigenic exposure model
              theta_cd8 = 1/0.02, # average exposure rate in antigenic exposure model
              e_cd4 = 1/0.1,
              e_cd8 = 1/0.1)


###################################################################################################
# Specify simulation time (units = days)
times <- seq(0,365, by = 1)


###################################################################################################
# Initial V1 Runs
start <- Sys.time() #Runtime ~ 0.5 seconds for compiled C++

v1_res <- ode(y = init,
              time = times, 
              func = Fraser_cpp_model,
              parms = parms, 
              method = "lsoda" ) %>% as.data.frame() %>% `colnames<-` (c('time', names(init)))

duration <- Sys.time()-start
duration 


###################################################################################################
# Output to File