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


###################################################################################################
##### Compile C++ Code #####
sourceCpp("./scripts/withinhost_fraser.cpp")


###################################################################################################
#### Set seed ####
set.seed(4472)


###################################################################################################
#### Set initial states ####
init <- list()

###################################################################################################
#### Set parameters ####
parms <- list(lambda_cd4 = , #daily thymic production of new cd4
              lambda_cd8 = parms["lambda_cd8"], #daily thymic production of new cd8
              a_0 = parms["a_0"], #average rate of T cell activation per antigenic exposure
              mu = parms["mu"], #daily rate of non-antigen-driven homeostaatic T cell division
              x_s = parms["x_s"], #relative t cell pool size, below which t cell activation fails due to exhaustion
              mu_a = parms["mu_a"], # activated t cell division rate
              p_a = parms["p_a"], # average probability of an activated t cell successfully dividing in an HIV negative control
              beta = parms["beta"], # average per virion infection rate of an activated CD4+ T cell
              alpha = parms["alpha"], #death rate of productively infected cell (in the absence of CTL)
              alpha_latent = parms["alpha_latent"], #rate of reactivation of latent infected cells
              frac_latent = parms["frac_latent"], #proportion of successful infections that result in latency
              p = parms["p"], #maximum proliferation rate of anti-HIV CTLs
              d = parms["d"], #death rate of anti-HIV CTL
              z_0 = parms["z_0"], #pre-infection frequency of anti-HIV CTLs
              infectious_0 = parms["infectious_0"],
              bc_ratio = parms["b"], # ratio of viral production rate in productively infected cells and viral lifetime
              sigma = parms["sigma"], # max rat of CTL killing of HIV infected cells
              # N_PB = parms["N_PB"],
              theta_cd4 = parms["theta_cd4"], # average clearance rate in antigenic exposure model
              theta_cd8 = parms["theta_cd8"], # average exposure rate in antigenic exposure model
              e_cd4 = parms["e_cd4"],
              e_cd8 = parms["e_cd8"])


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