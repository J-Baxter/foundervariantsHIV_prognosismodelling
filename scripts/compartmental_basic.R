###################################################################################################
#Basic deterministic models, single variant infection
###################################################################################################

# States (init)
# T <- target cells
# I <- actively infected
# V <- free virus
# M <- long-lived infected cells
# L <- latently infected cells

# Parameters (parms)
# gamma <- production rate target cells per microL per day
# d_t <- natural death target cells per day

# d_i <- death of infected cells per day
# beta <- rate constant describing efficacy of free virus particle infectivity per virion per day
# k_i <- rate of production of virus from productively infected cells

# p_m <- proportion of infected cells that become long lived (immediately upon infection)
# d_m <- death of long lived infected cells
# k_m <- rate of production of virus from long-lived infected cells

# a_q <- activation of latent cells
# p_q <- proportion of infected cells that become latent (immediately upon infection)
# d_q <- death of latent cells

# c <- virus clearance from blood plasma (virions/ml/day)


###################################################################################################
# Dependencies
library(deSolve)
library(tidyverse)


###################################################################################################
# Model Functions

# Basic Model: productively infected cells only
Basic_Model <- function(times, state, parms){
  with(as.list(c(parms,state)), {
    
    #ODEs
    dT <- gamma - d_t*T - beta*V*T
    
    dI <- beta*V*T - d_i*I
    
    dV <- k_i*I - c*V
    
    list(c(dT, dI, dV))
    
  })
}


# Latent Model: including additional state for latently infected cells
Latent_Model <- function(times, state, parms){
  with(as.list(c(parms,state)), {
    
    #ODEs
    dT <- gamma - d_t*T - beta*V*T
    
    dI <- (1-p_q)*beta*V*T + a_q*L - d_i*I
    
    dL <- p_q*beta*V*T - d_q*L - a_q*L
    
    dV <- k_i*I - c*V
    
    list(c(dT, dI, dL, dV))
    
  })
}


# Adding a long-lived infected cell state with a slower death & release of virions
MultiPhasic_Model <- function(times, state, parms){
  with(as.list(c(parms,state)), {
    
    #ODEs
    dT <- gamma - d_t*T - beta*V*T
    
    dI <- (1-(p_q+p_m))*beta*V*T + a_q*L - d_i*I
    
    dM <- p_m*beta*V*T - d_m*L 
    
    dL <- p_q*beta*V*T - d_q*L - a_q*L
    
    dV <- k_i*I + k_m*M - c*V
    
    list(c(dT, dI, dM, dL, dV))
    
    
  })
}


# Calculate R0
# Currently only for basic model
R0 <- function(parms, modeltype = 'basic'){
  if (modeltype == 'basic'){
    out <- (parms['beta']*parms['gamma']*parms['p']) / (parms['d_i']*parms['d_t']*parms['c'])
  }else if (modeltype == 'latent'){
    #out <- (parms['beta']*parms['gamma']*parms['p']) / (parms['d_i']*parms['d_t']*parms['c']) %>% unname()
  }else if (modeltype == 'multiphasic'){
    # out <- (parms['beta']*parms['gamma']*parms['p']) / (parms['d_i']*parms['d_t']*parms['c']) %>% unname()
  }
  
  return(unname(out))
}


###################################################################################################
# Specify initial conditions and parameters
parms <- c(gamma = 10**5,
           d_t = 0.1,
           
           d_i = 0.5,
           beta = 2*10**(-7),
           k_i = 100, 
           
           p_m = ,
           d_m = ,
           k_m = ,
           
           p_q = 0.1,
           a_q = 3.6*10**(-2),
           d_q = 1.36*10**(-3),
           
           c = 5) 


times <- seq(0,50, length.out = 500)


T_0 <- unname(parms['gamma']/parms['d_t']) # target cells at equilibrium at t0.
I_0 <- 0
L_0 <- 0
M_0 <- 0
V_0 <- 100


init <- c(T = T_0,
          I = I_0,
          L = L_0,
          M = M_0,
          V = V_0)


###################################################################################################
# Sanity check whether infection will spread
R0(parameters)
stopifnot(R0(parameters)>1)


###################################################################################################
# Run simulation
output <- ode(y = init, time = times, func = Basic_Model, parms = parms, method = "lsoda") %>% as.data.frame()


###################################################################################################
# Basic plot, faceted by compartment

PlotODE <- function(results){
  require(ggplot2)
  
  results_wide <- results
  results_long <- gather(results, key = "state" , value = "value", 2:4) 
  plt <- ggplot(results_long) +
    geom_line(aes(x = time, y = value, colour = state)) +
    theme_classic()+
    facet_wrap(.~state, scales = 'free')
  
  return(plt)
}

PlotODE(output)