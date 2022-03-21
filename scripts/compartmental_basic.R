#Basic deterministic model - discrete time



library(deSolve)
library(tidyverse)


Basic_Model <- function(times, state, parms){
  T <- state['T'] # target cells
  I <- state['I'] # actively infected
  V <- state['V'] # free virus
  N <- T + I + V 
  
  beta <- parms['beta'] #rate constant describing efficacy og free virus particle infectivity
  gamma <- parms['gamma'] # production rate target cells
  d_t <- parms['d_t'] # natural death target cells
  d_i <- parms['d_i'] # death of infected cells
  p <- parms['p'] # rate of production of virus from infected cells
  c <- parms['c'] # rate of clearance of free virus
  
  # define differential equations
  dT <- gamma - d_t*T - (beta*V*T/N)
  
  dI <- (beta*V*T/N) - d_i*I
  
  dV <- p*I - c*V
  
  res <- list(c(dT, dI, dV))
  
  return(res)
}


R0_basic <- function(parms){
  out <- (parms['beta']*parms['gamma']*parms['p']) / (parms['d_i']*parms['d_t']*parms['c']) %>% unname()
  return(out)
}

# As above, but including additional state for latently infected cells
Latent_Model <- function(times, state, parms){
  T <- state['T'] # target cells
  I <- state['I'] # actively infected
  L <- state['L'] # latently infected
  V <- state['V'] # free virus
  
  beta <- parms['beta'] #rate constant describing efficacy og free virus particle infectivity
  gamma <- parms['gamma'] # production rate target cells
  d_t <- parms['d_t'] # natural death target cells
  d_i <- parms['d_i'] # death of infected cells
  a_q <- parms['a_q'] # activation of latent cells
  p_q <- parms['p_q'] # proportion of infected cells that become latent (immediately upon infection)
  d_q <- parms['d_q'] # death of latent cells (currently == d_t)
  p <- parms['p'] # rate of production of virus from infected cells
  c <- parms['c'] # rate of clearance of free virus
  
  # define differential equations
  dT <- gamma - d_t*T - beta*V*T
  
  dI <- (1-p_q)*beta*V*T + a_q*L - d_i*I
  
  dL <- p_q*beta*V*T - d_q*L - a_q*L
  
  dV <- p*I - c*V
  
  res <- list(c(dT, dI, dL, dV))
  
  return(res)
}


# Adding a long-lived infected cell state with a slower death & release of virions

MultiPhasic_Model <- function(times, state, parms){
  T <- state['T'] # target cells
  I <- state['I'] # productively infected
  M <- state['M'] # long-lived infected
  L <- state['L'] # latently infected
  V <- state['V'] # free virus
  
  
  
  gamma <- parms['gamma'] # production rate target cells
  d_t <- parms['d_t'] # natural death target cells
  
  d_i <- parms['d_i'] # death of infected cells
  beta <- parms['beta'] #rate constant describing efficacy of free virus particle infectivity
  k_i <- parms['k_i'] # rate of production of virus from productively infected cells
  
  p_m <- parms['p_m'] # proportion of infected cells that become long lived (immediately upon infection)
  d_m <- parms['d_m'] # death of long lived infected cells
  k_m <- parms['k_m'] # rate of production of virus from long-lived infected cells
  
  a_q <- parms['a_q'] # activation of latent cells
  p_q <- parms['p_q'] # proportion of infected cells that become latent (immediately upon infection)
  d_q <- parms['d_q'] # death of latent cells
  
  
  c <- parms['c'] # rate of clearance of free virus
  
  # define differential equations
  dT <- gamma - d_t*T - beta*V*T
  
  dI <- (1-(p_q+p_m))*beta*V*T + a_q*L - d_i*I
  
  dM <- p_m*beta*V*T - d_m*L 
  
  dL <- p_q*beta*V*T - d_q*L - a_q*L
  
  dV <- k_i*I + k_m*M - c*V
  
  res <- list(c(dT, dI, dM, dL, dV))
  
  return(res)
}


# Define Parameters
parameters <- c(beta = 2*10**(-7), #per virion per day
                gamma = 10**5, #per microL per day
                d_t = 0.1,# per day
                d_i = 0.5, # per day
                #d_q = 1.36*10**(-3),
                #a_q = 3.6*10**(-2),
                #p_q = 0.1, 
                p = 100, #per cell per day
                c = 5) # per day

# Initialise timeframe and conditions
times <- seq(from = 0, to = 100 , by = 1)

T_0 <- 10000
I_0 <- 0
#L_0 <- 0
#M_0 <- 0
V_0 <- 100

init_states <- c(T = T_0, I = I_0, V = V_0)

#check whether infection will spread
R0_basic(parameters)
stopifnot(R0_basic(parameters)>1)

# Run simulation
output <- ode(y = init_states , time = times, func = Basic_Model, parms = parameters, method = "rk4") %>% as.data.frame()

output 



#function to plot data direct from model output

PlotODE <- function(results){
  require(ggplot2)
  
  results_wide <- results
  results_long <- gather(results, key = "state" , value = "value", 2:4) 
  plt <- ggplot(results_long) +
    geom_line(aes(x = time, y = value, colour = state)) +
    theme_classic()+
    facet_wrap(.~state)
  
  return(plt)
}

PlotODE(output)