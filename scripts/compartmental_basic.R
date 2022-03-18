#PBasic deterministic model

# only productively infected cells

library(deSolve)
library(tidyverse)

Basic_Model <- function(times, state, parms){
  X <- state['X']
  Y <- state['Y']
  V <- state['V']
  
  beta <- parms['beta'] #
  gamma <- parms['gamma'] # production rate target cells
  d_x <- parms['d_x'] # natural death target cells
  d_y <- parms['d_y'] # death of infected cells
  p <- parms['p'] # rate of production of virus from infected cells
  c <- parms['c'] # rate of clearance of free virus
  
  dX <- gamma - d_x*X - beta*V*X
  dY <- beta*V*X - d_y*Y
  dV <- p*Y - c*V
  
  res <- list(c(dX, dY, dV))
  
  return(res)
}

parameters <- c(beta = 2*10**-7,
                gamma = 10**5,
                d_x = 0.1,
                d_y = 0.5,
                p = 100,
                c = 5)

times <- seq(from = 0, to = 50 , by = 1)

X_0 <- parameters['gamma']/parameters['d_x']
Y_0 <- 0
V_0 <- 5

init_states <- c(X = X_0, Y = Y_0, V = V_0)

output <- ode(y = init_states , time = times, func = Basic_Model, parms = parameters, method = "rk4") %>% as.data.frame()

