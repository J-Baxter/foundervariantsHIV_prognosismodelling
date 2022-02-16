# Stochastic Implementation of Perelson Model (Single Founder)
# States:
# - V: free virus
# - T: Uninfected cells
# - I: Productively infected cells
# - L: Latently infected cells
# - M: Long-lived infected cells


# Parameters (per day):
# - d_T: natural death rate (per cell)
# - beta: per cell infection prob
# - gamma: production rate of target cells
# - d_I: death of actively infected cells due to virus (per cell) (delta > d_T to reflect shorter lifespan)
# - d_L: death of latently infected cells due to virus (per cell)
# - alpha: activation of latently infected cells
# - p: contsant production rate (per cell) of free virus by infected cell
# - c: clearance rate (per virus) of free virus from circulation


# Dependencies
library(tidyverse)
library(ggplot2)


# Event function
# -init_state: the initial state of compartmental model
#  (a named vector containing T, I, L, M, V)
# -params: the model parameters
#   (a named vector containing d_T, d_L, d_I, alpha, gamma, beta, p & c)
#  -t: the end time

Perelson_gillespie <- function(init_state, params, tf){
  
  time <- 0
  
  # assign parameters to easy access variables
  alpha <- params['alpha']
  beta <- params["beta"]
  gamma <- params["gamma"]
  d_T <- params["d_T"]
  d_I <- params["d_I"]
  d_L <- params["d_L"]
  p <- params["p"]
  c <- params["c"]
  
  # assign states to easy access variables
  
  T <- init_state["T"]
  I <- init_state["I"]
  L <- init_state["L"]
  M <- init_state["M"]
  V <- init_state["V"] # virus included as a state?
  
  N <- T + I + L + M + V
  
  # create results data frame
  results_df <- data.frame(time = 0, t(init_state))
  
  # run loop until end time is reached
  
  while (time < tf){
    
    # update current rates (Gillespie direct algorithm takes a set of possible events and a rate 
    # at which each outcome occurs)
    
    rates <- c()
    rates['susceptible_birth'] <- gamma * T/N #check for these whether denominator is necessary
    rates['susceptible_death'] <- d_T * T/N
    rates['infection_active'] <- beta * I * L / N
    rates['infection_latent'] <- beta * T * L / N
    rates['infection_longlived'] <- beta * I * M / N
    rates['activation'] <- alpha * L/N
    rates['death_active'] <- d_I * I/N
    rates['death_latent'] <- d_L * L/N
    rates['death_longlived'] <- d_M * M/N
    rates['viral_budding'] <-  p * M * I / N
    rates['viral_clearance'] <- c * V
    
    if (sum(rates) > 0) {
      
      time <- time + rexp(n = 1, rate = sum(rates))
      
      if (time <= tf){
        cumulative_rates <- cumsum(rates)
        
        event_type <- runif(n = 1, min = 0, max = sum(rates))
        
        if(type < cumulative_rates['susceptible_birth']){
          T <- T + 1
        }
        
        else if(type < cumulative_rates['susceptible_death']){
          T <- T - 1
        }
        
        else if(type < cumulative_rates['infection_active']){
          T <- T - 1
          I <- I + 1
        }
        
        else if(type < cumulative_rates['infection_latent']){
          T <- T - 1
          L <- L + 1
        }
        
        else if(type < cumulative_rates['infection_longlived']){
          T <- T - 1
          M <- M + 1
        }
        
        else if(type < cumulative_rates['activation']){
          L <- L - 1
          I <- I + 1
        }
        
        else if(type < cumulative_rates['death_active']){
          I <- I - 1
        }
        
        else if(type < cumulative_rates['death_latent']){
          L <- L - 1
        }
        
        else if(type < cumulative_rates['death_longlived']){
          M <- M - 1
        }
        
        else if(type < cumulative_rates['viral_budding']){
          V <- V + 1
        }
        
        else if(type < cumulative_rates['viral_clearance']){
          V <- V - 1
        }
      }
      else {
        time <- tf
      }
    }
    else {
      time <- tf
    }
    results_df <- rbind(results_df, c(time = time, T=T. I=I, L=L, M=M, V=V))
  }
  return(results_df)
}


# Set initial states
init_values <- c(T = ,
                 I = ,
                 L = ,
                 M = ,
                 V = )

# Set parameters
params <- c(alpha = , 
            beta = , 
            gamma = , 
            d_T = , 
            d_I = , 
            d_L = , 
            p = , 
            c = )

tmax <- 100

r <- Perelson_gillespie(init_state = init_values, params = params, tf = tmax)