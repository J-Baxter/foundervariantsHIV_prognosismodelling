# Function to allocate viral load to transmitter distribuion according to different
# methods:
# Distribution of transmitter and recipient viral loads inferred from Bonhoeffer Eqns,
# (Bonhoeffer et al. 2015) with a starting 'carrier' distribution approximated by a 
# Normal distribution N(4.35, 0.478) fitted to the Amsterdam Seroconverters 
# study (Fraser et al. 2007).

# The function, ClassifyVL, should be applied 'rowwise' and generates a new column with 
# the number corresponding to the allocated transmitter SpVL. The following methods are 
# implemented:
#   1) Random allocation
#   2) Allocate transmitter role to maximum viral load within each pair
#   3) Allocate transmitter role to the viral load with the greatest likelihood of belonging to
#     the transmitter distribution.
#   4) Calculate the posterior probability of a viral load being a transmitter/recipient and 
#   simulate allocations.

##################################### Dependencies #####################################
source('./scripts/BonhoefferEqns.R')


################################################################################################### 
AllocateTransmitter <- function(VL1, VL2, transmitter, recipient, method){
  
  # Transmitter and recipient mut be named vectors with mean an variance
  stopifnot(all(is.numeric(transmitter) && is.numeric(recipient)))
  stopifnot(all(c(names(transmitter), names(recipient)) %in% c('mean', 'var')))
  stopifnot(all(length(transmitter) && length(recipient) == 2))
  
  # Assign parameter names
  mu_t <- transmitter[['mean']]
  sd_t <- sqrt(transmitter[['var']])
  
  mu_r <- recipient[['mean']]
  sd_r <- sqrt(recipient[['var']])
  
  VL_pair <- c(VL1, VL2)
  
  if(method == 'random'){
    # Random allocation of transmitter VL
    transmitter <- sample(c(1,2), 1)
    
  } else if (method == 'max'){
    # Transmitter is the maximum value of VL1, VL2
    transmitter <- which.max(VL_pair)
    
  }else if (method == 'ml'){
    # Transmitter is the VL with highest likelihood of being in the transmitter distribution
    likVL1 <- dnorm(VL_pair[[1]], mu_t, sd_t)
    likVL2 <- dnorm(VL_pair[[2]], mu_t, sd_t) 
    
    transmitter <- which.max(c(likVL1, likVL2))
    
  }else if (method == 'bayes'){
    #warning('Global sensitivity analysis advised.')
    
    # Transmitter is asigned according to the probability that VL is transmitter 
    lik_VL1isT <- dnorm(VL_pair[[1]], mu_t, sd_t) * dnorm(VL_pair[[2]], mu_r, sd_r)
    lik_VL2isT <- dnorm(VL_pair[[2]], mu_t, sd_t) * dnorm(VL_pair[[1]], mu_r, sd_r)
    
    # Posterior probability 
    Q <- lik_VL1isT/(lik_VL2isT + lik_VL1isT)
    
    r <- runif(1)
    
    transmitter <- ifelse(r<Q, 1, 2) # if R  < Q then allocate transmitter identity to A otherwise, allocate it to B.
    
    }
  
  return(transmitter)
}


##################################### Bounds of ML Classifier #####################################
# Determine bounds of ML Classifier by solving quadratic equation
# a <- 1/v_r - 1/v_t
# b <- 2*mu_t/v_t - 2*mu_r/v_r
# c <- (mu_r**2)/v_r - (mu_t**2)/v_t + 2*log(v_r/v_t)

#ub <- (-b + sqrt(b**2 - 4*a*c))/(2*a)
#lb <- (-b - sqrt(b**2 - 4*a*c))/(2*a)


##################################### Testing #####################################
#mc <- 4.35
#vc <- 0.478**2

#t <- CalcTransmitter(mc,vc)
#r <- CalcRecipient(mc,vc)

#vlpairs_test <- shcs_data %>%
#  rowwise() %>%
#  mutate(transmitter = ClassifyVL(log10_SpVL_1, log10_SpVL_2, transmitter = t, recipient = r, method = 'bayes'))

##################################### END #####################################