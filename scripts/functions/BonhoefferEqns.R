# Linear equations from Bonhoeffer et al. 2015


CalcTransmitter <- function(MC,VC){
  # Variables
  # MC <- Carrier phenotype
  # VC <- Carrier variance
  mu_0 <- 4.46  # Gaussian approximation of Transmission Potential
  v_0 <- 0.96  # Gaussian approximation of Transmission Potential

 
  #Infer transmitter phenotype
  MT <- (MC*v_0 + mu_0*VC)/(v_0 + VC)
  VT <- (VC*v_0)/(v_0 + VC)
  
  out <- c('mean' = MT,
           'var' = VT)
  
  return(out)
}


CalcRecipient <- function(MC, VC){
  
  # Variables
  # MC <- Carrier phenotype
  # VC <- Carrier variance
  h2 <- 0.33
  mu_0 <- 4.46  # Gaussian approximation of Transmission Potential
  v_0 <- 0.96  # Gaussian approximation of Transmission Potential
  mu_e <- 3 # Distribution of environmental effects (NB This parameter value does not impact MR)
  v_t <- 0.15 # variance asssociated with sampling at transmission
  
  # Assume populations are at equilibrium to calculate environmental variance
  v_e <-  VC - VC*h2  
  
  # Subtract environmental effects to calculate genotypic distribution
  m_c <- MC - mu_e
  v_c <- VC - v_e #problem is here because v_e > VC therefore v_c < 0
  
  # Apply transmission potential to infer donor genotype
  m_d <- (m_c*(v_e + v_0) + (mu_0 - mu_e)*v_c)/(v_c + v_e + v_0)
  v_d <- (v_c*(v_e + v_0))/(v_c + v_e + v_0)

  # sample from donor genotypes for recipient genotypes
  m_r <- m_d 
  v_r <- v_d + v_t
  
  # Re-distribute environmental variance for recipient phenotype1
  MR <- m_r + mu_e
  VR <- v_r + v_e
  
  out <- c('mean' = MR,
           'var' = VR)
  
  return(out)
}


## END ##