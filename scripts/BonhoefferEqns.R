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
  mu_0 <- 4.46  # Gaussian approximation of Transmission Potential
  v_0 <- 0.96  # Gaussian approximation of Transmission Potential
  mu_e <- 3  # Distribution of environmental effects 
  v_e <- 0.5  # Distribution of environmental effects 
  v_t <- 0.15 # variance asssociated with sampling at transmission
  
  # Subtract environmental effects to calculate genotypic distribution
  m_c <- MC - mu_e
  v_c <- VC - v_e #problem is here because v_e > VC therefore v_c < 0
  
  # Apply transmission potential to infer donor genotype
  m_d <- (m_c*(v_e + v_0) + (mu_0 - mu_e)*v_c)/(v_c + v_e + v_0)
  v_d <- (v_c*(v_e + v_0))/(v_c + v_e + v_0)
  print(m_d)
  print(v_d)
  # sample from donor genotypes for recipient genotypes
  m_r <- m_d 
  v_r <- v_d + v_t
  
  # Re-distribute environmental variance for recipient phenotype
  MR <- m_r + mu_e
  VR <- v_r + v_e
  
  out <- c('mean' = MR,
           'var' = VR)
  
  return(out)
}


SpVLClassifier <- function(VL, transmitter, recipient){
  
  # Transmitter and recipient mut be named vectors with mean an variance
  stopifnot(all(is.numeric(transmitter) && is.numeric(recipient)))
  stopifnot(all(c(names(transmitter), names(recipient)) %in% c('mean', 'var')))
  stopifnot(all(length(transmitter) && length(recipient) == 2))
  
  # Assign parameter names
  mu_t <- transmitter[['mean']]
  v_t <- transmitter[['var']]
  
  mu_r <- recipient[['mean']]
  v_r <- recipient[['var']]
  
  # Determine bounds of ML Classifier
  a <- 1/v_r - 1/v_t
  print(a)
  b <- 2*mu_t/v_t - 2*mu_r/v_r
  print(b)
  c <- (mu_r**2)/v_r - (mu_t**2)/v_t + 2*log(v_r/v_t)
  print(c)
  
  ub <- (-b + sqrt(b**2 - 4*a*c))/(2*a)
  lb <- (-b - sqrt(b**2 - 4*a*c))/(2*a)
  
  print(ub)
  print(lb)
  
  # Classify VLs
  res <- matrix(nrow  = length(VL), ncol = 2)
  res[,1] <- VL
  
  res[,2] <- ifelse(lb > res[,1] & res[,1] > ub,  'Transmitter', 'Recipient')
  
  return(res)
}


# TEST
mc <- 4.35
vc <- 0.478**2


t <- CalcTransmitter(mc,vc)

r <- CalcRecipient(mc,vc)

v <- rnorm(100, 4.5, sd = 2)

test <- SpVLClassifier(VL = v, transmitter = t, recipient = r)

df <- tibble(transmitter = rnorm(1000, mean = t[['mean']], sd = sqrt(t[['var']])),
             recipient = rnorm(1000, mean = r[['mean']], sd = sqrt(r[['var']]))) %>%
  pivot_longer(cols = everything(), values_to = 'SpVL', names_to = 'Position')

ggplot(df)+
  geom_histogram(aes(x = SpVL, fill = Position))+
  geom_vline(xintercept = 3.739681) + 
  geom_vline(xintercept = 5.132123) + 
  facet_wrap(.~Position)
