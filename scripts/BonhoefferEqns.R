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


# Allocate SpVL to two distributions according to maximum likelihood discriminant rule
# Deprecated due to overlap in transmitter/recipient distributions.
# KEEP CODE FOR FUTURE USE!

#SpVLClassifier <- function(VL_pair, transmitter, recipient){
  
  ## Transmitter and recipient mut be named vectors with mean an variance
  #stopifnot(all(is.numeric(transmitter) && is.numeric(recipient)))
  #stopifnot(all(c(names(transmitter), names(recipient)) %in% c('mean', 'var')))
  #stopifnot(all(length(transmitter) && length(recipient) == 2))
  
  # Assign parameter names
  #mu_t <- transmitter[['mean']]
  #v_t <- transmitter[['var']]
  
  #mu_r <- recipient[['mean']]
  #v_r <- recipient[['var']]
  
  # Determine bounds of ML Classifier
  #a <- 1/v_r - 1/v_t
 # b <- 2*mu_t/v_t - 2*mu_r/v_r
  #c <- (mu_r**2)/v_r - (mu_t**2)/v_t + 2*log(v_r/v_t)
  
  #ub <- (-b + sqrt(b**2 - 4*a*c))/(2*a)
  #lb <- (-b - sqrt(b**2 - 4*a*c))/(2*a)
  
  #print(ub)
  #print(lb)
  
  # Likelihood of being a transmitter
  
  
  # Classify VLs
  
  #res <- 
 
  #res[,2] <- ifelse( > res[,1] & res[,1] > ub,  'Transmitter', 'Recipient')
  
  #return(res)
#}

# Tidy pipe for the abov
#vlpairs_alt <- shcs_data_long %>% mutate(transmitter = case_when(
  #dnorm(log10_SpVL, mean = 4.37, sd = sqrt(0.1845)) > dnorm(log10_SpVL, mean = 4.356, sd = sqrt(0.3737)) ~ 1 ,
 # .default = 0)) %>%
 # mutate(recipient= case_when(
  #  dnorm(log10_SpVL, mean = 4.37, sd = sqrt(0.1845)) < dnorm(log10_SpVL, mean = 4.356, sd = sqrt(0.3737)) ~ 1 ,
  #  .default = 0))


## THIS ONE ##
vlpairs <- shcs_data %>% mutate(transmitter = case_when(
  dnorm(log10_SpVL_1, mean = 4.37, sd = sqrt(0.1845)) > dnorm(log10_SpVL_2, mean = 4.37, sd = sqrt(0.1845)) ~ 1 ,
  dnorm(log10_SpVL_1, mean = 4.37, sd = sqrt(0.1845)) < dnorm(log10_SpVL_2, mean = 4.37, sd = sqrt(0.1845)) ~ 2 ,
  .default = NA
))



vl_pair_long <- vlpairs %>% 
  select(c(log10_SpVL_1, log10_SpVL_2, transmitter, allocatetransmitter)) %>% pivot_longer(cols = starts_with('log10_SpVL'), values_to = 'log10SpVL') %>%
  mutate(role_ml = case_when(name == 'log10_SpVL_1' &  transmitter == 1 ~ 't',
                                 name == 'log10_SpVL_2' & transmitter == 2 ~ 't',
                                 .default = 'r'
                                 )) %>%
  mutate(role_max = case_when(name == 'log10_SpVL_1' &  allocatetransmitter == 1 ~ 't',
                             name == 'log10_SpVL_2' & allocatetransmitter == 2 ~ 't',
                             .default = 'r'
  ))


a <- ggplot(vl_pair_long)+
  geom_histogram(aes(x = log10SpVL, fill = role_ml))
  
b <- ggplot(vl_pair_long)+
  geom_histogram(aes(x = log10SpVL, fill = role_max))

c <- ggplot(vl_pair_long)+
  geom_histogram(aes(x = log10SpVL, fill = role_ml))

plot_grid(a,b)
# TEST
#mc <- 4.35
#vc <- 0.478**2

#t <- CalcTransmitter(mc,vc)

#r <- CalcRecipient(mc,vc)

#v <- rnorm(100, 4.5, sd = 2)

#test <- SpVLClassifier(VL = v, transmitter = t, recipient = r)

#df <- tibble(transmitter = rnorm(1000, mean = t[['mean']], sd = sqrt(t[['var']])),
#             recipient = rnorm(1000, mean = r[['mean']], sd = sqrt(r[['var']]))) %>%
#  pivot_longer(cols = everything(), values_to = 'SpVL', names_to = 'Position')

#ggplot(df)+
 # geom_histogram(aes(x = SpVL, fill = Position))+
#  geom_vline(xintercept = 3.739681) + 
#  geom_vline(xintercept = 5.132123) + 
 # facet_wrap(.~Position)
## END ##