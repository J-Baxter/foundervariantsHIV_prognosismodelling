# CD4 decline model
# from Disentangling Human Tolerance and Resistance Against HIV
# Regoes et al. PLoS Biology 2014
# CD4 decline is 'change in CD4+ T cells per microlitre per day'
# Incorporates terms for age and recipient sex

ToleranceModel <- function(log10SPVL, age, sex){
  
  male <- ifelse(sex == 'M', 1, 0)
  ##################################################
  # parameters from Regoes et al. PLoS Biology 2014
  ##################################################
  
  #alpha_0 = -5.6e-3
  #c = -1.6e-4
  #deltaCD4 <- (alpha_0 + c*age)*(log10SPVL)**2
  
  alpha_0_F = -5.866e-3 # tolerance at birth for females
  eta_0_M = 4.91e-4 # sex difference of tolerance at birth
  c_F = - 1.67e-4 # tolerance change per life years for females
  z_M = -1e-6 # sex difference between the change of tolerance per life year
  
  deltaCD4 <- (alpha_0_F + male*eta_0_M + age * (c_F + male*z_M)) * (log10SPVL)**2
  
  return(deltaCD4)
} 

## END ##