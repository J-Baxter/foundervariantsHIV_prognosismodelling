# CD4 decline model
# from Disentangling Human Tolerance and Resistance Against HIV
# Regoes et al. PLoS Biology 2014
# CD4 decline is 'change in CD4+ T cells per microlitre per day'



CD4Decline <- function(log10SPVL,age,sex){
  
  ##################################################
  # parameters from Regoes et al. PLoS Biology 2014
  ##################################################
  
  alpha_f = -0.005866
  eta_m = 0.000491
  c_f = -0.000167
  z_m = 0.000001
  
  # Participant sex as dummy variable
  sex <- ifelse(sex =='M', 1, 0)
  
  deltaCD4 <- (alpha_f + eta_m*sex + (c_f + z_m*sex)*age)*(log10SPVL)**2
  
  return(deltaCD4)
} 

## END ##