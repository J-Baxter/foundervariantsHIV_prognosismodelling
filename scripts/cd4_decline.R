# CD4 decline model
# from Disentangling Human Tolerance and Resistance Against HIV
# Regoes et al. PLoS Biology 2014
# CD4 decline is 'change in CD4+ T cells per microlitre per day'


cd4_data <- read_delim('./data/pbio.1001951.s006.tsv', delim = '\t' )

# -0.1 + 0.051*log10(V) - 0.017*(log10(V))**2
cd4_model <- lm(CD4.decline ~ spVL + I(spVL**2), data = cd4_data )
#why does this fit differently to Poly(spVL, 2)?
summary(cd4_model)
confint(cd4_model, level = 0.95)


CD4Decline <- function(log10SPVL,age){
  
  ##################################################
  # parameters from Regoes et al. PLoS Biology 2014
  ##################################################
  
  alpha_0 = -5.6e-3
  c = -1.6e-4
  
  deltaCD4 <- (alpha_0 + c*age)*(log10SPVL)**2
  
  return(deltaCD4)
} 

## END ##