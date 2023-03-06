# CD4 decline model
# from Disentangling Human Tolerance and Resistance Against HIV
# Regoes et al. PLoS Biology 2014
# CD4 decline is 'change in CD4+ T cells per microlitre per day'

cd4_data <- read_delim('./data/pbio.1001951.s006.tsv', delim = '\t' )


cd4_plt <- ggplot(cd4_data, aes(x = spVL, y=CD4.decline))+
  geom_point(colour = '#ef654a', shape= 4) +
  scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**1)),  #
                     expand = c(0,0),
                     limits = c(-2,2))+
  scale_x_continuous(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                     expand = c(0,0),
                     )+
  
  geom_smooth(method = 'lm', formula = 'y ~ x + I(x**2)', colour = 'black') + 
  my_theme + 
  annotation_logticks(sides = 'b') 

setEPS()
postscript( paste(figs_dir, 'cd4_plt.eps', sep = '/'), width = 10, height = 10)
cd4_plt 
dev.off()


# -0.1 + 0.051*log10(V) - 0.017*(log10(V))**2
cd4_model <- lm(CD4.decline ~ spVL + I(spVL**2), data = cd4_data )
summary(cd4_model)
confint(cd4_model, level = 0.95)

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