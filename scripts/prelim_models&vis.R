###################################################################################################
###################################################################################################
# Modelling the association between the number of variants and viral load
# Script utilises methodology from Thompson et al 2019, https://github.com/robin-thompson/MultiplicityOfInfection
# implemented in R by Katie Atkins https://github.com/katiito/HIV_eclipse_models

# N Scenarios are presented:
#1. Null hypothesis
#2.
#3.
#4.


###################################################################################################
###################################################################################################
# Dependencies
library(tidyr)
library(dplyr)
library(ggplot2)

source("./scripts/populationdata_models.R")

GetRecipVL <- function(donor_loads, h2){
  n = length(donor_loads)
  ssr = sum((donor_loads-mean(donor_loads))^2) # sum of squared residuals
  e <- rnorm(n)
  e <- resid(lm(e ~ donor_loads))
  e <- e*sqrt((1-h2)/h2*ssr/(sum(e^2)))
  recip_loads <- donor_loads + e
  
  return(recip_loads)
}

###################################################################################################
# Data generation
# Normal distibution of log spvl from Hollingsworth Cohort (Link). mean = 4.39, sd = 0.84
# https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1000876
NPAIRS <- 100000
donor_spvl <- rnorm(NPAIRS, mean = 10^4.39, sd = 0.84)

# Heritability estimate
R2 <- 0.35

# Generate recipient data
recip_spvl <- GetRecipVL(donor_spvl, 0.35)

# Check calcualted R2  == input (to run in test script)
#summary(lm(recip_spvl ~ donor_spvl))$r.squared

# Data frame of paired viral loads
spvl_df <- cbind.data.frame(donor_spvl, recip_spvl)


###################################################################################################
# Probability that recipient infection is initiated by multiple founder variants
# applies populationmodel_fixedVL_Environment function written by Katie Atkins

prob_recip_multiple <- sapply(donor_spvl[1:5], populationmodel_fixedVL_Environment) # Check 0.01 > p > 0.0001


###################################################################################################
# Visualise probability that recipient infection is intitiated by multiple founders



ggplot() +
  geom_density(aes(x=diff), color = 'blue')

ggplot() +
  geom_density(aes(x=donor_spvl), color = 'blue') + 
  geom_density(aes(x=recip_spvl), color = 'red') + 
  geom_density(aes(x=recip_vl_mid), color = 'green') +
  #geom_density(aes(x=recip_vl_high), color = 'black') +
  theme_bw()+
  scale_y_continuous(limits = c(0,0.5), expand = c(0,0))+
  scale_x_continuous(limits = c(0,12), expand = c(0,0))+
  xlab('Log10  SPVL')

diff <-  rnorm(NPAIRS, mean = 0.84, sd = 1.28)


paired_spvl <- cbind.data.frame(donor = donor_spvl, recipient= donor_spvl+diff)

library(ggsci)
test_1 <-  ggplot(tidyr::gather(paired_spvl)) +
  stat_density(aes(x= value, colour = key), geom = 'line', position="identity") + 
  scale_colour_npg()+
  theme_bw()+scale_y_continuous(limits = c(0,0.5), expand = c(0,0))+
  scale_x_continuous(limits = c(0,12), expand = c(0,0))+
  xlab('Log10 SPVL') + theme(legend.position = c(0.8, 0.8), legend.background = element_blank())


test_2 <- ggscatter(paired_spvl, x='donor', y = 'recipient', add = "reg.line", shape = 4) +
  theme_bw()+
  xlab('Log10 Donor SPVL') +
  ylab('Log10 Recipient SPVL') +
  theme(legend.position = "none") + 
  scale_y_continuous(limits = c(0,12), expand = c(0,0))+
  scale_x_continuous(limits = c(0,12), expand = c(0,0))+
  stat_cor(aes(label = paste(..rr.label..)), label.x = 8, label.y = 11, digits = 3)

test_3 <- ggplot() +
  stat_density(aes(x=diff), geom = 'line', position="identity",  color = 'darkgreen') + 
  theme_bw()+scale_y_continuous(limits = c(0,0.5), expand = c(0,0))+
  xlab('difference')

test_panel <- cowplot::plot_grid(test_1, test_3, test_2, ncol = 3, labels = 'AUTO')
test_panel


fit <- lm(recip_vl_mid ~ donor_spvl, vl_df) %>% summary()

ggplot()+
  geom_point(aes(x = donor_spvl,y = recip_spvl), data = vl_df) 
  
 