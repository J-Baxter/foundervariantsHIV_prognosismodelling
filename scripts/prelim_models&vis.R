###################################################################################################
###################################################################################################
# Preliminary models and visualisation of published data concerning the relationship between founder
# variant multiplicity, the decline in CD4+ T-Cells and set point viral load

###################################################################################################
###################################################################################################
# Dependencies
library(tidyr)
library(lme4)
library(dplyr)
library(metafor)
library(ggplot2)
library(reshape2)
library(ggplot2)


# Donor data generation
# Normal distibution of log spvl from Hollingsworth Cohort (Link). Mean = 4.35, SD = 0.47 
DONORS <- 100000

donor_spvl <- rnorm(DONORS, mean = 4.39, sd = 0.84) 

ggplot() + geom_density(aes(x=donor_spvl), color = 'blue') + theme_bw()

#Inclusion of other covariates
#Sex
sex <- c('m', "f")
donor_sex <- sample(sex, DONORS, replace = T)

#Age
age <- c(NA, '15-24', '25-29', '30-39', '40-64')
donor_age <- sample(age, DONORS, replace = T)

#GUD
gud <- c(NA, 'present', "absent")
donor_gud <- sample(gud, DONORS, replace = T)

#Subtype
subtype <- c(NA, 'recombinant', 'a', 'c', 'd')
donor_subtype <- sample(subtype, DONORS, replace = T)

donor_data <- cbind.data.frame(donor_spvl, donor_sex, donor_age, donor_gud, donor_subtype)

# Null model (Hollingsworth)
# Predict recipient SPVL from donor SPVL using GLM fitted in Hollingsworth et al. to generate
# Null distribution of recipient SPLV


