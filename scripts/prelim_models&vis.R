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
NPAIRS <- 100000

donor_spvl <- rnorm(NPAIRS, mean = 4.39, sd = 0.84) 

ggplot() + geom_density(aes(x=donor_spvl), color = 'blue') + theme_bw()

# Recipient covariates required for model
# Sex
# 1 = male, 0 = female
recipient_sex <- rbinom(NPAIRS, 1, 0.5) %>% cbind.data.frame() #altered to represent bias in transmission direction?
colnames(recipient_sex) <- 'sex'

# Age
age <- c(NA, '15-24', '25-29', '30-39', '40-64') %>% as.factor()
recipient_age <- sample(age, NPAIRS, replace = T) %>% data.table() %>% one_hot(sparsifyNAs = TRUE)
colnames(recipient_age) <- c('a1', 'a2', 'a3', 'a4')

# GUD
gud <- c(NA, 'absent', 'present') %>% as.factor()
recipient_gud <- sample(gud, NPAIRS, replace = T) %>% data.table() %>% one_hot(sparsifyNAs = TRUE)
colnames(recipient_gud) <- c('gud_n', 'gud_p')

# Subtype
subtype <- c(NA, 'a', 'c', 'd','recombinant') %>% as.factor()
recipient_subtype <- sample(subtype, NPAIRS, replace = T) %>% data.table() %>% one_hot(sparsifyNAs = TRUE)
colnames(recipient_subtype) <- c('sub_a', 'sub_c', 'sub_d', 'sub_r')

# Position in chain
chain <- c('index', 'second') %>% as.factor()
recipient_chain <- sample(chain, NPAIRS, replace = T) %>% data.table() %>% one_hot(sparsifyNAs = TRUE)
colnames(recipient_chain) <- c('chain_index', 'chain_second')

recipient_covar <- cbind.data.frame(recipient_sex, recipient_age, recipient_gud, recipient_subtype, recipient_chain)

# Null model (Hollingsworth)
# Predict recipient SPVL from donor SPVL using GLM fitted in Hollingsworth et al. to generate
# Null distribution of recipient SPLV

NullModel <- function(donor, recipient){
  
  recip <- donor - 0.15*recipient$sex - 0.15*recipient$a1 - 0.32*recipient$a2 + 0.35*recipient$a3 + 0.73*recipient$a4 +
    0.42*recipient$gud_n + 0.62*recipient$gud_p + 0.068*recipient$sub_a + 1.85*recipient$sub_c + 0.39*recipient$sub_d +
    0.78*recipient$sub_r + 0.42*recipient$chain_index + 0.00*recipient$chain_second
  
  return(recip)
}



#recip_spvl = i - 0.15*sex - 0.15*a1 - 0.32*a2 + 0.35*a3 + 0.73*a4 + 0.42*gud_n + 0.62*gud_p + 0.068*sub_a + 1.85*sub_c + 0.39*sub_D + 0.78*sub_R + 0.42*chain_index + 0.00*chain_second



# Calculate difference between donor-recipient pairs


# Null Panel
# A) Distribution of donor SPVL
# B) Distribution of recipient SPVL
# C) Donor-recipient paired difference in SPVL



