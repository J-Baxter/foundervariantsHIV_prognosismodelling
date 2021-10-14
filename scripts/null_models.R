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
library(scales)
library(ggpmisc)
library(parallel)

source("./scripts/populationdata_models.R")
source("./scripts/misc_functions.R")

# Calculate recipient set point viral loads from a given population of donor set point viral loads.
# Back calculation from a set correlation coefficient (r2) and sum of squared residuals
# 'Parent ~ Offspring' regression R2 analogous to heritability of viral load. 
GetRecipVL <- function(donorload, h2 = 0.33){
  n <- length(donorload)
  ssr <- sum((donorloads-mean(donorload))^2) # sum of squared residuals
  e <- rnorm(n)
  e <- resid(lm(e ~ donorload))
  e <- e*sqrt((1-h2)/h2*ssr/(sum(e^2)))
  recipload <- donorload + e
  
  return(recipload)
}


# Infer weighting for a given viral load 
WeightPDF <- function(viralload){
  alpha = -3.55
  sigma <- 0.78/(sqrt(1 - ((2*alpha^2)/(pi*(1 + alpha^2)))))
  mu <-  4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha^2)))
  weight <-  (2/sigma)*dnorm((log10(viralload) - mu)/sigma)*pnorm(alpha*(log10(viralload) - mu)/sigma)
  return(weight)
} 


###################################################################################################
#Set seed
set.seed(4472)

# Data generation
# Normal distibution of log spvl from Amsterdam Cohort (Link). mean = 4.39, sd = 0.84
# https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1000876
NPAIRS <- 1000
donor_logspvl <- rnorm(NPAIRS,  mean = 4.39, sd = 0.84)

# Heritability estimate
R2 <- 0.33

# Generate recipient data
recip_logspvl <- GetRecipVL(donor_logspvl, R2)
stopifnot(min(recip_logspvl)>0)

# Check calcualted R2  == input (to run in test script)
#summary(lm(recip_spvl ~ donor_spvl))$r.squared

donor_spvl <- 10^donor_logspvl
recip_spvl <- 10^recip_logspvl

# Data frame of paired viral loads
spvl_df <- cbind.data.frame(donor_spvl = donor_spvl, 
                            recip_spvl = recip_spvl) 


###################################################################################################
# Probability that recipient infection is initiated by multiple founder variants
# applies populationmodel_fixedVL_Environment function written by Katie Atkins

prob_recip_multiple <- RunParallel(populationmodel_fixedVL_Environment, (donor_spvl[1:50])) %>%
  do.call(rbind.data.frame, .)
  

combined_data <- cbind.data.frame(spvl_df[1:50,], prob_recip_multiple)
head(combined_data)

###################################################################################################
# Probability that recipient infection is multiple founder, given a certain viral load


test_prob <- sapply(recip_spvl, function(x) WeightPDF(x)/sum(WeightPDF(recip_spvl)))
  
  
#sample args:
# x = either a vector of one or more elements from which to choose
# n = a positive number, the number of items to choose from
# size = a non-negative integer giving the number of items to choose
# prob = a vector of probability weights for obtaining the elements of the vector being sampled.
set_spvl <- 10^(1:9)
sample(set_spvl, size = 100, prob = test_prob, replace = T)



###################################################################################################
# Visualise probability that recipient infection is intitiated by multiple founders

# Fig 1a
fig_1a <- ggplot(spvl_df, aes(x = donor_spvl, recip_spvl)) +
  geom_point() +
  scale_x_log10(limits = c(1, 10^10),
                expand = c(0,0),
                name = 'Donor SPVL (log10)',
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(name = 'Recipient SPVL (log10)',
                limits = c(1, 10^10),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  geom_smooth(method = lm) + 
  stat_poly_eq(formula = y ~ x)
  
  
# Fig 1b
fig_1b <- ggplot(combined_data, aes(x = donor_spvl, multiple_founder_proportion))+
  geom_point()+
  scale_x_log10(name = 'Donor SPVL (log10)',
                limits = c(1, 10^10),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,0.6),
                     breaks = seq(0, 0.6, by = 0.2)) +
  theme_bw() +
  geom_smooth(method = lm)
  

# Fig 1c

