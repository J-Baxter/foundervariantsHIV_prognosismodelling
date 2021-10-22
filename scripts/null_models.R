###################################################################################################
###################################################################################################
# Modelling the association between the number of variants and viral load
# Script utilises methodology from Thompson et al 2019, https://github.com/robin-thompson/MultiplicityOfInfection
# implemented in R by Katie Atkins https://github.com/katiito/HIV_eclipse_models

###################################################################################################
###################################################################################################
# Dependencies
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggpmisc)
library(parallel)
library(cowplot)

source("./scripts/populationdata_models.R")
source("./scripts/misc_functions.R")

# Calculate recipient set point viral loads from a given population of donor set point viral loads.
# Back calculation from a set correlation coefficient (r2) and sum of squared residuals
# 'Parent ~ Offspring' regression R2 analogous to heritability of viral load. 
GetRecipVL <- function(donorload, h2 = 0.33){
  n <- length(donorload)
  ssr <- sum((donorload-mean(donorload))^2) # sum of squared residuals
  e <- rnorm(n)
  e <- resid(lm(e ~ donorload))
  e <- e*sqrt((1-h2)/h2*ssr/(sum(e^2)))
  recipload <- donorload + e
  
  return(recipload)
}


# Initialise donor~recipient pairs from a given distribution of donor logspvls and 
# regression coefficient corresponding to the heritability of viral load
InitPop <- function(popsize, donor_mean, donor_sd, herit = 0.33){
  init_donor <- rnorm(popsize,  mean = donor_mean, sd = donor_sd)
  init_recip <- -1
  while(min(init_recip)<0){
    init_recip <- GetRecipVL(init_donor, herit) 
  }
  out <- list(donor = init_donor, recip = init_recip)
  return(out)
}


# Infer weighting for a given viral load 
# Weighting 'g' from Thompson et al 2019
WeightPDF <- function(viralload){
  alpha = -3.55
  sigma <- 0.78/(sqrt(1 - ((2*alpha^2)/(pi*(1 + alpha^2)))))
  mu <-  4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha^2)))
  weight <-  (2/sigma)*dnorm((log10(viralload) - mu)/sigma)*pnorm(alpha*(log10(viralload) - mu)/sigma)
  return(weight)
} 


# Initialise simulation donor population from a specified min and max log spvl, generalised with
# a sampling distribution from WeightPDF (function 'g' from Thompson 2019 et al)
InitSimDonor <- function(pop_size, donor_min, donor_max, sample_prob){
  sim_range <- seq(donor_min, donor_max, length.out = pop_size)
  sim_logspvl <- sample(sim_range, size = pop_size, prob = sample_prob, replace = T)
  return(sim_logspvl)
  
}


# Infer simulation recipient population using donor~recipient regresssion model eqn from 
# initial donor~recipient pairs. takes coefficients, intercept and error term (sampled from
# normal dist with mean of residual standard error of model) for y = a + bx + e eqation.
InitSimRecip <- function(model, donor){
  out <- -1
  while (min(out) < 0) {
    if(class(model)=='lm'){
      out <- model$coefficients[1] + model$coefficients[2] * donor + rnorm(summary(model)$sigma)
    }else{
      stop("Function requires linear model as input.")
    }
  }
  
  return(out)
}


###################################################################################################
# Set seed
set.seed(4472)

# Data generation
# Normal distibution of log spvl from Rakkai Cohort, Hollingsworth et al. mean = 4.39, sd = 0.84
# https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1000876

NPAIRS <- 1000 # lower for test runs on low cpu machines
R2 <- 0.33 # Heritability estimate
RAKKAI_MEAN <- 4.39
RAKKAI_SD <- 0.84

model_population <- InitPop(NPAIRS, 
                            donor_mean = RAKKAI_MEAN, 
                            donor_sd = RAKKAI_SD, 
                            herit = R2)

donor_spvl <- 10^model_population$donor
recip_spvl <- 10^model_population$recip

# Data frame of paired viral loads
spvl_df <- cbind.data.frame(donor_spvl = donor_spvl, 
                            recip_spvl = recip_spvl) 

h2_model <- lm(recip ~ donor, data = model_population)
###################################################################################################
# Probability that recipient infection is initiated by multiple founder variants
# applies populationmodel_fixedVL_Environment function written by Katie Atkins
prob_recip_multiple <- RunParallel(populationmodel_fixedVL_Environment, donor_spvl) %>%
  do.call(cbind.data.frame, .) %>% t() # Time difference of 44.60442 mins on Macbook

combined_data <- cbind.data.frame(spvl_df, prob_recip_multiple)
head(combined_data)


###################################################################################################
# Probability that recipient infection is multiple founder, given a certain viral load

# Generate weightings using function g from Thompson et al
# Function generates a probability distribution (lognormal) which is used to sample from our range
# of simulated donor viral loads to generalise over a population
donor_prob <- sapply(log10(donor_spvl), function(x) WeightPDF(x)/sum(WeightPDF(log10(donor_spvl))))
hist(donor_prob) #visual check - should look normalish as already log

sim_donor_logspvl <- InitSimDonor(pop_size = NPAIRS,
                                  donor_min = 0.5, 
                                  donor_max = 8, 
                                  sample_prob = donor_prob)
hist(sim_donor_logspvl)

# Infer recipient viral loads of from simulated population
# predict.lm(h2_model, cbind.data.frame(donor_logspvl = sim_donor_logspvl)) (identical to donor)
# GetRecipVL(sim_donor_logspvl, R2) (too varied) suspect problem is attempting to pace a normal dist over uniform?
# Implementation below uses coefficients from recipient ~ donor model and normal dist error term
sim_recip_logspvl <- InitSimRecip(h2_model, sim_donor_logspvl) 

sim_donor_spvl <- 10^sim_donor_logspvl
sim_recip_spvl <- 10^sim_recip_logspvl

# Calculate probability of mulitple founder infection in recipient
sim_prob_multiple <- RunParallel(populationmodel_fixedVL_Environment, sim_donor_spvl)  %>%
  do.call(cbind.data.frame, .) %>% t()

sim_combined_data <- cbind.data.frame(sim_donor_spvl, sim_recip_spvl, sim_prob_multiple)
head(sim_combined_data)


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
fig_1b <- ggplot(combined_data, 
                 aes(x = donor_spvl, 
                     y = 1 - variant_distribution.V1))+
  geom_point()+
  scale_x_log10(name = 'Donor SPVL (log10)',
                limits = c(1, 10^10),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0,0),
                     limits = c(0,0.5),
                     breaks = seq(0, 0.6, by = 0.2)) +
  theme_bw() 

# Fig 1c
fig_1c <- ggplot(sim_combined_data, 
                 aes(x = sim_recip_spvl,
                     y = 1 - variant_distribution.V1))+
  geom_point()+
  scale_x_log10(name = 'Recipient SPVL (log10)',
                limits = c(1, 10^10),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0,0),
                     limits = c(0,0.75),
                     breaks = seq(0, 0.75, by = 0.25)) +
  theme_bw() 

# Panel 1
panel1 <- plot_grid(fig_1a, fig_1b, fig_1c, labels = 'AUTO', align = 'hv', ncol = 3)

panel1 


###################################################################################################
# Visualise probability distributions of the number of variants that initiate recipient infection ~
# recipient set point viral load 
# AKA panel 2

#store data long form & categorise recipient viral loads
sim_combined_data_long <- sim_combined_data %>%
  cbind(., ref = 1:nrow(sim_combined_data)) %>%
  gather(key = "variant_no",value = 'variant_prob',variant_distribution.V1:variant_distribution.V33) %>%
  mutate(spvl_cat = cut(sim_recip_spvl, breaks = 10^(0:9), labels = sapply(scientific(10^(1:9)), paste0, ' copies/ML')))

# replace categories of number of variants with integers
sim_combined_data_long$variant_no <- rep(1:33, each = 1000)

panel_2 <- ggplot(sim_combined_data_long, 
                  aes(x = variant_no, 
                      y =variant_prob/100)) + 
  geom_bar(stat = 'identity') + 
  scale_y_continuous(name = 'Probability', 
                     expand = c(0,0),
                     limits = c(0,1),
                     breaks = seq(0,1,0.2))+
  scale_x_continuous(name = 'Number of Variants', 
                     expand = c(0,0.6),
                     #limits = c(1,12),
                     breaks = seq(0,12,1))+
  facet_wrap(~spvl_cat) +
  theme_bw() +
  theme(panel.spacing = unit(1, "lines")) +
  coord_cartesian(xlim = c(1,8))


###################################################################################################
#                                                                                                 #
#                                               END                                               #
#                                                                                                 #
###################################################################################################