###################################################################################################
###################################################################################################
# Modelling the association between the number of variants and viral load
# Script utilises methodology from Thompson et al 2019, https://github.com/robin-thompson/MultiplicityOfInfection
# implemented in R by Katie Atkins https://github.com/katiito/HIV_eclipse_models
# Testing hypothesis 1: Stage of infection in donor

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
source("./scripts/null_models.R")
load("./nulldata.RData")

