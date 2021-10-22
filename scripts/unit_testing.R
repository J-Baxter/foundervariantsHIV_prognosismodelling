###################################################################################################
###################################################################################################
# Test script accompanying probabilistic models of the association between the number of variants
# and viral load
# Mostly unit testing using testthat
###################################################################################################
###################################################################################################
# Dependencies
library(testthat)
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


###################################################################################################
# Set seed
set.seed(4472)

# GetRecipVL

# InitPop

# WeightPDF

# InitSimDonor

# Recipient VL check

