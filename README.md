# RECONCILING FOUNDER VARIANT MULTIPLICITY OF HIV-1 INFECTION WITH THE RATE OF CD4+ DECLINE

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Results](#results)
- [License](./LICENSE)

# Overview
In the absence of effective antiretroviral therapy, almost all cases of HIV-1 infection will result in depletion of CD4+ T cells, such that the body can no longer formulate an adequate defence against opportunistic infection. However, the rate at which CD4+ T cells decline varies significantly across people living with HIV with unsuppressed viral replication, and the onset of AIDS can occur as little as two months to more than two decades post infection. While some studies report that infections that are initiated by multiple variants are associated with a faster CD4+ T cell decline, data are inconsistent and the mechanism through which this might occur is not understood. Here, we consider how variation in the number of variants initiating infection could explain variation in CD4+ T cell decline. Specifically, we use mathematical and statistical models that integrate clinical and epidemiological data to test the hypothesis that we should observe an association between greater genetic diversity at infection and a worse HIV prognosis.


Our multi-model framework encapsulates our present day characterisation of key relationships encompassing HIV transmission. We consider the statistical relationship between transmitter and recipient set point viral load (SpVL), the statistical relationship between recipient SpVL and rate of CD4+ T Cell decline; and the mechanistic relationship between transmitter SpVL and the probability that infection is initiated by multiple variants. We assume that the observed effect of multiple founder variants in the recipient on the rate of CD4 T Cell decline is mediated through SpVL.

# Repo Contents
- [scripts](./scripts): `R` code.
- [results](./results): raw results presented in the manuscript.
- [data](./data): joint posterior distributions for transmission model and paired SpVL data.

# System Requirements

## Hardware Requirements

To run our analyses requires only a standard computer. Results from the 'transmission model'
will be vastly accelerated with access to multiple cores. Using an iMac with 3.1 GHz 6-Core Intel i5 and 16 GB RAM, the pipeline with the original data takes ~30 hours. On a Linux server with 3.1 GHz Intel Xeon (32 cores) and 128GB RAM, the pipeline will run in ~2 hours. 


## Software Requirements

### OS Requirements

This code was tested on *Linux* and *MacOSX* operating systems:

Linux: Ubuntu 22.04 LTS
Mac OSX: Venutra 13.2.1

Users should have `R` version 4.1.0 or higher, and several packages set up from CRAN.


### Package dependencies
Users should install the following packages:

```
install.packages(c('tidyverse', 'parallel', 'brms', 'magrittr', 'performance', 'tmvtnorm','epitools', 'tidybayes', 'ggmcmc', 'ggpmisc', 'ggpubr', 'ggprism', 'cowplot', 'scales','extrafont', 'ggExtra', 'ggdag', 'RColorBrewer', 'stringi', 'matrixcalc','bayesplot', 'pdftools'))
```


# Description
Executing the 'main.R' script will run the analyses to test the hypothesis that we should observe an association between greater genetic diversity at infection and a worse HIV prognosis. The component models are sourced from 'HeritabilityModel.R', 'TransmissionModel.R' and 'ToleranceModel.R'.

The script 'observing_sigeffect' calculates how likely we are to observe a significant difference between the SpVLs of single and multiple variant infections, given a true effect size.


# Results
We found that we are unlikely to record a significant difference in SpVL between multiple and single variant infections, at frequencies of multiple variant infections consistent with empirical observations. Next, we found that we would not expect multiple variant infections to lead to higher SpVL or faster CD4+ T cell decline without a causal mechanism. Specifically, the probability that infection is initiated by multiple variants is greatest at the highest transmitter SpVLs, yet the relationship between transmitter and recipient SpVL is relatively weak. This finding supports the hypothesis that a within-patient causal mechanism is required to explain the association of multiple variant infection with higher viral loads and faster CD4+ T cell decline. Further investigation into events happening during and just after transmission are required to enhance our understanding of this association.
