# RECONCILING FOUNDER VARIANT MULTIPLICITY OF HIV-1 INFECTION WITH THE RATE OF CD4+ DECLINE

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Results](#results)
- [License](./LICENSE)

# Overview
IHIV-1 transmission precipitates a stringent genetic bottleneck, with 75% of new infections initiated by a single genetic variant. Where multiple variants initiate infection, recipient set point viral load (SpVL) and the rate of CD4+ T cell decline may be elevated, but these findings remain inconsistent. Here, we summarised the evidence for this phenomenon, then tested whether previous studies possessed sufficient statistical power to reliably identify a true effect of multiple variant infection associating with higher SpVL. Next, we combined models of HIV-1 transmission, heritability and disease progression to understand whether available data suggest a faster CD4+ T cell decline would be expected to associate with multiple variant infection, without an explicit dependency between the two. We found that most studies had insufficient power to identify a true significant difference, prompting an explanation for previous inconsistencies. Next, our model framework revealed we would not expect to observe a positive association between multiple variant infections and faster CD4+ T cell decline, in the absence of an explicit dependency. Consequently, while empirical evidence may be consistent with a positive association between multiple variant infection and faster CD4+ T cell decline, further investigation is required to establish a causal basis for this association

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
Executing the 'main.R' script will run the analyses to test the hypothesis that we should observe an association between greater genetic diversity at infection and a worse HIV prognosis, in the absence of an
explicit dependency between the two. The component models are sourced from 'HeritabilityModel.R', 'TransmissionModel.R' and 'ToleranceModel.R'.

The script 'observing_sigeffect' calculates how likely we are to observe a significant difference between the SpVLs of single and multiple variant infections, given a true effect size (ie. a power test of previously
conducted observational studies)


# Results
We found that most studies had insufficient power to identify a true significant difference, prompting an explanation for previous inconsistencies. Next, our model framework revealed we would not expect to observe a positive association between multiple variant infections and faster CD4+ T cell decline, in the absence of an explicit dependency. Specifically, the probability that infection is initiated by multiple variants is greatest at the highest transmitter SpVLs, yet the relationship between transmitter and recipient SpVL is relatively weak. Consequently, while empirical evidence may be consistent with a positive association between multiple variant infection and faster CD4+ T cell decline, further investigation is required to establish a causal basis for this association
