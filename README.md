# RECONCILING FOUNDER VARIANT MULTIPLICITY OF HIV-1 INFECTION WITH THE RATE OF CD4+ DECLINE

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Results](#results)
- [License](./LICENSE)

<<<<<<< HEAD
# Overview
In the absence of effective antiretroviral therapy, almost all cases of HIV-1 infection will result in depletion of CD4+ T cells, such that the body can no longer formulate an adequate defence against opportunistic infection. However, the rate at which CD4+ T cells decline varies significantly across people living with HIV with unsuppressed viral replication, and the onset of AIDS can occur as little as two months to more than two decades post infection. While some studies report that infections that are initiated by multiple variants are associated with a faster CD4+ T cell decline, data are inconsistent and the mechanism through which this might occur is not understood. Here, we consider how variation in the number of variants initiating infection could explain variation in CD4+ T cell decline. Specifically, we use mathematical and statistical models that integrate clinical and epidemiological data to test the hypothesis that we should observe an association between greater genetic diversity at infection and a worse HIV prognosis.


Our multi-model framework encapsulates our present day characterisation of key relationships encompassing HIV transmission. We consider the statistical relationship between transmitter and recipient set point viral load (SpVL), the statistical relationship between recipient SpVL and rate of CD4+ T Cell decline; and the mechanistic relationship between transmitter SpVL and the probability that infection is initiated by multiple variants. We assume that the observed effect of multiple founder variants in the recipient on the rate of CD4 T Cell decline is mediated through SpVL.

# Repo Contents
- [scripts](./scripts): `R` code.
- [results](./results): raw results presented in the manuscript.
- [data](./data): example data.

# System Requirements

## Hardware Requirements

To run our analyses requires only a standard computer. Results from the 'transmission model'
will be vastly accelerated with access to multiple cores. On a 32GB RAM, 8 core iMac, the pipeline with the original data takes ~30 hours. On a AWS instance of 128GB RAM and 64 cores, the pipeline will run in 2 hours. 


## Software Requirements

### OS Requirements

This code was tested on *Linux* and *MacOSX* operating systems:

Linux: Ubuntu 22.04 LTS
Mac OSX:  

Users should have `R` version 4.1.0 or higher, and several packages set up from CRAN.


### Package dependencies
Users should install the following packages:

```
install.packages(c('ggplot2', 'abind', 'irlba', 'knitr', 'rmarkdown', 'latex2exp', 'MASS', 'randomForest'))
```

If you are having an issue that you believe to be tied to software versioning issues, please drop us an [Issue](https://github.com/neurodata/lol/issues). 

### Package Installation

From an `R` session, type:

```
require(devtools)
install_github('neurodata/lol', build_vignettes=TRUE, force=TRUE)  # install lol with the vignettes
require(lolR)
vignette("lol", package="lolR")  # view one of the basic vignettes
```

# Results

In this [benchmark comparison](http://docs.neurodata.io/lol/lol-paper/figures/real_data.html), we show that LOL does better than all linear embedding techniques in supervised HDLSS settings when dimensionality is high (d > 100, ntrain <= d) on 20 benchmark problems from the [UCI](https://archive.ics.uci.edu/ml/index.php) and [PMLB](https://github.com/EpistasisLab/penn-ml-benchmarks) datasets. LOL provides a good tradeoff between maintaining the class conditional difference (good misclassification rate) in a small number of dimensions (low number of embedding dimensions).
=======
Three quarters of new HIV-1 infections are initiated by a single genetic variant, and infections initiated by multiple variants are linked with higher recipient viral loads and a faster rate CD4+ T cell decline. However, these findings have not been universally replicated and a mechanism through which multiple variants might lead to higher viral loads and a faster rate CD4+ T cell decline is yet to be elucidated. Furthermore, the role of the transmitter’s infection in determining both the number of founder variants initiating infection and the virulence of the infecting virus, is yet to be investigated. In this study, we first summarised the existing evidence for this phenomenon. We then quantified how likely we are to observe a true difference in set point viral load in an observational study, between multiple and single variant infections. We found the probability of observing a significantly greater viral load varied between studies. Specifically, at frequencies of multiple variant infection consistent with empirical observations, we are unlikely to record a significant difference in SpVL between multiple and single variant infections.  Finally, we combined data-driven models of transmission, heritability and HIV-1 disease progression to test whether we expect to observe an association between multiple variant infection and clinical progression, in the absence of an explicitly incorporated causal effect. Our model predicts that we would not expect multiple variant infections to lead to higher SpVL or faster CD4+ T cell decline without a causal mechanism. Although the probability that infection is initiated in a recipient is greatest at the highest transmitter SpVLs, the relatively low contribution of the virus genotype to recipient SpVL precludes a major role in determining the recipient prognosis. This finding supports the hypothesis that a within-patient causal mechanism is required to explain the association of multiple variant infection with higher viral loads and faster CD4+ T cell decline. Further investigation into events happening during and just after transmission are required to enhance our understanding of this association. 

./scripts/main will run the analyses. 
There is the option to export figures to file at the end of this script. 
>>>>>>> 6a169223f1bf7e2cec48d0e0cb54628a6eac8272

