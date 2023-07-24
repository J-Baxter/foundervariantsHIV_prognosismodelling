# for fixed spVL, per virion prob of infection and proportion of exposure with potential for infection
# this function calculates i) the probability of infection, ii) the probability of multiple founder variants

TransmissionModel2 <- function(sp_ViralLoad = 10^6, PerVirionProbability = 4.715*1e-8, PropExposuresInfective = 0.029){
        require(doFuture)
        require(dplyr)
        require(fitdistrplus)
        
        
        
        # The VLs during primary stage
        
        # Fiebig et al 2003 profiles. Stages I-V only.
        # fiebig_vl <- data.table::fread('reconciliation/data/fiebig_vl.csv')
        # fiebig_vl <- round(10^fiebig_vl[1:5,c(6,7,2,8,9)])
        # npVals <- as.vector(as.matrix(fiebig_vl))
        
        npVals <- c(97, 42084, 64344, 3878, 1714, 421, 92159, 136380, 27524, 7214, 2742, 262131, 491492,
                    88507, 25176, 22626, 967867, 1216484, 427765, 105021, 42084, 3574304, 4492440, 1341713, 471801)
        
        # The SPVLs of transmitters
        nSims <- 100
        
        # The different VLs during preAIDS (Mellors et al 1995)
        naVals <- round(c(10^4, 10^4.5, 10^5, 10^5.5, 10^6))
        
        # fraction of exposure that can lead to infection
        # rename for reducing code length
        f <- PropExposuresInfective 
        sp_ViralLoad <- round(sp_ViralLoad)
        ##################################################
        # Set weights for each VL
        # For SPVL parameters are from Thompson et al 2019
        ##################################################
        alpha = -3.55
        sigma = 0.78/(sqrt(1 - ((2*alpha^2)/(pi*(1 + alpha^2)))))
        mu = 4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha^2)))
        kappa = 1
        np = round(kappa*max(npVals))
        nc = round(kappa*max(sp_ViralLoad))
        na = round(kappa*max(naVals))
        Dmax = 25.4
        Dfifty = 3058
        Dk = 0.41
        
        #Duration of each stage
        fiebig_t <- c(5,5.3,3.2,5.6,69.5)/365 #Stage duration (Fiebig et al 2003)
        
        tauc = array(0, length(sp_ViralLoad))
        taup = sum(fiebig_t)
        taua = 0.75
        
        # The weight of each npVals is proportional to stage duration and 
        # the corresponding quantile
        
        fiebig_q <- c(0.1,0.15,0.5,0.15,0.1) #Quantiles 
        gp <- as.vector(outer(fiebig_t/sum(fiebig_t), fiebig_q))
        
        # For sp_ViralLoad is, SPVL informs i) duration of chronic infection and 
        # ii) weights for the viral loads
        tauc = Dmax*(Dfifty^Dk)/(sp_ViralLoad^Dk + (Dfifty^Dk))
        gc = (2/sigma) * dnorm((log10(sp_ViralLoad) - mu)/sigma) * pnorm(alpha*(log10(sp_ViralLoad) - mu)/sigma)
        
        maximumTime = max(taup + tauc + taua)
        
        # Now scale the distribution of donors with each SPVL, g(V), by the infectious periods
        # of the donors, since the expressions in Fraser et al. (2007) represent numbers of individuals of
        # each viral load at seroconversion, rather than absolute numbers of individuals in the population
        
        # total time of infectiousness across VL types
        totalOfTime = length(gc)*(taup + taua) + sum(tauc)
        
        # normalise g
        gc = gc * ((taup + tauc + taua)/totalOfTime)
        gc = gc/sum(gc)
        
        # The weight of each naVals is given by the density function of
        # data from Mellors et al 1995 fitted to a lnorm
        
        # latefit <- fitdist(log(c(20190, 62530, 103900, 188000, 1235000),10), distr = "lnorm", method = "mle")
        # lower_bounds <- seq(3.75, 5.75, by=0.5) #based on naVals
        # upper_bounds <- seq(4.25, 6.25, by=0.5)
        # ranges <- lapply(1:length(lower_bounds), function(i) c(10^lower_bounds[i], 10^upper_bounds[i]))
        # ga <- sapply(ranges, function(x) diff(plnorm(log(x, 10), 
        #                                             meanlog = latefit$estimate[1], 
        #                                             sdlog = latefit$estimate[2])))
        # ga <- ga/sum(ga)
        
        ga <- c(0.06053491, 0.23398873, 0.34935005, 0.25220102, 0.10392529)
        
        ##################################################
        # Probability of acquisition
        ##################################################
        
        # FIRST CALCULATE THE number of max Virions Considered 
        chronic_prob_fulldist = dbinom(0:sp_ViralLoad,
                                       sp_ViralLoad,
                                       PerVirionProbability)
        VirionsConsidered <- which.max(chronic_prob_fulldist)
        maxVirionsConsidered <- min(which(chronic_prob_fulldist[VirionsConsidered:length(chronic_prob_fulldist)] < 1e-15))
        maxVirionsConsidered <- VirionsConsidered+maxVirionsConsidered-2
        
        # This reduces computational time. Variant transmission does not change 
        # much for higher virion transmission
        maxVirionsConsidered <- ifelse(maxVirionsConsidered>500, 500, maxVirionsConsidered)
        
        #  NOW CALCULATE FULL DISTRIBUTION OF VIRIONS NUMBER
        
        primary_prob_fulldist = purrr::map(.x = npVals,
                                           ~dbinom(0:maxVirionsConsidered, 
                                                   .x, 
                                                   PerVirionProbability))
        
        chronic_prob_fulldist = purrr::map(.x = sp_ViralLoad,
                                           ~dbinom(0:maxVirionsConsidered, 
                                                   .x, 
                                                   PerVirionProbability))
        
        
        preaids_prob_fulldist = purrr::map(.x = naVals,
                                           ~dbinom(0:maxVirionsConsidered, 
                                                   .x, 
                                                   PerVirionProbability))
        
        #if maxVirionsConsidered>500 edit last particle value so each x_prob_fulldist sums to 1
        if(any(round(sapply(primary_prob_fulldist, sum),9)!=1)){
                for(i in 1:length(primary_prob_fulldist)){
                        if(round(sum(primary_prob_fulldist[[i]]),9)!=1){
                                primary_prob_fulldist[[i]][length(primary_prob_fulldist[[i]])] <- primary_prob_fulldist[[i]][length(primary_prob_fulldist[[i]])] + 1-sum(primary_prob_fulldist[[i]])
                        }
                }
        }
        if(any(round(sapply(chronic_prob_fulldist, sum),9)!=1)){
                for(i in 1:length(chronic_prob_fulldist)){
                        if(round(sum(chronic_prob_fulldist[[i]]),9)!=1){
                                chronic_prob_fulldist[[i]][length(chronic_prob_fulldist[[i]])] <- chronic_prob_fulldist[[i]][length(chronic_prob_fulldist[[i]])] + 1-sum(chronic_prob_fulldist[[i]])
                        }
                }
        }
        if(any(round(sapply(preaids_prob_fulldist, sum),9)!=1)){
                for(i in 1:length(preaids_prob_fulldist)){
                        if(round(sum(preaids_prob_fulldist[[i]]),9)!=1){
                                preaids_prob_fulldist[[i]][length(preaids_prob_fulldist[[i]])] <- preaids_prob_fulldist[[i]][length(preaids_prob_fulldist[[i]])] + 1-sum(preaids_prob_fulldist[[i]])
                        }
                }
        }
        
        probNVirionsTransmittedPerSexAct <- array(0, maxVirionsConsidered+1)
        probNVirionsTransmittedPerSexAct_PRIMARY <- array(0, maxVirionsConsidered+1)
        probNVirionsTransmittedPerSexAct_CHRONIC <- array(0, maxVirionsConsidered+1)
        probNVirionsTransmittedPerSexAct_PREAIDS <- array(0, maxVirionsConsidered+1)
        
        #The principle here is:
        # probNVirionsTransmittedPerSexAct_XXX[1] = (1-f) + f*XXX_prob_fulldist[1]
        # probNVirionsTransmittedPerSexAct_XXX[-1] = f*XXX_prob_fulldist[-1]
        
        myfun1 <- function(g0, distr0){
                return(g0 * ((1-f) + f*distr0[1]))
        }
        
        myfun2 <- function(g0, distr0){
                return(g0 * f * distr0[-1])
        }
        
        myfun3 <- function(g0, tau0, distr0){
                return(g0 * ((1-f) + f*((taup/(taup + tau0 + taua))*sum(unlist(purrr::pmap(.l = list(as.list(gp), primary_prob_fulldist),
                                                                                           .f = myfun1))) +
                                                (tau0/(taup + tau0 + taua))*distr0[1] +
                                                (taua/(taup + tau0 + taua))*sum(unlist(purrr::pmap(.l = list(as.list(ga), preaids_prob_fulldist),
                                                                                                   .f = myfun1))))))
        }
        
        myfun4 <- function(g0, tau0, distr0){
                return(g0 * (f*((taup/(taup + tau0 + taua))*colSums(matrix(unlist(purrr::pmap(.l = list(as.list(gp), primary_prob_fulldist),
                                                                                              .f = myfun2)), 
                                                                           byrow= TRUE, nrow = length(primary_prob_fulldist))) +
                                        (tau0/(taup + tau0 + taua))*distr0[-1] +
                                        (taua/(taup + tau0 + taua))*colSums(matrix(unlist(purrr::pmap(.l = list(as.list(ga), preaids_prob_fulldist),
                                                                                                      .f = myfun2)), 
                                                                                   byrow= TRUE, nrow = length(preaids_prob_fulldist))))))
        }
        
        probNVirionsTransmittedPerSexAct[1] <- sum(unlist(purrr::pmap(.l = list(as.list(gc),as.list(tauc), chronic_prob_fulldist),
                                                                      .f = myfun3)))
        
        probNVirionsTransmittedPerSexAct[-1] <- colSums(matrix(unlist(purrr::pmap(.l = list(as.list(gc),as.list(tauc), chronic_prob_fulldist),
                                                                                  .f = myfun4)), 
                                                               byrow= TRUE, nrow = length(chronic_prob_fulldist)))
        
        
        probNVirionsTransmittedPerSexAct_PRIMARY[1] <- sum(unlist(purrr::pmap(.l = list(as.list(gp), primary_prob_fulldist),
                                                                              .f = myfun1)))
        
        probNVirionsTransmittedPerSexAct_PRIMARY[-1] <- colSums(matrix(unlist(purrr::pmap(.l = list(as.list(gp), primary_prob_fulldist),
                                                                                          .f = myfun2)), 
                                                                       byrow= TRUE, nrow = length(primary_prob_fulldist)))
        
        probNVirionsTransmittedPerSexAct_CHRONIC[1] <- sum(unlist(purrr::pmap(.l = list(as.list(gc), chronic_prob_fulldist),
                                                                              .f = myfun1)))
        probNVirionsTransmittedPerSexAct_CHRONIC[-1] <- colSums(matrix(unlist(purrr::pmap(.l = list(as.list(gc), chronic_prob_fulldist),
                                                                                          .f = myfun2)),
                                                                       byrow= TRUE, nrow = length(chronic_prob_fulldist)))
        
        probNVirionsTransmittedPerSexAct_PREAIDS[1] <- sum(unlist(purrr::pmap(.l = list(as.list(ga), preaids_prob_fulldist),
                                                                              .f = myfun1)))
        
        probNVirionsTransmittedPerSexAct_PREAIDS[-1] <- colSums(matrix(unlist(purrr::pmap(.l = list(as.list(ga), preaids_prob_fulldist),
                                                                                          .f = myfun2)), 
                                                                       byrow= TRUE, nrow = length(preaids_prob_fulldist)))
        
        probTransmissionPerSexAct = 1 - probNVirionsTransmittedPerSexAct[1]
        probTransmissionPerSexAct_PRIMARY = 1 - probNVirionsTransmittedPerSexAct_PRIMARY[1] 
        probTransmissionPerSexAct_CHRONIC = 1 - probNVirionsTransmittedPerSexAct_CHRONIC[1]
        probTransmissionPerSexAct_PREAIDS = 1 - probNVirionsTransmittedPerSexAct_PREAIDS[1]
        
        ##################################################
        # Probability of multiplicity
        ##################################################
        
        # calculate the ID of variants according to gamma distribution depending on time since infection 
        # these distributions are parameterised as hard-coded based on the haplotype distribution calculations in 
        # Thompson et al. Virus Evolution 2018
        
        # numberFounderStrainDistribution <- data.table::fread(paste0('reconciliation/data/tables/', maxVirionsConsidered, '.csv'))
        
        # calculate across all the times of transmission
        nTimeSteps = 100
        timeWindowEdges = seq(0, maximumTime, by = maximumTime/nTimeSteps)
        timeVals = timeWindowEdges[-length(timeWindowEdges)] + diff(timeWindowEdges)/2
        
        # calculate the ID of variants according to gamma distribution depending on time since infection 
        # these distributions are parameterised as hard-coded based on the haplotype distribution calculations in 
        # Thompson et al. Virus Evolution 2018
        
        
        # functions to sample IDs of variants and then calculate distributions of unique variants
        
        # conditioning on infection having occurred, calculate variant transmission for different virion transmission
        sample_strain_ID <- function(nparticles, maxVirionsConsidered, TimeSinceInfection){
                
                # if only one particle then only one variant 
                if (nparticles < 2){
                        numberFounderStrains <- 1;
                }else{
                        
                        HTprobDist_fn <- function(TimeSinceInfection, maxVirionsConsidered){
                                HTprobDist <- dgamma(x = 1:maxVirionsConsidered, 
                                                     shape = 0.417, 
                                                     scale = TimeSinceInfection/0.563)
                                HTprobDist <- HTprobDist/sum(HTprobDist)
                                return(HTprobDist)
                        }
                        
                        HTprobDist <- HTprobDist_fn(TimeSinceInfection, maxVirionsConsidered)
                        # for each of the nparticles transmitted, sample its ID using haplotype distribution
                        ids <- sample(x = 1:maxVirionsConsidered,
                                      size = nparticles,
                                      prob = HTprobDist,
                                      replace = TRUE)
                        
                        numberFounderStrains <- length(unique(ids))
                }
                return(numberFounderStrains)
        }
        
        distribution_founder_strains <- function(numberFounderStrains_samples, maxVirionsConsidered)
        {
                return(hist(numberFounderStrains_samples, 
                            breaks = 0:maxVirionsConsidered,
                            plot = FALSE)$density)
        }
        
        
        
        
        # calculate number of transmitted strains across nparticles - repeat nSims times
        parameters <- tidyr::expand_grid(nparticles = 1:maxVirionsConsidered, 
                                         tvals = timeVals,           
                                         sim = 1:nSims) # fastest variable is nsims
        
        nstrains  <- as.list(as.data.frame(matrix(unlist(purrr::map2(.x = parameters$nparticles,
                                                                     .y = parameters$tvals,
                                                                     .f = ~sample_strain_ID(.x,
                                                                                            maxVirionsConsidered = maxVirionsConsidered,
                                                                                            .y))),
                                                  ncol = maxVirionsConsidered * length(timeVals),
                                                  byrow = FALSE)))
        
        
        # calculate strain number distribution across increasing nparticles from simulations
        numberFounderStrainDistribution <- as.data.frame(matrix(unlist(purrr::map(.x = nstrains,
                                                                                  .f = ~distribution_founder_strains(.x, maxVirionsConsidered = maxVirionsConsidered))),
                                                                ncol = maxVirionsConsidered,
                                                                byrow = TRUE))
        parameters_nt <- expand_grid(nparticles = 1:maxVirionsConsidered, 
                                     tvals = timeVals) 
        numberFounderStrainDistribution <- cbind(parameters_nt, numberFounderStrainDistribution)
        
        
        ## FROM HERE IT NEEDS TO BE WEIGHTED BY VL density 
        
        # repeat founder strain distribution length(VL) times so we can use as input
        listofnumberFounderStrainDistribution = rep(list(numberFounderStrainDistribution), length(sp_ViralLoad))
        
        # calculate the number of variants for each particle number, weighted by the prob of particle number
        weight_numberFounderStrainDistribution <- purrr::map2(.x = listofnumberFounderStrainDistribution,
                                                              .y = as.list(tauc),
                                                              .f = ~(mutate(., prob_nparticles =
                                                                                    case_when(
                                                                                            tvals <= taup ~ probNVirionsTransmittedPerSexAct_PRIMARY[nparticles+1],
                                                                                            (tvals > taup) & (tvals <= (taup + .y)) ~ probNVirionsTransmittedPerSexAct_CHRONIC[nparticles+1],
                                                                                            tvals > (taup + .y) & (tvals <= (taup + .y + taua)) ~ probNVirionsTransmittedPerSexAct_PREAIDS[nparticles+1],
                                                                                            (tvals > (taup + .y + taua)) ~ 0)) %>%
                                                                             mutate(
                                                                                     across(starts_with("V"), ~(.* prob_nparticles)))))
        rm(listofnumberFounderStrainDistribution)
        
        
        variant_distribution <-  as.data.frame(weight_numberFounderStrainDistribution) %>%
                mutate(across(dplyr::starts_with("V"), ~(. * gc * (1/( taup + tauc + taua))))) %>%
                group_by(nparticles) %>%
                summarise(across(.cols = dplyr::starts_with(c("V" , 'prob')), .fns = sum)) %>%
                ungroup()
                
                
        #variant_distribution <- as.data.frame(matrix(
                #unlist(purrr::map(.x = weight_numberFounderStrainDistribution,
                #                  .f = ~(dplyr::select(., dplyr::starts_with("V")) %>%
                #                                 colSums(.)))),
               # ncol = maxVirionsConsidered,
                #byrow = TRUE))  %>%
                #mutate(across(everything(), ~(. * gc * (1/( taup + tauc + taua))))) %>%
                #group_by(nparticles) %>%
                #summarise(across(.cols = dplyr::starts_with(c("V" , 'prob')), .fns = sum)) %>%
                #ungroup()
        
        variant_distribution <- variant_distribution  %>% 
                mutate(across(.cols = dplyr::starts_with("V"), .fns = ~ .x / sum(variant_distribution %>% dplyr::select(dplyr::starts_with('V'))))) 
        
        #multiple_founder_proportion <- 1 - as.numeric(variant_distribution[1])
        #rm(variant_distribution)
        
        # sum up across all times in primary infection for each spvl person and then average across all of them - weighting by duration of infection and freq in population                                                       
        
        #variant_distribution_primary <- as.data.frame(matrix(
        # unlist(purrr::map2(.x = weight_numberFounderStrainDistribution,
        #    .y = as.list(tauc),
        #  .f = ~( filter(., tvals <= taup) %>%
        #                 dplyr::select(., starts_with("V")) %>%
        #                 colSums(.)))),
        # ncol = ncol(numberFounderStrainDistribution)-2,
        #byrow = TRUE)) %>%
        #mutate(across(everything(), ~(. * gc * (1/(taup + tauc + taua))))) %>%
        #colSums()
        
        # variant_distribution_primary <- variant_distribution_primary / sum(variant_distribution_primary)
        #multiple_founder_proportion_primary <- 1 - as.numeric(variant_distribution_primary[1])
        
        # sum up across all times in chronic infection for each spvl person and then average across all of them - weighting by duration of infection and freq in population                                                       
        # variant_distribution_chronic <- as.data.frame(matrix(
        # unlist(purrr::map2(.x = weight_numberFounderStrainDistribution,
        #                   .y = as.list(tauc),
        #                   .f = ~(filter(., tvals > taup & tvals <= (taup + .y)) %>%
        #                                  dplyr::select(., starts_with("V")) %>%
        #                                  colSums(.)))),
        # ncol = ncol(numberFounderStrainDistribution)-2,
        # byrow = TRUE)) %>%
        #mutate(across(everything(), ~(. * gc * (1/(taup + tauc + taua))))) %>%
        # colSums()
        # variant_distribution_chronic <- variant_distribution_chronic / sum(variant_distribution_chronic)
        #multiple_founder_proportion_chronic <- 1 - as.numeric(variant_distribution_chronic[1])
        
        # sum up across all times in preaids infection for each spvl person and then average across all of them - weighting by duration of infection and freq in population                                                       
        #variant_distribution_preaids <- as.data.frame(matrix(
        #   unlist(purrr::map2(.x = weight_numberFounderStrainDistribution,
        #                     .y = as.list(tauc),
        #                    .f = ~(filter(., tvals > (taup + .y) & (tvals <= (taup + .y + taua))) %>%
        #                                  dplyr::select(., starts_with("V")) %>%
        #                                 colSums(.)))),
        #ncol = ncol(numberFounderStrainDistribution)-2,
        #byrow = TRUE)) %>%
        # mutate(across(everything(), ~(. * gc * (1/(taup + tauc + taua))))) %>%
        #colSums()
        #variant_distribution_preaids <- variant_distribution_preaids / sum(variant_distribution_preaids)
        #multiple_founder_proportion_preaids <- 1 - as.numeric(variant_distribution_preaids[1])
        
        # output the prob of transmission across all infectious period and the chance of multiple lineages
        output <- list(variant_distribution,
                       probTransmissionPerSexAct,
                       sp_ViralLoad)
        #list(probTransmissionPerSexAct = probTransmissionPerSexAct,
        # probTransmissionPerSexAct_primary = probTransmissionPerSexAct_PRIMARY,
        # probTransmissionPerSexAct_chronic = probTransmissionPerSexAct_CHRONIC,
        #probTransmissionPerSexAct_preaids = probTransmissionPerSexAct_PREAIDS,
        #multiple_founder_proportion = multiple_founder_proportion,
        # multiple_founder_proportion_primary = multiple_founder_proportion_primary,
        # multiple_founder_proportion_chronic = multiple_founder_proportion_chronic,
        #multiple_founder_proportion_preaids = multiple_founder_proportion_preaids)
        
        return(output)
        
        
        
}
