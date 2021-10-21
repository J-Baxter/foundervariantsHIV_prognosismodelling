# for fixed spVL, per virion prob of infection and proportion of exposure with potential for infection
# this function calculates i) the probability of infection, ii) the probability of multiple founder variants
# JB edit 21/10/21 - returns probability distribution instead of P(multiple)
populationmodel_fixedVL_Environment <- function(sp_ViralLoad = 10^5, PerVirionProbability = 4.715*1e-8, PropExposuresInfective = 0.029){

  # this is to speed up the code, but can be updated based on phylogenetic results
  
  nSims <- 100
  
  # rename for reducing code length
  f <- PropExposuresInfective # fraction of exposure that can lead to infection
  
  ##################################################
  # parameters from the Thomson paper for ....
  ##################################################
  alpha = -3.55
  sigma = 0.78/(sqrt(1 - ((2*alpha^2)/(pi*(1 + alpha^2)))))
  mu = 4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha^2)))
  Vp = 87173000
  Va = 24004000
  kappa = 1
  np = round(kappa*Vp)
  na = round(kappa*Va)
  Dmax = 25.4
  Dfifty = 3058
  Dk = 0.41
  nc = round(sp_ViralLoad)
  # g = zeros(length(ncVals), 1)
  # tauc = zeros(length(ncVals), 1)
  taup = 0.24;
  taua = 0.75;
  # don't need to weight by the fraction in a particular viral load because only using specific VLs
  # g <- (2/sigma)*normpdf((log10(ViralLoad) - mu)/sigma)*normcdf(alpha*(log10(ViralLoad) - mu)/sigma)
  tauc <- Dmax*(Dfifty^Dk)/(nc^Dk + (Dfifty^Dk))
  maximumTime <- taup + taua + tauc
  # g = g./sum(g);

        # Now scale the distribution of donors with each SPVL, g(V), by the infectious periods
        # of the donors, since the expressions in Fraser et al. (2007) represent numbers of individuals of
        # each viral load at seroconversion, rather than absolute numbers of individuals in the population
        
        
        totalOfTimeVals = taup + tauc + taua
        
        #  NOW CALCULATE FULL DISTRIBUTION OF VIRIONS NUMBER
        m = np*PerVirionProbability
        s = np*PerVirionProbability*(1-PerVirionProbability)
        iter = 1;
        
        primary_prob_fulldist = dbinom(0:ceiling(m+iter*s),
                                       np,
                                       PerVirionProbability)
        
        while (primary_prob_fulldist[length(primary_prob_fulldist)] > 1e-15){
              iter = iter + 1
              primary_prob_fulldist = dbinom(0:ceiling((m+iter*s)),np,PerVirionProbability)
        }
        
        maxVirionsConsidered = ceiling((m+iter*s))
        
        chronic_prob_fulldist = dbinom(0:maxVirionsConsidered, 
                                       nc, 
                                       PerVirionProbability)
        
        preaids_prob_fulldist = dbinom(0:maxVirionsConsidered, 
                                       na, 
                                       PerVirionProbability)
        
        probNVirionsTransmittedPerSexAct <- array(0, maxVirionsConsidered+1)
        probNVirionsTransmittedPerSexAct_PRIMARY <- array(0, maxVirionsConsidered+1)
        probNVirionsTransmittedPerSexAct_CHRONIC <- array(0, maxVirionsConsidered+1)
        probNVirionsTransmittedPerSexAct_PREAIDS <- array(0, maxVirionsConsidered+1)
        
        probNVirionsTransmittedPerSexAct[1] = (1-f) + f*((taup/(taup + tauc + taua))*primary_prob_fulldist[1] +
                  (tauc/(taup + tauc + taua))*chronic_prob_fulldist[1] +
                  (taua/(taup + tauc + taua))*preaids_prob_fulldist[1])
        
        probNVirionsTransmittedPerSexAct[-1] = f*((taup/(taup + tauc + taua))*primary_prob_fulldist[-1] +
                 (tauc/(taup + tauc + taua))*chronic_prob_fulldist[-1] +
                 (taua/(taup + tauc + taua))*preaids_prob_fulldist[-1])
        
        probNVirionsTransmittedPerSexAct_PRIMARY[1] = (1-f) + f*primary_prob_fulldist[1]
        
        probNVirionsTransmittedPerSexAct_PRIMARY[-1] = f*primary_prob_fulldist[-1]
        
        probNVirionsTransmittedPerSexAct_CHRONIC[1] = (1-f) + f*chronic_prob_fulldist[1]
        
        probNVirionsTransmittedPerSexAct_CHRONIC[-1] = f*chronic_prob_fulldist[-1]
        
        probNVirionsTransmittedPerSexAct_PREAIDS[1] = (1-f) + f*preaids_prob_fulldist[1]
        
        probNVirionsTransmittedPerSexAct_PREAIDS[-1] = f*preaids_prob_fulldist[-1]
        
        # end
        
        probTransmissionPerSexAct = 1 - probNVirionsTransmittedPerSexAct[1]
        
      
        
        # calculate across all the times of transmission
        nTimeSteps = 100;
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
        parameters <- expand_grid(nparticles = 1:maxVirionsConsidered, 
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
                                  tvals = timeVals) # fastest variable is nsims
        
        numberFounderStrainDistribution <- cbind(parameters_nt, numberFounderStrainDistribution)

        
        
        
        
        weight_numberFounderStrainDistribution <- numberFounderStrainDistribution %>% 
                    # add the prob of n particles depending on the duration of infection
                    mutate(prob_nparticles = 
                             case_when(
                               tvals < taua ~ probNVirionsTransmittedPerSexAct_PRIMARY[nparticles+1],
                               (tvals > taua) & (tvals <= (taua+tauc)) ~ probNVirionsTransmittedPerSexAct_CHRONIC[nparticles+1],
                               tvals > (taua+tauc) & (tvals <= (taua+tauc+taup)) ~ probNVirionsTransmittedPerSexAct_PREAIDS[nparticles+1],
                               tvals > (taua+tauc+taup) ~ 0
                               )) %>%      
                    # weight the variant distirbution by the chance of the number of particles
                      mutate(across(starts_with("V"), ~(.x * prob_nparticles))) 
                      
          variant_distribution <- weight_numberFounderStrainDistribution %>%
                select(starts_with("V")) %>%
                colSums(.)
          variant_distribution <- variant_distribution / sum(variant_distribution)
          multiple_founder_proportion <- 1 - as.numeric(variant_distribution[1])
        
          
          # output the prob of transmission across all infectious period and the chance of multiple lineages
          output <-  c(probTransmissionPerSexAct = probTransmissionPerSexAct, variant_distribution = variant_distribution)
                   
          return(output)
        

}