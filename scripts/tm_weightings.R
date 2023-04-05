
PerVirionProbability = 4.715*1e-8

#SPVL = virions available for transmission
# fixed ppp 
# higher spvl = lower tau_c (due to increased rate of decline)
# 


sp_ViralLoad = 10**6
PerVirionProbability = 4.715*1e-8
f= 0.029
alpha = -3.55
sigma = 0.78/(sqrt(1 - ((2*alpha^2)/(pi*(1 + alpha^2)))))
mu = 4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha^2)))
Vp = 87173000
Va = 24004000
kappa = 1
np = round(kappa*Vp)
na = round(kappa*Va)
sp_ViralLoad = round(kappa*sp_ViralLoad)
Dmax = 25.4
Dfifty = 3058
Dk = 0.41
tauc = array(0, length(sp_ViralLoad))
taup = 0.24
taua = 0.75
taucrit = 2  
  
# spVL specific quantities: i) weights for the viral loads, ii) duration of chronic infection

tauc = Dmax*(Dfifty^Dk)/(sp_ViralLoad^Dk + (Dfifty^Dk))
g = (2/sigma) * dnorm((log10(sp_ViralLoad) - mu)/sigma) * pnorm(alpha*(log10(sp_ViralLoad) - mu)/sigma)

# Now scale the distribution of donors with each SPVL, g(V), by the infectious periods
# of the donors, since the expressions in Fraser et al. (2007) represent numbers of individuals of
# each viral load at seroconversion, rather than absolute numbers of individuals in the population
numIndividualTypes = length(g)

# total time of infectiousness across VL types
totalOfTimeVals = numIndividualTypes*(taup + taua) + sum(tauc)

# normalise g
g = g * ((taup + tauc + taua)/totalOfTimeVals)
g = g/sum(g)

maximumTime = max(taup + tauc + taua);

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

# calculate distribution for each spVL
chronic_prob_fulldist = purrr::map(.x = sp_ViralLoad,
                                   ~dbinom(0:maxVirionsConsidered, 
                                           .x, 
                                           PerVirionProbability))

preaids_prob_fulldist = dbinom(0:maxVirionsConsidered, 
                               na, 
                               PerVirionProbability)

probNVirionsTransmittedPerSexAct <- array(0, maxVirionsConsidered+1)
probNVirionsTransmittedPerSexAct_PRIMARY <- array(0, maxVirionsConsidered+1)
probNVirionsTransmittedPerSexAct_CHRONIC <- array(0, maxVirionsConsidered+1)
probNVirionsTransmittedPerSexAct_PREAIDS <- array(0, maxVirionsConsidered+1)

p_0_primary = g * ((1-f) + f*((taup/(taup + tauc + taua))*primary_prob_fulldist[1])) 
p_0_chronic = g * ((1-f) + f*((tauc/(taup + tauc + taua))*chronic_prob_fulldist[[1]][1]))
p_0_preaids = g * ((1-f) + f*((taua/(taup + tauc + taua))*preaids_prob_fulldist[1]))


p_0 = g * ((1-f) + f*((taup/(taup + tauc + taua))*primary_prob_fulldist[1] + (tauc/(taup + tauc + taua))*chronic_prob_fulldist[[1]][1] + (taua/(taup + tauc + taua))*preaids_prob_fulldist[1]))
c(taup, tauc, taua)
c(1-p_0_primary, 1-p_0_chronic, 1-p_0_preaids)
1 -(p_0_primary + p_0_chronic + p_0_preaids - 2*g*(1-f))
df <- cbind.data.frame(tau = c('p', 'c', 'a'), value= c(taup, tauc, taua), vl= c(87173000, 1e+06, 2400400), p = c(1-p_0_primary, 1-p_0_chronic, 1-p_0_preaids))


df$right <- cumsum(df$value) 
df$left <- c(0,  cumsum(df$value)[-3])


plot_acq <- ggplot(df) + 
  geom_rect(aes(xmin = left, xmax = right, ymax = p, ymin = 0, fill = tau))+
  #geom_hline(aes(yintercept = 1 -(p_0_primary + p_0_chronic + p_0_preaids - 2*(1-f)))) + 
  scale_x_continuous(expand = c(0,0), name = 'Time', limits = c(0,6)) + 
  scale_y_continuous(expand = c(0,0), name = 'Probability of Acquisition') + 
  scale_fill_brewer(palette = 'OrRd')+
  annotate(geom="text", x=0.12, y=0.004, label = expression(tau['p'] ), color="black", size = 10) + 
  annotate(geom="text", x=2.69, y=0.004, label = expression(tau['c'] ), color="black", size = 10) + 
  annotate(geom="text", x=5.5, y=0.004, label = expression(tau['a'] ), color="black", size = 10) + 
  my_theme+
  theme(legend.position = 'none')

var_t <- read_csv('./data/variant_distribution_time.csv') %>% mutate(x_var = as.factor(x_var)) %>% rename(time = T)

plot_var <- ggplot(var_t, aes(x = time, y = P, fill = x_var)) + 
  geom_area() + 
  scale_x_continuous(name = 'Time', expand = c(0,0))+
  scale_y_continuous(name = 'Variant Proportion', expand = c(0,0))+
  scale_fill_brewer(palette = 'OrRd', na.value = 'black')+
  my_theme+
  theme(legend.position = 'none')

tm_grid <- cowplot::plot_grid( plot_var,plot_acq, nrow = 1, align = 'hv', labels = 'AUTO')

ggsave("acquisition_time.eps", device =  cairo_ps , plot =tm_grid , width = 21, height = 12)


# calculate distribution for each spVL


myfun1 <- function(g0, tau0, distr0){
  return(g0 * ((1-f) + f*((taup/(taup + tau0 + taua))*primary_prob_fulldist[1] +
                            (tau0/(taup + tau0 + taua))*distr0[1] +
                            (taua/(taup + tau0 + taua))*preaids_prob_fulldist[1])))
}

myfun2 <- function(g0, tau0, distr0){
  return(g0 * (f*((taup/(taup + tau0 + taua))*primary_prob_fulldist[-1] +
                    (tau0/(taup + tau0 + taua))*distr0[-1] +
                    (taua/(taup + tau0 + taua))*preaids_prob_fulldist[-1])))
}

myfun3 <- function(g0, distr0){
  return(g0 * ((1-f) + f*distr0[1]))
}

myfun4 <- function(g0, distr0){
  return(g0 * f * distr0[-1])
}

probNVirionsTransmittedPerSexAct[1] <- sum(unlist(purrr::pmap(.l = list(as.list(g),as.list(tauc), chronic_prob_fulldist),
                                                              .f = myfun1)))
probNVirionsTransmittedPerSexAct[-1] <- colSums(matrix(unlist(purrr::pmap(.l = list(as.list(g),as.list(tauc), chronic_prob_fulldist),
                                                                          .f = myfun2)), byrow= TRUE, nrow = length(chronic_prob_fulldist)))


probNVirionsTransmittedPerSexAct_CHRONIC[1] <- sum(unlist(purrr::pmap(.l = list(as.list(g), chronic_prob_fulldist),
                                                                      .f = myfun3)))   
probNVirionsTransmittedPerSexAct_CHRONIC[-1] <- colSums(matrix(unlist(purrr::pmap(.l = list(as.list(g), chronic_prob_fulldist),
                                                                                  .f = myfun4)), byrow= TRUE, nrow = length(chronic_prob_fulldist)))


probNVirionsTransmittedPerSexAct_PRIMARY[1] <- sum(unlist(purrr::map(.x = g,
                                                                     .f = ~(.x*((1-f) + f*primary_prob_fulldist[1])) )))  
probNVirionsTransmittedPerSexAct_PRIMARY[-1] <- colSums(matrix(unlist(purrr::map(.x = g,
                                                                                 .f = ~(.x*f*primary_prob_fulldist[-1]))), byrow= TRUE, nrow = length(chronic_prob_fulldist)))


probNVirionsTransmittedPerSexAct_PREAIDS[1] = sum(unlist(purrr::map(.x = g,
                                                                    .f = ~(.x * ((1-f) + f*preaids_prob_fulldist[1])))))
probNVirionsTransmittedPerSexAct_PREAIDS[-1] <- colSums(matrix(unlist(purrr::map(.x = g,
                                                                                 .f = ~(.x*f*preaids_prob_fulldist[-1]))), byrow= TRUE, nrow = length(chronic_prob_fulldist)))


probTransmissionPerSexAct = 1 - probNVirionsTransmittedPerSexAct[1]



###################################################################################################
# For a given t value, weighting is given by:
nSims = 100
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


# infectiousperiod <- taup + tauc + taua


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
                                                                 across(dplyr::starts_with("V"), ~(.* prob_nparticles)))))



#variant_distribution <- as.data.frame(matrix(
# unlist(purrr::map(.x = weight_numberFounderStrainDistribution,
#                  .f = ~(dplyr::select(., starts_with("V")) %>%
#                            colSums(.)))),
# ncol = maxVirionsConsidered,
# byrow = TRUE))  %>%
# mutate(across(everything(), ~(. * g * (1/( taup + tauc + taua))))) %>%
# colSums()
#variant_distribution <- variant_distribution / sum(variant_distribution)
# multiple_founder_proportion <- 1 - as.numeric(variant_distribution[1])

w = 1
variant_distribution_unweighted <- as.data.frame(weight_numberFounderStrainDistribution) %>%
  mutate(case_when(tvals <= taucrit & taucrit < (taup + tauc + taua) ~ (dplyr::select(., starts_with("V")) * g * (w/((w-1)*taucrit + taup + tauc + taua))), 
                   tvals > taucrit ~ (dplyr::select(., starts_with("V")) * g * (1/( (w-1)*taucrit +taup + tauc + taua))))) %>%
  group_by(nparticles) %>%
  summarise(across(.cols = starts_with(c("V" , 'prob')), .fns = sum)) %>%
  ungroup()%>%
  summarise(across(.cols = starts_with('V'), .fns = sum))

variant_distribution_unweighted  <- variant_distribution_unweighted   %>% 
  mutate(across(.cols = starts_with("V"), .fns = ~ .x / sum(variant_distribution_unweighted  %>% select(starts_with('V'))))) %>%
  mutate(w = 1)

w = 10
variant_distribution_weighted<- as.data.frame(weight_numberFounderStrainDistribution) %>%
  mutate(case_when(tvals <= taucrit & taucrit < (taup + tauc + taua) ~ (dplyr::select(., starts_with("V")) * g * (w/((w-1)*taucrit + taup + tauc + taua))), 
                   tvals > taucrit ~ (dplyr::select(., starts_with("V")) * g * (1/( (w-1)*taucrit +taup + tauc + taua))))) %>%
  group_by(nparticles) %>%
  summarise(across(.cols = starts_with(c("V" , 'prob')), .fns = sum)) %>%
  ungroup() %>%
  summarise(across(.cols = starts_with('V'), .fns = sum))

variant_distribution_weighted  <- variant_distribution_weighted   %>% 
  mutate(across(.cols = starts_with("V"), .fns = ~ .x / sum(variant_distribution_weighted  %>% select(starts_with('V'))))) %>%
  mutate(w = 10)

compare_weights = rbind.data.frame(variant_distribution_unweighted, variant_distribution_weighted) %>%
  pivot_longer(cols = starts_with('V'), 
               names_to = 'variants',
               values_to = 'p') %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric()) 

ggplot(compare_weights) + 
  geom_bar(aes(x =variants, y= p, fill = w), stat = 'identity', position = "dodge2")


