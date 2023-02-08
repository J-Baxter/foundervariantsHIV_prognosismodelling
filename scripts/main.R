# Main Analysis for "RECONCILING THE FOUNDER VARIANT MULTIPLICITY OF HIV-1 INFECTION WITH THE RATE
# OF CD4+ DECLINE"
# J Baxter

# This script combines the functions for the analysis of paired SPVL distributions, occuring in
# three main stages:

# 1. Initialisation of study population
# 2. Null hypothesis
# 3. Ecological hypotheses


# Dependencies
source('./scripts/dependencies.R')
source('./scripts/populationdata_acrossVL_models.R')
source('./scripts/init_population.R')
source('./scripts/ecologicalhypos.R')

# Set Seed
set.seed(4472)
options(scipen = 100) #options(scipen = 100, digits = 4)


###################################################################################################
# Initialise Population
#Distributions of SPVL for couples with strong support for transmission Hollingsworth et al, 2010

index <- c('mean'= 4.61, 'sd' = 0.63) 
secondary <- c('mean' = 4.60 , 'sd' = 0.85)
NPAIRS <- 500 #500

pop <- InitPop(N = NPAIRS, 
               H2 = 0.67, #0.33
               donor_vl = index, 
               recipient_vl = secondary)


# Fit mixed model to simulated data to re-estimate/confirm heritability estimate
h2_model <- lm(recipient ~ transmitter, data = select(pop, c('transmitter', 'recipient')))

###################################################################################################
# Probability that recipient infection is initiated by multiple founder variants
# applies populationmodel_fixedVL_Environment function written by Katie Atkins

prob_recip_multiple <- RunParallel(populationmodel_acrossVL_Environment, pop$transmitter) %>%
  do.call(cbind.data.frame, .) %>% t()

transmitterspvl_variantdist <- cbind.data.frame(pop$transmitter, prob_recip_multiple)


###################################################################################################
# Probability that recipient infection is multiple founder, given a certain viral load
# Generate weightings using function g from Thompson et al
# This function generates a probability distribution (lognormal) which is used to sample from our range
# of simulated transmitter viral loads to generalise over a population

WeightPDF <- function(viralload){
  alpha = -3.55
  sigma <- 0.78/(sqrt(1 - ((2*alpha**2)/(pi*(1 + alpha**2)))))
  mu <-  4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha**2)))
  weight <-  (2/sigma)*dnorm((log10(viralload) - mu)/sigma)*pnorm(alpha*(log10(viralload) - mu)/sigma)
  return(weight)
} 


transmitter_prob <- sapply(pop$transmitter_log10SPVL, function(x) WeightPDF(x)/sum(WeightPDF(pop$transmitter_log10SPVL)))
hist(transmitter_prob) #visual check - should look normal(ish) as already log


# Initialise simulation donor population from a specified min and max log spvl
InitSimTransmitter <- function(pop_size, transmitter_min, transmitter_max, sample_prob){
  sim_range <- seq(transmitter_min, transmitter_max, length.out = pop_size)
  sim_logspvl <- sample(sim_range, size = pop_size, prob = sample_prob, replace = T) %>%
    cbind.data.frame() %>%
    `colnames<-`('transmitter')
  
  return(sim_logspvl)
}


sim_transmitter_log10SPVL <- InitSimTransmitter(pop_size = NPAIRS,
                                              transmitter_min = 1, 
                                              transmitter_max = 7, 
                                              sample_prob = transmitter_prob)
hist(sim_transmitter_log10SPVL$transmitter) #visual check - should look uniform


# Leverage heritability model to estimate recipient SPVL
sim_spvl <- predict(h2_model, newdata=sim_transmitter_log10SPVL) %>% 
  cbind.data.frame(recipientlog10SPVL = ., donor = sim_transmitter_log10SPVL) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SPVL')}"))


# Calculate probability of mulitple founder infection in recipient
sim_prob_multiple <- RunParallel(populationmodel_acrossVL_Environment, sim_spvl$transmitter)  %>%
  do.call(cbind.data.frame, .) %>% t()

sim_recipientspvl_variantdist <- cbind.data.frame(sim_spvl$recipient, sim_prob_multiple)


###################################################################################################
# Write to file
# population only
write_csv(pop, file = paste(results_dir, 'base_population.csv', sep = '/'))

# transmitter SPVL & variant distribution
write_csv(transmitterspvl_variantdist, file = paste(results_dir, 'transmitterspvl_variantdist', sep = '/'))

# recipient SPVL & variant distribution
write_csv(sim_recipientspvl_variantdist, file = paste(results_dir, 'recipientspvl_variantdist.csv', sep = '/'))


###################################################################################################
# Alternative

# Generate a matrix of transmitter/founding particles and associated viral loads
FoundingVars <- function(spvl, sd, p_dists){
  
  #Round Log10 SPVL to nearest .5
  spvl_disc <- round(spvl/0.5)*0.5
  
  #Identify appropriate probability distribution
  lookup <- which(spvl_disc == seq(2,7,by=0.5))
  prob_var <- p_dists[,lookup]
  
  n_var <- sample(x = 1:33, 1, replace = T, prob = prob_var)
  
  within_host_genotypes <- rnorm(n = n_var, mean = spvl, sd = sd)
  
  m <- matrix(NA, nrow = 1, ncol = 33)
  
  m[1:length(within_host_genotypes)] <- within_host_genotypes
  
  
  return(m)
  
}


test <- sapply(pop$transmitter_log10SPVL, FoundingVars, sd = 0.5, p_dists = prob_dists) %>% t()
variants <- apply(test, 1, function(x) sum(!is.na(x)))

recip_max <- apply(test, 1, max, na.rm = TRUE)

recip_mean <-  apply(test, 1, mean, na.rm = TRUE)

v1 <- cbind.data.frame(transmitter_log10SPVL = pop$transmitter_log10SPVL, test)



pop1 <- cbind.data.frame(recipient_log10SPVL = pop$recipient_log10SPVL, transmitter_log10SPVL = pop$transmitter_log10SPVL, hypothesis = 'null', v = variants)
pop2 <- cbind.data.frame(recipient_log10SPVL = 10**recip_max, transmitter_log10SPVL = pop$transmitter_log10SPVL, hypothesis = 'max', v = variants)
pop3 <- cbind.data.frame(recipient_log10SPVL = 10**recip_mean, transmitter_log10SPVL = pop$transmitter_log10SPVL, hypothesis = 'mean', v = variants)

pops_com <- rbind.data.frame(pop1, pop2, pop3)



# CD4 Extension

pops_decline <- pops_com %>% 
  mutate(CD4_decline = (0.0111*log10(recipient_log10SPVL))**2) 

ggplot(pops_decline, aes(x = recipient_log10SPVL, y = decline, colour = hypothesis, shape = as.character(v))) +
  geom_point(alpha = 0.5) +
  scale_x_log10(limits = c(1, 10**10),
                expand = c(0,0),
                name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = "CD4 Decline", expand = c(0,0))+
  #scale_y_log10(name = expression(paste("CD4 Decline")),
                #limits = c(1, 10**10),
               # expand = c(0,0)#,
                #breaks = trans_breaks("log10", function(x) 10**x),
                #labels = trans_format("log10", math_format(.x))
                #) +
  theme_classic()+
  theme(
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )+
  annotation_logticks(sides = 'b') 



ggplot(pops_decline, aes(x = v, y = decline, colour = hypothesis, shape = as.character(v))) +
  geom_boxplot() +
  theme_classic()+
  theme(
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )



