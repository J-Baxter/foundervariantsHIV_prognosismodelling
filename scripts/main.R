# Main Analysis for "RECONCILING THE FOUNDER VARIANT MULTIPLICITY OF HIV-1 INFECTION WITH THE RATE
# OF CD4+ DECLINE"
# J Baxter

# This script combines the functions for the analysis of paired SPVL distributions, occuring in
# three main stages:

# 1. Initialisation of study population
# 2. Null hypothesis
# 3. Ecological hypotheses


################################### Dependencies ###################################
source('./scripts/dependencies.R')
source('./scripts/populationdata_acrossVL_models.R')
source('./scripts/init_population.R')
source('./scripts/ecologicalhypos.R')
source('./scripts/cd4_delcine.R')

# Set Seed
set.seed(4472)
options(scipen = 100) #options(scipen = 100, digits = 4)


################################### Initialise Population ###################################
# Distributions of SPVL for couples with strong support for transmission Hollingsworth et al, 2010
# Dummy population until we incorporate Rakai data

index <- c('mean'= 4.61, 'sd' = 0.63) 
secondary <- c('mean' = 4.60 , 'sd' = 0.85)

NPAIRS <- 50 #500

pop <- InitPop(N = NPAIRS, 
               H2 = 0.33,
               donor_vl = index, 
               recipient_vl = secondary) %>%
  cbind.data.frame(sex = sample(c('M', 'F'), nrow(.), replace = T)) %>%
  cbind.data.frame(age = sample(18:50, nrow(.), replace = T))


################################### Estimate Heritability ###################################
h2_model <- lm(recipient_log10SPVL ~ transmitter_log10SPVL, data = test_data) # Building this H2 model is a key step BRMS?

h2_priors <- prior(normal(1, 2), nlpar = "b1") + prior(normal(0, 1), nlpar = "b2") + prior(normal(0, 2), nlpar = "b3") 

h2_fit <- brm(bf(recipient_log10SPVL ~ b1*transmitter_log10SPVL + b2*sex + b3*age, b1 + b2 + b3 ~ 1, nl = TRUE),
              data = pop, prior = h2_priors )

summary(h2_fit) #Check output makes sense + ESS
plot(h2_fit) # Check convergence
bayes_R2(h2_fit) #Estimate R2

pred.data = expand.grid(transmitter_log10SPVL = seq(min(pop$transmitter_log10SPVL), max(pop$transmitter_log10SPVL), length=20),
                        sex = c('M', 'F'),
                        age = seq(min(pop$age), max(pop$age), length=20))

pred.data <- pred.data %>% cbind.data.frame(., predict(h2_fit, newdata=pred.data))

pop%>%
  add_predicted_draws(h2_fit) %>%  # adding the posterior distribution
  ggplot(aes(x = transmitter, y = recipient)) +  
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                  alpha = 0.5, colour = "black") +
  geom_point(data = France, colour = "darkseagreen4", size = 3)

ggplot()+
  geom_point(data = pop, aes(x = transmitter, recipient), #'#CB6015' #'#66c2a4','#2ca25f','#006d2c'
    colour = '#ef654a',
    shape = 4, size = 3) +
  scale_x_log10(limits = c(1, 10**10),
                expand = c(0,0),
                name = expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_log10(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**10),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  ggplot(pred.data, aes(y = Estimate, x = transmitter_log10SPVL))+
  geom_line()
  
  annotation_logticks() +
  my_theme

################################### Initialise Simulated Populations ###################################
transmitter_prob <- sapply(pop$transmitter_log10SPVL, function(x) WeightPDF(x)/sum(WeightPDF(pop$transmitter_log10SPVL)))

#hist(transmitter_prob) #visual check - should look normal(ish) as already log

sim_transmitter_log10SPVL <- InitSimTransmitter(pop_size = NPAIRS, transmitters = pop
                                                #transmitter_min = 1, 
                                                #transmitter_max = 7, 
                                                sample_prob = transmitter_prob)

#hist(sim_transmitter_log10SPVL) #visual check - should look uniform

sim_spvl <- predict(h2_model, newdata= data.frame(transmitter_log10SPVL = sim_transmitter_log10SPVL)) %>% 
  cbind.data.frame(sim_recipient_log10SPVL = ., sim_transmitter_log10SPVL = sim_transmitter_log10SPVL) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SPVL')}")) %>%
  `colnames<-` (str_remove(colnames(.), 'sim_'))


################################### P(Multiple Variants) | Transmitter SPVL (H0) ###################################
# Probability that recipient infection is initiated by multiple founder variants
# applies populationmodel_acrossVL_Environment function written by Katie Atkins
# extracts estimates for 1) Chronic transmitters and 2) Integrated over the course of transmitter infection

transmitterspvl_variantdist <- RunParallel(populationmodel_acrossVL_Environment, pop$transmitter) %>%
  
  # Post-Processing 
  do.call(cbind.data.frame, .) %>% 
  t() %>%
  cbind.data.frame(transmitter = pop$transmitter, .) %>%
  rename(probTransmissionPerSexAct_average = probTransmissionPerSexAct) %>%
  rename_with(function(x) gsub('\\.(?<=\\V)', '\\1_average.\\2',x, perl= TRUE), 
              .cols = !contains('chronic') & contains('V')) %>%
  
  # Long Format
  pivot_longer(cols = c(contains('chronic'),contains('average')), 
               names_to = c('type', 'stage', 'variants'), 
               names_pattern = "(^[^_]+)(_[^.]*)(.*)", 
               values_to = 'p') %>% 
  mutate(stage = gsub('^.*_', '', stage)) %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric()) 


################################### P(Multiple Variants) | Recipient SPVL (H0) ###################################
# Generate weightings using function g from Thompson et al
# This function generates a probability distribution (lognormal) which is used to sample from our range
# of simulated transmitter viral loads to generalise over a population


# Calculate probability of multiple founder infection in recipient
sim_recipientspvl_variantdist <- RunParallel(populationmodel_acrossVL_Environment, sim_spvl$transmitter)  %>%
  do.call(cbind.data.frame, .) %>% 
  t() %>%
  cbind.data.frame(recipient = sim_spvl$recipient, .) %>%
  rename(probTransmissionPerSexAct_average = probTransmissionPerSexAct) %>%
  rename_with(function(x) gsub('\\.(?<=\\V)', '\\1_average.\\2',x, perl= TRUE), 
              .cols = !contains('chronic') & contains('V')) %>%
  
  # Long Format
  pivot_longer(cols = c(contains('chronic'),contains('average')), 
               names_to = c('type', 'stage', 'variants'), 
               names_pattern = "(^[^_]+)(_[^.]*)(.*)", 
               values_to = 'p') %>% 
  mutate(stage = gsub('^.*_', '', stage)) %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric()) %>%
  
  # Incorporate CD4 decline
  mutate(CD4_decline = CD4Decline(recipient_log10SPVL, age, sex))


################################### Write to file ###################################
# population only
write_csv(pop, file = paste(results_dir, 'base_population.csv', sep = '/'))

# transmitter SPVL & variant distribution
write_csv(transmitterspvl_variantdist, file = paste(results_dir, 'transmitterspvl_variantdist.csv', sep = '/'))

# recipient SPVL & variant distribution
write_csv(sim_recipientspvl_variantdist, file = paste(results_dir, 'recipientspvl_variantdist.csv', sep = '/'))


###################################################################################################
################################### Non-Linear Models (No Segregation) ###################################
# Four groups of nonlinear relationships trialled: Polynomial, Concave/Convex Curves, Sigmoidal, Curves with Max/Min

quad_model <- lm(recipient_log10SPVL ~ transmitter_log10SPVL+ I(transmitter_log10SPVL**2) + age + sex , data = pop) 

exp_model <- 

asymp_model <-

logistic_model <- 

  
################################### Non-Linear Models (Probability MV weights) ###################################
# Four groups of nonlinear relationships trialled: Polynomial, Concave/Convex Curves, Sigmoidal, Curves with Max/Min

# Calculate weights based on transmitter SPVL


poly_model <- 
  
concave_model <- 
  
convex_model <- 
  
sigmoid_model <- 
  
max_model <- 
  



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

spvl <- seq(2,7,by=0.5)
prob_dists <- lapply(spvl, populationmodel_acrossVL_Environment) %>% sapply(., function(x) x[2:length(x)])

variant_matrix <- sapply(pop$transmitter_log10SPVL, FoundingVars, sd = 0.5, p_dists = prob_dists) %>% t()

# Sanity check
apply(variant_matrix, 1, function(x) sum(!is.na(x))) %>%
  all() %>%
  stopifnot()


################################### Exclusion ###################################
recip_log10SPVL_exclusion <- apply(variant_matrix, 1, Exclusion, na.rm = TRUE)


################################### Additive ###################################
recip_log10SPVL_additive <- apply(variant_matrix, 1, Additive, na.rm = TRUE)


################################### Interaction ###################################
recip_log10SPVL_interaction <- apply(variant_matrix, 1, Interaction, na.rm = TRUE)


################################### Post-Processing ###################################
v1 <- cbind.data.frame(transmitter_log10SPVL = pop$transmitter_log10SPVL, test)


pop1 <- cbind.data.frame(recipient_log10SPVL = pop$recipient_log10SPVL, transmitter_log10SPVL = pop$transmitter_log10SPVL, hypothesis = 'null', v = variants)
pop2 <- cbind.data.frame(recipient_log10SPVL = 10**recip_max, transmitter_log10SPVL = pop$transmitter_log10SPVL, hypothesis = 'max', v = variants)
pop3 <- cbind.data.frame(recipient_log10SPVL = 10**recip_mean, transmitter_log10SPVL = pop$transmitter_log10SPVL, hypothesis = 'mean', v = variants)

pops_com <- rbind.data.frame(pop1, pop2, pop3)



################################### CD4 decline ###################################

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



