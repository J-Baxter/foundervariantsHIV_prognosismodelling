###################################################################################################
###################################################################################################
# Modelling the association between the number of variants and viral load
# Script utilises methodology from Thompson et al 2019, https://github.com/robin-thompson/MultiplicityOfInfection
# implemented in R by Katie Atkins https://github.com/katiito/HIV_eclipse_models

###################################################################################################
###################################################################################################
# Dependencies
source("./scripts/populationdata_models.R")
source("./scripts/dependencies.R")

# Calculate recipient set point viral loads from a given population of donor set point viral loads.
# Back calculation from a set correlation coefficient (r2) and sum of squared residuals
# 'Parent ~ Offspring' regression R2 analogous to heritability of viral load. 
InitDonor <- function(carrier, TP){
  #stopifnot() #stop if TP does not contain mean and variance
  
  carrier_mean <- carrier[['vl_mean']]
  carrier_var <- carrier[['vl_sd']] ** 2
  carrier_h2 <- carrier[['h2']]
  TP_mean <- TP[['tp_mean']]
  TP_var <- TP[['tp_sd']] ** 2
  
  # Calculate donor population summary statistics from B8 (Eq 7) Bonhoeffer
  donor_mean <- (carrier_mean*TP_var + TP_mean*carrier_var) / (carrier_var + TP_var)
  donor_sd <- ((carrier_var*TP_var) / (carrier_var + TP_var)) %>% sqrt()
  
  return(c(vl_mean = donor_mean, vl_sd = donor_sd, h2 = carrier_h2))
  
  
}


GetRecipVL <- function(donorload, h2 = 0.33){
  n <- length(donorload)
  ssr <- sum((donorload-mean(donorload))**2) # sum of squared residuals
  e <- rnorm(n)
  e <- resid(lm(e ~ donorload))
  e <- e*sqrt((1-h2)/h2*ssr/(sum(e**2)))
  recipload <- donorload + e
  
  return(recipload)
}



# Initialise donor~recipient pairs from a given distribution of donor logspvls and 
# regression coefficient corresponding to the heritability of viral load
InitDonorRecip <- function(popsize, donor_population){
  #must contain xyz
  
  donor_mean <- donor_population[['vl_mean']]
  donor_sd <- donor_population[['vl_sd']]
  donor_h2 <- donor_population[['h2']]
  
  init_donor <- rnorm(popsize,  mean = donor_mean, sd = donor_sd)
  init_recip <- -1
  
  while(min(init_recip)<0){
    init_recip <- GetRecipVL(init_donor, donor_h2) 
  }
  
  out <- list('donor' = cbind.data.frame(vl=init_donor) , 'recipient' = cbind.data.frame(vl = init_recip))
  return(out)
}



# Infer weighting for a given viral load 
# Weighting 'g' from Thompson et al 2019
WeightPDF <- function(viralload){
  alpha = -3.55
  sigma <- 0.78/(sqrt(1 - ((2*alpha**2)/(pi*(1 + alpha**2)))))
  mu <-  4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha**2)))
  weight <-  (2/sigma)*dnorm((log10(viralload) - mu)/sigma)*pnorm(alpha*(log10(viralload) - mu)/sigma)
  return(weight)
} 


# Initialise simulation donor population from a specified min and max log spvl, generalised with
# a sampling distribution from WeightPDF (function 'g' from Thompson 2019 et al)
InitSimDonor <- function(pop_size, donor_min, donor_max, sample_prob){
  sim_range <- seq(donor_min, donor_max, length.out = pop_size)
  sim_logspvl <- sample(sim_range, size = pop_size, prob = sample_prob, replace = T) %>%
    cbind.data.frame() %>%
    `colnames<-`('donor')
  return(sim_logspvl)
  
}



###################################################################################################
# Set seed
set.seed(4472)

# Data generation
# Normal distibution of log spvl from Rakkai Cohort, Hollingsworth et al. mean = 4.39, sd = 0.84
# https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1000876


NPAIRS <- 500 # lower for test runs on low cpu machines


# Initial distributions describe carrier population
zambia_carrier <- c(vl_mean = 4.74, vl_sd = 0.61, h2 = NA)
amsterdam_carrier <- c(vl_mean = 4.35, vl_sd = 0.68, h2 = NA)
rakkai_carrier <- c(vl_mean = 4.51, vl_sd = 0.79, h2 = 0.33)
TP <- c(tp_mean = 4.64, tp_sd = 0.96)

zambia_donor <- InitDonor(zambia_carrier, TP)
amsterdam_donor <- InitDonor(amsterdam_carrier, TP)
rakkai_donor <- InitDonor(rakkai_carrier, TP) #Must check whether TP can be generally applicable or requires recalculation
rakkai_donor <- rakkai_carrier

rakkai_donorrecip <- InitDonorRecip(NPAIRS, rakkai_donor)


rakkai_cohort <- c(list('carrier' = cbind.data.frame(vl = rnorm(NPAIRS, rakkai_carrier[['vl_mean']], rakkai_carrier[['vl_sd']]))), 
                   rakkai_donorrecip) %>% 
  bind_rows(.id = 'population') %>%
  mutate(vl = 10**vl)

ggplot(rakkai_cohort, aes(x = vl, fill = population)) +
  geom_density(
    alpha = 0.5) +
  scale_x_log10(limits = c(1, 10**10),
                expand = c(0,0),
                name = expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(expand = c(0,0))+
  theme_classic() 


rakkai_pairs <- rakkai_cohort %>% filter(population %in% c('donor', 'recipient')) %>%  
  group_by(population) %>%
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = population, values_from = vl) %>%
  select(-row)



###################################################################################################
# Probability that recipient infection is initiated by multiple founder variants
# applies populationmodel_fixedVL_Environment function written by Katie Atkins
prob_recip_multiple <- RunParallel(populationmodel_fixedVL_Environment, rakkai_pairs$donor) %>%
  do.call(cbind.data.frame, .) %>% t() # Time difference of 44.60442 mins on Macbook

combined_data <- cbind.data.frame(rakkai_pairs, prob_recip_multiple)
head(combined_data)


###################################################################################################
# Probability that recipient infection is multiple founder, given a certain viral load

# Generate weightings using function g from Thompson et al
# Function generates a probability distribution (lognormal) which is used to sample from our range
# of simulated donor viral loads to generalise over a population
donor_prob <- sapply(log10(rakkai_pairs$donor), function(x) WeightPDF(x)/sum(WeightPDF(log10(rakkai_pairs$donor))))
hist(donor_prob) #visual check - should look normalish as already log

sim_donor_logspvl <- InitSimDonor(pop_size = NPAIRS,
                                  donor_min = 1, 
                                  donor_max = 7, 
                                  sample_prob = donor_prob)

hist(sim_donor_logspvl$donor)

# Infer recipient viral loads of from simulated population
# predict.lm(h2_model, cbind.data.frame(donor_logspvl = sim_donor_logspvl)) (identical to donor)
# GetRecipVL(sim_donor_logspvl, R2) (too varied) suspect problem is attempting to pace a normal dist over uniform?
# Implementation below uses coefficients from recipient ~ donor model and normal dist error term
h2_model <- lm(recipient ~ donor, data = log10(rakkai_pairs))
sim_spvl <- predict(h2_model, newdata=sim_donor_logspvl) %>% 
  cbind.data.frame(recipient = ., donor = sim_donor_logspvl) %>%
  mutate(across(.cols = everything()), 10 **. )


# Calculate probability of mulitple founder infection in recipient
sim_prob_multiple <- RunParallel(populationmodel_fixedVL_Environment, sim_spvl$donor)  %>%
  do.call(cbind.data.frame, .) %>% t()

sim_combined_data <- cbind.data.frame(sim_spvl, sim_prob_multiple)
head(sim_combined_data)


###################################################################################################
# Visualise probability that recipient infection is intitiated by multiple founders

# Fig 1a

font_add_google("Questrial", "Questrial")
showtext_auto()

fig_1a <- 
  ggplot(rakkai_pairs, aes(x = donor, recipient)) +
  geom_point(colour = '#ef6548', #'#CB6015' #'#66c2a4','#2ca25f','#006d2c'
             alpha = 0.5) +
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
  theme_classic(base_family = "Questrial")+
  theme(
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )+
  annotation_logticks() + 
  geom_smooth(method = lm, colour = 'black', se = F) + 
  stat_poly_eq(formula = y ~ x)

setEPS()
<<<<<<< HEAD
#postscript("donor-recipient.eps", width = 10, height = 10)
ggsave(plot = fig_1a, 'donorrecip.eps', width = 12, height = 12, device = cairo_ps)
fig_1a
dev.off()

=======
postscript("donor-recipient.eps", width = 10, height = 10)
fig_1a
dev.off()
  
>>>>>>> 5045211448c88d1384ee7f31d06b6a88b1ff2c86
# Fig 1b
fig_1b <- ggplot(combined_data, 
                 aes(x = donor, 
                     y = 1 - variant_distribution.V1))+
  geom_point(colour = '#2ca25f',  #usher colours?
             alpha = 0.5)+
  scale_x_log10(name = expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0,0),
                     limits = c(0,0.6),
                     breaks = seq(0, 0.6, by = 0.2)) +
  theme_classic() 

# Fig 1c
fig_1c <- ggplot(sim_combined_data, 
                 aes(x = recipient,
                     y = 1 - variant_distribution.V1))+
  geom_point(colour = '#2ca25f',  #usher colours?
             alpha = 0.5)+
  scale_x_log10(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0,0),
                     limits = c(0,0.6),
                     breaks = seq(0, 0.6, by = 0.2)) +
  theme_classic() 

# Panel 1
panel1 <- plot_grid(fig_1a, fig_1b, fig_1c, labels = 'AUTO', align = 'hv', ncol = 3)


###################################################################################################
# Visualise probability distributions of the number of variants that initiate recipient infection ~
# recipient set point viral load 
# AKA panel 2
lab <- expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml**-1, ')'))

#store data long form & categorise recipient viral loads
sim_combined_data_long <- sim_combined_data %>%
  cbind(., ref = 1:nrow(sim_combined_data)) %>%
  gather(key = "variant_no",value = 'variant_prob',variant_distribution.V1:variant_distribution.V33) %>%
  mutate(spvl_cat = cut(recipient, 
                        breaks = 10**(0:9),
                        labels = 1:9))

# replace categories of number of variants with integers
sim_combined_data_long$variant_no <- rep(1:33, each = 1000)

# within each spvl cat. calculate median probability for number of variants
median_probs <- sim_combined_data_long %>%
  split.data.frame(sim_combined_data_long$spvl_cat) %>%
  lapply(., function(x) {group_by(x, variant_no) %>%
      summarise(variant_prob = median(variant_prob))}) %>%
  bind_rows(.,.id = 'source')

colnames(median_probs)[1] <- 'spvl_cat'

# alt panel 2 taking probability distributions from single specified donor spvls
panel2 <- cbind.data.frame(donor = 1:9, recipient = predict(h2_model, cbind.data.frame(donor = 1:9))) %>% # this isn't really working
  mutate(across(.cols = everything()), )
panel2_donors <- 10 ** panel2_donors
panel2_recips <- 10 ** panel2_recips

panel2_probs <- RunParallel(populationmodel_fixedVL_Environment, panel2$donor)  %>%
  do.call(cbind.data.frame, .) %>% t()
panel2_data <- cbind.data.frame(panel2_donors, panel2_recips, panel2_probs)

panel2_data_long <- panel2_data %>%
  cbind(., ref = 1:nrow(panel2_data)) %>%
  gather(key = "variant_no",value = 'variant_prob',variant_distribution.V1:variant_distribution.V33) %>%
  mutate(spvl_cat = cut(panel2_donors, 
                        breaks = 10**(0:9),
                        labels = 1:9))

# replace categories of number of variants with integers
panel2_data_long$variant_no <- rep(1:33, each = 9)

panel2_labeller <- as_labeller(sapply(1:9, function(x) paste(x, "~log[10]~copies~ml**-1")) %>% `names<-` (1:9),
                               default = label_parsed)

panel2 <- ggplot(median_probs, 
                 aes(x = variant_no, 
                     y =variant_prob)) + 
  geom_bar(stat = 'identity', fill = '#2ca25f') + 
  scale_y_continuous(name = 'Probability', 
                     expand = c(0,0),
                     limits = c(0,0.8),
                     breaks = seq(0,1,0.2))+
  scale_x_continuous(name = 'Number of Variants', 
                     expand = c(0,0.6),
                     #limits = c(1,12),
                     breaks = seq(0,12,1))+
  facet_wrap(~spvl_cat,
             labeller = panel2_labeller) +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines")) +
  coord_cartesian(xlim = c(1,8))

# Output to file
paste(figs_dir, 'panel1.jpeg', sep = '/')

jpeg(filename = paste(figs_dir, 'panel1.jpeg', sep = '/'), width = 23, height = 9, unit = 'cm', res = 350)
panel1 
dev.off()

setEPS()
ggsave(paste(figs_dir, 'panel1.eps', sep = '/'), width = 23, height = 9, device = cairo_ps)
panel1 
dev.off()


jpeg(filename = paste(figs_dir, 'panel2.jpeg', sep = '/'), width = 18, height = 18, unit = 'cm', res = 350 )
panel2 
dev.off()

setEPS()
postscript(paste(figs_dir, 'panel2.eps', sep = '/'), width = 18, height = 18, unit = 'cm')
panel2 
dev.off()
###################################################################################################
#                                                                                                 #
#                                               END                                               #
#                                                                                                 #
###################################################################################################