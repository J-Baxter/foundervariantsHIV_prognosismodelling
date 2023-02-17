source('./scripts/dependencies.R')


# Read in and format data
julian_tp_data <- read_csv("data/julian_tp_data.csv", 
                           col_types = cols(
                             sequences.n.pair = col_number(), 
                             sequences.n.source = col_number(), 
                             sequences.n.recipient = col_number(), 
                             HBX2.coordinates.ini = col_number(), 
                             HBX2.coordinates.end = col_number(), 
                             sequence.length = col_number(), 
                             VL.transmitter = col_number(), 
                             VL.recipient = col_number(), 
                             VL.transmitter.other = col_number()))


aba5443_data_s2 <- read_csv("data/aba5443_data_s2.csv", 
                            col_types = cols(particle.model.best.fit = col_number(), 
                                             variants.mean = col_number(), variants.mode = col_number(), 
                                             variants.max = col_number()))


test_data <- left_join(julian_tp_data, aba5443_data_s2) %>%
  filter(!is.na(VL.transmitter)) %>%
  filter(!is.na(VL.recipient)) %>%
  mutate(variants.multiple = ifelse(variants.mean == 1, 'single', 'multiple')) %>%
  mutate(across(contains('VL'), .fns = log10, .names = 'log10{.col}')) %>%
  #mutate(log.VL.transmitter = log10(VL.transmitter)) %>%
  #mutate(log.VL.recipient = log10(VL.recipient)) %>%
  separate(type, c('transmitter', 'recipient'),  -1) %>%
  dplyr::select(transmitter, recipient, log10VL.transmitter, log10VL.recipient, VL.transmitter, VL.recipient,subtype, variants.mean, variants.multiple) %>%
  `colnames<-` (str_remove(colnames(.), '[:punct:]')) %>%
  mutate(couple = seq(1,nrow(.)))
  
# Null
prob_recip_multiple_test <- RunParallel(populationmodel_acrossVL_Environment, test_data$VLtransmitter) %>%
  do.call(cbind.data.frame, .) %>% t()

transmitterspvl_variantdist_test <- cbind.data.frame(transmitter = test_data$VLtransmitter, prob_recip_multiple_test)

ggplot(transmitterspvl_variantdist_test, 
       aes(x = transmitter, 
           y = 1 - variant_distribution.V1))+
  geom_point(colour = '#ef654a',  
             shape = 4, size = 3)+
  scale_x_log10(name = expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0,0),
                     limits = c(0,0.6),
                     breaks = seq(0, 0.6, by = 0.1)) +
  annotation_logticks(sides = 'b') + my_theme

################################### P(Multiple Variants) | Recipient SPVL (H0) ###################################
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


transmitter_prob_test <- sapply(test_data$log10VLtransmitter, function(x) WeightPDF(x)/sum(WeightPDF(test_data$log10VLtransmitter)))
hist(transmitter_prob_test) #visual check - should look normal(ish) as already log


# Initialise simulation donor population from a specified min and max log spvl
InitSimTransmitter <- function(pop_size, transmitter_min, transmitter_max, sample_prob){
  sim_range <- seq(transmitter_min, transmitter_max, length.out = pop_size)
  sim_logspvl <- sample(sim_range, size = pop_size, prob = sample_prob, replace = T)
  
  return(sim_logspvl)
}


sim_transmitter_log10SPVL_test <- InitSimTransmitter(pop_size = 47,
                                                transmitter_min = 1, 
                                                transmitter_max = 7, 
                                                sample_prob = transmitter_prob_test)
hist(sim_transmitter_log10SPVL_test ) #visual check - should look uniform


# Leverage heritability model to estimate recipient SPVL
sim_spvl <- predict(h2_model, newdata= data.frame(transmitter_log10SPVL = sim_transmitter_log10SPVL)) %>% 
  cbind.data.frame(sim_recipient_log10SPVL = ., sim_transmitter_log10SPVL = sim_transmitter_log10SPVL) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SPVL')}")) %>%
  `colnames<-` (str_remove(colnames(.), 'sim_'))


# Calculate probability of mulitple founder infection in recipient
sim_prob_multiple <- RunParallel(populationmodel_acrossVL_Environment, sim_spvl$transmitter)  %>%
  do.call(cbind.data.frame, .) %>% t()

sim_recipientspvl_variantdist <- cbind.data.frame(recipient = sim_spvl$recipient, sim_prob_multiple)

##############3
sim_prob_multiple_test <- RunParallel()  %>%
  do.call(cbind.data.frame, .) %>% t()

sim_recipientspvl_variantdist_test <- cbind.data.frame(transmitter = test_data$VLtransmitter, sim_prob_multiple_test)





# Plot recip_vl ~ donor_vl

lm <- lmer(log10VLrecipient ~ log10VLtransmitter + subtype + recipient + (log10VLtransmitter|couple), data = test_data)

#Hollingsworth
test_data <- left_join(julian_tp_data, aba5443_data_s2) %>%
  filter(!is.na(VL.transmitter)) %>%
  filter(!is.na(VL.recipient)) %>%
  mutate(variants.multiple = ifelse(variants.mean == 1, 'single', 'multiple')) %>%
  mutate(across(contains('VL'), .fns = log10, .names = 'log10{.col}')) %>%
  #mutate(log.VL.transmitter = log10(VL.transmitter)) %>%
  #mutate(log.VL.recipient = log10(VL.recipient)) %>%
  separate(type, c('transmitter', 'recipient'),  -1) %>%
  dplyr::select(transmitter, recipient, log10VL.transmitter, log10VL.recipient, subtype, variants.mean, variants.multiple) %>%
  `colnames<-` (str_remove(colnames(.), '[:punct:]')) %>%
  mutate(couple = seq(1,nrow(.))) %>%
  pivot_longer(cols = contains('log10VL'), values_to = 'Log10VL', names_to = 'position') %>%
  mutate(position = str_remove(position, 'log10VL'))

lm <- lmer(Log10VL ~ position + subtype + recipient + (1|couple), data = test_data)
lmer(log10VL.recipient ~ log10VL.transmitter + subtype + recipient + (1|couple), data = test_data)

nlform <- bf(log10VLrecipient ~ log10VLtransmitter,
             log10VLrecipient ~ 1)

nlprior <- c(
             prior(normal(0,1), nlpar = "log10VLrecipient"))

fit_loss1 <- brm(formula = nlform, data = test_data, family = gaussian(),
                 prior = nlprior, control = list(adapt_delta = 0.9))


ggplot(aes(x = log10.VL.transmitter, y= log10.VL.recipient, colour = variants.multiple, shape = recipient), data = test_data) + geom_point(size = 3) +
  scale_shape(labels = c('Female-to-Male', 'Male-to-Female', 'MSM'), name = 'Exposure')+
  scale_colour_discrete(labels = c('Multiple', 'Single'), name = 'Founder Variants')+
  theme_classic()+
  scale_x_continuous(name = 'Transmitter Log10 SPVL', expand = c(0,0), limits = c(0, 9))+
  scale_y_continuous(name = 'Recipient Log10 SPVL', expand = c(0,0), limits = c(0, 9)) + 
  geom_smooth(method = 'lm', aes(group = variants.multiple))


nlform <- bf(log10.VL.recipient ~ Additive(variants) + recipient + E,
             recipient~ 1, 
             E ~ 1,
             variants ~ rnorm(N, log10VLrecipient, sd),
             log10VLrecipient ~ 1,
             N ~ 1,
             sd ~ 1,
             nl = TRUE)

nlprior <- c(prior(normal(0,1), nlpar = "recipient"),
             prior(normal(0,1), nlpar = "E"),
             prior(normal(0,1), nlpar = "sd"),
             prior(constant(1), nlpar = "N"),
             prior(constant(1), nlpar = "log10VLrecipient"))

fit_loss1 <- brm(formula = nlform, data = test_data, family = gaussian(),
                 prior = nlprior, control = list(adapt_delta = 0.9))



nlform <- bf(CD4_decline ~ alpha * SPVL_log10^2 + sex,
             SPVL_log10 ~ Additive(variants) + E,
             sex ~ 1, 
             alpha ~ 1,
             E ~ 1,
             variants ~ rnorm(theta, T, x),
             x ~ 1,
             nl = TRUE)






CD4_decline = alpha * SPVL_log10^2 + age + sex
SPVL_log10 = func(variants)
variants = & normal distirbution with mean (transmitter_vl)

transmitter_vl (fixed)