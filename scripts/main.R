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
source('./scripts/init_population.R')
source('./scripts/populationdata_models.R')
source('./scripts/populationdata_acrossVL_models.R')


# Null
prob_recip_multiple <- RunParallel(populationmodel_fixedVL_Environment, pop$transmitter) %>%
  do.call(cbind.data.frame, .) %>% t()

# Alternative

# Generate a matrix of transmitter/founding particles and associated viral loads
test <- sapply(pop$transmitter, FoundingVars, sd = 0.5, p_dists = prob_dists) %>% t()
variants <- apply(test, 1, function(x) sum(!is.na(x)))

recip_max <- apply(test, 1, max, na.rm = TRUE)

recip_mean <-  apply(test, 1, mean, na.rm = TRUE)

pop1 <- cbind.data.frame(recipient_log10SPVL = pop$recipient_log10SPVL, transmitter_log10SPVL = pop$transmitter_log10SPVL, hypothesis = 'null', v = variants)
pop2 <- cbind.data.frame(recipient_log10SPVL = 10**recip_max, transmitter_log10SPVL = pop$transmitter_log10SPVL, hypothesis = 'max', v = variants)
pop3 <- cbind.data.frame(recipient_log10SPVL = 10**recip_mean, transmitter_log10SPVL = pop$transmitter_log10SPVL, hypothesis = 'mean', v = variants)

pops_com <- rbind.data.frame(pop1, pop2, pop3)

ggplot(pops_com, aes(x = transmitter_log10SPVL, recipient_log10SPVL, colour = hypothesis)) +
  geom_point(alpha = 0.5) +
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
  facet_wrap(.~v)+
  theme_classic()+
  theme(
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )+
  annotation_logticks() 

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
