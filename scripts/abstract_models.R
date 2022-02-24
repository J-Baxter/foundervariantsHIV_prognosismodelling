#Considering VL a quantitative trait


# Hypotheses:
# 1) each virus particle transmitted is associated with a viral load (within N(donor, 1)) and the 
#    ‘fittest’ (read highest) particle establishes infection. Hence, in multi-variant infection
#    there is a greater probability of observing higher viral loads.

# 2) When two particles are transmitted, there is an additive effect between their assigned viral
#    loads, such that the viral load in the recipient is equal to the sum of the viral loads of 
#    the transmitted particles 


# 3) When two particles are transmitted, there is an multiplicative effect between their assigned 
#    viral loads, such that the viral load in the recipient is greater than the sum of the viral 
#    loads of the transmitted particles

# 4) When two particles are transmitted, there is a stabilising effect between their assigned viral 
#   loads, such that the viral load in the recipient is the median of the the viral loads of the
#   transmitted particles (ie one transmitted particle is very high, one low, and the recipient vl
#   is moderately high).


#Dependencies
library(tidyverse)
library(ggpmisc)
library(ggpubr)
library(ggprism)
library(cowplot)


# Simulate population of donor infections from a known distribution of donor SPVLs
InitDonor <- function(popsize, donor_mean, donor_sd, n_particles = 10000){
  init_donor <- rnorm(popsize,  mean = donor_mean, sd = donor_sd)
  
  donor_particles <-  sapply(init_donor, rnorm, n = n_particles, sd = 0.1, simplify = FALSE)
  
  out <- setNames(donor_particles, init_donor)
  
  return(out)
}


# A given number of founders are sampled from donor particles
# The founder particle with the highest assigned vl determines the new vl in the host

# PPP <- per particle probability
# 
H1Recipient <- function(number_of_founders, donor, ppp){
  recip_particles <- sapply(donor, sample, size = number_of_founders, simplify = FALSE)
  recip_vl <- sapply(recip_particles, max, simplify = FALSE) %>% 
    unlist() %>% 
    cbind.data.frame()
  
  return(recip_vl)
}
  

# Run H1 analysis
Main <- function(input, hypothesis){
  
  npairs <- input['npairs']
  spvl_mean <- input['spvl_mean']
  spvl_sd <- input['spvl_sd']
  perpart_prob <- input['perpart_prob']
  max_founders <- input['max_founders']
  
  donor_particles <- InitDonor(npairs, spvl_mean, spvl_sd)
  
  paired_df <- data.frame('donor' = names(donor_particles) %>% as.numeric())
  
  ######################
  if (hypothesis == 'h1'){
    for (i in 1:max_founders){
      paired_df  <- cbind.data.frame(paired_df, H1Recipient(i, donor_particles))
    }
    colnames(paired_df) <- c('donor', paste(1:max_founders, sep = "")) 
    out <- paired_df
    
  }
  
  ######################
  else if (hypothesis == 'h2'){
    
  }
  
  ######################
  else if (hypothesis == 'h3'){
    
  }
  
  ######################
  else if (hypothesis == 'h4'){
    
  }
  
  return(out)
}


FormatTukey <- function(anova){
  stat <- anova %>%
    TukeyHSD()
  stat_summary <- round(stat$multiplicity,  digits = 3)
  
  xpos <- rownames(stat_summary) %>% 
    strsplit('-') %>% 
    do.call(rbind.data.frame,.) %>%
    `colnames<-` (c('group1', 'group2'))
  
  out <- cbind.data.frame(xpos, label = stat_summary[,4], y.position = 8)
  rownames(out) <- NULL
  
  return(out)
}
###################################################################################################
# Set seed
set.seed(4472)

# Data generation
# Normal distibution of log spvl from Rakkai Cohort, Hollingsworth et al. mean = 4.39, sd = 0.84
# https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1000876

NPAIRS <- 100
RAKKAI_MEAN <- 4.39
RAKKAI_SD <- 0.84
PERPART_PROB <- 4.715*1e-8
MAX_FOUNDER_VARS <- 3

input <- c('npairs' = NPAIRS, 
           'spvl_mean' = RAKKAI_MEAN,
           'spvl_sd' = RAKKAI_SD,
           'perpart_prob' = PERPART_PROB,
            'max_founders' = MAX_FOUNDER_VARS)

out <- Main(input)


###################################################################################################
# Visualisation

# Transform wide <- long format
paired_df_long <- gather(out, multiplicity, recipient, '1':as.character(max_founders)) %>%
  dplyr::mutate(multiplicity = factor(multiplicity))


# scatter plot donor ~ recipient spvl
plot_a <- ggplot() +
  geom_point(data = paired_df_long , aes(x = 10^donor, y = 10^recipient, colour = multiplicity), shape = 4) +
  theme_classic() + 
  scale_colour_manual(values = c('#66c2a4','#2ca25f','#006d2c'), 
                      labels = c('1' = 'Single Variant',
                                 '2'= 'Two Variants',
                                 '3'='Three Variants'),
                      name = 'Founder Multiplicity') +
  scale_x_log10(name = expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml^-1, ')')),
                limits = c(1, 10^10),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_log10(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml^-1, ')')),
                limits = c(1, 10^10),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(.x))) +
  annotation_logticks()



# Statistical analysis of recipient_spvl ~ founder variant multiplicity
# 1. ANOVA
# 2. Tukey pairwise comparisons
stat_anova <- aov(recipient ~ multiplicity, data = paired_df_long)
stat_anova_pval <- summary(stat_anova)[[1]][['Pr(>F)']][1] %>% round(digits = 3)

tukey_results <- FormatTukey(stat_anova)
#add_pvalue(tukey_results, xmin = 'group1', xmax = 'group2', label = 'label', y.position = 'y.position')

my_comparisons <- list(c('1', '2'), c('1', '3'))

# box plot recipient spvl ~ number of founder variants 
# incorporates additional analyses to clarify associations
plot_b <- ggplot(data = paired_df_long, aes(x = multiplicity, y = 10^recipient, colour = multiplicity)) +
  geom_boxplot() +
  theme_classic() + 
  scale_colour_manual(values = c('#66c2a4','#2ca25f','#006d2c')) +
  scale_y_log10(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml^-1, ')')),
                limits = c(1, 10^10),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(.x))) + 
  scale_x_discrete(labels = c('1' = 'Single Variant',
                              '2'= 'Two Variants',
                              '3'='Three Variants'),
                   name = 'Founder Multiplicity') +
  stat_compare_means(method = 'anova', label.x = 0.7, label.y = 9) +
  stat_compare_means(comparisons = my_comparisons, method = 't.test')+
  annotation_logticks(sides = 'l')


plot <- plot_grid(plot_a + theme(legend.position= c(0.85,0.857)), 
                  plot_b + theme(legend.position = 'none'), 
                  labels = 'AUTO', align = 'hv')

setEPS()
postscript("./figures/allelseequal_primitive.eps", width = 16, height = 10)
plot
dev.off()

  