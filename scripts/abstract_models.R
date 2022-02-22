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

# Set seed
set.seed(4472)

# Data generation
# Normal distibution of log spvl from Rakkai Cohort, Hollingsworth et al. mean = 4.39, sd = 0.84
# https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1000876

NPAIRS <- 1000 # lower for test runs on low cpu machines
RAKKAI_MEAN <- 4.39
RAKKAI_SD <- 0.84

# H1
InitDonor <- function(popsize, donor_mean, donor_sd, n_particles = 100){
  init_donor <- rnorm(popsize,  mean = donor_mean, sd = donor_sd)
  
  donor_particles <-  sapply(init_donor, rnorm, n = n_particles, sd = 0.1, simplify = FALSE)
  
  out <- setNames(donor_particles, init_donor)
  
  return(out)
}


H1Recipient <- function(number_of_founders, donor){
  recip_particles <- sapply(donor, sample, size = number_of_founders, simplify = FALSE)
  recip_vl <- sapply(recip_particles, max, simplify = FALSE) %>% 
    unlist() %>% 
    cbind.data.frame()
  
  return(recip_vl)
}
  

donor_particles <- InitDonor(NPAIRS,  RAKKAI_MEAN, RAKKAI_SD, n_particles = 100)

single_recip_vl <- H1Recipient(1, donor_particles) 
double_recip_vl <- H1Recipient(2, donor_particles) 
triple_recip_vl <- H1Recipient(3, donor_particles) 

paired_df <- data.frame('donor' = names(donor_particles) %>% as.numeric(), 
                        'single' = single_recip_vl[,1],  
                        'double' = double_recip_vl[,1],
                        'triple' = triple_recip_vl[,1])
paired_df_long <- gather(paired_df, multiplicity, recipient, single:triple)

paired_df_long$multiplicity <- factor(paired_df_long$multiplicity, levels = c('single', 'double', 'triple')) 

plot_a <- ggplot() +
  geom_point(data = paired_df_long , aes(x = 10^donor, y = 10^recipient, colour = multiplicity), shape = 4) +
  theme_classic() + 
  scale_colour_manual(values = c('#e5f5f9','#99d8c9','#2ca25f')) +
  scale_x_log10(name = expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml^-1, ')')),
                limits = c(1, 10^10),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_log10(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml^-1, ')')),
                limits = c(1, 10^10),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(.x)))


stat_anova <- aov(recipient ~ multiplicity, data = paired_df_long)
stat_anova_pval <- summary(stat_anova)[[1]][['Pr(>F)']][1] %>% round(digits = 3)

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

tukey_results <- FormatTukey(stat_anova)
#add_pvalue(tukey_results, xmin = 'group1', xmax = 'group2', label = 'label', y.position = 'y.position')

my_comparisons <- list(c('single', 'double'), c('single', 'triple'))

plot_b <- ggplot(data = paired_df_long, aes(x = multiplicity, y = 10^recipient, colour = multiplicity)) +
  geom_boxplot() +
  theme_classic() + 
  scale_colour_manual(values = c('#e5f5f9','#99d8c9','#2ca25f')) +
  scale_y_log10(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml^-1, ')')),
                limits = c(1, 10^10),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(.x))) + 
  stat_compare_means(method = 'anova', label.x = 0.7, label.y = 9) +
  stat_compare_means(comparisons = my_comparisons, method = 't.test')

library(cowplot)
plot <- plot_grid(plot_a, plot_b, labels = 'AUTO', align = 'hv')

setEPS()
postscript("./figures/allelseequal_primitive.eps", width = 16, height = 10)
plot
dev.off()

  