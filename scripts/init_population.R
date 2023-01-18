# Initialise study population
# SPVL distributions from Hollingsworth et al., PLoS 2010 https://doi.org/10.1371/journal.ppat.1000876

# Dependencies
source('./scripts/depedencies.R')


# Set Seed
set.seed(4472)


# Initialisation function
InitPop <- function(N, H2, donor_vl, recipient_vl){
  require(MASS)
  require(tidyverse)
  
  matrix <- matrix(c(1, 1-(1/0.33),
                     1-(1/0.33), 1), ncol = 2)
  
  d <-diag(c(donor_vl[["sd"]], recipient_vl[["sd"]]))
  
  S <- d * matrix * d
  
  paired_spvl <- MASS::mvrnorm(n = sample_size, 
                               mu = c('transmitter' = donor_vl[["mean"]],
                                      'recipient' = recipient_vl[["mean"]]), 
                               Sigma = S) %>% 
    data.frame() %>%
    mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{.col}_log10SPVL"))
  
  return(paired_spvl)
}


# Distributions of SPVL for couples with strong support for transmission
index <- c('mean'= 4.61, 'sd' = 0.63) 
secondary <- c('mean' = 4.60 , 'sd' = 0.85)


# Run Initialisation
InitPop(N = 500, 
        H2 = 0.33, 
        donor_vl = index, 
        recipient_vl = secondary)


# Plot 
ggplot(dist, aes(x = transmitter_log10SPVL, recipient_log10SPVL)) +
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
  theme_classic()+
  theme(
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )+
  annotation_logticks() 

## END ##
