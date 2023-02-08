library(tidyverse)
library(MASS)
library(scales)

sample_size <- 100

sample_meanvector <- c('transmitter'= 4.61, 'recipient' = 4.60)

sample_sdvector <- c('transmitter' = 0.63 , 'recipient' = 0.85)

test <- 1-(1/0.33)

sample_correlation_matrix <- matrix(c(1,0.33,0.33,1), ncol = 2)

d <-diag(sample_sdvector)

S <- d * sample_correlation_matrix * d

dist = mvrnorm(n = sample_size, mu = sample_meanvector, Sigma = S) %>% data.frame() %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{.col}_log10SPVL"))

lm(recipient ~ transmitter, data = dist) %>% summary()


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




prob_recip_multiple <- RunParallel(populationmodel_fixedVL_Environment, dist$transmitter_log10SPVL) %>%
  do.call(cbind.data.frame, .) %>% t() # Time difference of 44.60442 mins on Macbook

combined_data <- cbind.data.frame(dist, prob_recip_multiple)
head(combined_data)

fig_1b <- ggplot(
  combined_data , 
                 aes(x = transmitter_log10SPVL, 
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



