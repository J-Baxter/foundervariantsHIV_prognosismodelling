# Probability of Acquisition plots
library(tidyverse)
library(showtext)


p_aq <- data.frame(route = c('FM', 'MF', 'MSM (IA)', 'MSM (RA)', 'PWID', 'MTC'),
                   p_aq.ci.lower = c(1/10000, 6/10000, 4/10000, 102/10000, 41/10000, 1700/10000),
                    p_aq = c(4/10000, 8/10000, 11/10000, 138/10000, 63/10000, 2260/10000),
                   p_aq.ci.upper = c(14/10000, 11/10000, 28/10000, 186/10000, 92/10000, 2900/10000)) %>% 
  mutate(q = 1-p_aq) %>%
  pivot_longer(!route, names_to = 'type', values_to ='probability') %>%
  mutate(route = factor(route, levels = c('FM', 'MF', 'MSM (IA)', 'PWID', 'MSM (RA)','MTC')))

font_add_google("Questrial", "Questrial")
showtext_auto()

plt1 <- ggplot() + 
  geom_bar(aes(fill = type, y = probability, x = route), position = 'fill', stat = 'identity', data = p_aq %>% filter(type %in% c('p', 'q')) %>% mutate(type = factor(type, levels = c('q', 'p')))) +
  scale_x_discrete(name = 'Exposure', expand = c(0.1,0.1)) + 
  scale_y_continuous(name = 'Probability of Aquisition', breaks = seq(0,1,by = 0.1), expand = c(0.01,0.01)) + 
  scale_fill_discrete(type = c('#fee6ce', '#e84c21'))+
  coord_flip()+
  theme_minimal(base_family = "Questrial")+
  theme(
    panel.grid = element_blank(),
    legend.position = 'none',
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))
  )
  

setEPS()
postscript("probability_acqusition.eps", width = 16, height = 9)
plt1 
dev.off()

#1-novariants transmitted
#zerotruncated Poisso
# Calculate p(MV) for a given p(acquisition)

#lambda = p(acq)

# frequency of transmissions of one, two, or more than two variants
# of 0.1 (Poisson mean, 0.1), 0.25 (Poisson mean, 0.5), and 0.4 (Poisson mean, 0.5). 

library(VGAM)
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


data <- left_join(julian_tp_data, aba5443_data_s2) 

m1 <- vglm(variants.mode ~ 1, family = pospoisson(), data = data)

abrams <- data.frame(variants.mode = c(3,3,5,2,3,2,2,3,3,3,3,3,3,2,2, rep(1,54)))

pppos <- function(k,gamma){
  prob <-  (gamma**k)/((exp(gamma)-1)*factorial(k))
  return(prob)
}

p_1 <- c(pppos(1, 0.001), pppos(1, 0.01), pppos(1, 0.1), pppos(1, 0.2), pppos(1, 0.4))
p_2 <- c(pppos(2, 0.001), pppos(2, 0.01), pppos(2, 0.1), pppos(2, 0.2), pppos(2, 0.4))
p_grtr2 <- rep(1,5) - c(p_1 + p_2)

p_mv_df <- data.frame('p_mv' = c(p_1, p_2, p_grtr2), 'p_acq' = rep(c('0.001', '0.01','0.1', '0.2', '0.4'),3), variants = factor(c(rep('1', 5), rep('2', 5), rep('>2', 5)), levels = c('1', '2', '>2')) )

plt2 <- ggplot() + 
  geom_blank(data=p_mv_df, aes(x = variants, y = p_mv)) +
  geom_bar(data=p_mv_df, aes(x = variants, y = p_mv, fill = p_acq), stat = 'identity', position = position_dodge()) + 
  geom_ribbon(aes(ymin = 0.21, ymax = 0.29, x = 0:4), fill = 'black', alpha = 0.1) +
  geom_abline(intercept = 0.25, slope = 0, linetype = 'longdash', size = 0.7)+
  scale_x_discrete(name = 'Number of Variants', expand = c(0,0)) + 
  scale_y_continuous(name = 'Frequency of Infection', expand = c(0,0))+
  scale_fill_brewer(type = seq, palette = 'OrRd', name = 'P(Acquisition)') + 
  theme_classic(base_family = "Questrial")+
  theme(
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )

setEPS()
postscript("probability_numbervariants.eps", width = 16, height = 9)
plt2
dev.off()


# Probabilities MV

p_mv_sysreview <- data.frame('exposure' = factor(c('HSX:mf', 'HSX:fm', 'HSX:unknown', 'MSM', 'MTC:pre', 'MTC:intra', 'MTC:post', 'MTC:unknown', 'PWID'), 
                                                 levels = c('HSX:mf', 'HSX:fm', 'HSX:unknown', 'MSM', 'MTC:pre', 'MTC:intra', 'MTC:post', 'MTC:unknown', 'PWID')),
                             p_mv.ci.lower = c(0.14, 0.08, 0.12, 0.22, 0.08, 0.14, 0.03, 0.08, 0.24),
                             p_mv = c(0.21, 0.13, 0.32, 0.30, 0.17, 0.27, 0.18, 0.18, 0.37),
                             p_mv.ci.upper = c(0.31, 0.21, 0.61, 0.4, 0.33, 0.45, 0.57, 0.35, 0.53))


aq_mv_df <- filter(p_aq, type != 'q') %>%
  pivot_wider(names_from = type, values_from = probability) %>%
  mutate(exposure = c('HSX:fm', 'HSX:mf', 'MSM', 'MSM', 'PWID', 'MTC:unknown')) %>%
  right_join(p_mv_sysreview, by = 'exposure') %>%
  filter(!is.na(route)) 


plt3 <- ggplot(aq_mv_df) +
  geom_point(aes(x = p_aq, y = p_mv, color = exposure))+
  geom_linerange(aes(xmin = p_aq.ci.lower, xmax = p_aq.ci.upper, y = p_mv, color = exposure ))+
  geom_linerange(aes(ymin = p_mv.ci.lower, ymax = p_mv.ci.upper, x = p_aq, color = exposure)) +
  scale_x_log10(name = 'Probability of Acquisition', expand = c(0,0), limits = c(0.0001, 1), breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", label_math())) + 
  scale_y_continuous(name = 'Probability of Multiple Founders', expand = c(0,0), limits = c(0,.6)) +
  #scale_color_brewer(name = 'Exposure', palette = 'OrRd')+
  scale_colour_manual(values = c('#fdd49e','#fdbb84','#fc8d59','#e34a33','#b30000')) +
  annotation_logticks()+
  theme_classic(base_family = "Questrial")+
  theme(legend.position = 'right',
        text = element_text(size=20),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )


setEPS()
postscript("probability_comparison.eps", width = 14, height = 9)
plt3
dev.off()
