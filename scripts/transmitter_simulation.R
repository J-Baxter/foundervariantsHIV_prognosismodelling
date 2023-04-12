# Simulating Model Pop (Transmitters)
source('./scripts/SimPop.R')

################## Simulate transmitter population stage-by-stage ################## 
n = 1000
SpVL_vec <- SpVL_vec <- seq(2, 7, length.out = 10000) %>% raise_to_power(10,.)

p_SpVL <- sapply(SpVL_vec, function(x) WeightSerocons(x)/sum(WeightSerocons(SpVL_vec)))

p_obs <- Obs(SpVL_vec)

sim_SpVL <- sample(SpVL_vec, size = n, prob = p_SpVL*p_obs, replace = F) 


############### Plot probability distirbution at each stage ###################
cols <- c("Seroconverter" = "#fdd49e", 'Transmission Potential' = '#fc8d59', "P(Observed)" = "#b30000")
a <- ggplot()+
  geom_line(aes(x = SpVL_vec, y = p_SpVL, colour = 'Seroconverter'), linejoin = "round",  size = 2)+
  scale_colour_manual(values = cols)+
  scale_y_continuous(name = 'Proportion Seroconverter') + 
  scale_x_log10(name = expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme+ theme(legend.position = 'bottom')


b <- ggplot()+
  geom_line(aes(x = SpVL_vec , y = p_obs, colour = 'P(Observed)'), linejoin = "round",  size = 2)+
  scale_y_continuous(name = 'P(Observed)') + 
  scale_colour_manual(values = cols)+
  scale_x_log10(name = expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme + theme(legend.position = 'bottom')

c <- ggplot()+
  geom_line(aes(x = SpVL_vec , y = B(SpVL_vec)*D(SpVL_vec), colour = 'Transmission Potential'), linejoin = "round",  size = 2)+
  scale_y_continuous(name = 'Transmission Potential') + 
  scale_colour_manual(values = cols)+
  scale_x_log10(name = expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme + theme(legend.position = 'bottom')

panel_s2a <- plot_grid(a+theme(legend.position = 'none'), c+theme(legend.position = 'none'),b+theme(legend.position = 'none'), nrow = 1, align = 'h')


############### Plot intermediate and final histograms ####################
sim_serocons <- sample(SpVL_vec, size = 5000, prob = p_SpVL, replace = F) 
sim_observed <- sample(SpVL_vec, size = 5000, prob = p_obs, replace = F) 
sim_both <- sample(SpVL_vec, size = 5000, prob = p_SpVL*p_obs, replace = F) 

v <- ggplot()+
  geom_histogram(aes(x = SpVL_vec, y = after_stat(count / sum(count))))+
  scale_x_log10(name = expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme + ggtitle('Initial SpVL Vector')

x <- ggplot()+
  geom_histogram(aes(x = sim_serocons , y = after_stat(count / sum(count))), fill = '#ef654a')+
  scale_y_continuous(expand = c(0,0), limit = c(0, 0.07), breaks = seq(0, 0.07, by = 0.01), name = 'Frequency')+
  scale_x_log10(name = expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme + ggtitle(str_wrap('Weighted by Proportion of Seroconverters at each SpVL (G)', 30))

y <- ggplot()+
  geom_histogram(aes(x = sim_observed , y = after_stat(count / sum(count))), fill = '#ef654a')+
  scale_y_continuous(expand = c(0,0), limit = c(0, 0.07), breaks = seq(0, 0.07, by = 0.01), name = 'Frequency')+
  scale_x_log10(name = expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme + ggtitle(str_wrap('Weighted by Probability of Observing Transmission', 30))

z <- ggplot()+
  geom_histogram(aes(x = sim_both, y = after_stat(count / sum(count))), fill = '#ef654a')+
  scale_y_continuous(expand = c(0,0), limit = c(0, 0.07), breaks = seq(0, 0.07, by = 0.01), name = 'Frequency')+
  scale_x_log10(name = expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme + ggtitle(str_wrap('Weighted by product of both distributions',30))

panel_s2b <- plot_grid(plot_grid(x,y,z,nrow = 1))
panel_s2 <- plot_grid(panel_s2a,NULL, panel_s2b, nrow = 3, rel_heights= c(1,0.05,1))
## END ##