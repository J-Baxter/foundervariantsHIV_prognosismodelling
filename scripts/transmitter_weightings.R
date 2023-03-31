# Simulating Model Pop (Transmitters)

# Duration of Infection
D <- function(viralload){
  Dmax = 25.4
  Dfifty = 3058
  Dk = 0.41
  
  out <- Dmax*(Dfifty**Dk)/(viralload**Dk + Dfifty**Dk)
  
  
  return(out)
  
}


# The fraction of individuals at seroconversion (the point at which HIV-1 antibodies develop and become
# detectable) with each log(SPVL), vc, in the population at any given time is described by 'g'
WeightSerocons <- function(viralload){
  alpha = -3.55
  sigma <- 0.78/(sqrt(1 - ((2*alpha**2)/(pi*(1 + alpha**2)))))
  mu <-  4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha**2)))
  weight <-  (2/sigma)*dnorm((log10(viralload) - mu)/sigma)*pnorm(alpha*(log10(viralload) - mu)/sigma)
  return(weight)
} 


# proportion of all infected individuals that have each SPVL when they are in the chronic phase
# of infection (including individuals currently in the primary or pre-AIDS phases), g(V). This distribution is
# obtained by weighting gv(vc) by the relative lengths of the infectious periods of individuals with each SPVL 
# and renormalising the resulting distribution.
WeightInfecteds <- function(viralload){
  taup = 0.24
  taua = 0.75
  tauc = D(viralload)
  
  g = WeightSerocons(viralload)
  
  numIndividualTypes = length(g)
  totalOfTimeVals = numIndividualTypes*(taup + taua) + sum(tauc)
  
  g = g * ((taup + tauc + taua)/totalOfTimeVals)
  g = g/sum(g)
  
  return(g)
}


# Infectiousness
B<- function(viralload){
  Bmax = 0.317
  Bfifty = 13938
  Bk = 1.02
  
  out <- (Bmax*viralload**Bk)/(viralload**Bk + Bfifty**Bk)
  
  return(out)
}


# Probability of observing at least one infection for a given viral load
Obs <- function(viralload){
  tauc = D(viralload)
  
  tp <- B(viralload) * tauc
  
  p_0 <- ((tp*d)**0*exp(-tp*d))/factorial(tp)
  
  return(1-p_0)
  
}

# Infer probs

# Starting Vector 
spvl_vec <- seq(2, 7, length.out = 10000) %>% raise_to_power(10,.)

p_spvl <- sapply(spvl_vec, function(x) WeightSerocons(x)/sum(WeightSerocons(spvl_vec)))

p_obs <- Obs(spvl_vec)

# Visualise

cols <- c("Seroconverter" = "#fdd49e", 'Transmission Potential' = '#fc8d59', "P(Observed)" = "#b30000")
a <- ggplot()+
  geom_line(aes(x = spvl_vec, y = p_spvl, colour = 'Seroconverter'), linejoin = "round",  size = 2)+
  scale_colour_manual(values = cols)+
  scale_y_continuous(name = 'Proportion Seroconverter') + 
  scale_x_log10(name = expression(paste("SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme+ theme(legend.position = 'bottom')


b <- ggplot()+
  geom_line(aes(x = spvl_vec , y = p_obs, colour = 'P(Observed)'), linejoin = "round",  size = 2)+
  scale_y_continuous(name = 'P(Observed)') + 
  scale_colour_manual(values = cols)+
  scale_x_log10(name = expression(paste("SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme + theme(legend.position = 'bottom')

c <- ggplot()+
  geom_line(aes(x = spvl_vec , y = B(spvl_vec)*D(spvl_vec), colour = 'Transmission Potential'), linejoin = "round",  size = 2)+
  scale_y_continuous(name = 'Transmission Potential') + 
  scale_colour_manual(values = cols)+
  scale_x_log10(name = expression(paste("SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme + theme(legend.position = 'bottom')

panel_S2a <- plot_grid(a+theme(legend.position = 'none'), c+theme(legend.position = 'none'),b+theme(legend.position = 'none'), nrow = 1, align = 'h')



sim_serocons <- sample(spvl_vec, size = 5000, prob = p_spvl, replace = F) 
sim_observed <- sample(spvl_vec, size = 5000, prob = p_obs, replace = F) 
sim_both <- sample(spvl_vec, size = 5000, prob = p_spvl*p_obs, replace = F) 

v <- ggplot()+
  geom_histogram(aes(x = spvl_vec, y = after_stat(count / sum(count))))+
  scale_x_log10(name = expression(paste("SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme + ggtitle('Initial SPVL Vector')

x <- ggplot()+
  geom_histogram(aes(x = sim_serocons , y = after_stat(count / sum(count))), fill = '#ef654a')+
  scale_y_continuous(expand = c(0,0), limit = c(0, 0.07), breaks = seq(0, 0.07, by = 0.01))+
  scale_x_log10(name = expression(paste("SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme + ggtitle(str_wrap('Weighted by Proportion of Seroconverters at each SPVL (G)', 30))

y <- ggplot()+
  geom_histogram(aes(x = sim_observed , y = after_stat(count / sum(count))), fill = '#ef654a')+
  scale_y_continuous(expand = c(0,0), limit = c(0, 0.07), breaks = seq(0, 0.07, by = 0.01))+
  scale_x_log10(name = expression(paste("SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme + ggtitle(str_wrap('Weighted by Probability of Observing Transmission', 30))

z <- ggplot()+
  geom_histogram(aes(x = sim_both, y = after_stat(count / sum(count))), fill = '#ef654a')+
  scale_y_continuous(expand = c(0,0), limit = c(0, 0.07), breaks = seq(0, 0.07, by = 0.01))+
  scale_x_log10(name = expression(paste("SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme + ggtitle(str_wrap('Weighted by product of both distributions',30))



panel_S2b <- plot_grid(plot_grid(x,y,z,nrow = 1))
panel_S2 <- plot_grid(panel_S2a,NULL, panel_S2b, nrow = 3, rel_heights= c(1,0.05,1))
ggsave("panel_S2.jpeg", device =  jpeg , plot =panel_S2 , width = 16, height = 16)