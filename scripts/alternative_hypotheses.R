# For a given TP (paired SPVL) 
# Probability distributions of multiple variants at fixed viral loads
# PerVirionProbability = 4.715*1e-8, PropExposuresInfective = 0.029

spvl <- 10**seq(2,7,by=0.5)
prob_dists <- lapply(spvl, populationmodel_acrossVL_Environment) %>% sapply(., function(x) x[2:length(x)])
colnames(prob_dists) <- spvl
rownames(prob_dists) <- 1:33

prob_dists_plot <- cbind.data.frame(prob_dists) %>% 
  pivot_longer(cols = everything(), names_to = 'spvl', values_to = 'p') %>% 
  mutate(spvl = as.numeric(spvl))%>%
  cbind.data.frame('variants' = rep(1:33, each = 11))

my_breaks = c(-1, -2, -4, -6, -8, -10, -12)

ggplot(prob_dists_plot[prob_dists_plot$variants<8,]) +
  geom_tile(aes(x = spvl, y = variants, fill = p))+
 
  #
  scale_fill_steps(high = '#e34a33',  low = '#FFFFFF', breaks = seq(0,0.7,by=0.05))+
  #scale_fill_gradient(high = '#e34a33'  , low = '#FFFFFF', na.value = '#FFFFFF', breaks = seq(0,0.8,by=0.1))+
  scale_x_log10(limits = c(2, 10**8),
                expand = c(0,0),
                name = expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
 theme


ggplot(prob_dists_plot[prob_dists_plot$variants<8,]) +
  geom_point(aes(x = spvl, y = variants, size = p))+
  scale_radius( range = c(0.1,15)) +
  scale_x_log10(limits = c(2, 10**8),
                expand = c(0,0),
                name = expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  scale_y_continuous(name = 'Variants', breaks = 1:12)+
  #scale_fill_gradient2(high= '#990000' , mid = '#e34a33',low = '#FFFFFF', trans = "log",  na.value = '#FFFFFF',
  #breaks = my_breaks, labels = my_breaks)+
  theme_classic()+
  theme(
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )

ggplot(aes(x = ))
