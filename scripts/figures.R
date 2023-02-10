# Figures

# Dependencies
source('./scripts/dependencies.R')
font_add_google("Questrial", "Questrial")
showtext_auto()


# Import data
pop <- read.csv('./results/08Feb23/base_population.csv') # Make sure you change the date!!!
transmitterspvl_variantdist <- read.csv('./results/08Feb23/transmitterspvl_variantdist.csv')
recipientspvl_variantdist <- read.csv('./results/08Feb23/recipientspvl_variantdist.csv')


# Set theme
my_theme <- theme_classic(base_family = "Questrial")+
  theme(
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )


# Transmitter SPVL ~ Recipient SPVL (include distributions on axis?)
plt1 <- ggplot(pop, aes(x = transmitter, recipient)) +
  geom_point( #'#CB6015' #'#66c2a4','#2ca25f','#006d2c'
    colour = '#ef654a',
    shape = 4, size = 3) +
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
  annotation_logticks() +
  theme

setEPS()
postscript( paste(figs_dir, 'population.eps', sep = '/'), width = 10, height = 10)
plt1 
dev.off()


# Null Scenario (Probability of multiple variants ~ transmitter SPVL)
fig_2a <- ggplot(transmitterspvl_variantdist, 
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
  annotation_logticks(sides = 'b') +
  theme


fig_2b <- ggplot(sim_recipientspvl_variantdist, 
                 aes(x = recipient,
                     y = 1 - variant_distribution.V1))+
  geom_point(colour = '#ef654a',  
             shape = 4, size = 3)+
  scale_x_log10(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0,0),
                     limits = c(0,0.6),
                     breaks = seq(0, 0.6, by = 0.1)) +
  annotation_logticks(sides = 'b') +
  theme


# Null Scenario (Variant distribution ~ transmitter SPVL)
transmitterspvl_variantdist_plot <- transmitterspvl_variantdist %>% 
  pivot_longer(cols = contains('V'), names_to = 'variants', values_to = 'p') %>% 
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric())


# geom_density_ridges(alpha=0.6, stat="binline", bins=20)
fig_2c <- ggplot(transmitterspvl_variantdist_plot, 
                 aes(x = transmitter, 
                     y = variants, height = p))+
  geom_ridgeline()
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
  annotation_logticks(sides = 'b') +
  theme




# Alternative
# Toy Plots
  
data <- expand_grid(v1 = seq(4,6,by=0.2), v2 = seq(4,6,by=0.2)) %>%
  rowwise() %>% 
  mutate(Additive = Additive(5, c(v1, v2))) %>%
  mutate(Interaction = Interaction(5, c(v1, v2))) %>%
  mutate(Exclusion = Exclusion(5, c(v1, v2))) %>%
  pivot_longer(cols = c(Additive, Interaction, Exclusion), names_to = 'func', values_to = 'recip')
  

toy_models <- ggplot(data) +
  geom_raster(aes(x = v1, y = v2, fill= recip))+
  scale_fill_viridis_c('Recipient SPVL', option = 'inferno') +
  theme_classic(base_family = "Questrial")+
  facet_wrap(~func, ncol = 1)+
  theme(
    legend.position = 'bottom',
    text = element_text(size=14)
  )
setEPS()
postscript( paste(figs_dir, 'toy_models.eps', sep = '/'), width = 8, height = 16)
toy_models
dev.off()


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