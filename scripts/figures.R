# Figures

# Dependencies
source('./scripts/dependencies.R')
font_add_google("Questrial", "Questrial")
showtext_auto()


# Import data
pop <- read.csv('./results/08/Feb23/base_population.csv') # Make sure you change the date!!!
transmitterspvl_variantdist <- read.csv('./results/08/Feb23/transmitterspvl_variantdist.csv')
recipientspvl_variantdist <- read.csv('./results/08/Feb23/recipientspvl_variantdist.csv')

# Transmitter SPVL ~ Recipient SPVL (include distributions on axis?)
plt1 <- ggplot(pop, aes(x = transmitter, recipient)) +
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


# Null
fig_1b <- ggplot(transmitterspvl_variantdist, 
                 aes(x = transmitter, 
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


fig_1c <- ggplot(sim_recipientspvl_variantdist, 
                 aes(x = recipient,
                     y = 1 - variant_distribution.V1))+
  geom_point(colour = '#2ca25f',  #usher colours?
             alpha = 0.5)+
  scale_x_log10(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0,0),
                     limits = c(0,0.6),
                     breaks = seq(0, 0.6, by = 0.2)) +
  theme_classic() 



# Alternative
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