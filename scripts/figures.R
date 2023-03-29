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

#################### Panel 1 is graphical abstract, so starting from panel 2 ####################

############################################## Panel 2 ##############################################
#3 model outputs (tolerance, herit, jp transmission)

############################################## Panel 3 ##############################################
#Confounder

############################################## Panel 4 ##############################################
#Timing

############################################## Panel 5 ##############################################
#Non-Linear



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
   + 
  annotation_logticks() +
  my_theme + theme(legend.position = 'none')

setEPS()
postscript( paste(figs_dir, 'population.eps', sep = '/'), width = 10, height = 10)
plt1 
dev.off()


# Null Scenario (Probability of multiple variants ~ transmitter SPVL)
fig2_data <- transmitterspvl_variantdist %>% 
  filter(variants == 1) 
fig_2a <- ggplot(transmitterspvl_variantdist, 
                 aes(x = transmitter, 
                     y = 1 - p,
                     colour = stage))+
  geom_point(shape = 3, size = 4) +
  scale_colour_manual(values = c('#e34a33','#fdbb84')) + 
  scale_x_log10(name = expression(paste("Transmitter SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.6),
                     breaks = seq(0, 0.6, by = 0.1)) +
  annotation_logticks(sides = 'b') +
  my_theme + theme(legend.position = 'none')

setEPS()
postscript( paste(figs_dir, 'fig_2a .eps', sep = '/'), width = 10, height = 10)
fig_2a 
dev.off()

fig2b_data <- sim_recipientspvl_variantdist %>% 
  filter(variants == 1) 
fig_2b <- ggplot(fig2b_data , 
                 aes(x = recipient,
                     y = 1 - p,
                     colour = stage))+
  geom_point(shape = 3, size = 4) +
  scale_colour_manual(values = c('#e34a33','#fdbb84'))+
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
  my_theme+ theme(
                  legend.position = 'none')


fig2c_data <- sim_recipientspvl_variantdist %>% 
  filter(variants == 1) 
fig_2c <- ggplot(fig2c_data , 
                 aes(y = CD4_decline,
                     x = 1-p,
                     colour = stage,
                     group = stage))+
  geom_point(shape = 3, size = 4) +
  scale_colour_manual(values = c('#e34a33','#fdbb84'))+
  scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**1)),  #
                     expand = c(0,0),
                     limits = c(-2,2)) +
  scale_x_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0,0),
                     limits = c(0,0.6),
                     breaks = seq(0, 0.6, by = 0.1))+
  geom_smooth(method = 'lm', formula = 'y~x') +
  
  my_theme + theme(axis.title.y = element_text(family = 'sans'), 
                   legend.position = 'none')


grid_1 <- cowplot::plot_grid(fig_2a, fig_2b, fig_2c, align = 'hv', nrow = 1, labels = 'AUTO')
legend <- get_legend(
  fig_2a +theme(legend.position = "bottom")
)
grid_2 <- cowplot::plot_grid(grid_1,legend , nrow = 2, rel_heights = c(1, .1) )

setEPS()
postscript( paste(figs_dir, 'grid_2.eps', sep = '/'), width = 20, height = 10)
grid_2 
dev.off()


# Null Scenario (Variant distribution ~ transmitter SPVL)
transmitterspvl_variantdist_plot <- transmitterspvl_variantdist %>% 
  pivot_longer(cols = contains('V'), names_to = 'variants', values_to = 'p') %>% 
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric())


# geom_density_ridges(alpha=0.6, stat="binline", bins=20)
fig_2c <- ggplot(transmitterspvl_variantdist_plot, 
                 aes(x = transmitter, 
                     y = variants, height = p))+
  geom_ridgeline() +
  geom_point(colour = '#ef654a',  
             shape = 4, size = 3) +
  scale_x_log10(name = expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(1, 10**8),
                expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = 'Variants',
                     expand = c(0,0),
                     limits = c(0,0.6),
                     breaks = seq(0, 0.6, by = 0.1)) +
  annotation_logticks(sides = 'b') +
  my_theme





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

fig_3a= ggplot(nls_df) + 
  geom_density(aes(x = recipient_log10SPVL, fill = sim), alpha = 0.4)+
  my_theme +
  scale_fill_brewer(palette = 'OrRd')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                     limits = c(0,8), expand = c(0,0))+ theme(
                       legend.position = 'none')
fig3_data <- combined %>% 
  filter(variants == 1) %>%
  filter (stage == 'average')

fig_3b <- ggplot(fig3_data , 
                 aes(x = recipient,
                     y = 1 - p,
                     colour = model))+
  geom_point(shape = 3, size = 4) +
  scale_colour_brewer(palette = 'OrRd')+
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
  my_theme+ theme(
    legend.position = 'none')



fig_3c <- ggplot(fig3_data , 
                 aes(y = CD4_decline,
                     x = 1-p,
                     colour = model,
                     group = model))+
  geom_point(shape = 3, size = 4) +
  scale_colour_brewer(palette = 'OrRd')+
  scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)),  #
                     expand = c(0,0),
                     limits = c(-2,2)) +
  scale_x_continuous(name = 'P(Multiple Founder Recipient)',
                     expand = c(0,0),
                     limits = c(0,0.6),
                     breaks = seq(0, 0.6, by = 0.1))+
  geom_smooth(method = 'lm', formula = 'y~x') +
  
  my_theme + theme(axis.title.y = element_text(family = 'sans'), 
                   legend.position = 'none')



grid_3 <- cowplot::plot_grid(fig_3a, fig_3b, fig_3c, align = 'hv', nrow = 1, labels = 'AUTO')

legend <- get_legend(
  fig_3a +theme(legend.position = "bottom")
)
grid_4 <- cowplot::plot_grid(grid_3,legend , nrow = 2, rel_heights = c(1, .1) )

ggsave(paste0(figs_dir, '/test.eps'), device = cairo_ps , plot = grid_4, width = 16, height = 10)


donor_diversity <- cbind.data.frame(time = 1:300, v = rgamma() )