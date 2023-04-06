# This script plots extracted estimates/data from the Thompson transmission model to produce 
# plots describing the within-host processes encoded within the model

require(tidyverse)
require(RColorBrewer)

my_palette <- RColorBrewer::brewer.pal(name="OrRd",n=9)[4:9]

# For a given SpVL, plot the probabilities of transmission over time since infection
spvl_infdur <- read_csv('./data/spvl_infectionduration.csv')
plt_tda <- 
  ggplot(spvl_infdur)+
  geom_line(aes(x = t, y = p, colour = as.factor(spvl)), size = 1.5)+
  scale_color_manual(values = my_palette, name = 'SpVL')+
  my_theme+
  scale_x_continuous(expand = c(0,0), limits = c(0,17), 'Time') +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.04), 'P(Acquisition)')+
  theme(legend.position = c(0.92,0.88))


# Plot the proportion of the xth most common variants over time since infection
var_t <- read_csv('./data/variant_distribution_time.csv') %>% 
  mutate(x_var = as.factor(x_var)) %>% 
  rename(time = T)

plot_var <- ggplot(var_t, aes(x = time, y = P, fill = x_var)) + 
  geom_area() + 
  scale_x_continuous(name = 'Time', expand = c(0,0))+
  scale_y_continuous(name = 'Proportion of Xth Most Common Variant', expand = c(0,0))+
  scale_fill_brewer(palette = 'OrRd', na.value = 'black')+
  my_theme+
  theme(legend.position = 'none')

panel_s1 <- cowplot::plot_grid( plot_var,plt_tda, nrow = 1, align = 'hv', labels = 'AUTO')

## END ##