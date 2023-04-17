# Plotting CD4 ~ Grouped VL

mellors_data <- read_csv('./data/mellors_data.csv') %>% 
  mutate(spvl = factor(spvl, levels = c('≤500', '501-3000', '3001-10000', '10001-30000', '>30000')))

mellors_plot <-  ggplot(mellors_data, aes(x = spvl, y = est))+
  geom_point(size = 3, colour = '#ef654a')+
  geom_linerange(aes(ymin = ci.lb, ymax = ci.ub, x= spvl), colour = '#ef654a') + 
  scale_x_discrete('Binned SpVL')+
  scale_y_continuous(limits = c(-90,-20), expand = c(0,0), name =expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', year**-1)))+
  my_theme +
  theme(axis.title.y = element_text(family = 'sans'))

ggsave(plot = mellors_plot, filename = paste(figs_dir,sep = '/', "mellors_plot.jpeg"), device = jpeg, width = 10, height = 10)


x <- ggplot()+
  geom_histogram(aes(x = sim_serocons , y = after_stat(count / sum(count))), fill = '#ef654a')+
  scale_y_continuous(expand = c(0,0), limit = c(0, 0.07), breaks = seq(0, 0.07, by = 0.01), name = 'Frequency')+
  scale_x_log10(name = expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.02,0.02),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x)))+
  my_theme 

ppt_panelx <- plot_grid(x, mellors_plot,  labels = 'AUTO', align = 'hv', nrow = 1)

ggsave(plot = ppt_panelx , filename = paste(figs_dir,sep = '/', "ppt_panelx .jpeg"), device = jpeg, width = 18, height = 10)
