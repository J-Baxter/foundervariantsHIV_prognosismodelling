# Graphical Abstract/Cartoon Plots

# Alternative Hypotheses
data <- expand_grid(v1 = seq(4,6,by=0.2), v2 = seq(4,6,by=0.2)) %>%
  rowwise() %>% 
  mutate(Additive = Additive(5, c(v1, v2))) %>%
  mutate(Interaction = Interaction(5, c(v1, v2))) %>%
  mutate(Exclusion = Exclusion(5, c(v1, v2))) %>%
  pivot_longer(cols = c(Additive, Interaction, Exclusion), names_to = 'func', values_to = 'recip')


toy_models <- ggplot(data) +
  geom_raster(aes(x = v1, y = v2, fill= recip))+
  scale_fill_distiller('Recipient SPVL', palette = 'OrRd') +
  #scale_fill_viridis_c('Recipient SPVL', option = 'inferno') +
  theme_classic(base_family = "Questrial")+
  facet_wrap(~func, ncol = 1)+
  theme(
    legend.position = 'none',
    text = element_text(size=20)
  )
setEPS()
postscript( paste(figs_dir, 'toy_models.eps', sep = '/'), width = 8, height = 16)
toy_models
dev.off()


# Probability of the number of variants | donor viral load

spvl <- 10**seq(2,7,by=0.5)
prob_dists <- lapply(spvl, populationmodel_acrossVL_Environment) %>% sapply(., function(x) x[2:length(x)])
colnames(prob_dists) <- spvl
rownames(prob_dists) <- 1:33

prob_dists_plot <- cbind.data.frame(prob_dists) %>% 
  pivot_longer(cols = everything(), names_to = 'spvl', values_to = 'p') %>% 
  mutate(spvl = as.numeric(spvl))%>%
  cbind.data.frame('variants' = rep(1:33, each = 11))