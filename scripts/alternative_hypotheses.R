# For a given TP (paired SPVL) 
# Probability distributions of multiple variants at fixed viral loads
# PerVirionProbability = 4.715*1e-8, PropExposuresInfective = 0.029

spvl <- seq(2,7,by=0.5)
prob_dists <- lapply(spvl, populationmodel_fixedVL_Environment) %>% sapply(., function(x) x[2:length(x)])
colnames(prob_dists) <- spvl
rownames(prob_dists) <- 1:33

prob_dists_plot <- cbind.data.frame(prob_dists) %>% 
  pivot_longer(cols = everything(), names_to = 'spvl', values_to = 'p') %>% 
  cbind.data.frame('variants' = rep(1:33, each = 11))

my_breaks = c(-1, -2, -4, -6, -8, -10, -12)

ggplot(prob_dists_plot) +
  geom_tile(aes(x = spvl, y = variants, fill = p))+
  scale_fill_gradient(low = '#990000' , high = '#fff7ec', trans = "log",  na.value = '#7f0000',
                      breaks = my_breaks, labels = my_breaks)+
  theme_classic()+
  theme(
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )

FoundingVars <- function(spvl, sd, p_dists){
  
  #Round Log10 SPVL to nearest .5
  spvl_disc <- round(spvl/0.5)*0.5
  
  #Identify appropriate probability distribution
  lookup <- which(spvl_disc == seq(2,7,by=0.5))
  prob_var <- p_dists[,lookup]
  
  n_var <- sample(x = 1:33, 1, replace = T, prob = prob_var)
  
  within_host_genotypes <- rnorm(n = n_var, mean = spvl, sd = sd)
  
  m <- matrix(NA, nrow = 1, ncol = 33)
  
  m[1:length(within_host_genotypes)] <- within_host_genotypes
  
  
  return(m)
  
}

