# Under what conditions might we expect to observe an association between multiple variant infection and SpVL
#
# This script follows the linear equations set out by Bonhoeffer et al. 2015 that describes the 
# changes of the distribution of set point viral load in the HIV carrier population over a full 
# cycle of transmission
# https://doi.org/10.1371/journal.ppat.1004634
#
# We extend the linear equations to decompose the recipient distribution of SpVL into distributions for
# single and multple variant infections, leveraging the effect sizes and means calculated by Janes et al. 
# 2015 
# https://doi.org/10.1038/nm.3932 


################################### Dependencies ###################################

BonhoefferEqns <- function(MC, VC){
  
  # Variables
  # MC <- Carrier phenotype
  # VC <- Carrier variance
  u_0 <- 4.46  # Gaussian approximation of Transmission Potential
  v_0 <- 0.96  # Gaussian approximation of Transmission Potential
  mu_e <- 3  # Distribution of environmental effects 
  v_e <- 1  # Distribution of environmental effects 
  v_t <- 0.15 # variance asssociated with sampling at transmission
  
  # Subtract environmental effects to calculate genotypic distribution
  m_c <- MC - mu_e
  v_c <- VC - v_e
  
  # Apply transmission potential to infer donor genotype
  m_d <- (m_c*(v_e + v_0) + (mu_0 - mu_e)*v_c)/(v_c + v_e + v_0)
  v_d <- (v_c*(v_e + v_0))/(v_c + v_e + v_0)
  
  # Infer donor phenotype
  MD <- (MC*v_0 + mu_0*VC)/(v_0 + VC)
  VD <- (VC*v_0)/(v_0 + VC)
  
  # sample from donor genotypes for recipient genotypes
  m_r <- m_d 
  v_r <- v_d + v_t
  
  # Re-distribute environmental variance for recipient phenotype
  MR <- m_r + mu_e
  VR <- v_r + v_e
  
  out <- c('mean' = MR,
           'variance' = VR)
  
  return(out)
}


################################### Calculate Recipient Dist ###################################
zambia_variance <- 0.61 
zambia_mean <- 4.74 

recipient_dist <- BonhoefferEqns(zambia_mean, 
                                 zambia_variance)


################################### Decompose Recipient Distribution ###################################




################################### Simulate proportion of significant observations ###################################




# Simulate for different effect sizes and p(MVs)
effect_size <- seq(0.01, 1, by = 0.01)
p_mv <- seq(0.01, 0.6, by = 0.01)

test <- matrix(data = NA, nrow = 100, ncol = 60)


for (i in 1:length(effect_size)){
  
  MR_sv <- m_d + mu_e
  MR_mv <- m_d + effect_size[i] + mu_e 
  
  for(j in 1:length(p_mv)){
    n <- 200
    sig.count <- 0

    for (z in 1:100){
      sv <- rnorm(n*(1-p_mv[j]), mean = MR_sv, sd = sqrt(VR))
      mv <- rnorm(n*p_mv[j], mean = MR_mv, sd = sqrt(VR))
      
      if(t.test(sv, mv)$p.value <= 0.05){
        sig.count = sig.count + 1
      }
      
      test[i,j] <- sig.count/100
    }
   
  }
}

test_tibble <- test %>%
  reshape2::melt() %>%
  mutate(p_mv = rep(p_mv, each = 100)) %>%
  mutate(effect_size = rep(effect_size,  60)) %>%
  select(-contains('Var'))


################################### Export to file ###################################

plt1a <- ggplot(test_tibble) +
  geom_tile(aes(x = p_mv, y = effect_size, fill = value))+    
  geom_point(x = 0.32, y = 0.293, shape = 4, colour = 'white') + 
  annotate("text", label = "Janes et al. RV144",
    x = 0.4, y = 0.32, size = 8, colour = "white"
  )+
  geom_point(x = 0.25, y = 0.372, shape = 4, colour = 'white') + 
  annotate("text", label = "Janes et al. STEP",
           x = 0.32, y = 0.4, size = 8, colour = "white"
  )+
  geom_vline(xintercept = 0.21, colour = "white", linetype = 2) + 
  annotate("text", label = "MF",
           x = 0.23, y = 0.95, size = 8, colour = "white"
  )+
  geom_vline(xintercept = 0.13, colour = "white", linetype = 2)+
  annotate("text", label = "FM",
           x = 0.15, y = 0.95, size = 8, colour = "white"
  )+
  geom_vline(xintercept = 0.3, colour = "white", linetype = 2) + 
  annotate("text", label = "MM",
           x = 0.32, y = 0.95, size = 8, colour = "white"
  )+
  scale_y_continuous('SpVL Effect Size of Multiple Variants', expand= c(0,0)) +
  scale_x_continuous('Cohort P(Multiple Variants)', expand= c(0,0))+
  scale_fill_distiller(palette = 'OrRd', 'P(SpVL ~ P(MV) is Significant)', direction = 1) +
  my_theme 


vls <- tibble(recipient_mean = rnorm(10000, MR, sqrt(VR)), 
              recipient_mv = rnorm(10000, MR+0.3, sqrt(VR)),
              recipient_sv = rnorm(10000, MR-0.7, sqrt(VR))) %>%
  pivot_longer(cols = everything(), names_to = 'stage', values_to = 'vl')

plt1b <- ggplot(vls ) + 
  geom_density(aes(x = vl ,fill = stage), alpha = 0.5)+
  scale_x_continuous(expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')), expand= c(0,0))+
  scale_y_continuous('Density', expand= c(0,0))+
  scale_fill_brewer(palette = 'OrRd')+
  geom_vline(xintercept = MR+0.3,  linetype = 2)+
  geom_vline(xintercept = MR-0.7,  linetype = 2)+
  geom_vline(xintercept = MR, linetype = 2)+
  facet_wrap(.~stage, nrow = 3)+
  my_theme

cowplot::plot_grid(plt1b, plt1a, align = 'hv', nrow = 1, labels = 'AUTO')
