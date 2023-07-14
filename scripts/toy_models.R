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


# Assume recipient phenotype distribution is a mixture of two gaussian components where
# the difference in means is set by an empirically determined effect size
DecomposeRecipientSpVL <- function(recipient_mean, recipient_var, p_mv, effect_size){
  
  # recipient_mean = (1 - p_mv)*x + p_mv*(x + effect_size)
  # Because weights of finite mixture must sum to 1; rearranges to:
  sv_mean <- recipient_mean - p_mv*effect_size  
  mv_mean <- sv_mean + effect_size
  
  
  #recipient_var <- p_m*(y + mv_mean^2) + (1-p_mv)*(y + sv_mean^2) - recipient_mean^2
  # Rearranges for the same reasons
  # Assumption of equal variances justified by Janes et al.
  decomposed_var <- recipient_var + recipient_mean^2 - p_mv * mv_mean^2 - (1 - p_mv) * sv_mean^2
  
  out <- list('sv' = c('mean' = sv_mean, 'var' = decomposed_var),
              'mv' = c('mean' = mv_mean, 'var' = decomposed_var))
  
  return(out)
}



SimEffectSizePMV <- function(n, e = effect_size, p = p_mv){
  # Change P_mv to fraction of n
  m <- matrix(data = NA, nrow = 100, ncol = 60)
  for (i in 1:length(e)){
    for(j in 1:length(p)){
      spvls <- DecomposeRecipientSpVL(effect_size = e[i],
                                      p_mv = p[j],
                                      recipient_mean = 4.74,
                                      recipient_var = 0.61)
      sig.count <- 0
      for (z in 1:100){
        sv <- rnorm(n*(1-p[j]), mean = spvls$sv['mean'], sd = sqrt(spvls$sv['var']))
        mv <- rnorm(n*p[j], mean = spvls$mv['mean'], sd = sqrt(spvls$mv['var']))
        if(t.test(sv, mv)$p.value <= 0.05){
          sig.count = sig.count + 1
        }
        m[i,j] <- sig.count/100
      }
    }
  }
  out <- m %>%
    reshape2::melt() %>%
    mutate(p_mv = rep(p_mv, each = 100)) %>%
    mutate(effect_size = rep(effect_size,  60)) %>%
    select(-contains('Var')) %>%
    mutate(sample_size = n)
  return(out)
}


################################### Calculate Recipient Dist ###################################
zambia_variance <- 0.61 
zambia_mean <- 4.74 

recipient_dist <- BonhoefferEqns(zambia_mean, 
                                 zambia_variance)


################################### Decompose Recipient Distribution ###################################
decomp_vl <- DecomposeRecipientSpVL(recipient_mean = 4.74, 
                                    recipient_var = 0.61,
                                    p_mv =  0.3, 
                                    effect_size = 0.3)


vls <- tibble(recipient_mean = rnorm(100000, 4.74, sqrt(0.61)), 
              recipient_mv = rnorm(100000, decomp_vl$mv['mean'], sqrt(decomp_vl$mv['var'])),
              recipient_sv = rnorm(100000, decomp_vl$sv['mean'], sqrt(decomp_vl$sv['var']))) %>%
  pivot_longer(cols = everything(), names_to = 'stage', values_to = 'vl')

facet.labs <- c("Recipient Population", "Multiple Variant Recipients", "Single Variant Recipients")
names(facet.labs) <- c('recipient_mean', 'recipient_mv', 'recipient_sv')

plt1b <- ggplot(vls ) + 
  geom_density(aes(x = vl ,fill = stage), alpha = 0.5, colour = 'white')+
  scale_x_continuous(expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')')), expand= c(0,0))+
  scale_y_continuous('Density', expand= c(0,0))+
  scale_fill_brewer(palette = 'OrRd')+
  geom_vline(xintercept = decomp_vl$mv['mean'],  linetype = 2, colour = '#fdbb84')+
  geom_vline(xintercept = 4.74 ,  linetype = 2, colour = '#fee8c8')+
  geom_vline(xintercept = decomp_vl$sv['mean'], linetype = 2, colour = '#e34a33')+
  facet_wrap(.~stage, nrow = 3,labeller = labeller(stage = facet.labs))+
  my_theme






################################### Simulate proportion of significant observations ###################################

# Vectors of effect size and probability that infection is initiated by multiple variants
effect_size <- seq(0.01, 1, by = 0.01)
p_mv <- seq(0.02, 0.6, by = 0.01)

cl <- detectCores() %>% `-` (3) 

simsignificances <- mclapply(c(25,50,100,200), SimEffectSizePMV, 
                             mc.cores = cl,
                             mc.set.seed = FALSE) %>%
  bind_rows()

simsignificances <- SimEffectSizePMV(50)
################################### Export to file ###################################

plt1a <- ggplot(simsignificances ) +
  geom_tile(aes(x = p_mv, y = effect_size, fill = value))+    
  geom_point(x = 0.32, y = 0.293, shape = 2, colour = 'white') + 
  
  geom_point(x = 0.25, y = 0.372, shape = 5, colour = 'white') + 
 
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
  facet_wrap(.~sample_size) + 
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
