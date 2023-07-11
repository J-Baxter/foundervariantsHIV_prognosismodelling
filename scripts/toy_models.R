# Simulating Model Pop (Transmitters)
source('./scripts/SimPop.R')

##########################
# Gaussian approximation of Transmission Potential
u_0 <- 4.46 #TP 
v_0 <- 0.96 #TP

# Distribution of environmental effects 
mu_e <- 3 #env
v_e <- 1 #env

# Change of the viral genotype due to intrahost evolution 
mu_i <- 0.2 # genetic change in the virus that affects log spVL across all patients in the same wa
v_i <- 0.15 # genetic changes that affect log spVL in a manner that is specific to the patient

# variance asssociated with sampling at transmission
v_t <- 0.15

# SpVL distribution from Zambia Cohort
VC <- 0.61 
MC <- 4.74 

# Subtract environmental effects
m_c <- MC - mu_e
v_c <- VC - v_e

# Infer donor genotype
m_d <- (m_c*(v_e + v_0) + (mu_0 - mu_e)*v_c)/(v_c + v_e + v_0)
v_d <- (v_c*(v_e + v_0))/(v_c + v_e + v_0)

# Infer donor phenotype
MD <- (MC*v_0 + mu_0*VC)/(v_0 + VC)
VD <- (VC*v_0)/(v_0 + VC)

# sample from donor genotypes for recipient genotypes
m_r <- m_d #(m_d)*0.7 + (m_d+0.3)*0.3
v_r <- v_d + v_t

# Re-distribute environmental variance for recipient phenotype
MR <- m_r + mu_e
VR <- v_r + v_e       

n <- 200
vls <- tibble(donor = rnorm(n, MD, sqrt(VD)), 
              carrier = rnorm(n, MC, VC),
              recipient = rnorm(n, MR, VR)) %>%
  pivot_longer(cols = everything(), names_to = 'stage', values_to = 'vl')

ggplot(vls ) + 
  geom_density(aes(x = vl ,fill = stage), alpha = 0.5)



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

ggplot(test_tibble) +
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
  my_theme + theme(legend.position = 'bottom')

