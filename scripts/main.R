# Main Analysis for "RECONCILING THE FOUNDER VARIANT MULTIPLICITY OF HIV-1 INFECTION WITH THE RATE
# OF CD4+ DECLINE"
# J Baxter

# This script combines the functions for the analysis of paired SPVL distributions, occuring in
# three main stages:

# 1. Initialisation of study population
# 2. Null hypothesis
# 3. Ecological hypotheses


################################### Dependencies ###################################
source('./scripts/dependencies.R')
source('./scripts/populationdata_acrossVL_models.R')
source('./scripts/init_population.R')
source('./scripts/ecologicalhypos.R')
source('./scripts/cd4_delcine.R')

# Set Seed
set.seed(4472)
options(scipen = 100) #options(scipen = 100, digits = 4)


################################### Initialise Population ###################################
# Distributions of SPVL for couples with strong support for transmission Hollingsworth et al, 2010
# Dummy population until we incorporate Rakai data

index <- c('mean'= 4.61, 'sd' = 0.63) 
secondary <- c('mean' = 4.60 , 'sd' = 0.85)

NPAIRS <- 50 #500

pop <- InitPop(N = NPAIRS, 
               H2 = 0.33,
               donor_vl = index, 
               recipient_vl = secondary) %>%
  cbind.data.frame(sex = sample(c('M', 'F'), nrow(.), replace = T)) %>%
  cbind.data.frame(age = sample(18:70, nrow(.), replace = T))


################################### Estimate Heritability ###################################
h2_model <- lm(recipient_log10SPVL ~ transmitter_log10SPVL, data = pop) # Building this H2 model is a key step BRMS?

h2_priors <- prior(normal(1, 2), nlpar = "b1") + prior(normal(0, 1), nlpar = "b2") + prior(normal(0, 2), nlpar = "b3") 

h2_fit <- brm(recipient_log10SPVL ~ transmitter_log10SPVL + sex + age,
              data = pop)

summary(h2_fit) #Check output makes sense + ESS
plot(h2_fit) # Check convergence
bayes_R2(h2_fit) #Estimate R2

model_fit <- pop %>%
    add_predicted_draws(h2_fit) %>%  # adding the posterior distribution
    ggplot(aes(x = transmitter_log10SPVL , y = recipient_log10SPVL )) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = pop, colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("Calidris canutus abundance\n") +  # latin name for red knot
    xlab("\nYear") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85))



pred.data = expand.grid(transmitter_log10SPVL = seq(min(pop$transmitter_log10SPVL), max(pop$transmitter_log10SPVL), length=20),
                        sex = c('M', 'F'),
                        age = seq(min(pop$age), max(pop$age), length=20))

pred.data <- pred.data %>% cbind.data.frame(., predict(h2_fit, newdata=pred.data))

pop%>%
  add_predicted_draws(h2_fit) %>%  # adding the posterior distribution
  ggplot(aes(x = transmitter, y = recipient)) +  
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                  alpha = 0.5, colour = "black") +
  geom_point(data = France, colour = "darkseagreen4", size = 3)

ggplot()+
  geom_point(data = pop, aes(x = transmitter, recipient), #'#CB6015' #'#66c2a4','#2ca25f','#006d2c'
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
  ggplot(pred.data, aes(y = Estimate, x = transmitter_log10SPVL))+
  geom_line()
  
  annotation_logticks() +
  my_theme

################################### Initialise Simulated Populations ###################################
transmitter_prob <- sapply(pop$transmitter_log10SPVL, function(x) WeightPDF(x)/sum(WeightPDF(pop$transmitter_log10SPVL)))

hist(transmitter_prob) #visual check - should look normal(ish) as already log

sim_transmitter_log10SPVL <- InitSimTransmitter(pop_size = NPAIRS,
                                                transmitter_min = 1, 
                                                transmitter_max = 7, 
                                                sample_prob = transmitter_prob)

hist(sim_transmitter_log10SPVL) #visual check - should look uniform

sim_spvl <- predict(h2_model, newdata= data.frame(transmitter_log10SPVL = sim_transmitter_log10SPVL)) %>% 
  cbind.data.frame(sim_recipient_log10SPVL = ., sim_transmitter_log10SPVL = sim_transmitter_log10SPVL) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SPVL')}")) %>%
  `colnames<-` (str_remove(colnames(.), 'sim_')) %>% mutate(sim = 'linear')


################################### P(Multiple Variants) | Transmitter SPVL (H0) ###################################
# Probability that recipient infection is initiated by multiple founder variants
# applies populationmodel_acrossVL_Environment function written by Katie Atkins
# extracts estimates for 1) Chronic transmitters and 2) Integrated over the course of transmitter infection

transmitterspvl_variantdist <- RunParallel(populationmodel_acrossVL_Environment, pop$transmitter) %>%
  
  # Post-Processing 
  do.call(cbind.data.frame, .) %>% 
  t() %>%
  cbind.data.frame(transmitter = pop$transmitter, .) %>%
  rename(probTransmissionPerSexAct_average = probTransmissionPerSexAct) %>%
  rename_with(function(x) gsub('\\.(?<=\\V)', '\\1_average.\\2',x, perl= TRUE), 
              .cols = !contains('chronic') & contains('V')) %>%
  
  # Long Format
  pivot_longer(cols = c(contains('chronic'),contains('average')), 
               names_to = c('type', 'stage', 'variants'), 
               names_pattern = "(^[^_]+)(_[^.]*)(.*)", 
               values_to = 'p') %>% 
  mutate(stage = gsub('^.*_', '', stage)) %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric()) 


################################### P(Multiple Variants) | Recipient SPVL (H0) ###################################
# Generate weightings using function g from Thompson et al
# This function generates a probability distribution (lognormal) which is used to sample from our range
# of simulated transmitter viral loads to generalise over a population


# Calculate probability of multiple founder infection in recipient
sim_recipientspvl_variantdist <- RunParallel(populationmodel_acrossVL_Environment, sim_spvl$transmitter)  %>%
  do.call(cbind.data.frame, .) %>% 
  t() %>%
  cbind.data.frame(recipient = sim_spvl$recipient, .) %>%
  rename(probTransmissionPerSexAct_average = probTransmissionPerSexAct) %>%
  rename_with(function(x) gsub('\\.(?<=\\V)', '\\1_average.\\2',x, perl= TRUE), 
              .cols = !contains('chronic') & contains('V')) %>%
  
  # Long Format
  pivot_longer(cols = c(contains('chronic'),contains('average')), 
               names_to = c('type', 'stage', 'variants'), 
               names_pattern = "(^[^_]+)(_[^.]*)(.*)", 
               values_to = 'p') %>% 
  mutate(stage = gsub('^.*_', '', stage)) %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric()) %>%
  mutate(model = 'linear')

sim_recipientspvl_variantdist$CD4_decline <- predict(cd4_model, newdata = data.frame(spVL = log10(sim_recipientspvl_variantdist$recipient)))
################################### Write to file ###################################
# population only
write_csv(pop, file = paste(results_dir, 'base_population.csv', sep = '/'))

# transmitter SPVL & variant distribution
write_csv(transmitterspvl_variantdist, file = paste(results_dir, 'transmitterspvl_variantdist.csv', sep = '/'))

# recipient SPVL & variant distribution
write_csv(sim_recipientspvl_variantdist, file = paste(results_dir, 'recipientspvl_variantdist.csv', sep = '/'))


###################################################################################################
################################### Non-Linear Models (No Segregation) ###################################
# Four groups of nonlinear relationships trialled: Polynomial, Concave/Convex Curves, Sigmoidal, Curves with Max/Min

quad_model <- lm(recipient_log10SPVL ~ transmitter_log10SPVL+ I(transmitter_log10SPVL**2), data = pop) 
sim_spvl_quad <- predict(quad_model, newdata= data.frame(transmitter_log10SPVL = sim_transmitter_log10SPVL)) %>% 
  cbind.data.frame(sim_recipient_log10SPVL = ., sim_transmitter_log10SPVL = sim_transmitter_log10SPVL) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SPVL')}")) %>%
  `colnames<-` (str_remove(colnames(.), 'sim_'))  %>% mutate(sim = 'quad')

exp_model <- lm(recipient_log10SPVL ~ exp(transmitter_log10SPVL), data = pop) 
sim_spvl_exp <- predict(exp_model, newdata= data.frame(transmitter_log10SPVL = sim_transmitter_log10SPVL)) %>% 
  cbind.data.frame(sim_recipient_log10SPVL = ., sim_transmitter_log10SPVL = sim_transmitter_log10SPVL) %>%
  mutate(across(.cols = everything(), .fns = ~ 10**.x, .names = "{str_remove(col, '_log10SPVL')}")) %>%
  `colnames<-` (str_remove(colnames(.), 'sim_')) %>% mutate(sim = 'exp')

nls_df <- rbind.data.frame(sim_spvl_quad , sim_spvl_exp, sim_spvl )

fig_3a= ggplot(nls_df) + 
  geom_density(aes(x = recipient_log10SPVL, fill = sim), alpha = 0.4)+
  my_theme +
  scale_fill_brewer(palette = 'OrRd')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                     limits = c(0,8), expand = c(0,0))+ theme(
                       legend.position = 'none')



quad_recipientspvl_variantdist <- RunParallel(populationmodel_acrossVL_Environment, sim_spvl_quad$transmitter)  %>%
  do.call(cbind.data.frame, .) %>% 
  t() %>%
  cbind.data.frame(recipient = sim_spvl_quad$recipient, .) %>%
  rename(probTransmissionPerSexAct_average = probTransmissionPerSexAct) %>%
  rename_with(function(x) gsub('\\.(?<=\\V)', '\\1_average.\\2',x, perl= TRUE), 
              .cols = !contains('chronic') & contains('V')) %>%
  
  # Long Format
  pivot_longer(cols = c(contains('chronic'),contains('average')), 
               names_to = c('type', 'stage', 'variants'), 
               names_pattern = "(^[^_]+)(_[^.]*)(.*)", 
               values_to = 'p') %>% 
  mutate(stage = gsub('^.*_', '', stage)) %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric()) %>%
  mutate(model = 'quad')

quad_recipientspvl_variantdist$CD4_decline <- predict(cd4_model, newdata = data.frame(spVL = log10(quad_recipientspvl_variantdist$recipient)))

exp_recipientspvl_variantdist <- RunParallel(populationmodel_acrossVL_Environment, sim_spvl_exp$transmitter)  %>%
  do.call(cbind.data.frame, .) %>% 
  t() %>%
  cbind.data.frame(recipient = sim_spvl_exp$recipient, .) %>%
  rename(probTransmissionPerSexAct_average = probTransmissionPerSexAct) %>%
  rename_with(function(x) gsub('\\.(?<=\\V)', '\\1_average.\\2',x, perl= TRUE), 
              .cols = !contains('chronic') & contains('V')) %>%
  
  # Long Format
  pivot_longer(cols = c(contains('chronic'),contains('average')), 
               names_to = c('type', 'stage', 'variants'), 
               names_pattern = "(^[^_]+)(_[^.]*)(.*)", 
               values_to = 'p') %>% 
  mutate(stage = gsub('^.*_', '', stage)) %>%
  mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric()) %>%
  mutate(model = 'exp')

exp_recipientspvl_variantdist$CD4_decline <- predict(cd4_model, newdata = data.frame(spVL = log10(exp_recipientspvl_variantdist$recipient)))


combined <- rbind.data.frame(exp_recipientspvl_variantdist, quad_recipientspvl_variantdist, sim_recipientspvl_variantdist)

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


setEPS()
postscript( paste(figs_dir, 'grid_2.eps', sep = '/'), width = 20, height = 10)
grid_2 
dev.off()

################################### Non-Linear Models (Probability MV weights) ###################################
# Four groups of nonlinear relationships trialled: Polynomial, Concave/Convex Curves, Sigmoidal, Curves with Max/Min

# Calculate weights based on transmitter SPVL


poly_model <- 
  
concave_model <- 
  
convex_model <- 
  
sigmoid_model <- 
  
max_model <- 
  



# Generate a matrix of transmitter/founding particles and associated viral loads
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

spvl <- seq(2,7,by=0.5)
prob_dists <- lapply(spvl, populationmodel_acrossVL_Environment) %>% sapply(., function(x) x[2:length(x)])

variant_matrix <- sapply(pop$transmitter_log10SPVL, FoundingVars, sd = 0.5, p_dists = prob_dists) %>% t()

# Sanity check
apply(variant_matrix, 1, function(x) sum(!is.na(x))) %>%
  all() %>%
  stopifnot()


################################### Exclusion ###################################
recip_log10SPVL_exclusion <- apply(variant_matrix, 1, Exclusion, na.rm = TRUE)


################################### Additive ###################################
recip_log10SPVL_additive <- apply(variant_matrix, 1, Additive, na.rm = TRUE)


################################### Interaction ###################################
recip_log10SPVL_interaction <- apply(variant_matrix, 1, Interaction, na.rm = TRUE)


################################### Post-Processing ###################################
v1 <- cbind.data.frame(transmitter_log10SPVL = pop$transmitter_log10SPVL, test)


pop1 <- cbind.data.frame(recipient_log10SPVL = pop$recipient_log10SPVL, transmitter_log10SPVL = pop$transmitter_log10SPVL, hypothesis = 'null', v = variants)
pop2 <- cbind.data.frame(recipient_log10SPVL = 10**recip_max, transmitter_log10SPVL = pop$transmitter_log10SPVL, hypothesis = 'max', v = variants)
pop3 <- cbind.data.frame(recipient_log10SPVL = 10**recip_mean, transmitter_log10SPVL = pop$transmitter_log10SPVL, hypothesis = 'mean', v = variants)

pops_com <- rbind.data.frame(pop1, pop2, pop3)



################################### CD4 decline ###################################

pops_decline <- pops_com %>% 
  mutate(CD4_decline = (0.0111*log10(recipient_log10SPVL))**2) 

ggplot(pops_decline, aes(x = recipient_log10SPVL, y = decline, colour = hypothesis, shape = as.character(v))) +
  geom_point(alpha = 0.5) +
  scale_x_log10(limits = c(1, 10**10),
                expand = c(0,0),
                name = expression(paste("Recipient SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(name = "CD4 Decline", expand = c(0,0))+
  #scale_y_log10(name = expression(paste("CD4 Decline")),
                #limits = c(1, 10**10),
               # expand = c(0,0)#,
                #breaks = trans_breaks("log10", function(x) 10**x),
                #labels = trans_format("log10", math_format(.x))
                #) +
  theme_classic()+
  theme(
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )+
  annotation_logticks(sides = 'b') 



ggplot(pops_decline, aes(x = v, y = decline, colour = hypothesis, shape = as.character(v))) +
  geom_boxplot() +
  theme_classic()+
  theme(
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )



