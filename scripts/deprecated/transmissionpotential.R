NPAIRS <- 100 # lower for test runs on low cpu machines


# Initial distributions describe carrier population
zambia_carrier <- c(vl_mean = 4.74, vl_sd = 0.61, h2 = NA)
amsterdam_carrier <- c(vl_mean = 4.35, vl_sd = 0.68, h2 = NA)
rakkai_carrier <- c(vl_mean = 4.39, vl_sd = 0.84, h2 = 0.33)
TP <- c(tp_mean = 4.64, tp_sd = 0.96)


InitDonor <- function(carrier, TP){
  #stopifnot() #stop if TP does not contain mean and variance
  
  carrier_mean <- carrier[['vl_mean']]
  carrier_var <- carrier[['vl_sd']] ** 2
  carrier_h2 <- carrier[['h2']]
  TP_mean <- TP[['tp_mean']]
  TP_var <- TP[['tp_sd']] ** 2
  
  # Calculate donor population summary statistics from B8 (Eq 7) Bonhoeffer
  donor_mean <- (carrier_mean*TP_var + TP_mean*carrier_var) / (carrier_var + TP_var)
  donor_sd <- ((carrier_var*TP_var) / (carrier_var + TP_var)) %>% sqrt()
  
  return(c(vl_mean = donor_mean, vl_sd = donor_sd, h2 = carrier_h2))

  
}

zambia_donor <- InitDonor(zambia_carrier, TP)
amsterdam_donor <- InitDonor(amsterdam_carrier, TP)
rakkai_donor <- InitDonor(rakkai_carrier, TP) #Must check whether TP can be generally applicable or requires recalculation


# Calculate recipient set point viral loads from a given population of donor set point viral loads.
# Back calculation from a set correlation coefficient (r2) and sum of squared residuals
# 'Parent ~ Offspring' regression R2 analogous to heritability of viral load. 
GetRecipVL <- function(donorload, h2 = 0.33){
  n <- length(donorload)
  ssr <- sum((donorload-mean(donorload))**2) # sum of squared residuals
  e <- rnorm(n)
  e <- resid(lm(e ~ donorload))
  e <- e*sqrt((1-h2)/h2*ssr/(sum(e**2)))
  recipload <- donorload + e
  
  return(recipload)
}


# Initialise donor~recipient pairs from a given distribution of donor logspvls and 
# regression coefficient corresponding to the heritability of viral load
InitDonorRecip <- function(popsize, donor_population){
  #must contain xyz
  
  donor_mean <- donor_population[['vl_mean']]
  donor_sd <- donor_population[['vl_sd']]
  donor_h2 <- donor_population[['h2']]
  
  init_donor <- rnorm(popsize,  mean = donor_mean, sd = donor_sd)
  init_recip <- -1
  
  while(min(init_recip)<0){
    init_recip <- GetRecipVL(init_donor, donor_h2) 
  }
  
  out <- list('donor' = cbind.data.frame(vl=init_donor) , 'recipient' = cbind.data.frame(vl = init_recip))
  return(out)
}


rakkai_donorrecip <- InitDonorRecip(NPAIRS, rakkai_donor)


rakkai_cohort <- c(list('carrier' = cbind.data.frame(vl = rnorm(NPAIRS, rakkai_carrier[['vl_mean']], rakkai_carrier[['vl_sd']]))), 
                     rakkai_donorrecip) %>% 
  bind_rows(.id = 'population') %>%
  mutate(vl = 10**vl)


ggplot(rakkai_cohort, aes(x = vl, fill = population)) +
  geom_density(
    alpha = 0.5) +
  scale_x_log10(limits = c(1, 10**10),
                expand = c(0,0),
                name = expression(paste(" SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_continuous(expand = c(0,0))+
  theme_classic() 


dotplot_data <- rakkai_cohort %>% filter(population %in% c('donor', 'recipient')) %>%  
  group_by(population) %>%
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = population, values_from = vl) %>%
  select(-row)

ggplot(dotplot_data, aes(x = donor, recipient)) +
  geom_point(colour = '#2ca25f', #'#CB6015' #'#66c2a4','#2ca25f','#006d2c'
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
  theme_classic() +
  geom_smooth(method = lm, colour = 'black', se = F) + 
  stat_poly_eq(formula = y ~ x)


