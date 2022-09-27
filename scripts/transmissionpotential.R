NPAIRS <- 100 # lower for test runs on low cpu machines
R2 <- 0.33 # Heritability estimate
RAKKAI_MEAN <- 4.39
RAKKAI_SD <- 0.84


# Initial distributions describe carrier population
zambia_carrier <- c(vl_mean = 4.74, vl_sd = 0.61)
amsterdam_carrier <- c(vl_mean = 4.35, vl_sd = 0.68)
TP <- c(tp_mean = 4.64, tp_sd = 0.96)


InitDonor <- function(carrier, TP){
  #stopifnot() #stop if TP does not contain mean and variance
  
  carrier_mean <- carrier[['vl_mean']]
  carrier_var <- carrier[['vl_sd']] ** 2
  TP_mean <- TP[['tp_mean']]
  TP_var <- TP[['tp_sd']] ** 2
  
  # Calculate donor population summary statistics from B8 (Eq 7) Bonhoeffer
  donor_mean <- (carrier_mean*TP_var + TP_mean*carrier_var) / (carrier_var + TP_var)
  donor_sd <- ((carrier_var*TP_var) / (carrier_var + TP_var)) %>% sqrt()
  
  return(c(vl_mean = donor_mean, vl_sd = donor_sd))

  
}

zambia_donor <- InitDonor(zambia_carrier, TP)
amsterdam_donor <- InitDonor(amsterdam_carrier, TP)




donor_mean <- (zambia[['vl_mean']]*0.96+4.64*zambia[['vl_sd']]**2)/(zambia[['vl_sd']]**2+0.96)
donor_sd <- (zambia[['vl_sd']]**2*0.96)/(zambia[['vl_sd']]**2+0.96) %>% sqrt()

zambia_donor <- rnorm(1000, donor_mean, donor_sd)
zambia_carrier <- rnorm(1000, 4.74, 0.61)


zambia_cohort <- list('carrier' = cbind.data.frame(vl = zambia_carrier), 'donor' = cbind.data.frame(vl = zambia_donor)) %>% bind_rows(.id = 'population')


ggplot(zambia_cohort, aes(x = vl, fill = population)) +
  geom_density(
    alpha = 0.5) +
  theme_classic() +
  scale_x_log10(limits = c(1, 10**10),
                expand = c(0,0),
                name = expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  theme_classic() 

ggplot()
#from bonhoeffer: 
(1 / sqrt(2*pi*0.96)) * exp( -( ((zambia_carrier-4.64)**2) / (2*0.96)) )
transmission_prob <- (1/(zambia_carrier*0.96*sqrt(2*pi))*exp(-(zambia_carrier-4.64)**2)/(2*0.96**2)) 
transmission_potential <- rnorm(1000, 4.6, sqrt(0.96)) %>% hist()

plot(zambia_carrier, transmission_potential )

zambia_population <- InitPop(1000, 
                            donor_mean = zambia[['mean']], 
                            donor_sd = zambia[['sd']], 
                            herit = R2)

amsterdam_population <- InitPop(1000, 
                             donor_mean = amsterdam[['mean']], 
                             donor_sd = amsterdam[['sd']], 
                             herit = R2)

populations <- lapply(list(zambia, amsterdam), function(x) {
  model_population <- InitPop(NPAIRS, donor_mean = x[['mean']], donor_sd = x[['sd']], herit = R2)
  donor_spvl <- 10**model_population$donor
  recip_spvl <- 10**model_population$recip
  spvl_df <- cbind.data.frame(donor_spvl = donor_spvl, 
                              recip_spvl = recip_spvl)})

names(populations) <- c('zambia', 'amsterdam')

population_df <- bind_rows(populations, .id = 'cohort')





ggplot(population_df, aes(x = donor_spvl, fill = cohort)) +
  geom_histogram(#usher colours?
    alpha = 0.5) +
  scale_x_log10(limits = c(1, 10**10),
                expand = c(0,0),
                name = expression(paste("Donor SPVL", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  theme_classic() 



ggplot(population_df, aes(x = donor_spvl, recip_spvl,colour = cohort)) +
  geom_point(#usher colours?
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