NPAIRS <- 100 # lower for test runs on low cpu machines
R2 <- 0.33 # Heritability estimate
RAKKAI_MEAN <- 4.39
RAKKAI_SD <- 0.84

zambia <- c(vl_mean = 4.74, 
            vl_sd = 0.61)
zambia_carrier <- rnorm(1000, 4.74, 0.61)
amsterdam <- c(vl_mean = 4.35, 
               vl_sd = 0.47)

transmission_prob <- (1/(zambia_carrier*0.96*sqrt(2*pi))*exp(-(zambia_carrier-4.64)**2)/(2*0.96**2)) 
transmission_potential <- norm(1000, 4.6, sqrt(0.96)) %>% hist()

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