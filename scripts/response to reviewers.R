################################### Dependencies ###################################
require(tidyverse)
source('./scripts/BonhoefferEqns.R')

DecomposeRecipientSpVL <- function(recipient_mean, recipient_var, p_mv, effectsize){
  
  # recipient_mean = (1 - p_mv)*x + p_mv*(x + effect_size)
  # Because weights of finite mixture must sum to 1; rearranges to:
  sv_mean <- recipient_mean - p_mv*effectsize  
  mv_mean <- sv_mean + effectsize
  
  
  #recipient_var <- p_m*(y + mv_mean^2) + (1-p_mv)*(y + sv_mean^2) - recipient_mean^2
  # Rearranges for the same reasons
  # Assumption of equal variances justified by Janes et al.
  decomposed_var <- recipient_var + recipient_mean^2 - p_mv * mv_mean^2 - (1 - p_mv) * sv_mean^2
  
  out <- list('sv' = c('mean' = sv_mean,
                       'var' = decomposed_var),
              'mv' = c('mean' = mv_mean, 
                       'var' = decomposed_var))
  
  return(out)
}


zambia_variance <- 0.78**2 #0.78 is the standard deviation

zambia_mean <- 4.74 

recipient_dist <- CalcRecipient(zambia_mean, 
                                zambia_variance)



################################### Plot component distributions ###################################
effect_size <- seq(0.1, 1, by = 0.1)

expanded_df <- expand_grid(recipient_dist[['mean']], 
                           recipient_dist[['var']],
                           p_mv = 0.25,
                           effect_size) %>%
  mutate(test = pmap(., ~DecomposeRecipientSpVL(..1, ..2, ..3, ..4))) %>% 
  unnest_wider(test) %>%
  unnest_wider(c(sv,mv), names_sep = '_') %>%
  pivot_longer(starts_with(c('mv', 'sv')),
               names_to = c("multiplicity", ".value"),
               names_sep = "_") %>%
  mutate(sd = sqrt(var)) %>%
  mutate(p_mv = as_factor(p_mv),
         effect_size = as_factor(effect_size))%>%
  select(-c(starts_with('recipient'), var, -p_mv)) %>%
  crossing(x = seq(1, 8, length.out = 100)) %>%  # Create combinations of each row with x values
  mutate(y = dnorm(x, mean, sd))
  
ggplot(expanded_df, aes(x = x, y = y, colour = multiplicity)) +
  geom_line() +
  facet_wrap(.~effect_size, nrow = 2) +
  theme_minimal() +
  scale_x_continuous(expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')'))) +
  scale_y_continuous('Probability Density') +
  scale_colour_discrete('Founder Variant Multiplicity', labels = as_labeller(c('mv' = 'Multiple', 'sv' = 'Single'))) +
  theme(legend.position = 'bottom')
 
################################### Plot component qq plots ###################################

effect_size <- seq(0.1, 1, by = 0.1)
p_mv <- seq(0.1, 0.5, by = 0.1)

expanded_df_2 <- expand_grid(recipient_dist[['mean']], 
                             recipient_dist[['var']],
                           p_mv,
                           effect_size)  %>%
  mutate(test = pmap(., ~DecomposeRecipientSpVL(..1, ..2, ..3, ..4))) %>%
  unnest_wider(test) %>%
  unnest_wider(c(sv,mv), names_sep = '_') %>%
  mutate(mv_samples = pmap(., ~rnorm(5000, mean = ..7, sd = sqrt(..8))),
         sv_samples = pmap(., ~rnorm(5000, mean = ..5, sd = sqrt(..6)))) %>%
  unnest_longer(c(mv_samples, sv_samples), indices_include = F) %>%
  pivot_longer(c(mv_samples, sv_samples), values_to = 'viral_load', names_to = 'multiplicity') %>%
  mutate(p_mv = as_factor(p_mv),
         effect_size = as_factor(effect_size))


ggplot(expanded_df_2, aes(sample = viral_load)) +
  geom_qq(size = 0.5) +
  geom_qq_line() +
  theme_minimal() +
  facet_grid(cols = vars(p_mv), rows = vars(effect_size)) +
  labs(x = "Theoretical Quantiles",
       y = "Sample Quantiles")



################################### Plot meta analysis distributions ###################################

# first ensure only one measurement is present
meta_data_vls <- meta_data %>% 
  group_by(cohort, participant.id) %>% 
  slice_sample(n = 1) %>% 
  ungroup() %>% 
  mutate(viral.load = as.numeric(viral.load)) %>% 
  filter(!is.na(viral.load)) %>% 
  filter(!is.na(multiple.founders))

meta_data_vls_summary <- meta_data_vls %>%
  summarise(m = mean(log10(viral.load)), .by = multiple.founders)

ggplot(meta_data_vls) + 
  geom_density(aes(x = log10(viral.load), fill = multiple.founders, colour = multiple.founders), alpha = 0.5) + 
  geom_vline(aes(xintercept = m,colour = multiple.founders), data = meta_data_vls_summary, ) + 
  theme_minimal() +
  scale_x_continuous(expression(paste("SpVL", ' (', Log[10], " copies ", ml**-1, ')'))) +
  scale_y_continuous('Probability Density') +
  scale_fill_discrete('Founder Variant Multiplicity', labels = as_labeller(c('multiple' = 'Multiple', 'single' = 'Single'))) +
  scale_colour_discrete('Founder Variant Multiplicity', labels = as_labeller(c('multiple' = 'Multiple', 'single' = 'Single'))) +
  theme(legend.position = 'bottom')

