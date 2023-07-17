# Validating the composition of simulated cohorts in "RECONCILING THE FOUNDER VARIANT MULTIPLICITY 
# OF HIV-1 INFECTION WITH THE RATE OF CD4+ DECLINE"

# This script calculates normalised measures of bias and prepares tidy data frames to be plotted
# as confusion matrices in the supplementary materials. 

require('./scripts/simulate_cohort_funcs.R')
# MUST HAVE RUN ./scripts/main.R UP TO LINE 85

# Sanity Check
stopifnot('stratified_data' %in% ls(envir = .GlobalEnv))


# Infer numerical repressentation of simulated data to facilitate calculation of covariance
# matrices and mean to compare with empirical data
sim_data_int_list <- stratified_data %>%
  pivot_wider(names_from = partner, values_from = c(sex, age.inf)) %>%
  mutate(across(starts_with('sex'), .fns = ~ match(.x, c('M', 'F')))) %>% 
  mutate(across(starts_with('riskgroup'), .fns = ~ match(.x, c('HET', 'MSM', 'PWID')))) %>%
  relocate(ends_with('couplemean'), .after = last_col())  %>%
  group_by(riskgroup) %>%
  group_split() %>%
  lapply(., function(x) x %>% select(where(~n_distinct(.) > 1))) %>%
  setNames(c('HET', 'MSM', 'PWID'))


# Estimate cov-var matrices
sim_data_covmats <- lapply(sim_data_int_list, cov)
shcs_data_covmats <- lapply(shcs_data_int_list, cov)


# Calculate normalise
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

sim_data_means <- lapply(sim_data_int_list, function(x) x %>% mutate(across(everything(), .fns = ~min_max_norm(.x)))) %>%
  lapply(., colMeans)

shcs_data_means <- lapply(shcs_data_int_list, function(x) x %>% mutate(across(everything(), .fns = ~min_max_norm(.x)))) %>%
  lapply(., colMeans)


# Calculate bias (x_sim - x)
bias_means <- mapply(function(x_sim, x_emp) (x_sim - x_emp), x_sim = sim_data_means, x_emp = shcs_data_means, SIMPLIFY = F) %>%
  lapply(., bind_rows) %>%
  lapply(., function(x) add_column(x, !!!c(sex_1 = 0, sex_2 = 0)[setdiff(c('sex_1', 'sex_2'), colnames(x))])) %>%
  bind_rows() %>%
  add_column(riskgroup = c('HET', 'MSM', 'PWID')) %>%
  pivot_longer(cols = -riskgroup, names_to = 'variable', values_to = 'bias')

bias_covmats <- mapply(function(x_sim, x_emp) x_sim - x_emp, x_sim = sim_data_covmats, x_emp = shcs_data_covmats, SIMPLIFY = F) %>%
  reshape2::melt() %>% 
  tibble() %>%
  rename(riskgroup = L1)

