table(shcs_data$riskgroup_1, shcs_data$riskgroup_2) %>% 
  as.matrix() %>%
  reshape2::melt() %>% 
  tibble() %>%
  rename(riskgroup_1 = Var1, riskgroup_2 = Var2) %>%
  filter(!if_any(contains('riskgroup'), ~. %in% c('OTHER', 'UNKNOWN'))) %>%
  ggplot(aes(x = riskgroup_1, y = riskgroup_2, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  geom_text(aes(label = value), color = "black", size = 16) +
  scale_fill_distiller(palette = 'OrRd') + 
  my_theme

shcs_data %>%
  rowwise() %>%
  filter(riskgroup_1 == riskgroup_2) %>%
  ggplot(aes(y=SpVL_couplemean, x = riskgroup_1)) + 
  geom_boxplot() + 
  #facet_wrap(. ~ riskgroup_1)+ 
  #scale_y_continuous(expand = c(0,0)) + 
  scale_y_log10(limits = c(10**1, 10**8),
                expand = c(0,0),
                name = expression(paste("SpVL Pair Mean", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  my_theme + 
  stat_compare_means(method = "kruskal.test",label.x = 1.7, label.y = 7) +
  geom_jitter()

# Encode categorical levels as integers, then split dataframe by riskgroup
shcs_data_int_list <- shcs_data %>%
  rowwise() %>%
  filter(riskgroup_1 == riskgroup_2) %>%
  filter(!if_any(starts_with('riskgroup'), ~ . %in% c('UNKNOWN', 'OTHER'))) %>%
  select(contains(c('sex', 'age', 'riskgroup', 'log10_SpVL_couplemean'))) %>%
  mutate(across(starts_with('sex'), .fns = ~ match(.x, c('M', 'F')))) %>% 
  mutate(across(starts_with('riskgroup'), .fns = ~ match(.x, c('HET', 'MSM', 'PWID')))) %>%
  group_by(riskgroup_1) %>%
  group_split() %>%
  lapply(., function(x) x %>% select(where(~n_distinct(.) > 1))) %>%
  setNames(c('HET', 'MSM', 'PWID'))

# infer cumulative probabilites of discrete variables (within each riskgroup)
cum_probs_list <- lapply(shcs_data_int_list, function(x) x %>% select(where(is.integer)) %>% apply(2, function(i) cumsum(table(i))/length(i), simplify = FALSE))

# Function to simulate new data from a multivariate normal distribution, parametersied by 
# the variance - covariance matrix for each group. categorical variables represented as integers
# and re-factored according to cumulative probability distributions calculated. 
TestFunc <- function(dataframe,probs) {
  lb <- ifelse(grepl("age.inf", names(dataframe)), 16, -Inf)
  
  out <- tmvtnorm::rtmvnorm(n = 200, 
                          mean = as.vector(colMeans(dataframe)), 
                          sigma = cov(dataframe),
                          lower = lb,
                          upper = rep(Inf, length = ncol(dataframe))) %>%
    as_tibble(name_repair = NULL) %>%
    setNames(colnames(dataframe)) %>%
    mutate(across(starts_with(c('sex', 'riskgroup')), 
                  .fns = ~cut(.x,
                              right = T,
                              include.lowest = T,
                              breaks = as.vector(c(min(.x), quantile(.x, probs[[cur_column()]]))),
                              labels = 1:length(probs[[cur_column()]])))) %>%
    add_column(., !!!c(sex_1 = 1, sex_2 = 1)[setdiff(c('sex_1', 'sex_2'), colnames(.))]) %>%
    mutate(across(starts_with('sex'), .fns = ~ case_when(.x == 1 ~ 'M',
                                                         .x == 2 ~ 'F'))) 

  return(out)}


sim_data <- mapply(TestFunc,
                   data = shcs_data_int_list,
                   probs = cum_probs_list,
                   SIMPLIFY = F) %>%
  setNames(., c('HET', 'MSM', 'PWID')) %>%
  bind_rows(., .id = "riskgroup") %>% 
  pivot_longer(cols = - c(contains('couplemean'), riskgroup), 
               names_to = c(".value", "partner"), 
               names_pattern  = "^(.*)_([0-9])$") %>% 
  mutate(partner = factor(partner, levels = c("1", "2"))) %>%
  relocate(c(partner, sex), .before = riskgroup) %>%
  relocate(ends_with('couplemean'), .after = last_col()) 

sim_data_int_list <- sim_data %>%
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


ggplot(bias_means) +
  geom_col(aes(x = variable, y = bias))+ 
  facet_wrap(.~riskgroup) +
  my_theme + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Plot bias of new data vs original swiss data
bias_covmats <- mapply(function(x_sim, x_emp) x_sim - x_emp, x_sim = sim_data_covmats, x_emp = shcs_data_covmats, SIMPLIFY = F) %>%
  reshape2::melt() %>% 
  tibble() %>%
  rename(riskgroup = L1)

ggplot(bias_covmats, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  geom_text(aes(label = signif(value, 3)), color = "black", size = 3) +
  scale_fill_distiller(palette = 'OrRd') + 
  facet_wrap(.~riskgroup) +
  my_theme+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

