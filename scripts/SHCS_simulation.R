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


shcs_data_int_list <- shcs_data %>%
  rowwise() %>%
  filter(riskgroup_1 == riskgroup_2) %>%
  filter(!if_any(starts_with('riskgroup'), ~ . %in% c('UNKNOWN', 'OTHER'))) %>%
  select(contains(c('sex', 'age', 'riskgroup', 'log10_SpVL_couplemean'))) %>%
  mutate(across(starts_with('sex'), .fns = ~ match(.x, c('M', 'F')))) %>% 
  mutate(across(starts_with('riskgroup'), .fns = ~ match(.x, c('HET', 'MSM', 'PWID')))) %>%
  group_by(riskgroup_1) %>%
  group_split() %>%
  lapply(., function(x) x %>% select(where(~n_distinct(.) > 1)))

  
cum_probs_list <- lapply(shcs_data_int_list, function(x) x %>% select(where(is.integer)) %>% apply(2, function(i) cumsum(table(i))/length(i), simplify = FALSE))

TestFunc <- function(dataframe,probs) {
  out <- mvtnorm::rmvnorm(n = 200, mean = as.vector(colMeans(dataframe)), sigma = cov(dataframe)) %>%
    as_tibble(name_repair = NULL) %>%
    setNames(colnames(dataframe)) %>%
    mutate(across(starts_with(c('sex', 'riskgroup')), 
                  .fns = ~cut(.x,
                              right = T,
                              include.lowest = T,
                              breaks = as.vector(c(min(.x), quantile(.x, probs[[cur_column()]]))),
                              labels = 1:length(probs[[cur_column()]])))) %>%
    mutate(across(starts_with('sex'), .fns = ~ case_when(.x == 1 ~ 'M',
                                                         .x == 2 ~ 'F'))) %>%
    mutate(across(starts_with('riskgroup'), .fns = ~ case_when(.x == 1 ~ "HET",
                                                               .x == 2 ~ "MSM",
                                                               .x == 3 ~ "PWID",
                                                               .x == 4 ~ "UNKNOWN",
                                                               .x == 5 ~ "OTHER")))
  return(out)}


sim_data_int_list <- mapply(TestFunc ,
  data = shcs_data_int_list,
  probs = cum_probs_list,
  SIMPLIFY = F)

# Rejoing other levels (transmission route, sex for msm)

# Plot covariance matrices

