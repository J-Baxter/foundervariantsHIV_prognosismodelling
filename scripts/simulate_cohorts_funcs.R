# Functions facilitating simualtion of patient characteristics, defined by variance-covariance
# matrices from SHCS data

################################### Dependencies ###################################
require(tidyverse)
require(tmvtnorm)

################################### EnumerateLevels ###################################
# Encode categorical levels as integers, then split dataframe by riskgroup

EnumerateLevels <- function(x){
  
  out <- x  %>%
    select(contains(c('sex', 'age', 'riskgroup', 'log10_SpVL_couplemean'))) %>%
    unite(sex ,c(sex_1, sex_2), sep = '') %>%
    mutate(across(starts_with('sex'), .fns = ~ match(.x, c('MF', 'FM', 'MM', 'FF')))) %>% 
    mutate(across(starts_with('riskgroup'), .fns = ~ match(.x, c('HET', 'MSM', 'PWID')))) %>%
    group_by(riskgroup_1) %>%
    group_split() %>%
    # Remove columns where there is only one variable
    lapply(., function(x) x %>% select(where(~n_distinct(.) > 1))) %>%
    setNames(c('HET', 'MSM', 'PWID'))
  
  return(out)
}


################################### SimCohorts ###################################
# Function to simulate new data from a multivariate normal distribution, parametersied by 
# the variance - covariance matrix for each group. categorical variables represented as integers
# and re-factored according to cumulative probability distributions calculated. 

SimCohorts <- function(dataframe,probs) {
  
  # Log10 transform age.inf
  dataframe <- dataframe %>%
    mutate(across(starts_with('age.inf'), .fns = ~ log10(.x)))
  
  # Set lower bound for age (log10(16) = 1.20412)
  # Set upper bound for age (log10(80) = 1.90309)
  lb <- ifelse(grepl("age.inf", names(dataframe)), 1.20412 , -Inf)
  ub <- ifelse(grepl("age.inf", names(dataframe)), 1.90309 , Inf)
  
  # Simulate from truncated multivariate normal
  out <- tmvtnorm::rtmvnorm(n = 200, 
                            mean = as.vector(colMeans(dataframe)), 
                            sigma = cov(dataframe),
                            lower = lb,
                            upper = ub,
                            algorithm = 'rejection' # Use rejection sampling instead of Gibbs
  ) %>%
    as_tibble(name_repair = NULL) %>%
    setNames(colnames(dataframe)) %>%
    # exponentiate age.inf 
    mutate(across(starts_with('age.inf'), .fns = ~ raise_to_power(10, .x))) %>%
  
    # Cut estimates to infer categorical variables
    mutate(across(starts_with(c('sex', 'riskgroup')), 
                  .fns = ~cut(.x,
                              right = T,
                              include.lowest = T,
                              breaks = as.vector(c(min(.x), quantile(.x, probs[[cur_column()]]))),
                              labels = 1:length(probs[[cur_column()]])))) %>%
    
    # add column with int = 3 if a column named sex does not exist (exclusively MSM)
    add_column(., !!!c(sex= 3)[setdiff('sex', colnames(.))]) %>%
    mutate(across(starts_with('sex'), .fns = ~ case_when(.x == 1 ~ 'MF', #c('MF', 'FM', 'MM', 'FF')
                                                         .x == 2 ~ 'FM',
                                                         .x == 3 ~ 'MM', 
                                                         .x == 4 ~'FF'))) 
  
  return(out)
  }


#END