# Import and process SHCS data for "RECONCILING THE FOUNDER VARIANT MULTIPLICITY OF HIV-1 
# INFECTION WITH THE RATE OF CD4+ DECLINE"
#
# Paired set point viral load (SpVL), age, riskgroup and sex data provided on request
# by the Swiss HIV Cohort Study. These data were used in the analysis of Bertels et. al
# 2018 https://academic.oup.com/mbe/article/35/1/27/4210012

# NB: partner number is allocated randomly and does not assign transmitter/recipient status

shcs_data <-read_csv("data/shcs_data.csv", 
                     col_types = cols(sex.1 = readr::col_factor(levels = c("M", "F")), 
                                      sex.2 = readr::col_factor(levels = c("M", "F")), 
                                      riskgroup.1 = readr::col_factor(levels = c("HET", "MSM",
                                                                                 "IDU", "OTHER",
                                                                                 "UNKNOWN")), 
                                      riskgroup.2 = readr::col_factor(levels = c("HET", "MSM", 
                                                                                 "IDU", "OTHER", 
                                                                                 "UNKNOWN"))))%>%
  select(-contains('ID')) %>%
  rowid_to_column( "ID.pair") %>%
  mutate(across(contains('riskgroup'), 
                .fns = ~ fct_recode(.x, PWID = "IDU"))) %>%
  rename_with( ~ stri_replace_last_fixed(.x, 'spVL', 'SpVL')) %>%
  mutate(across(contains('SpVL'), 
                .fns = ~ raise_to_power(10, .x))) %>%
  mutate(SpVL.couplemean = rowMeans(across(contains('SpVL')))) %>%
  mutate(across(contains('SpVL'),
                .fns = ~log10(.x), .names = "log10_{.col}")) %>%
  rowwise() %>%
  mutate(transmitter.allocaterandom = ClassifyVL(log10_SpVL.1, log10_SpVL.2, transmitter = t, recipient = r, method = 'random')) %>%
  mutate(transmitter.allocateML = ClassifyVL(log10_SpVL.1, log10_SpVL.2, transmitter = t, recipient = r, method = 'ml')) %>%
  rename_with( ~ stri_replace_last_fixed(.x, '.', '_')) 


# Transform to long-format data table
shcs_data_long_transmitterrandom <- shcs_data %>% 
  pivot_longer(cols = -c(ID_pair, starts_with('transmitter'), contains('couplemean')), 
               names_to = c(".value", "partner"),
               names_pattern  = "^(.*)_([0-9])$") %>% #matches('SpVL_[[:digit:]]')
  #mutate(across(ends_with('SpVL'),
               # .fns = ~ mean(.x), .names = "{.col}_cohortmean")) %>%
  mutate(partner = as.integer(partner)) %>%
  mutate(partner = case_when(partner == transmitter_allocaterandom ~ 'transmitter',
                             .default = 'recipient')) %>%
  mutate(partner = factor(partner, levels = c("transmitter", "recipient"))) %>%
  mutate(age.inf_category = cut(age.inf, 
                                breaks = c(15,24,29,39,80), 
                                labels = c('15-24', '25-29', '30-39','40-80'))) %>% # cut by default is exclusive of the lower bound
  select(-starts_with('transmitter')) %>%
  relocate(ID_pair) %>%
  relocate(age.inf, .after = sex) %>%
  relocate(ends_with('SpVL'), .after = riskgroup) %>%
  relocate(ends_with('couplemean'), .after = log10_SpVL)


# Long data with max(SpVL|Pair) as transmitter
shcs_data_long_transmitterML<- shcs_data %>% 
  pivot_longer(cols = -c(ID_pair, starts_with('transmitter'), contains('couplemean')), 
               names_to = c(".value", "partner"),
               names_pattern  = "^(.*)_([0-9])$") %>% 
  mutate(partner = as.integer(partner)) %>%
  mutate(partner = case_when(partner == transmitter_allocateML~ 'transmitter',
                             .default = 'recipient')) %>%
  mutate(partner = factor(partner, levels = c("transmitter", "recipient"))) %>%
  mutate(age.inf_category = cut(age.inf, 
                                breaks = c(15,24,29,39,80), 
                                labels = c('15-24', '25-29', '30-39','40-80'))) %>% # cut by default is exclusive of the lower bound
  select(-starts_with('transmitter')) %>%
  mutate(ID_pair = ID_pair + 196) %>% #So that ID pair numbers are contiguous
  relocate(ID_pair) %>%
  relocate(age.inf, .after = sex) %>%
  relocate(ends_with('SpVL'), .after = riskgroup) %>%
  relocate(ends_with('couplemean'), .after = log10_SpVL)

# END #