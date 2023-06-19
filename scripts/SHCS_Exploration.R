#Swiss Cohort Analysis
source('./scripts/dependencies.R')
library(performance)
library(stringi)
library(lme4)

shcs_data <-read_csv("data/shcs_data.csv", 
                     col_types = cols(sex.1 = readr::col_factor(levels = c("M", "F")), 
                                      sex.2 = readr::col_factor(levels = c("M", "F")), 
                                      riskgroup.1 = readr::col_factor(levels = c("HET", "MSM", "IDU", "OTHER", "UNKNOWN")), 
                                      riskgroup.2 = readr::col_factor(levels = c("HET", "MSM", "IDU", "OTHER", "UNKNOWN"))))%>%
  select(-contains('ID')) %>%
  rowid_to_column( "ID.pair") %>%
  mutate(across(contains('riskgroup'), .fns = ~ fct_recode(.x, PWID = "IDU"))) %>%
  rename_with( ~ stri_replace_last_fixed(.x, 'spVL', 'SpVL')) %>%
  mutate(across(contains('SpVL'), ~ raise_to_power(10,.x))) %>%
  mutate(SpVL.couplemean = rowMeans(across(contains('SpVL')))) %>%
  mutate(across(contains('SpVL'), .fns = ~log10(.x), .names = "log10_{.col}")) %>%
  rename_with( ~ stri_replace_last_fixed(.x, '.', '_')) 
  

# Transform to long-format data table
shcs_data_long <- shcs_data %>% 
  pivot_longer(cols = -c(ID_pair, contains('couplemean')), names_to = c(".value", "partner"), names_pattern  = "^(.*)_([0-9])$") %>%
  mutate(across(ends_with('SpVL'), .fns = ~ mean(.x), .names = "{.col}_cohortmean")) %>%
  mutate(partner = factor(partner, levels = c("1", "2"))) %>%
  relocate(ID_pair) %>%
  relocate(age.inf, .after = sex) %>%
  relocate(ends_with('SpVL'), .after = riskgroup) %>%
  relocate(ends_with('couplemean'), .after = log10_SpVL) 



shcs_plt1a <- ggplot(shcs_data , aes(x = SpVL_1, SpVL_2)) +
  geom_point( #'#CB6015' #'#66c2a4','#2ca25f','#006d2c'
    colour = '#ef654a',
    shape = 4, size = 3) +
  scale_x_log10(limits = c(10**2, 10**7),
                expand = c(0.05,0),
                name = expression(paste("SpVL Partner 1", ' (', Log[10], " copies ", ml**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_log10(name = expression(paste("SpVL Partner 2", ' (', Log[10], " copies ", ml**-1, ')')),
                limits = c(10**2, 10**7),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x),
                labels = trans_format("log10", math_format(.x))) +
   annotation_logticks() +
  my_theme

# Sex
shcs_plt1b <- ggplot(data_model , aes(x = sex, fill = partner)) +
  geom_bar() + 
  scale_y_continuous(expand = c(0,0), name = 'Count') + 
  scale_x_discrete(expand = c(0,0), name = 'Sex') +
  scale_fill_manual(values = c('#ef654a','#fdd49e'), name = 'Partner')+
  my_theme


#Risk Group
shcs_plt1b <- ggplot(data_model , aes(x =  riskgroup)) +
  geom_bar( fill = '#ef654a') + 
  scale_y_continuous(expand = c(0,0), name = 'Count') + 
  scale_x_discrete(expand = c(0,0), name = 'Risk Group') +
  #scale_fill_brewer(palette = 'OrRd', name = 'Sex')+
  my_theme


# Age
shcs_plt1c <-   #scale_fill_manual(values = c('#ef654a','#fdd49e'), name = 'Partner')+
  #scale_fill_manual(values = c('#ef654a','#fdd49e'), name = 'Partner')+
  my_theme

shcs_plt1c <- ggplot(data_model) + 
  geom_histogram(aes(x = age.inf, fill = partner)) + 
  scale_y_continuous(expand = c(0.01,0.01), name = 'Count') + 
  scale_x_continuous(expand = c(0.01,0.01), name = 'Age at Infection',
                     breaks = seq(0,12, by = 1)) +
  scale_fill_manual(values = c('#ef654a','#fdd49e'), name = 'Partner')+
  my_theme


shcs_panel <- plot_grid(shcs_plt1a, shcs_plt1b, shcs_plt1c, align = 'hv', labels = 'AUTO', nrow = 1)
ggsave(plot = shcs_panel, filename = paste(figs_dir,sep = '/', "shcs_panel.jpeg"), device = jpeg, width = 14, height = 6) #Within-host dynamics -functional

# Linear Models

# Null Model: one coefficient, the overall mean for all individuals
null_model <- lm(log10_SpVL ~ log10_SpVL_cohortmean, data = shcs_data_long)

# Unadjusted Model: each viral load is predicted by the coefficient for the couple
unadjusted_model <- lmer(log10_SpVL ~ (1|log10_SpVL_couplemean), data = shcs_data_long) 

# Adjusted Model
adjusted_model <- lmer(log10_SpVL ~ (1|log10_SpVL_couplemean) + partner + sex + age.inf + riskgroup, data = shcs_data_long) 

# Adjusted Model - No couple effect
adjustednocouple_model <- lm(log10_SpVL ~ 1 + partner + sex + age.inf + riskgroup, data = shcs_data_long) 

performance::r2(adjusted_model)

# Adjusted -  adjustednocouple
# Q: not coming up with adjusted/unadjustd R2 - suspect this is lmm vs lm. Not clear 
# how couple intercept acheived without this.

# Simulated data
shcs_data_int <- shcs_data_long %>%
  select(-starts_with('ID'), -starts_with('SpVL'), -contains('cohortmean'), -log10_SpVL) %>%
  mutate(across(where(is.factor), .fns = ~ match(.x, unique(.x))))

cum_probs <- shcs_data_int %>% 
  select(where(is.integer)) %>% 
  apply(.,2, function(x) cumsum(table(x))/length(x))

sim_data <- mvtnorm::rmvnorm(n = 100, mean = as.vector(colMeans(shcs_data_int)), sigma = cov(shcs_data_int)) %>%
  as_tibble(name_repair = NULL) %>%
  setNames(colnames(shcs_data_int)) %>%
  mutate(across(c(partner,riskgroup,sex), 
                .fns = ~cut(.x,
                            right = T,
                            include.lowest = T,
                            breaks = as.vector(c(min(.x), quantile(.x,cum_probs[[cur_column()]]))),
                            labels = 1:length(cum_probs[[cur_column()]])))) %>%
  mutate(sex = case_when(sex == 1 ~ 'M',
                         sex == 2 ~ 'F')) %>%
  mutate(riskgroup = case_when(riskgroup == 1 ~ "HET",
                               riskgroup == 2 ~ "MSM",
                               riskgroup == 3 ~ "PWID",
                               riskgroup == 4 ~ "UNKNOWN",
                               riskgroup == 5 ~ "OTHER")) %>%
  mutate(across(where(is.character), .fns = ~as.factor(.x))) %>%
  mutate(sex = relevel(sex, 'M')) %>%
  mutate(riskgroup = factor(riskgroup, levels = c('PWID', 'HET', 'MSM', 'OTHER', 'UNKNOWN'))) %>%
  filter(min(shcs_data_long$log10_SpVL_couplemean)<log10_SpVL_couplemean && log10_SpVL_couplemean<max(shcs_data_long$log10_SpVL_couplemean))

preds <-predict(adjusted_model, newdata = sim_data, re.form = ~(1|log10_SpVL_couplemean)) #newlevels detected in data - current suspicion is random effects?
summary(adjusted_model)$coefficients[, 1]


