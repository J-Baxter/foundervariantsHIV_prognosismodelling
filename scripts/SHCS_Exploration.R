#Swiss Cohort Analysis
source('./scripts/dependencies.R')

shcs_data <-read_csv("data/shcs_data.csv", 
                     col_types = cols(sex.1 = col_factor(levels = c("M", "F")), 
                                      sex.2 = col_factor(levels = c("M", "F")), 
                                      riskgroup.1 = col_factor(levels = c("HET", "MSM", "IDU", "OTHER", "UNKNOWN")), 
                                      riskgroup.2 = col_factor(levels = c("HET", "MSM", "IDU", "OTHER", "UNKNOWN"))))%>%
  rowid_to_column( "ID.pair") %>%
  mutate(across(contains('riskgroup'), .fns = ~ recode_factor(.x, 'IDU' = "PWID"))) %>%
  mutate(across(contains('spVL'), ~ raise_to_power(10,.x))) %>%
  rename_with( ~ stri_replace_last_fixed(.x, '.', '_')) %>%
  rename_with( ~ stri_replace_last_fixed(.x, 'spVL', 'SpVL')) 

# Transform to long-format data table
shcs_data_long <- shcs_data %>% 
  pivot_longer(cols = -c(ID_pair), names_to = c(".value", "partner"), names_sep = "_" ) %>%
  group_by(ID_pair) %>%
  mutate(SpVL_couplemean = mean(SpVL)) %>%
  ungroup() %>% 
  mutate(SpVL_cohortmean = mean(SpVL)) %>%
  mutate(across(contains('SpVL'), .fns = ~log10(.x), .names = "log10_{.col}")) %>%
  mutate(partner = factor(partner, levels = c("1", "2")))


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


# Adjusted Model (but one that might work within the context of our framework)
adjusted_direct_model <- lm(log10_SpVL_1 ~ log10_SpVL_2 + sex_1 + age.inf_1 + riskgroup_1, data = shcs_data) 
 
# Adjusted Model - No couple effect
adjustednocouple_model <- lm(log10_SpVL ~ 1 + partner + sex + age.inf + riskgroup, data = shcs_data_long) 

performance::r2(adjusted_model)

# Adjusted -  adjustednocouple
# Q: not coming up with adjusted/unadjustd R2 - suspect this is lmm vs lm. Not clear 
# how couple intercept acheived without this.

# Simulated data
t <- shcs_data_long %>%
  fastDummies::dummy_cols(remove_selected_columns = T) %>%
  select(-starts_with('ID'), -starts_with('SpVL'), -contains('cohortmean')) 

probs <- t %>% 
  select(where(is.integer)) %>% 
  colMeans()

sim_data_b <- mvtnorm::rmvnorm(n = 100, mean = as.vector(colMeans(t)), sigma = cov(t)) %>%
  as_tibble(name_repair = NULL) %>%
  setNames(colnames(t)) %>%
  mutate(across(partner_1:riskgroup_UNKNOWN, .fns = ~cut(.x,
                                                         right = T,
                                                         include.lowest = T,
                                                         breaks = as.vector(c(min(.x), quantile(.x,probs[[cur_column()]]))),
                                                         labels = 1,
                                                         dig.lab = 0) %>% #Appears to be some threshold problem for risk group
                  as.integer() %>%
                  replace_na(0)))  %>%
  pivot_longer(cols = partner_1:riskgroup_UNKNOWN, 
               names_to = c('variable', 'level'),
               names_sep = '_'
               ) %>% 

  filter(value == 1) %>%
  pivot_wider(names_from = variable, 
              values_from = level) 
