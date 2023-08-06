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
  pivot_longer(cols = -c(ID_pair, contains('couplemean')), names_to = c(".value", "partner"), names_pattern  = "^(.*)_([0-9])$") %>% #matches('SpVL_[[:digit:]]')
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

# SpVL ~ Riskgroup/Sex

shsc_plt2a <- ggplot(shcs_data_long,aes(x = riskgroup, y = SpVL_couplemean)) +
  geom_boxplot(fill = '#ef654a') + 
  scale_y_log10(name = expression(paste("SpVL Partner 2", ' (', Log[10], " copies ", ml**-1, ')')),
                #limits = c(10**2, 10**7),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x)) + 
  stat_compare_means(method = "kruskal.test",label.x = 2.5, label.y = 6) +
  theme_classic(base_family = "Questrial")+
  theme(
    text = element_text(size=16),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines")
  )

shsc_plt2b <- ggplot(shcs_data_long, aes(x = sex, y = SpVL_couplemean)) +
  geom_boxplot(fill = '#ef654a') + 
  scale_y_log10(name = expression(paste("SpVL Partner 2", ' (', Log[10], " copies ", ml**-1, ')')),
                #limits = c(10**2, 10**7),
                expand = c(0.05,0),
                breaks = trans_breaks("log10", function(x) 10**x)) + 
  stat_compare_means(method = "wilcox.test",   label.x = 1.5, label.y = 6) + 
  
  theme_classic(base_family = "Questrial")+
  theme(
    text = element_text(size=16),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines")
  )

shcs_panel2 <- plot_grid(shcs_plt2a, shcs_plt2b, align = 'hv', labels = 'AUTO', nrow = 1)
ggsave(plot = shcs_panel2, filename = paste(figs_dir,sep = '/', "shcs_panel2.jpeg"), device = jpeg, width = 14, height = 6) #Within-host dynamics -functional
