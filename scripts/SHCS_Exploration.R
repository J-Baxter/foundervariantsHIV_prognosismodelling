#Swiss Cohort Analysis
source('./scripts/dependencies.R')

shcs_data <- read_csv("data/shcs_data.csv") %>% 
  rowid_to_column( "ID.pair") %>%
  mutate(across(contains('spVL'), ~ raise_to_power(10,.x))) %>%
  rename_with( ~ stri_replace_last_fixed(.x, '.', '_')) %>%
  rename_with( ~ stri_replace_last_fixed(.x, 'spVL', 'SpVL'))


# Transform to long-format data table
data_model <- shcs_data %>% 
  pivot_longer(cols = -c(ID_pair), names_to = c(".value", "partner"), names_sep = "_" )


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
shcs_plt1c <- ggplot(data_model) + 
  geom_histogram(aes(x = age.inf, fill = partner)) + 
  scale_y_continuous(expand = c(0.01,0.01), name = 'Count') + 
  scale_x_continuous(expand = c(0.01,0.01), name = 'Age at Infection',
                     breaks = seq(0,12, by = 1)) +
  scale_fill_manual(values = c('#ef654a','#fdd49e'), name = 'Partner')+
  my_theme

shcs_panel <- plot_grid(shcs_plt1a, shcs_plt1b, shcs_plt1c, align = 'hv', labels = 'AUTO', nrow = 1)
ggsave(plot = shcs_panel, filename = paste(figs_dir,sep = '/', "shcs_panel.jpeg"), device = jpeg, width = 14, height = 6) #Within-host dynamics -functional
