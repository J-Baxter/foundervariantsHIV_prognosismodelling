#Swiss Cohort Analysis
source('./scripts/dependencies.R')

shcs_data <- read_csv("data/shcs_data.csv") %>% 
  rowid_to_column( "ID.pair") %>%
  mutate(across(contains('spVL'), ~ raise_to_power(10,.x)))

ggplot(shcs_data , aes(x = spVL.1, spVL.2)) +
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
  my_theme + 
  theme(legend.position = 'none')

## test models
lm(spVL.2 ~ spVL.1 + sex.2 + riskgroup.2 + age.inf.2 , data = shcs_data) %>% summary()