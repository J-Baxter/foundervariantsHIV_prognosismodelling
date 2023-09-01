library(sn)
library(fGarch)

# sn
data <- rsn(1000, xi=-5, omega=2, alpha=5)
# fGarch
data <- rsnorm(1000, mean = 0, sd = 1, xi = 1.5)


set.seed(123)
data <- rsn(1000, xi=-5, omega=2, alpha=5)
# Fit curve
fit <- selm(data ~ 1, family = "SN")
plot(fit)

fit_sn@param$dp

t <- shcs_data_long_transmitterML%>%  select(log10_SpVL) %>% unlist() %>% as.numeric() #filter(partner == 'transmitter') %>%
fit_sn <- selm(t ~ 1, family = "SN")
fit_n <- fitdistrplus::fitdist(t, 'norm')

summary(fit_n)

ggplot(data = data.frame(x = t), aes(x = t)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.25) +
  stat_function(fun = dsnorm, args = c(mean = 4.27463, sd =  0.6886 , xi = 0.5154 ), colour = 'red', linewidth = 1.5) + 
  stat_function(fun = dnorm,  args = c(mean = 4.274059, sd = 0.687981), colour = 'blue', linewidth = 1.5) +
  stat_function(fun = dnorm,  args = c(mean = 4.74, sd = 0.78), colour = 'darkgreen', linewidth = 1.5) +
  scale_x_continuous(limits = c(1,8), breaks = seq(1,8, by =1))+
  my_theme

# relationship between population distribution and distribution of means?
mu <- shcs_data_long_transmitterML %>%  select(log10_SpVL_couplemean) %>% unlist() %>% as.numeric() %>% unique() #filter(partner == 'transmitter') %>%

ggplot(data = data.frame(x = mu), aes(x = mu)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.25) +
  stat_function(fun = dsnorm, args = c(mean = 4.27463, sd =  0.6886 , xi = 0.5154 ), colour = 'red', linewidth = 1.5) + 
  stat_function(fun = dnorm,  args = c(mean = 4.274059, sd = 0.687981), colour = 'blue', linewidth = 1.5) +
  stat_function(fun = dnorm,  args = c(mean = 4.74, sd = 0.78), colour = 'darkgreen', linewidth = 1.5) +
  scale_x_continuous(limits = c(1,8), breaks = seq(1,8, by =1))+
  my_theme

# Can we assume that the distribution of vls follows the parms specified for Zambia?: adjust variance in cov matrix (diag) for couple SpVL
test_cov <- cov(shcs_data_int_list$HET) 
test_cov[4,4] <- 0.6084

