# fix the coefficient of determination
R2 <- 0.2


zambia_variance <- 0.61**2 #0.61 is the standard deviation
zambia_mean <- 4.74 


# simulate x and y variables
test <- function(n, mean, variance, h2){
  
  repeat {
    x <- rnorm(n, mean, sqrt(variance))
    ssr <- sum((x-mean(x))^2) # sum of squared residuals
    e <- rnorm(n)
    e <- resid(lm(e ~ x))
    e <- e*sqrt((1-h2)/h2*ssr/(sum(e^2)))
    y <- x + e
    
    if(all(y > 0 && y <8))
      break
  }

  return(cbind.data.frame(y,x,h2))
}



test <- function(n, mean, variance, h2){
  
  repeat {
    x <- rnorm(n, mean, sqrt(variance))
    ssr <- sum((x-mean(x))^2) # sum of squared residuals
    
    # Residuals
    e <- rnorm(n)
    e <- resid(lm(e ~ x))
    e <- e*sqrt((1-h2)/h2*ssr/(sum(e^2)))
    y <- x + e
    
    if(all(y > 0 && y <8))
      break
  }
  
  return(cbind.data.frame(y,x,h2))
}


varyh2 <- sapply(seq(0.1,0.5, by = 0.1),  test, n = 100, mean = zambia_mean, variance = zambia_variance, simplify = F) %>% 
  bind_rows() %>%
  rowid_to_column() %>%
  pivot_longer(cols = c('x', 'y'), values_to = 'log10_SpVL', names_to = 'role') %>%
  mutate(rowid = factor(rowid)) %>%
  group_by(rowid) %>%
  mutate(pair_log10_SpVL = mean(log10_SpVL)) %>%
  ungroup() %>% filter(h2 == 0.2)

lmm <- lme4::lmer(log10_SpVL ~ pair_log10_SpVL +role+ (1|rowid), data = varyh2) %>% summary()




h2_test <- test(0.22, n = 100, mean = zambia_mean, variance = zambia_variance)

cov(h2_test$x, h2_test$y)/var(h2_test$x)


cov(log10(shcs_data$SpVL_1), log10(shcs_data$SpVL_2))/var(log10(shcs_data$SpVL_1))

cov(log10(shcs_data$SpVL_1), log10(shcs_data$SpVL_2))

sigma <- matrix(c(0.45422050, 0.09996507,0.09996507, 0.45422050),2,2) #cov(cbind.data.frame(shcs_data$log10_SpVL_1, shcs_data$log10_SpVL_2))

test <- MASS::mvrnorm(10000, Sigma =sigma , mu = c(4.290355, 4.257764))

oops <- sapply(seq(-0.4,0.4, by = 0.01), function(x)  MASS::mvrnorm(10000, Sigma =matrix(c(0.45422050, x, x, 0.45422050),2,2) , 
                                                                     mu = c(4.290355, 4.257764)), simplify = F) %>%
  #lapply(., function(x) var(x[,1])) %>%
  lapply(., function(x) cov(x[,1],x[,2])/var(x[,1])) %>%
  unlist() %>% round(., digits = 2)
 
cov(test[,1],test[,2])/var(test[,1]) #Find covariances for these means that allow us to synthesise populations with h2 = 0.05-0.5



seq(0.05, 0.5, by = 0.05)*0.45422050
# Remaining Problems:
# 1) Adjusting the transmitter distribution (drop within-sample predictions; apply sampling procedure to simulated data)

# TASKS:
# 1) final outputs a) CD4 b) SpVL (?)
# 2) brute search 
cor(h2_test$x, h2_test$y)
# check if R2 has the right value

ggplot(varyh2)+
  geom_point(aes(x=x, y=y))+
  my_theme+
  facet_wrap(.~factor(h2))

varyh2_pmv <- parallel::mclapply(round(10**varyh2$x), populationmodel_acrossVL_Env, theta=c(1.765e-06, 0.013296), mc.cores = cl) 

varyh2$pmv <- unlist(varyh2_pmv)

ggplot(varyh2, aes(x=pmv, y=y))+
  geom_point()+
  my_theme+
  facet_wrap(.~factor(h2)) +
  geom_smooth(method='lm', formula= y~x)


ggplot(varyh2, aes(y=pmv, x=x))+
  geom_point()+
  my_theme+
  facet_wrap(.~factor(h2)) +
  geom_smooth(method='lm', formula= y~x)

# What population of SpVls do you need to recover the mean pmv?
