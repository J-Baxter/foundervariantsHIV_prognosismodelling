# Sim Empirical Results presentation

# a) Effect size (SpVL)


# b) Survival Analysis (CD4+)

SimCD4_decline <- function(rate, initcd4){
  
  cd4_count <- initcd4
  time <- 0
  
  while(cd4_count>350){
    cd4_count = cd4_count - abs(rate)*60
    time = time + 60
  }
  
  return(time)
}


IsItAIDSoClockYet <- function(AIDSONSET, timesteps = seq(100, 15000, by = 100)){
  
  out <- purrr::map(.x = AIDSONSET, function(x) as.numeric(x > timesteps)) %>% 
    unlist() %>% 
    matrix(nrow = length(timesteps)) %>% 
    t() %>% 
    
    apply(., 2, mean)
  
  return(out)
}


# Bootstrap which infections are multiple founder (1000 replicate samples)
Replicates <- function(){
  
}


test_results <-  bind_rows(results[3:5]) %>% 
  filter(transmitterallocation == 'ML') %>% 
  filter(riskgroup_recipient == 'MF') %>% 
  mutate(p_mv = 1 -p_variants_1 ) %>%
  rowwise() %>%
  mutate(AIDS_onset = SimCD4_decline(rate =delta_CD4_recipient, initcd4 = 1000 )) 



time_vals <- seq(100, 15000, by = 100)
t <- test_results$AIDS_onset
purrr::map(.x = t , function(x) x > time_vals)
#I(log10_SpVL_recipient**2) +



> test <- 
> plot(time_vals, test)





test_model <- lm(delta_CD4_recipient ~  I(log10_SpVL_recipient**2) + p_mv, data = test_results)

preds <- predict(test_model, data.frame(p_mv = c(0,1,0,1), log10_SpVL_recipient =c(4, 4, 6, 6)), type = 'response' ) %>% 
  as_tibble() %>% 
  mutate(multiplicity = c('S', 'M', 'S', 'M')) %>% 
  mutate(cd4 = 1000) %>% 
  mutate( log10_SpVL_recipient =c(4, 4, 6, 6))

geom_

ggplot(preds)+
  geom_abline(aes(intercept = cd4, slope = value*365,  colour = multiplicity)) + 
  xlim(c(0,20)) + 
  ylim(c(0,1500)) +
  facet_wrap(.~factor(log10_SpVL_recipient))+
  my_theme


V <- list(birth=1,death=-1)
beta <- 00
mu <- 0.195
N0 <- 1000
nevent <- 101
T <- numeric(nevent)
N <- numeric(nevent) 
N[1] <- N0
k <- 1
while (T[k] < Inf && k <= nevent) {
  alpha <- c(beta,mu)*N[k] 
  Lambda <- sum(alpha)
  if (Lambda > 0) {
    tau <- rexp(n=1,rate=Lambda)
    event <- sample.int(n=length(V),size=1,prob=alpha/Lambda)
    N[k+1] <- N[k]+V[[event]]
  } else {
    tau <- Inf
  }
  T[k+1] <- T[k]+tau
  k <- k+1 }
T <- head(T,k-1)
N <- head(N,k-1)

plot(T,N)
