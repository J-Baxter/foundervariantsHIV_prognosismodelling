# Sim Empirical Results presentation

# a) Effect size (SpVL) Maybe do as statistical model accounting for repeat observations
Get

# b) Survival Analysis (CD4+)

SimCD4_decline <- function(rate, initcd4){
  
  cd4_count <- initcd4
  time <- 0
  
  while(cd4_count>=350){
    cd4_count = cd4_count - abs(rate)*60
    time = time + 60
  }
  
  return(time)
}


IsItAIDSoClockYet <- function(dataframe, timesteps = seq(100, 15000, by = 100)){
  AIDSONSET <- dataframe$AIDS_onset
  
  out <- purrr::map(.x = AIDSONSET, function(x) as.numeric(x > timesteps)) %>% 
    unlist() %>% 
    matrix(nrow = length(timesteps)) %>% 
    t() %>% 
    apply(., 2, mean)
  
  return( t(out))
}


# Calculate proportion of individuals with a cd4 count greater than 350 over time
GetCD4Survival <- function(data, replicates = 100, quantiles = c(0.01, 0.5, 0.99)){
  
  data <- as_tibble(data)
  stopifnot(length(unique(data$riskgroup_recipient)) == 1, 
            'More than one riskgroup present.')
  stopifnot(length(unique(data$transmitter_allocation)) == 1, 
            'More than one transmitter allocation method.')
  
  if(!('p_mv' %in% colnames(data))){
    data$p_mv <- 1 - data$p_variants_1
  }
  
  # Generate founder variant multiplicity assignment for each bootstrap replicate
  replicate_classes <- t(replicate(replicates, rbinom(length(data$p_mv), 1, data$p_mv)))
  
  # Infer temporal change in CD4 count per-patient
  data_cd4 <- data %>%
    rowwise() %>%
    mutate(AIDS_onset = SimCD4_decline(rate =delta_CD4_recipient, initcd4 = 1000)) %>%
    as_tibble()
    
    # Combine with data and split by multiplicity
  data_replicates_split <- apply(replicate_classes, 1, function(x) data_cd4  %>%
                                   as_tibble() %>% 
                                   bind_cols(guide = x) %>%
                                   mutate(guide = factor(guide)) %>% 
                                   split(., .['guide'])) %>%
    
    setNames(., 1:length(.)) %>%
    lapply(., setNames, c('single', 'multiple')) %>%
    list_flatten() 
  
  # At time t, calculate the proportion of individuals with <350 CD4
  replicates_aidsonset <- lapply(data_replicates_split, IsItAIDSoClockYet) %>%
    do.call(rbind.data.frame,  .) %>%
    as_tibble() %>%
    mutate(replicate = names(data_replicates_split)) %>% 
    rename_with( ~ tolower(gsub("V", "", .x, fixed = TRUE)))
  
  # Format quantile function
  quantile_df <- function(x, probs = quantiles) {
    tibble(
      val = quantile(x, probs, na.rm = TRUE)
    )
  }
  
  # Summarise replicates (calculate quantiles)
  replicates_aidsonset_summary <- replicates_aidsonset %>%
    separate_wider_delim(replicate, '_', names = c('replicate', 'multiplicity')) %>%
    group_by(multiplicity) %>% 
    reframe(across(1:150, quantile_df, .unpack = T)) %>%
    mutate(quant = rep(quantiles, nrow(.)/length(quantiles))) %>%
    pivot_longer(contains('val'), names_to = c('time', 'type'), names_sep = '_', values_to = 'p') %>%
    select(-type) %>%
    pivot_wider(names_from = quant, values_from = p) %>%
    mutate(time = as.integer(time)*100) # return time to original scale
  
   return(replicates_aidsonset_summary)
}


test_function <-  bind_rows(results[3:5]) %>% 
  filter(transmitterallocation == 'ML') %>% 
  filter(riskgroup_recipient == 'MF') %>% 
  GetCD4Survival(.)
  
 
 ggplot(test_function) +
   geom_line(aes(x = time, y= `0.5`, colour = multiplicity)) +
   geom_ribbon(aes(x = time, ymin = `0.01`, ymax = `0.99`, fill = multiplicity), alpha = 0.5)+
   scale_x_continuous('Days Post Infection') + 
   scale_y_continuous('Proportion of Cohort with < 350 CD4 mm3') +
   coord_cartesian(xlim = c(0,365*10))+
   my_theme
 

## Pure death process (deprecated)
#V <- list(birth=1,death=-1)
#beta <- 00
#mu <- 0.195
#N0 <- 1000
#nevent <- 101
#T <- numeric(nevent)
#N <- numeric(nevent) 
#N[1] <- N0
#k <- 1
#while (T[k] < Inf && k <= nevent) {
##  alpha <- c(beta,mu)*N[k] 
#  Lambda <- sum(alpha)
#  if (Lambda > 0) {
#    tau <- rexp(n=1,rate=Lambda)
#    event <- sample.int(n=length(V),size=1,prob=alpha/Lambda)
#    N[k+1] <- N[k]+V[[event]]
#  } else {
#    tau <- Inf
#  }
#  T[k+1] <- T[k]+tau
#  k <- k+1 }
#T <- head(T,k-1)
#N <- head(N,k-1)

#plot(T,N)
