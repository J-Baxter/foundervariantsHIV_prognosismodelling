# Sim Empirical Results presentation

###################################################################################################
# Difference between SpVLs
GetSpVLDiff <- function(data, replicates = 1000){
  data <- as_tibble(data)
  
  stopifnot({
    length(unique(data$riskgroup_recipient)) == 1
    length(unique(data$transmitterallocation)) == 1
  })
  
  if(!('p_mv' %in% colnames(data))){ 
    data$p_mv <- 1 - data$p_variants_1
  }
  
  # Generate founder variant multiplicity assignment for each bootstrap replicate
  replicate_classes <- t(replicate(replicates, rbinom(length(data$p_mv), 1, data$p_mv)))
  
  # Combine with data and split by multiplicity
  SpVL_replicates_diff <- apply(replicate_classes, 1, function(x) data %>%
                                  as_tibble() %>% 
                                  bind_cols(guide = x) %>%
                                  mutate(guide = factor(guide)) %>% 
                                  split(., .['guide'])) %>%
    setNames(., 1:length(.)) %>%
    lapply(., setNames, c('single', 'multiple')) %>%
    list_flatten() %>%
    lapply(., dplyr::select, 'log10_SpVL_recipient') %>%
    
    # Formatted dataframe
    bind_rows(., .id = 'replicate') %>%
    separate_wider_delim(replicate, '_', names = c('replicate', 'multiplicity')) %>%
    
    # Grouped means
    summarise(mean_log10_SpVL = mean(log10_SpVL_recipient), .by = c('replicate', 'multiplicity')) %>%
    pivot_wider(., names_from = multiplicity, values_from = mean_log10_SpVL, names_prefix = 'log10_SpVL_') %>%
    
    # Calculate unadjusted difference
    mutate(dif = log10_SpVL_multiple - log10_SpVL_single) %>% 
    pivot_longer(cols = contains('SpVL'), names_to = 'multiplicity', values_to = 'log10_SpVL')
  
  
  return(SpVL_replicates_diff)
  
}

###################################################################################################
# Resample CD4+ T Cell decline according to P(Multiple Variants)
# Plot proportion of individuals with over 350 cells/mm from a fixed starting point

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
    t() 
  
  return(out)
}


# Calculate proportion of individuals with a cd4 count greater than 350 over time
GetCD4Survival <- function(data, replicates = 1000, quantiles = c(0.01, 0.5, 0.99)){
  
  data <- as_tibble(data)

  stopifnot({
    length(unique(data$riskgroup_recipient)) == 1
    length(unique(data$transmitterallocation)) == 1
  })
  
  if(!('p_mv' %in% colnames(data))){
    data$p_mv <- 1 - data$p_variants_1
  }
  
  # Generate founder variant multiplicity assignment for each bootstrap replicate
  replicate_classes <- t(replicate(replicates, rbinom(length(data$p_mv), 1, data$p_mv)))
  
  # Infer temporal change in CD4 count per-patient
  data_cd4 <- data %>%
    rowwise() %>%
    #mutate(AIDS_onset = case_when(recipient_sex == 'F'~ SimCD4_decline(rate = delta_CD4_recipient, initcd4 = rlnorm(1, meanlog = 6.633318, sd = 0.2847)),
                                 # recipient_sex == 'M'~ SimCD4_decline(rate = delta_CD4_recipient, initcd4 = rlnorm(1, meanlog = 6.476972, sd = 0.3389))))
    #rlnorm(1, meanlog = 6.709304, sd = 0.3)
    mutate(AIDS_onset = SimCD4_decline(rate = delta_CD4_recipient, initcd4 = 1000)) %>%
    as_tibble()
  
  #geometric mean ± 2sd (men = 650, 330-1280\\ women = 760, 430-1350) 
  # rlnorm(1000, meanlog = 6.476972, sd = 0.3389) M
  # rlnorm(1000, meanlog = 6.633318, sd = 0.2847) W
    
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
  aidsonset_observations <- lapply(data_replicates_split, IsItAIDSoClockYet)
  
  replicates_aidsonset  <- lapply(aidsonset_observations, function(x) t(apply(x, 2, mean))) %>%
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
      reframe(across(1:150, quantile_df, .unpack = TRUE)) %>%
      mutate(quant = rep(quantiles, nrow(.)/length(quantiles))) %>%
      pivot_longer(contains('val'), names_to = c('time', 'type'), names_sep = '_', values_to = 'p') %>%
      dplyr::select(-type) %>%
      pivot_wider(names_from = quant, values_from = p) %>%
      mutate(time = as.integer(time)*100) # return time to original scale
  
  
   return(replicates_aidsonset_summary)
}


# END #