data <- modelresults %>% filter(transmitterallocation == 'ML') %>% group_split(dataset_id)
SimCD4_decline <- function(rate, initcd4){
  
  cd4_count <- initcd4
  time <- 0
  
  while(cd4_count>=350){
    cd4_count = cd4_count - abs(rate)*60
    time = time + 60
  }
  
  return(time)
}
GetCD4Diff <- function(data, replicates = 1000){
  data <- as_tibble(data)
  rg <- data %>% 
    select(dataset_id) %>%
    unique() %>% 
    as.vector()
  
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
  cd4_replicates_diff <- apply(replicate_classes, 1, function(x) data %>%
                                 as_tibble() %>% 
                                 bind_cols(guide = x) %>%
                                 mutate(guide = factor(guide)) %>% 
                                 split(., .['guide'])) %>%
    setNames(., 1:length(.)) %>%
    lapply(., setNames, c('single', 'multiple')) %>%
    list_flatten() %>%
    lapply(., dplyr::select, 'delta_CD4_recipient') %>%
    
    # Formatted dataframe
    bind_rows(., .id = 'replicate') %>%
    separate_wider_delim(replicate, '_', names = c('replicate', 'multiplicity')) %>%
    rowwise() %>%
    mutate(aidsonset = SimCD4_decline(delta_CD4_recipient, initcd4 = 1000)) %>%
    as_tibble() %>%
    
    # Grouped means
    summarise(delta_CD4_recipient = mean(delta_CD4_recipient), aidsonset = mean(aidsonset), .by = c('replicate', 'multiplicity')) %>%
    pivot_wider(., names_from = multiplicity, values_from = c('delta_CD4_recipient', 'aidsonset')) %>%
    
    # Calculate unadjusted difference
    mutate(dif = delta_CD4_recipient_multiple - delta_CD4_recipient_single) %>% 
    rename(deltaCD4recipient_multiple = delta_CD4_recipient_multiple) %>%
    rename(deltaCD4recipient_single = delta_CD4_recipient_single) %>%
    pivot_longer(cols = contains(c('deltaCD4', 'aidsonset')), 
                 names_to = c(".value", "multiplicity"), 
                 names_sep = "_"
                 ) %>%
    mutate(riskgroup = rg)
  
  
  return(cd4_replicates_diff)
  
}


test <- lapply(data, GetCD4Diff) %>% bind_rows()

test %>%
  summarise(med_dif = median(dif),
            med_delta = median(deltaCD4recipient), 
            sd_delta = sd(deltaCD4recipient), 
            .by = c('riskgroup', 'multiplicity'))

cd4_tble <- test %>%
  summarise(med_dif = median(dif),
            med_delta = median(deltaCD4recipient), 
            med_time = median(aidsonset),
            sd_time = sd(aidsonset), 
            .by = c('riskgroup', 'multiplicity')) %>%
  mutate(med_time = med_time/365) %>%
  mutate(sd_time = sd_time/365) %>%
  select(multiplicity, sd_time, riskgroup) %>%
  pivot_wider(names_from = multiplicity, values_from = sd_time) %>%
  
  
  t.test(cd4_tble$multiple*365, cd4_tble$single*365, alternative = 'g', paired = T)
  





