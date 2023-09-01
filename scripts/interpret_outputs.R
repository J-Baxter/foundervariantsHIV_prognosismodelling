# Sim Empirical Results presentation

# a) Effect size (SpVL) Maybe do as statistical model accounting for repeat observations
#
GetSpVLDiff <- function(data, replicates = 100){
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
    lapply(., select, 'log10_SpVL_recipient') %>%
    
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

  
# b) Survival Analysis (CD4+)
# No point in stochastic analysis as only variation will be exponential time between state changes
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
GetCD4Survival <- function(data, replicates = 100, quantiles = c(0.01, 0.5, 0.99)){
  
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
  coord_cartesian(xlim = c(0,365*10))+ #cut at 10 years
  my_theme +
  theme(legend.position = 'right')

# Log rank or CPH - 
iskgroup_recipient)) == 1
length(unique(data$transmitterallocation)) == 1
z=data.frame( p_variants_1 = 1-rep(seq(0.25, 0.32, length.out = 27), 2), delta_CD4_recipient= c(seq(-0.31, -0.18, by = 0.005), seq(-0.3, -0.17, by = 0.005)),
            riskgroup_recipient = 'N',
            transmitterallocation = 'X')


z=data.frame( delta_CD4_recipient= seq(-0.31, -0.18, by = 0.005), 
              p_variants_1 = 1-c(0.25,0.255, 0.26, 0.265, 0.26, 0.27,0.275, 0.28, 0.275, 0.28, 0.30, 0.29, 
      0.305, 0.31, 0.315, 0.31,0.30, 0.315,0.31, 0.305,0.31, 0.315,0.30, 0.32,0.32, 0.32,0.32),
      riskgroup_recipient = 'N',
      transmitterallocation = 'X') %>%
  bind_rows(.,.,.) %>%
  mutate(p_variants_1 = jitter(p_variants_1, factor = 3)) %>%
  mutate(delta_CD4_recipient= jitter(delta_CD4_recipient, factor = 15))


z_func <- GetCD4Survival(z, replicates = 100)


z_plota <- ggplot(z) +
  geom_point(aes(x = 1-p_variants_1 , y= delta_CD4_recipient ),colour = '#ef654a' ,shape = 4) +
  scale_x_continuous(name = 'P(Multiple Variant Recipient)',
                     expand = c(0.02,0.02),
                     limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1))+
  scale_y_continuous(name = expression(paste(Delta, ' CD4+ ', mu, l**-1, ' ', day**-1)),  #
                     expand = c(0,0),
                     limits = c(-0.5,0),
                     axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
                     axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 8),
                     axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
                     strip.text  = element_text(size = 8),
                     legend.text = element_text(size = 7),
                     legend.title = element_text(size = 8),) +
  my_theme 

z_plotb <- ggplot(z_func) +
  geom_line(aes(x = time, y= `0.5`, colour = multiplicity)) +
  geom_ribbon(aes(x = time, ymin = `0.01`, ymax = `0.99`, fill = multiplicity), alpha = 0.5)+
  scale_x_continuous('Days Post Infection') + 
  scale_y_continuous('Proportion of Cohort with < 350 CD4 mm3') +
  scale_color_manual(values = c("#fdbb84","#ef6548")) +
  scale_fill_manual(values = c("#fdbb84","#ef6548")) +
  coord_cartesian(xlim = c(0,365*10))+ #cut at 10 years
  my_theme +
  theme(legend.position = c(0.9,0.7),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 8),
        axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
        strip.text  = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),)

ggsave(plot = cowplot::plot_grid(z_plota, z_plotb , labels = 'AUTO', nrow = 1, align = 'hv', label_size = 10) , filename =  "zplot.jpeg", 
       device = jpeg,  width = 180, height = 80,  units = 'mm')

ggplot(z) + geom_point(aes(x=x, y=y))



