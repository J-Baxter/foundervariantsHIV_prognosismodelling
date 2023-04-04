# Functions for result post-processing
# test: all probs sum to 1

# Returns marginal probabilities and expected number of variants initiating infection
VariantPP <- function(preds, pop){
  
  reps <- length(preds)/nrow(pop)

  out <- preds %>%
    # format to dataframes
    lapply(., cbind.data.frame) %>%
    do.call(rbind.data.frame,.) %>%
    rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
    dplyr::select(-contains('nparticles')) %>%

    # Sum for marginal probability (Variant) for any number of virions
    dplyr::summarise(across(starts_with('V'), .fns =sum),.by = c(transmitter, w)) %>%
    
    # For each SPVL, simulate (sample) the number of variants initiating infection according to 
    # estimated probability distribution
    #mutate(exp.var = apply(test_df [,grepl( "V" , names( test_df ) ) ],1, function(x) sample(1:max(length(x)), 1, replace = T, prob = unlist(x)))) %>%
    mutate(recipient = rep(pop[['recipient']], each = reps )) %>%

    #Pivot 
    pivot_longer(cols = starts_with('V'),
                 names_to = 'variants',
                 values_to = 'p') %>%
    mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% 
             as.numeric())  %>%
    
    #Filter highly unlikely events
    filter(variants <= 12) %>%
    
    # Get probabilities for binned vls
    mutate(recipient_rounded = 10**(round(log10(recipient)/0.5)*0.5))%>%
    mutate(transmitter_rounded = 10**(round(log10(transmitter)/0.5)*0.5)) %>%
    mutate(mean_p = mean(p), .by = c(recipient_rounded, variants)) 
  
  return(out)
}



# Returns marginal probabilities and expected number of virions initiating infection
VirionPP <- function(preds, pop){
  
  reps <- length(preds)/nrow(pop)
  maxvirions <- 33
  
  out <- preds %>%
    # format to dataframes
    lapply(., cbind.data.frame) %>%
    do.call(rbind.data.frame,.) %>%
    rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
    rowwise() %>%
    
    # Sum for marginal probability (Virions) for any number of variants
    mutate(p_virion = sum(c_across(starts_with('V')))) %>%
    select(transmitter, p_virion, nparticles, w) %>%
    group_by(transmitter, w) %>%
    
    # For each SPVL, simulate (sample) the number of virions initiating infection according to 
    # estimated probability distribution
    #mutate(exp.virions = sample(1:maxvirions, 1, replace = T, prob = unlist(p_virion))) %>%
    ungroup() %>%
    mutate(recipient = rep(pop[['recipient']], each = maxvirions*reps )) %>%
    
    # Get probabilities for binned vls
    mutate(recipient_rounded = 10**(round(log10(recipient)/0.5)*0.5))%>%
    mutate(transmitter_rounded = 10**(round(log10(transmitter)/0.5)*0.5)) %>%
    mutate(mean_p_virion = mean(p_virion), .by = c(recipient_rounded, nparticles)) 
    
  return(out)
}


