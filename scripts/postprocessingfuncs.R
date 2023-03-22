# Functions for result post-processing

# Returns marginal probabilites and expected number of variants initiating infection
VariantPP <- function(preds, pop){
  out <- preds %>%
    lapply(., cbind.data.frame) %>%
    do.call(rbind.data.frame,.) %>%
    rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
    dplyr::select(-contains('nparticles')) %>%
    dplyr::summarise(across(starts_with('V'), .fns =sum),.by = c(transmitter, w)) %>%
    select(starts_with('V')) %>%
    
    #Infer expected number of variants, given probability distribution for each transmitter
    #mutate(exp.var = apply(.,1, function(x) sample(1:max(length(x)), 1, replace = T, prob = unlist(x)))) %>%
    mutate(exp.var = apply(.,1, function(x) sample(1:max(length(x)), 1, replace = T, prob = unlist(x)))) %>%
    
    mutate(transmitter = pop[['transmitter']], recipient = pop[['recipient']]) %>%
    
    #Pivot 
    pivot_longer(cols = starts_with('V'),
                 names_to = 'variants',
                 values_to = 'p') %>%
    mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% 
             as.numeric())  %>%
    
    #Filter highly unlikely events
    filter(variants <= 12) %>%
    
    # Annotation
    mutate(model = 'linear') 
  
  return(out)
}



# Returns marginal probabilites and expected number of virions initiating infection
VirionPP <- function(preds, pop){
  
  out <- preds %>%
    lapply(., cbind.data.frame) %>%
    do.call(rbind.data.frame,.) %>%
    rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
    dplyr::select(-starts_with('V')) %>% 
    pivot_wider(names_from = nparticles, names_prefix = 'virions_', values_from = prob_nparticles) %>%
    select(contains('virions_')) %>%
    
    #Infer expected number of virions, given probability distribution for each transmitter
    mutate(exp.virions = apply(.,1, function(x) sample(1:max(length(x)), 1, replace = T, prob = unlist(x)))) %>%
    mutate(transmitter = pop[['transmitter']], recipient = pop[['recipient']]) %>%
    
    #Pivot 
    pivot_longer(cols = starts_with('v'),
                 names_to = 'virions',
                 values_to = 'p') %>%
    mutate(virions = str_remove_all(virions,'[:alpha:]|[:punct:]') %>% 
             as.numeric())  %>%
    
    # Annotation
    mutate(model = 'linear') %>%
    
  return(out)
}


