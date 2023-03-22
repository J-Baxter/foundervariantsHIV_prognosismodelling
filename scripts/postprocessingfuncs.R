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
  
  out <- transmitter_tm %>%
    lapply(., cbind.data.frame) %>%
    do.call(rbind.data.frame,.) %>%
    rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
    rowwise() %>%
    mutate(p_virion = sum(c_across(starts_with('V')))) %>%
    select(transmitter, p_virion, nparticles) %>%
    group_by(transmitter) %>%
    mutate(exp.virions = sample(1:33, 1, replace = T, prob = unlist(p_virion))) %>%
    ungroup() %>%
    mutate(recipient = rep(pop[['recipient']], each = 33))  %>%
    
    # Annotation
    mutate(model = 'linear') %>%
    
  return(out)
}


