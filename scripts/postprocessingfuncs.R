# Functions for result post-processing
# test: all probs sum to 1

# Returns marginal probabilities and expected number of variants initiating infection
VariantPP <- function(preds, pop){
  
  #Suggest refactoring as per VirionPP
  reps <- length(preds)/nrow(pop)

  out <- preds %>%
    lapply(., cbind.data.frame) %>%
    do.call(rbind.data.frame,.) %>%
    rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
    dplyr::select(-contains('nparticles')) %>%
    dplyr::summarise(across(starts_with('V'), .fns =sum),.by = c(transmitter, w)) %>%
    select(starts_with('V')) %>% #is this necessary
    
    #Infer expected number of variants, given probability distribution for each transmitter
    #mutate(exp.var = apply(.,1, function(x) sample(1:max(length(x)), 1, replace = T, prob = unlist(x)))) %>%
    mutate(exp.var = apply(.,1, function(x) sample(1:max(length(x)), 1, replace = T, prob = unlist(x)))) %>%
    
    mutate(transmitter = rep(pop[['transmitter']], reps) , recipient = rep(pop[['recipient']],reps), w = rep(c(1,5,10,20), each = nrow(pop))) %>%
    
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



# Returns marginal probabilities and expected number of virions initiating infection
VirionPP <- function(preds, pop){
  
  reps <- length(preds)/nrow(pop)
  maxvirions <- 33
  
  out <- transmitter_timing  %>%
    lapply(., cbind.data.frame) %>%
    do.call(rbind.data.frame,.) %>%
    rename_with(function(x) gsub('variant.distribution.', '', x)) %>%
    rowwise() %>%
    mutate(p_virion = sum(c_across(starts_with('V')))) %>%
    select(transmitter, p_virion, nparticles, w) %>%
    group_by(transmitter, w) %>%
    mutate(exp.virions = sample(1:maxvirions, 1, replace = T, prob = unlist(p_virion))) %>%
    ungroup() %>%
    mutate(recipient = rep(pop[['recipient']], each = maxvirions*reps )) %>%
    
    # Annotation
    mutate(model = 'linear') %>%
    
  return(out)
}


