# Function that implements the Transmission Model (from Thompson et al., 2019) on 
# simulated data, and subsequently predicts CD4 decline using the Tolerance Model
# (from Regoes et al., 2014).

RunTM4Sim <- function(data, modeltype = 'linear'){
  
  # Sanity check
  stopifnot('cd4_model' %in% ls(envir = .GlobalEnv))
  stopifnot('populationmodel_acrossVL_Environment' %in% ls(envir = .GlobalEnv))
  
  
  # Run transmission model in parallel on simulated data
  results <- RunParallel(populationmodel_acrossVL_Environment, data['transmitter']) %>%
    do.call(cbind.data.frame, .) %>% 
    t() %>%
    cbind.data.frame(recipient = data['recipient'], .) %>%
    rename(probTransmissionPerSexAct_average = probTransmissionPerSexAct) %>%
    rename_with(function(x) gsub('\\.(?<=\\V)', '\\1_average.\\2',x, perl= TRUE), 
                .cols = !contains('chronic') & contains('V')) %>%
    
    # Long Format
    pivot_longer(cols = c(contains('chronic'),contains('average')), 
                 names_to = c('type', 'stage', 'variants'), 
                 names_pattern = "(^[^_]+)(_[^.]*)(.*)", 
                 values_to = 'p') %>% 
    mutate(stage = gsub('^.*_', '', stage)) %>%
    mutate(variants = str_remove_all(variants,'[:alpha:]|[:punct:]') %>% as.numeric()) %>%
    mutate(model = modeltype) %>%
    
    # Predict CD4 decline using fitted CD4 model
    mutate(cd4_decline = predict(cd4_model, newdata = data.frame(spVL = log10(recipient)))) # May need to be changed according to model covariates and data preprocessing
  
  
  return(results)
}

