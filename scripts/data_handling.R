t <- transmitterspvl_variantdist %>% 
  rename(probTransmissionPerSexAct_average = probTransmissionPerSexAct) %>%
  rename_with(function(x) gsub('\\.(?<=\\V)', '\\1_average.\\2',x, perl= TRUE), .cols = !contains('chronic') & contains('V')) %>%
  pivot_longer(cols = c(contains('chronic'),contains('average')), names_to = c('type', 'stage', 'variants'), names_pattern = "()-()-()", values_to = 'p') #1) everything before last _ 2)after last _ and before last . 3) after last .