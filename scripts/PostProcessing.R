PostProcess <- function(tmresult, metadata){

  variant_distribtion <- tmresult[['variant_distribution']]%>%
    dplyr::select(-nparticles)
  
  max_particles <- 32
  
  p_particles <- rep(0,max_particles)
  p_variants <- rep(0,max_particles)
  
  for (i in 1:max(dim(variant_distribtion))){
    
    p_particles[i] <- sum(variant_distribtion[i,])
    p_variants[i ] <- sum(variant_distribtion[,i])
      
  }
  
  probabilities <- bind_cols(p_particles= p_particles, p_variants = p_variants, n = 1:max_particles) %>%
    pivot_wider(names_from = n, values_from = c(p_variants, p_particles), names_sep = '_', names_vary = 'fastest')%>%
    bind_rows() 
  
  out <- bind_cols(metadata, probabilities)
  
  return(out)
}
