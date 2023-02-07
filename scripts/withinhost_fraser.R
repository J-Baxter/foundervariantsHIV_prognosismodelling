model <- function(times, state, parms){
  with(as.list(c(parms,state)), {
    
    v = b*y/c
    a_4 = a_0*k_4*(x4/(x4+x_s))
    a_8 = a_0*k_8*(x8/(x8+x_s))
    epsilon = y/(y+y_t)
    
    #ODEs
    dx4 <- gamma_cd4 - a_4 + 2*p_a*mu_a*x_4a + mu*x4 - 
    
    dx8 <- gamma_cd8 - a_8 + 2*p_a*mu_a*x_8a + mu*x8
    
    dx_4a <- a_4 - mu_a*x_4a - beta*v*x_4a
    
    dx_8a <- a_8 - mu_a*x_8a
    
    dy <- (1-f_latent)*beta*v*x_4a + alpha_latent*y_l - (alpha + sigma*epsilon*z)*y
    
    dy_l <- f_l*beta*v*x_4a - alpha_latent*y_l
    
    dz <- d(z_0 - z) + p*epsilon*z
    
    list(c(dx4, dx8, dx4a, dx8a))
    
  })
}