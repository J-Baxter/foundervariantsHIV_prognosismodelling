# Under what conditions might we expect to observe an association between multiple variant infection and SpVL
#
# This script follows the linear equations set out by Bonhoeffer et al. 2015 that describes the 
# changes of the distribution of set point viral load in the HIV carrier population over a full 
# cycle of transmission
# https://doi.org/10.1371/journal.ppat.1004634
#
# We extend the linear equations to decompose the recipient distribution of SpVL into distributions for
# single and multple variant infections, leveraging the effect sizes and means calculated by Janes et al. 
# 2015 
# https://doi.org/10.1038/nm.3932 


################################### Dependencies ###################################
require(tidyverse)
require(parallel)
require(epitools)


# Set CPU cluster size
cl <- detectCores() %>% `-` (3)


# Set seed
set.seed(4472)


# Assume recipient phenotype distribution is a mixture of two gaussian components where
# the difference in means is set by an empirically determined effect size
DecomposeRecipientSpVL <- function(recipient_mean, recipient_var, p_mv, effectsize){
  
  # recipient_mean = (1 - p_mv)*x + p_mv*(x + effect_size)
  # Because weights of finite mixture must sum to 1; rearranges to:
  sv_mean <- recipient_mean - p_mv*effectsize  
  mv_mean <- sv_mean + effectsize
  
  
  #recipient_var <- p_m*(y + mv_mean^2) + (1-p_mv)*(y + sv_mean^2) - recipient_mean^2
  # Rearranges for the same reasons
  # Assumption of equal variances justified by Janes et al.
  decomposed_var <- recipient_var + recipient_mean^2 - p_mv * mv_mean^2 - (1 - p_mv) * sv_mean^2
  
  out <- list('sv' = c('mean' = sv_mean, 'var' = decomposed_var),
              'mv' = c('mean' = mv_mean, 'var' = decomposed_var))
  
  return(out)
}


# For a specified effect and sample size, calculate the expected proportion of times we would 
# observe a significant difference between single and multiple founder variant infections
SimEffectSizePMV <- function(n, e , specifyPMV = FALSE){
  
  if (all(is.logical(specifyPMV))){
    p_mv <- (2:(n-2))/n
    p <- p_mv[p_mv<=0.6]
  }else{
    p <- specifyPMV
  }
 
  m <- matrix(data = NA, nrow = length(e), ncol = length(p))
  
  for (i in 1:length(e)){
    
    for(j in 1:length(p)){
      spvls <- DecomposeRecipientSpVL(effectsize = e[i],
                                      p_mv = p[j],
                                      recipient_mean = 4.74,
                                      recipient_var = 0.6084)
      sig.count <- 0
      
      for (z in 1:1000){
        sv <- rnorm(n*(1-p[j]), mean = spvls$sv['mean'], sd = sqrt(spvls$sv['var']))
        mv <- rnorm(n*p[j], mean = spvls$mv['mean'], sd = sqrt(spvls$mv['var']))
        
        if(t.test(sv, mv, var.equal = T)$p.value <= 0.05){
          sig.count = sig.count + 1
        }
        m[i,j] <- sig.count/1000
      }
    }
  }
  
  out <- m %>%
    reshape2::melt() %>%
    mutate(p_mv = rep(p, each = length(e))) %>%
    mutate(effect_size = rep(e,  length(p))) %>%
    select(-contains('Var')) %>%
    mutate(sample_size = n)
  return(out)
}


################################### Calculate Recipient Dist ###################################
zambia_variance <- 0.78**2 #0.78 is the standard deviation
zambia_mean <- 4.74 

recipient_dist <- CalcRecipient(zambia_mean, 
                                zambia_variance)


################################### Decompose Recipient Distribution ###################################
decomp_vl <- DecomposeRecipientSpVL(recipient_mean = recipient_dist[['mean']], 
                                    recipient_var = recipient_dist[['var']],
                                    p_mv =  0.25, 
                                    effectsize = 0.3)


vls <- tibble(recipient_mean = rnorm(100000, 4.74, 0.6084), 
              recipient_mv = rnorm(100000, decomp_vl$mv['mean'], sqrt(decomp_vl$mv['var'])),
              recipient_sv = rnorm(100000, decomp_vl$sv['mean'], sqrt(decomp_vl$sv['var']))) %>%
  pivot_longer(cols = everything(), names_to = 'stage', values_to = 'vl')


################################### Simulate proportion of significant observations ###################################
# Vector of effect and sample sizes
effect_size <- seq(0.01, 1, by = 0.01)
sample_size <- c(25,50,100,200)

simsignificances <- mclapply(sample_size,
                             SimEffectSizePMV, 
                             e = effect_size,
                             mc.cores = cl,
                             mc.set.seed = FALSE) %>%
  bind_rows()


################################### Compare Riskgroups ###################################

p_mv <- c(0.21, 0.13, 0.30)

comparesignificances <- SimEffectSizePMV(66, 0.3, specifyPMV = p_mv) %>%
  rowwise() %>%
  mutate(ci.lower = prop.test(value*100000, 100000, conf.level = .95, correct = F)$conf.int[1])%>%
  mutate(ci.upper = prop.test(value*100000, 100000, conf.level = .95, correct = F)$conf.int[2])

or_matrix <- matrix(c(1-comparesignificances$p_mv,comparesignificances$p_mv), nrow = 3) *100000

outcome <- c('single', 'multiple')
riskgroup <- c('MF', 'FM', 'MM')
dimnames(or_matrix) <- list('riskgroup' = riskgroup, 'outcome' = outcome)

epitools::oddsratio(or_matrix)

# END #
