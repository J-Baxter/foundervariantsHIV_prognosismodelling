# Create free-scale network from power-law distribution
InitNetwork <- function(n, exponent, k_max){
  require(tidyverse)
  require(fastnet)
  require(igraph)
  
  exponent = as.integer(exponent)
  
  net <-  net.random.plc(n, k_max, exponent) %>% 
    graph_from_adj_list() %>%
    as.undirected()
  
  return(net)
}

# create networks assuming different epi characteristics
msm_net <- InitNetwork (10000, 250, 2) 
hsx_net <- InitNetwork (10000,  60, 3) 

# Simulate HIV epidemic over networks (split infectious period into primary, chronic and pre-aids)


# plot epidemic curve
test_msm <- net.random.plc(10000, 250, 2) %>%  graph_from_adj_list() %>% intergraph::asNetwork()
ergm(test_msm ~ edges)
test_msm_2 <- 

# distribution of infections by stage of transmitter(?)



