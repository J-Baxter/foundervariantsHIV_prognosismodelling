ShortenDagArrows <- function(tidy_dag, proportion){
  # Update underlying ggdag object
  tidy_dag$data <- dplyr::mutate(tidy_dag$data, 
                                 xend = (1-proportion/2)*(xend - x) + x, 
                                 yend = (1-proportion/2)*(yend - y) + y,
                                 xstart = (1-proportion/2)*(x - xend) + xend,
                                 ystart = (1-proportion/2)*(y-yend) + yend)
  return(tidy_dag)
}

coords <- list(
  x = c(SpVLr = 10, SpVLt = 2.5, deltaCD4 = 10, PMV = 2.5, PAcq = 1),
  y = c(SpVLr = 10, SpVLt = 10, deltaCD4 = 2, PMV = 2, PAcq = 6)) %>% 
  coords2df() %>%
  coords2list()

test_dag <- dagify(
  deltaCD4 ~ PMV + SpVLr,
  SpVLr ~ PMV + SpVLt,
  PMV ~ SpVLt + PAcq,
  PAcq ~ SpVLt) 

coordinates(test_dag) <- coords

model_cols <- c('Heritability' = '#fdbb84', 'Transmission' ='#d7301f', 'Tolerance' = '#7f0000')
node_labels <- c(expression(P[Acq]), expression(P[MV]), expression(SpVL[R]), expression(SpVL['T']),
                 expression(paste(Delta, ' CD4+')))

 
plt_1a <- test_dag %>%
  tidy_dagitty() %>%
  ShortenDagArrows(proportion = .09) %>%
  mutate(linetype = ifelse(name == "PMV", "dashed", "solid")) %>% 
  arrange(name) %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) + 
  
  geom_rect(xmin = 0.1, xmax = 3.5, ymin = 1, ymax = 11.5, aes(fill = 'Transmission'), colour=NA, alpha = 0.05, size = 2) +
  geom_rect(xmin = 8.5, xmax = 11.5, ymin = 1, ymax = 11.5, aes(fill = 'Tolerance'), colour=NA, alpha = 0.05, size = 2) +
  geom_rect(xmin = 1.5, xmax = 11, ymin = 8.75, ymax = 11.25, aes(fill = 'Heritability'), colour=NA, alpha = 0.05,  size = 2) +
  
  geom_dag_point(colour = '#ef654a', size = 30, shape = 'square') +
  geom_dag_edges(aes(x = xstart, y = ystart, xend = xend, yend = yend, edge_linetype =  linetype),edge_width = 1.5) +
 
  geom_dag_text(parse = TRUE, label =node_labels  , colour = 'white', size = 6) +
  
  scale_x_continuous(limits = c(0,12), expand= c(0,0), name = NA )+ 
  scale_y_continuous(limits = c(0,12), expand= c(0,0), name = NA )+ 
  
  scale_fill_manual(values = model_cols, 'Model')+ 
  theme_dag() + 
  theme(legend.position = 'bottom',
        text = element_text(size=20, family = "Questrial"))
  
