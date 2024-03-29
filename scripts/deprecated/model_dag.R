# DAG
require(ggdag)
require(dagitty)

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
  x = c(SpVLr = 10, SpVLt = 2.5, deltaCD4 = 10, PMV = 2.5),
  y = c(SpVLr = 8, SpVLt = 8, deltaCD4 = 2, PMV = 2)) %>% 
  coords2df() %>%
  coords2list()

test_dag <- dagify(
  deltaCD4 ~ SpVLr,
  SpVLr ~ SpVLt,
  PMV ~ SpVLt ) +
  edge_label("deltaCD4", "SpVLr", label = "Edge Label 1") +
  edge_label("SpVLr", "SpVLt", label = "Edge Label 2") +
  edge_label("PMV", "SpVLt", label = "Edge Label 3")


dagitty::coordinates(test_dag) <- coords
label(test_dag) <- c('deltaCD4' = expression(paste(Delta, ' CD4+ ')),
                              'SpVLr' = expression(paste(SpVL[R])),
                              'SpVLt' = expression(paste(SpVL['T'])),
                              'PMV' = 'P(MV)')
plt_2a <- test_dag %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(size = 14, col = "white", fill = "white") +
  geom_dag_edges() +
  geom_dag_text(col = "black", parse = T, size = 3, 
                label = c('P(MV)',
                          
                          expression(paste(SpVL[R])),
                          expression(paste(SpVL[T])),
                          expression(paste(Delta, ' CD4+ ')))) +
  annotate("text", x = 6.2, y = 8.3, label = "Heritability", family ="LM Sans 10", size = 3) + 
  annotate("text", x = 2, y = 5, label = "Transmission", angle = 90, family ="LM Sans 10", size = 3) + 
  annotate("text", x = 10.4, y = 5, label = "Tolerance", angle = 270, family ="LM Sans 10", size = 3) + 
  scale_x_continuous(limits = c(1.5,11), expand = c(0,0)) + 
  scale_y_continuous(limits = c(1,9), expand = c(0,0)) + 
  theme_dag(base_family = "LM Sans 10") +
  theme(
    text = element_text(size=10),
    #axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    #axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 8),
    #axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    strip.text  = element_text(size = 8),
    #legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank()
  )

ggdag(test_dag)

model_cols <- c('Heritability' = '#fdbb84', 'Transmission' ='#d7301f', 'Tolerance' = '#7f0000')
node_labels <- list( 'P[MV]', 'SpVL[R]', 'SpVL[T]' , 'DeltaCD4')