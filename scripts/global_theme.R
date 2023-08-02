# Set theme
require(ggplot2)
require(showtext)
require(showtextdb)

# Format fonts for output
font_add("lmsans10", 'lmsans10-regular.otf')
showtext_auto()

my_theme <- theme_classic(base_family = "lmsans10", base_size = 10)+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank()
  )
