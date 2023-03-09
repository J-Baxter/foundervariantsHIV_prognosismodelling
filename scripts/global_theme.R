# Set theme
require(ggplot2)
require(showtext)
require(showtextdb)

# Format fonts for output
font_add_google("Questrial", "Questrial")
showtext_auto()

my_theme <- theme_classic(base_family = "Questrial")+
  theme(
    text = element_text(size=20),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )