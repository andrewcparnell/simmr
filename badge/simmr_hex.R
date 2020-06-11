# Design a Bchron hex sticker

rm(list = ls())
library(simmr)
library(hexSticker)
library(ggplot2)
library(showtext)
library(ggimage)
library(simmr)

data("geese_data_day1")

simmr_1 = with(geese_data_day1, 
               simmr_load(mixtures=mixtures,
                          source_names=source_names[-2],
                          source_means=source_means[-2,],
                          source_sds=source_sds[-2,]))

p = plot(simmr_1, title = '') + 
  theme_void() + 
  theme(legend.position= 'None')

# ages1 = BchronCalibrate(ages=11553,
#                         ageSds=230,
#                         calCurves='intcal13',
#                         ids='Ox-123456')
# p = plot(ages1, dateLabels = FALSE, fillCol = "#F0AB00") + 
#   geom_line(data = as.data.frame(ages1$`Ox-123456`), aes(x = ageGrid, y = densities), col = "#822327", size = 2) +
#   geom_line(data = as.data.frame(ages1$`Ox-123456`), aes(x = ageGrid, y = 0), col = "#822327", size = 2) +
#   ggtitle('') + theme_void() + theme_transparent()


## Loading Google fonts (http://www.google.com/fonts)
fname = "Alegreya Sans" #sample(font_families_google(), 1)
font_add_google(fname)
## Automatically use showtext to render text for future devices
showtext_auto()

sticker(p, package="simmr", p_size=10, s_x=1, s_y=0.8, 
        s_width=1.3, s_height=1,p_family = fname,h_size = 2,
        p_color = viridis(3)[2],
        h_fill="#FFFFFF", h_color="#440154FF",
        filename="badge/simmr_badge.png")
