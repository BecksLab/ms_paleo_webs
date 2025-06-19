library(ggplot2)
library(extrafont)
library(systemfonts)
library(showtext)
library(ggtext)
library(patchwork)
library(ggfx)
library(ggforce)

##Fonts----

font_add_google("Noto Sans",
                "Noto")

font_paths()  
font_files()
font_families()

trace(grDevices::png, exit = quote({
  showtext::showtext_begin()
}), print = FALSE)

##General Theme for plotting----

figure_theme = 
  theme_classic() +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = colorspace::darken("#dddddd", 0.1),
                                  size = 0.3),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_blank(),
        text = element_text(color = "#5e5e5e"), #font change
        plot.margin = margin(10, 5, 5, 10),
        legend.margin = margin(1, 2, 1, 2)
  )


##Colour palete----
colors <- c("niche" = "#D760F6",
            "random" = "#B01DFF",
            "adbm" = "#F5BD63",
            "lmatrix" = "#9E8324",
            "bodymassratio" = "#FAE2A6",
            "pfim" = "#26ECC9")

# Join colors with categories
pal_df <- data.frame(c = unname(colors)[1:6], l = names(colors)[1:6])
