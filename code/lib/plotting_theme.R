library(ggplot2)
library(extrafont)
library(systemfonts)
library(showtext)
library(ggtext)
library(patchwork)
library(ggfx)
library(ggforce)

## 1. Typography Setup ----
# Loads "Noto Sans" from Google Fonts for a clean, modern look
font_add_google("Noto Sans", "Noto")

# These commands help R find and manage the system fonts
font_paths()  
font_files()
font_families()

# Automatically triggers showtext to ensure custom fonts 
# render correctly in high-res PNG outputs
trace(grDevices::png, exit = quote({
  showtext::showtext_begin()
}), print = FALSE)

## 2. General Figure Theme ----
# A "naked" classic theme with specific professional overrides
figure_theme = 
  theme_classic() +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"), # Transparent fill
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = colorspace::darken("#dddddd", 0.1),
                                  size = 0.3),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_blank(),
        text = element_text(color = "#5e5e5e"), # Soft gray text for better readability
        plot.margin = margin(10, 5, 5, 10),
        legend.margin = margin(1, 2, 1, 2)
  )

## 3. Model Color Palette ----
# Assigns specific hex codes to each modelling framework
# This ensures a model is always the same colour across all plots
colors <- c("niche" = "#D760F6",
            "random" = "#B01DFF",
            "ADBM" = "#F5BD63",
            "ATN" = "#9E8324",
            "log ratio" = "#FAE2A6",
            "PFIM" = "#26ECC9")         # Pink

# Converts the named vector into a dataframe for easier mapping in ggplot
pal_df <- data.frame(l = names(colors), c = colors)