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
                                  linewidth = 0.3),
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
# 1. DISCRETE MODEL COLORS (Paleo-Safe)
colours <- c("niche"           = "#B8A9C9", # Light Shale
             "random"          = "#45354F", # Deep Pyrite (Darker for contrast)
             "ADBM"            = "#D98A47", # Ironstone
             "ATN"             = "#A66128", # Mudstone
             "Body-size ratio" = "#E6C18B", # Sandstone
             "PFIM"            = "#507D77"  # Glauconite Teal
)

# 2. CONTINUOUS OCEAN RAMP
# Low: #F2E8CF (Sand) -> Mid: #6A994E (Algae) -> High: #154734 (Anoxic)

# 3. DIVERGING (Non-White Midpoint)
# Low:  #364F6B (Cool Deep Blue)
# Mid:  #E3D5B8 (Parchment)
# High: #B03A2E (Extinction Red)



# Converts the named vector into a dataframe for easier mapping in ggplot
pal_df <- data.frame(l = names(colors), c = colors)
