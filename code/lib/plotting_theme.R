library(ggplot2)
library(extrafont)
library(systemfonts)
library(showtext)
library(ggrepel)
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
        text = element_text(color = "#5e5e5e", family = "Noto"), # Soft gray text for better readability
        plot.margin = margin(10, 5, 5, 10),
        legend.margin = margin(1, 2, 1, 2)
  )

## 3. Model Color Palette ----
# Assigns specific hex codes to each modelling framework
# This ensures a model is always the same colour across all plots
# 1. DISCRETE MODEL COLORS (Paleo-Safe)
colours <- c(
  # Group 1: The Anchor
  "PFIM"            = "#1E7548", # Deep Teal
  "random"          = "#5F249F", # Deep Violet
  "niche"           = "#69B3E7", # Soft Violet
  "ATN"             = "#6F263D", # Rich Rust
  "ADBM"            = "#A9431E", # Amber
  "Body-size ratio" = "#B9975B"  # Soft Goldenrod
)

pal_df <- data.frame(l = names(colours), c = colours)

# 2. CONTINUOUS OCEAN RAMP (Tied to Group 1)
# High-vibrancy transition from ice-blue to PFIM Teal
col_cont <- c("#AEfED5", "#32BF76", "#1E7548")

# 3. DIVERGING (Tied to Group 2 & 3)
# Contrasts the 'ATN' warmth against the 'random' cool violet
col_div <- c("#6F263D", "#F4F1DE", "#5F249F")