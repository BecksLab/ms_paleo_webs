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

accent_col <- "#5e5e5e"
yotes_black <- "#292D32"

## 2. General Figure Theme ----
# A "naked" classic theme with specific professional overrides
figure_theme =
  theme_classic(base_size = 14) +
  theme(text = element_text(colour = accent_col,
                            family = "Noto"),
        plot.title = element_text(size = rel(1.5), face = "bold"),
        plot.subtitle = element_text(size = rel(1.2)),
        plot.caption = element_text(size = rel(0.9)),
        axis.title = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1.0)),
        strip.text = element_text(size = rel(1.1), face = "bold"),
        legend.title = element_text(size = rel(1.0), face = "bold"),
        legend.text = element_text(size = rel(0.95)),
        panel.border = element_rect(colour = yotes_black, fill = "#ffffff00"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = colorspace::darken("#dddddd", 0.1),
                                  linewidth = 0.3),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_blank(),
        plot.margin = margin(10, 5, 5, 10),
        legend.margin = margin(1, 2, 1, 2)
  )

## 3. Model Color Palette ----
# Assigns specific hex codes to each modelling framework
# This ensures a model is always the same colour across all plots
# 1. DISCRETE MODEL COLORS (Paleo-Safe)
colours <- c(
  "PFWIM"           = "#1E7548", # Deep Teal
  "Random"          = "#5F249F", # Deep Violet
  "Niche"           = "#69B3E7", # Soft Violet
  "ATN"             = "#6F263D", # Rich Rust
  "ADBM"            = "#A9431E", # Amber
  "Body-size ratio" = "#B9975B"  # Soft Goldenrod
)

pal_df <- data.frame(l = names(colours), c = colours)
# make named vector
names(pal_df$c) <- pal_df$l

# 2. CONTINUOUS OCEAN RAMP
col_cont <- c("#DED0BE", "#6F263D")

# 3. DIVERGING
col_div <- c("#6F263D", "#F4F1DE", "#5F249F")

# 4. LOOKUP TABLE (for consistent renaming)
metric_lookup <- tibble::tribble(
  ~statistic, ~stat_label, ~level,
  "connectance", "Connectance", "Macro",
  "trophic_level", "Max trophic level", "Macro",
  "generality", "Generality", "Micro",
  "vulnerability", "Vulnerability", "Micro",
  "NDTI", "Trophic index", "Meso",
  "NDCI", "Competition index", "Meso"
)
