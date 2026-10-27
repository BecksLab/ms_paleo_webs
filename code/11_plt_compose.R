
library(tidyverse)
library(patchwork)

# load plot objects
setwd(here("code"))

source("lib/plotting_theme.R")

# load plot objects
linear_all <- readRDS("../data/outputs/linear_time.rds")
anova_summary <- readRDS("../data/outputs/anova_summary.rds")
loadings_plot <- readRDS("../data/outputs/loadings_plot.rds")
manova_plot <- readRDS("../data/outputs/manova_plot.rds")


layout <- "
AB
CB
"

combined_plot <- manova_plot  + linear_all + anova_summary +
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = 'A')

ggsave("../figures/combined_fg1.png",
       combined_plot,
       width = 5500,
       height = 3500,
       units = "px",
       dpi = 300)
