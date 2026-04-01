# -----------------------------
# 1. Libraries
# -----------------------------
library(here)         # robust file paths
library(tidyverse)    # data wrangling & plotting
library(genzplyr)     # utility wrappers (vibe_check, glow_up, yeet)
library(ggplot2)      # plotting
library(ggtext)       # enhanced text in plots

#load script that determines plotting aesthetics
source("lib/plotting_theme.R")

setwd(here("code"))

# -----------------------------
# 2. Load and prepare data
# -----------------------------
df <- read_csv(here("data/processed/pfim_downsampling_similarity.csv"))

df %>%
  squad_up(downsample_y) %>%
  no_cap(mean = mean(similarity))


ggplot(df,
       aes(x = downsample_y,
           y = similarity)) +
  geom_boxplot(aes(group = downsample_y),
               fill = colours["PFIM"]) +
  labs(x = "Downsampling parameter",
       y = "Jaccard similarity") +
  ylim(0,1) +
  figure_theme

ggsave(here("notebooks/figures/pfim_downsample.png"), 
       width = 10, height = 6, 
       dpi = 600)
