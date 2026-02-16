# libraries

library(genzplyr)
library(here) 
library(insight)
library(lme4)
library(performance)
library(tidyverse)
library(vegan)

# =========================
# 1. Data Loading & Cleaning
# =========================
setwd(here("code"))

source("lib/plotting_theme.R")

df <- read_csv("../data/processed/topology.csv") %>%
  vibe_check(-c(richness, distance, n_rep, redundancy, complexity, diameter)) %>%
  yeet(model != "pfim_metaweb") %>%
  glow_up(model = case_when(
    model == "pfim_downsample" ~ "PFIM",
    model == "bodymassratio" ~ "Body-size ratio",
    model == "adbm" ~ "ADBM",
    model == "lmatrix" ~ "ATN",
    TRUE ~ as.character(model)
  ), trophic_level = round(trophic_level, 0)) %>%
  na.omit()

dist_mat <- dist(scale(df[3:ncol(df)]))

adonis2(dist_mat ~ model * time, data = df, by = "margin")

# make factor
df$model <- factor(df$model)
df$time  <- factor(df$time)

ad <- adonis2(dist_mat ~ model * time, data = df)
