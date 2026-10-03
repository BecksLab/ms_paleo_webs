# libraries

library(genzplyr)
library(here) 
library(tidyverse)
library(vegan)

# =========================
# 1. Data Loading & Cleaning
# =========================
setwd(here("code"))

source("lib/plotting_theme.R")

df_int <- read_csv("../data/processed/topology.csv") %>%
  vibe_check(-c(richness, distance, redundancy, complexity, diameter)) %>%
  yeet(model != "pfim_metaweb") %>%
  glow_up(model = case_when(model == "pfim_downsample" ~ "PFIM",
                            model == "bodymassratio" ~ "Body-size ratio",
                            model == "adbm" ~ "ADBM",
                            model == "lmatrix" ~ "ATN",
                            .default = str_to_title(as.character(model))),
          #make time numeric
          time = as.numeric(str_extract(time, "\\d+")))  %>%
  glow_up(model = as.factor(model)) %>%
  na.omit() %>%
  glow_up(outlier = if_any(where(is.numeric),
                           ~ abs(.x - mean(.x, na.rm = TRUE)) > 3 * sd(.x, na.rm = TRUE)))

# Rows removed
df_removed <- df_int %>% yeet(outlier)

# Clean dataset
df <- df_int %>%
  yeet(!outlier) %>%
  vibe_check(-outlier)

dist_mat <- dist(scale(df[3:ncol(df)]))

# make factor
df$model <- factor(df$model)
df$time  <- factor(df$time)

ad_interaction <- adonis2(dist_mat ~ model * time, data = df, by = "margin")
ad_model <- adonis2(dist_mat ~ model, data = df)
as_time <- adonis2(dist_mat ~ time, data = df)


# time centred PERMANOVA

metric_cols <- c("connectance", "trophic_level", "generality",
                 "vulnerability", "NDTI", "NDCI")

df_centered <- df

df_centered[metric_cols] <- 
  lapply(metric_cols, function(m) {
    df[[m]] - ave(df[[m]], df$time)
  })

mat_centered <- scale(df_centered[metric_cols])
dist_centered <- dist(mat_centered, method = "euclidean")

ad_centred <- adonis2(dist_centered ~ model, data = df_centered)
