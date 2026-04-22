# =============================================================================
# Script:  param_effects.R
# Purpose: Quantify the effect of different parameter specifications for the ADBM
#          BRM models
# Output:  TODO
# =============================================================================

# -----------------------------
# 1. Libraries
# -----------------------------
library(effectsize)   # for eta-squared calculations
library(here)         # robust file paths
library(tidyverse)    # data wrangling & plotting
library(genzplyr)     # utility wrappers (vibe_check, glow_up, yeet)
library(ggplot2)      # plotting
library(ggtext)       # enhanced text in plots
library(rstatix)
library(effectsize)

# set path to code sub dir
setwd(here("code"))

#load script that determines plotting aesthetics
source("lib/plotting_theme.R")

# -----------------------------
# 6. BMR parms
# -----------------------------

df <- read_csv("../data/processed/topology_bmr_params.csv") %>%
  vibe_check(param_set, time,
             richness, connectance, diameter,
             trophic_level, generality, vulnerability,
             S1, S2, S4, S5) %>%
  na.omit()

Y <- df %>%
  vibe_check(param_set, time,
             richness, connectance, diameter,
             trophic_level, generality, vulnerability,
             S1, S2, S4, S5)

df %>%
  squad_up(param_set, time) %>%
  no_cap(across(richness:S5, ~ mean(.x, na.rm = TRUE)))

# -----------------------------
# 7. PERMANOVA
# -----------------------------
library(vegan)

Y_mat <- scale(Y %>% vibe_check(connectance, trophic_level, generality, vulnerability))

adonis2(Y_mat ~ param_set, data = Y, method = "euclidean")

rda_mod <- rda(Y_mat ~ param_set, data = Y)

anova(rda_mod, permutations = 999)
anova(rda_mod, by = "axis", permutations = 999)

RsquareAdj(rda_mod)

dispersion <- betadisper(dist(Y_mat), df$param_set)
anova(dispersion)

# RDA
rda_mod <- rda(Y_mat ~ param_set, data = df)
scores <- as.data.frame(scores(rda_mod, display = "sites"))
scores$param_set <- df$param_set

loadings <- scores(rda_mod, display = "species")

loadings_RDA1 <- sort(abs(loadings[,1]), decreasing = TRUE)

centroids <- scores %>%
  group_by(param_set) %>%
  summarise(RDA1 = mean(RDA1),
            RDA2 = mean(RDA2))

# convex hulls
hulls <- scores %>%
  group_by(param_set) %>%
  slice(chull(RDA1, RDA2))

ggplot(scores, 
       aes(RDA1, 
           RDA2, 
           color = param_set)) +
  geom_point(alpha = 0.35, 
             size = 1.2) +
  geom_polygon(data = hulls,
               aes(fill = param_set),
               alpha = 0.15,
               color = NA) +
  geom_point(data = centroids,
             aes(RDA1, RDA2, color = param_set),
             size = 4) +
  scale_colour_manual(values = c("#BA975B", "#5B6BBA", "#5BA7BA", "#BA785B", "#655F53", "#535665", "#536265"))+
  scale_fill_manual(values = c("#BA975B", "#5B6BBA", "#5BA7BA", "#BA785B", "#655F53", "#535665", "#536265")) +
  labs(x = "RDA1 (structural gradient)",
       y = "RDA2 (residual variation)") +
  figure_theme

ggsave(here("notebooks/figures/parm_vals_BMR.png"), 
       width = 10, height = 6, 
       dpi = 600)
