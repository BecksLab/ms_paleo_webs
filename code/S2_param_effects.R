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
library(vegan)

# set path to code sub dir
setwd(here("code"))

#load script that determines plotting aesthetics
source("lib/plotting_theme.R")

models <- c("BMR", "ADBM")

results <- list()

for (m in models) {

  df_m <- read_csv("../data/processed/topology_param_effects.csv") %>%
    yeet(model == m) %>%
    vibe_check(param_set, time,
               richness, connectance, diameter,
               trophic_level, generality, vulnerability,
               S1, S2, S4, S5) %>%
    na.omit()

  Y <- df_m %>%
    vibe_check(param_set, time,
               connectance, trophic_level,
               generality, vulnerability)

  Y_mat <- scale(Y %>%
    vibe_check(connectance, trophic_level,
               generality, vulnerability))

  # PERMANOVA
  ad <- adonis2(Y_mat ~ param_set, data = Y, method = "euclidean")

  # RDA
  rda_mod <- rda(Y_mat ~ param_set, data = df_m)

  results[[m]] <- list(
    permanova = ad,
    rda = rda_mod,
    r2 = RsquareAdj(rda_mod)
  )
}

df_summary <- tibble(
  model = names(results),
  R2_permanova = map_dbl(results, ~ .x$permanova$R2[1]),
  R2_rda = map_dbl(results, ~ .x$r2$adj.r.squared)
)

ggplot(df_summary, aes(x = model, y = R2_permanova, fill = model)) +
  geom_col() +
  labs(y = "Variance explained by parameters",
       x = "Model") +
  figure_theme

# -----------------------------
# 7. BMR
# -----------------------------

df <- read_csv("../data/processed/topology_param_effects.csv") %>%
  yeet(model == "BMR") %>%
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

Y_mat <- scale(Y %>% 
                 vibe_check(connectance, trophic_level, generality, vulnerability))

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
