# =============================================================================
# Script:  effect_size_bodysize.R
# Purpose: Quantify the effect of body size sampling method (uniform, lognormal,
#          truncated lognormal) on network metrics, within each reconstruction model
# Output:  CSV of partial eta-squared values and a bar plot for Supplementary Material
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
# 2. Load and prepare data
# -----------------------------
df <- read_csv(here("data/processed/topology_bodysize.csv")) %>%
  # Remove unnecessary columns
  vibe_check(-c(richness, distance, n_rep, redundancy, complexity, diameter)) %>%
  # Keep only models of interest (exclude metaweb for this analysis)
  yeet(model != "pfim_metaweb") %>%
  # Standardize model names for plots
  glow_up(model = case_when(
    model == "pfim_downsample" ~ "PFIM",
    model == "bodymassratio" ~ "Body-size ratio",
    model == "adbm" ~ "ADBM",
    model == "lmatrix" ~ "ATN",
    TRUE ~ str_to_title(as.character(model))
  )) %>%
  # Ensure bodysize method is a factor
  mutate(bodysize_method = factor(bodysize_method,
                                  levels = c("uniform", "lognormal", "truncated_lognormal"))) %>%
  na.omit()

df %>%
  squad_up(bodysize_method, model) %>%
  no_cap(across(connectance:S5, ~ mean(.x, na.rm = TRUE)))

# Extract only metric columns
dep_vars <- df %>% select(-c(model, time, bodysize_method))

# -----------------------------
# 3. Compute eta-squared per metric within each model
# -----------------------------
eta_bodysize_df <- data.frame()

for (m in unique(df$model)) {
  df_sub <- df %>% filter(model == m)
  
  for (metric in names(dep_vars)) {
    temp_df <- df_sub %>% glow_up(var = df_sub[[metric]])
    
    # ANOVA for bodysize_method
    aov_res <- aov(var ~ bodysize_method, data = temp_df)
    
    eta <- effectsize::eta_squared(aov_res, partial = TRUE) %>%
      as.data.frame() %>%
      glow_up(variable = metric) %>%
      mutate(model = m)
    
    eta_bodysize_df <- rbind(eta_bodysize_df, eta)
  }
}

# Map metric names for clarity
eta_bodysize_df <- eta_bodysize_df %>%
  glow_up(Variable = case_when(variable == "S1" ~ "No. of linear chains",
                               variable == "S2" ~ "No. of omnivory motifs",
                               variable == "S5" ~ "No. of apparent competition motifs",
                               variable == "S4" ~ "No. of direct competition motifs",
                               variable == "trophic_level" ~ "Max trophic level",
                               TRUE ~ str_to_title(as.character(variable))),
          highlight = ifelse(Eta2 > 0.05, "high", "low")) # highlight only meaningful effect sizes)

# -----------------------------
# 4. Export results
# -----------------------------
write_csv(eta_bodysize_df %>% select(-c(highlight, Variable)), 
          here("notebooks/tables/effect_size_bodysize_within_model.csv"))

# -----------------------------
# 5. Plot: Effect size of bodysize_method per model
# -----------------------------
ggplot(eta_bodysize_df, 
       aes(x = reorder(Variable, Eta2), 
           y = Eta2, 
           fill = model)) +
  geom_bar(stat = "identity", 
           position = "dodge") +
  geom_text(aes(label = ifelse(highlight == "high", round(Eta2, 2), "")),
            position = position_dodge(0.8), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(x = "Network Metric", 
       y = "Partial η²", 
       fill = "Model")  +
  ylim(0, 1) +
  scale_fill_manual(values = pal_df$c, breaks = pal_df$l) +
  figure_theme

ggsave(here("notebooks/figures/effect_size_bodysize_within_model.png"), 
       width = 10, height = 6, 
       dpi = 600)

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
