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
    TRUE ~ as.character(model)
  )) %>%
  # Ensure bodysize method is a factor
  mutate(bodysize_method = factor(bodysize_method,
                                  levels = c("uniform", "lognormal", "truncated_lognormal"))) %>%
  na.omit()

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
