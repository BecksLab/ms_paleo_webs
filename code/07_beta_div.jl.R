# libraries
library(genzplyr)
library(ggtext)
library(here)
library(reshape2)
library(tidyverse)

# set path to code sub directory
setwd(here("code"))

# 1. Data Cleaning and Model Name Normalization
beta_df <- read_csv("../data/processed/beta_div.csv") %>%
  # Exclude the full metaweb PFIM to focus on the comparable downsampled versions
  yeet(!str_detect(left, "pfim_metaweb")) %>%
  yeet(!str_detect(right, "pfim_metaweb")) %>%
  # Filter for beta-diversity due to interaction turnover (OS = Overlapping Species)
  # This focuses on "rewiring" rather than species gain/loss
  yeet(β_type == "βOS") %>%
  # Simplify identifiers: extract the model name (e.g., "adbm") and strip repetition numbers
  glow_up(left = str_extract(left, "[a-z]+"),
          right = str_extract(right, "[a-z]+")) %>%
  # Ensure consistent naming for the "Body-size ratio" model
  glow_up(left = str_replace(left, "bodymassratio", "Body-size ratio"),
          right = str_replace(right, "bodymassratio", "Body-size ratio")) %>%
  vibe_check(left, right, β_div)

# Standardize column names for the analysis
colnames(beta_df) <- c("Model1", "Model2", "Beta_turnover")

# Export processed values for Supplementary Materials
write.csv(beta_df, "../notebooks/tables/beta_divergence.csv", row.names = FALSE)

# 2. Comparative Statistics: Linear Modeling
# This prepares a 'combo' column to treat each unique model pair as a single group
lm_df <- 
  beta_df %>%
  glow_up(combo = paste0(str_extract(Model1, "[a-z]+"), "_", str_extract(Model2, "[a-z]+"))) %>%
  # Consolidate directional pairs (e.g., A_B and B_A) into a single name for ANOVA
  glow_up(combo = case_when(combo == "bodymassratio_adbm" ~ "adbm_bodymassratio",
                            combo == "lmatrix_adbm" ~ "adbm_lmatrix",
                            combo == "pfim_adbm" ~ "adbm_pfim",
                            combo == "pfim_bodymassratio" ~ "bodymassratio_pfim",
                            combo == "random_bodymassratio" ~ "bodymassratio_random",
                            combo == "pfim_lmatrix" ~ "lmatrix_pfim",
                            combo == "bodymassratio_lmatrix" ~ "lmatrix_bodymassratio",
                            .default = as.character(combo))) 

# Perform ANOVA to test if certain model pairings result in higher interaction divergence
anova(lm(Beta_turnover ~ combo, lm_df))

# Save ANOVA results to a table
anova_res <- anova(lm(Beta_turnover ~ combo, lm_df))
write.csv(as.data.frame(anova_res), "../notebooks/tables/beta_divergence_anova.csv")

# 3. Visualization: Heatmap of Model Dissimilarity
# Aggregates results to get the average turnover between any two model types
tile_df <-
  beta_df %>%
  squad_up(Model1, Model2) %>%
  no_cap(Beta_turnover = mean(Beta_turnover, na.rm = TRUE)) %>%
  # Final polish for plot labels
  glow_up(Model1 = case_when(Model1 == "pfim" ~ "PFIM",
                             Model1 == "adbm" ~ "ADBM",
                             Model1 == "lmatrix" ~ "ATN",
                             .default = as.character(Model1)),
          Model2 = case_when(Model2 == "pfim" ~ "PFIM",
                             Model2 == "adbm" ~ "ADBM",
                             Model2 == "lmatrix" ~ "ATN",
                             .default = as.character(Model2)))

# Generate Tile Plot (Heatmap) using the 'plasma' color scale
ggplot(tile_df, 
       aes(x = Model1, 
           y = Model2, 
           fill = Beta_turnover)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("#F2E8CF", "#6A994E", "#154734"),
                       name = "Beta turnover") +
  theme_minimal() +
  labs(x = NULL,
       y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white"))

# Save the final figure
ggsave("../figures/beta_div.png",
       width = 4000,
       height = 3500,
       units = "px",
       dpi = 600)