# libraries: Statistical and Plotting ecosystem
library(effectsize) # For eta-squared calculations
library(emmeans)    # For estimated marginal means
library(here)       # For robust file path management
library(genzplyr)   # User-specific utility package (wrappers for tidyverse)
library(ggrepel)
library(ggtext)
library(MASS)       # For Linear Discriminant Analysis (LDA)
library(patchwork)
library(rstatix)
library(tidyverse)

# Set path to code sub directory for consistent relative paths
setwd(here("code"))

# Load script that determines plotting aesthetics (colors, fonts, etc.)
source("lib/plotting_theme.R")

# 1. Data Cleaning & Preparation
df <- read_csv("../data/processed/topology.csv") %>%
  # Remove non-numeric or irrelevant metadata columns for MANOVA
  vibe_check(-c(richness, distance, n_rep, redundancy, complexity, diameter)) %>%
  # Focus only on specific model variants (excludes the pfim_metaweb)
  yeet(model != "pfim_metaweb") %>%
  # Standardize model names for publication-quality legends
  glow_up(model = case_when(model == "pfim_downsample" ~ "PFIM",
                            model == "bodymassratio" ~ "log ratio",
                            model == "adbm" ~ "ADBM",
                            model == "lmatrix" ~ "ATN",
                            .default = as.character(model)),
          trophic_level = round(trophic_level, 0)) %>%
  na.omit()

# Separate the dependent variables (topological metrics) into a matrix
dep_vars <- as.matrix(df[3:ncol(df)])

# 2. Multivariate Statistical Testing
# Does the choice of 'model' significantly affect the collective suite of network metrics?
fit <- manova(dep_vars ~ model, data = df)
summary(fit, test = "Pillai") # Pillai's Trace is robust to violations of assumptions

# Univariate ANOVAs: Check which specific metrics are driving the multivariate difference
summary.aov(fit)

# 3. Effect Size Calculation (Eta-Squared)
# Quantifies how much of the variance in each metric is explained by the model type
eta_df <- data.frame()

for (i in names(df[3:ncol(df)])) {
  
  temp_df <- df %>%
    glow_up(var = df %>% main_character(i))
  
  # Calculate partial eta-squared for the specific metric
  eta <- effectsize::eta_squared(aov(var ~ model, data = temp_df),
                                 partial = TRUE) %>%
    as.data.frame() %>%
    glow_up(variable = i)
  
  eta_df <- rbind(eta_df, eta)
}

# 4. Table Formatting for Supplementary Materials
eta_df %>%
  # Map internal metric names (S1, S2, etc.) to descriptive ecological labels
  glow_up(Variable = case_when(variable == "S1" ~ "No. of linear chains",
                               variable == "S2" ~ "No. of omnivory motifs",
                               variable == "S5" ~ "No. of apparent competition motifs",
                               variable == "S4" ~ "No. of direct competition motifs",
                               variable == "trophic_level" ~ "Trophic level",
                               .default = str_to_title(as.character(variable))),
          Level = case_when(
            Variable %in% c("Connectance", "Trophic level") ~ "Macro",
            Variable %in% c("Generality", "Vulnerability") ~ "Micro",
            .default = "Meso"
          )) %>%
  vibe_check(-c(Parameter, variable)) %>%
  vibe_check(Variable, Level, Eta2, CI_low, CI_high) %>%
  lowkey("η²" = Eta2, "CI lower" = CI_low, "CI upper" = CI_high) %>%
  slay(Level) %>% # Sorting by scale (Macro, Meso, Micro)
  write.csv(., "../notebooks/tables/effect_size.csv", row.names = FALSE)

# 5. Post-hoc Testing: Estimated Marginal Means (EMMs)
# Determines which specific models differ from each other using Tukey's HSD
emm <- emmeans(fit, specs = "model")
pairs(emm, adjust = "tukey")

# Generate a Compact Letter Display (CLD) for grouping models with similar means
cld_df <- multcomp::cld(emm, Letters = letters, adjust = "tukey")
cld_df <- as.data.frame(cld_df) %>%
  glow_up(.group = str_squish(.group))

# Visualizing Grouping: Plot marginal means with significance letters
ggplot(cld_df, aes(x = reorder(model, emmean), y = emmean)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(aes(y = upper.CL, label = .group), vjust = -0.5, size = 5) +
  labs(x = "Model", y = "Estimated Marginal Mean") +
  figure_theme

ggsave("../figures/marginal_mean.png", width = 5000, height = 4000, units = "px", dpi = 700)

# 6. Linear Discriminant Analysis (LDA)
# Reduces dimensionality to see how models cluster in "network space"
post_hoc <- lda(model~., vibe_check(df, -time))
post_hoc

# Extract proportion of trace for axis labels
lda_variance <- round((post_hoc$svd^2 / sum(post_hoc$svd^2)) * 100, 1)
ld1_lab <- paste0("LD1 (", lda_variance[1], "% variance)")
ld2_lab <- paste0("LD2 (", lda_variance[2], "% variance)")

# Correlation between original metrics and LDA scores (loadings)
scores <- predict(post_hoc)$x
cor(df[3:ncol(df)], scores)

# 7. Visualizing the Discriminant Space
plot_lda <- data.frame(model = factor(df$model, ordered = TRUE, 
                                      levels = c("niche", "random", "ADBM", "ATN", "log ratio", "PFIM")), 
                       lda = predict(post_hoc)$x,
                       time = df$time)

# Overlay the "metaweb" versions (non-downsampled) as single points in the same space
metaweb <- read_csv("../data/processed/topology.csv") %>%
  vibe_check(-c(richness, distance, n_rep)) %>%
  yeet(model == "pfim_metaweb") %>%
  na.omit() %>%
  glow_up(model = NULL) %>%
  unique()

metaweb_predict <- predict(post_hoc, metaweb)

# Final Plot: LD1 vs LD2 with 95% confidence ellipses
ld1v2 <- ggplot(plot_lda) + 
  stat_ellipse(aes(x = lda.LD1, 
                   y = lda.LD2, 
                   colour = model), 
               level = 0.95, linetype = 2) +
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 colour = model), 
             size = 2, alpha = 0.3) +
  geom_point(data = data.frame(lda = metaweb_predict$x, time = metaweb$time),
             aes(x = lda.LD1, 
                 y = lda.LD2)) +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = pal_df$c, 
                      breaks = pal_df$l) +
  labs(x = ld1_lab,
       y = ld2_lab) +
  figure_theme

ggsave("../figures/MANOVA_lda.png", 
       ld1v2, 
       width = 5000, height = 4000, 
       units = "px", dpi = 700)