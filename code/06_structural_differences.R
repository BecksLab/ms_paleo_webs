# =========================
# NETWORK ANALYSIS WORKFLOW
# =========================

# 0. Libraries
library(effectsize)   # eta-squared
library(here)         # file paths
library(genzplyr)
library(ggpubr)
library(ggrepel)
library(ggtext)
library(MASS)         # LDA
library(MVN)
library(patchwork)
library(rstatix)      # stats utilities
library(tidyverse)
library(heplots)      # MANOVA assumption checks
library(biotools)     # Box's M test
library(candisc)      # canonical discriminant analysis

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

# Dependent variable matrix for MANOVA/CDA
dep_vars <- as.matrix(df[3:ncol(df)])

# =========================
# 2. MANOVA + Assumption Checks
# =========================
fit <- manova(dep_vars ~ model, data = df)

# --- Multivariate test
summary(fit, test = "Pillai")

# --- Univariate ANOVAs
summary.aov(fit)

# --- MANOVA Assumption Checks

## 2a. Multivariate Normality (Henze-Zirkler)
result_hz <- hz(data = dep_vars)

## 2b. Homogeneity of Covariance (Box's M test)
boxM(dep_vars, df$model) # p < 0.05 indicates heterogeneity

## 2c. Multicollinearity
cor_matrix <- cor(dep_vars)
high_cor <- which(abs(cor_matrix) > 0.9, arr.ind = TRUE) # flag very high correlations


# =========================
# 3. Canonical Discriminant Analysis (CDA)
# =========================
cda <- candisc(fit) # integrates MANOVA with canonical variates
summary(cda)

# --- Extract canonical loadings (which metrics drive variance)
round(cda$structure[, 1:3], 2)

# --- Canonical means (model positions in CV space)
round(cda$means[, 1:3], 2)

# =========================
# 4. LDA for visualization
# =========================
lda_fit <- lda(model ~ ., data = df %>% vibe_check(-time))
lda_scores <- predict(lda_fit)$x

lda_variance <- round((lda_fit$svd^2 / sum(lda_fit$svd^2)) * 100, 1)
ld1_lab <- paste0("LD1 (", lda_variance[1], "% variance)")
ld2_lab <- paste0("LD2 (", lda_variance[2], "% variance)")

# Correlation between LDA axes and original metrics
cor(df[3:ncol(df)], lda_scores)

# =========================
# 5. LDA Visualization
# =========================
plot_lda <- data.frame(
  model = factor(df$model, levels = c("niche","random","ADBM","ATN","Body-size ratio","PFIM")),
  lda = lda_scores,
  time = df$time
)

ggplot(plot_lda, aes(x = lda.LD1, 
                     y = lda.LD2, 
                     colour = model, 
                     fill = model,
                     shape = model)) + 
  stat_ellipse(aes(x = lda.LD1, 
                   y = lda.LD2, 
                   colour = model), 
               level = 0.95, linetype = 2)  +
  geom_point(size = 2, alpha = 0.4) +
  labs(x = ld1_lab, y = ld2_lab) +
  scale_colour_manual(values = pal_df$c, breaks = pal_df$l) +
  scale_fill_manual(values = pal_df$c, breaks = pal_df$l) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha = 1, size = 3))) +
  figure_theme
ggsave("../figures/MANOVA_lda.png", width = 5000, height = 4000, units = "px", dpi = 700)

# =========================
# 6. Canonical Loadings Table
# =========================

# Extract canonical loadings
loadings <- round(cda$structure[, 1:3], 2)
loadings_df <- as.data.frame(loadings)
loadings_df$Metric <- rownames(loadings_df)

# Reorder columns: Metric first
loadings_df <- loadings_df[, c("Metric", "Can1", "Can2", "Can3")]

# Save table for supplementary materials
write.csv(loadings_df, "../notebooks/tables/canonical_loadings.csv",
          row.names = FALSE)

loadings_df <- as.data.frame(cda$structure[, 1:2]) %>%
  rownames_to_column("Metric") %>%
  glow_up(CV1 = Can1,
          CV2 = Can2)  %>%
  mutate(
    abs_loading = pmax(abs(CV1), abs(CV2)),
    keep = abs_loading >= 0.4
  ) %>%
  mutate(Level = case_when(
    Metric %in% c("richness", "connectance", "trophic_level") ~ "Macro",
    Metric %in% c("S1", "S2", "S4", "S5") ~ "Meso",
    Metric %in% c("generality", "vulnerability") ~ "Micro",
    TRUE ~ "Other"
  ))

# =========================
# 7. Canonical Loadings Plot
# =========================

ggplot(loadings_df, aes(x = CV1, y = CV2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_segment(
    aes(x = 0, y = 0, xend = CV1, yend = CV2),
    arrow = arrow(length = unit(0.02, "npc")),
    linewidth = 0.8,
    alpha = 0.8
  ) +
  geom_segment(
    aes(x = 0, y = 0, xend = CV1, yend = CV2, colour = Level),
    arrow = arrow(length = unit(0.02, "npc")),
    linewidth = 0.9
  ) +
  scale_colour_manual(values = c(
    Macro = "#1b9e77",
    Meso  = "#d95f02",
    Micro = "#7570b3"
  )) +
  geom_text_repel(
    data = subset(loadings_df, keep),
    aes(label = Metric),
    size = 4,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    x = "CV1",
    y = "CV2"
  ) +
  figure_theme

# Save figure
ggsave("../figures/canonical_loadings_plot.png", width = 8, height = 5, dpi = 300)
