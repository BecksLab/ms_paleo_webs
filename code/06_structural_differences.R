# =========================
# NETWORK ANALYSIS WORKFLOW
# =========================

# 0. Libraries
library(effectsize)
library(emmeans)
library(here)         # file paths
library(genzplyr)
library(ggpubr)
library(ggrepel)
library(ggtext)
library(lme4)
library(lmerTest)
library(MASS)         # LDA
library(multcomp)
library(MVN)
library(patchwork)
library(purrr)
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

# =========================
# 6. Pairwise - for supp matt
# =========================

run_model_comparisons <- function(data,
                                  metrics = c("connectance",
                                              "trophic_level",
                                              "generality",
                                              "vulnerability",
                                              "S1",
                                              "S2",
                                              "S4",
                                              "S5"),
                                  random_effect = "time",
                                  p_adjust = "sidak") {
  
  # Ensure factors
  data$model <- as.factor(data$model)
  data[[random_effect]] <- as.factor(data[[random_effect]])
  
  all_emm <- list()
  all_pairs <- list()
  all_models <- list()
  
  for (metric in metrics) {
    
    # Build model formula
    formula <- as.formula(
      paste(metric, "~ model + (1|", random_effect, ")")
    )
    
    # Fit mixed model
    model_fit <- lmer(formula, data = data)
    
    # Estimated marginal means (explicit spec avoids family confusion)
    emm <- emmeans(model_fit, specs = "model")
    
    # Explicit pairwise contrasts
    pairwise_res <- contrast(emm,
                             method = "pairwise",
                             adjust = p_adjust)
    
    # Compact letter display
    cld_res <- multcomp::cld(emm,
                             adjust = p_adjust,
                             Letters = letters)
    
    # Clean dataframe
    cld_df <- as.data.frame(cld_res) %>%
      mutate(metric = metric) %>%
      mutate(.group = gsub(" ", "", .group))
    
    # Store results
    all_emm[[metric]] <- cld_df
    all_pairs[[metric]] <- as.data.frame(pairwise_res)
    all_models[[metric]] <- model_fit
  }
  
  combined_emm <- bind_rows(all_emm)
  
  return(list(
    emmeans = combined_emm,
    pairwise = all_pairs,
    models = all_models,
    adjust_method = p_adjust
  ))
}

emm_results <- run_model_comparisons(df)

# clean for plotting
emm_df <- emm_results$emmeans %>%
  glow_up(stat = case_when(metric == "connectance" ~ "Connectance",
                           metric == "trophic_level" ~ "Max trophic level",
                           metric == "generality" ~ "Generality",
                           metric == "vulnerability" ~ "Vulnerability",
                           metric == "S1" ~ "No. of linear chains",
                           metric == "S2" ~ "No. of omnivory motifs",
                           metric == "S5" ~ "No. of apparent competition motifs",
                           metric == "S4" ~ "No. of direct competition motifs"),
          level = case_when(stat %in% c("Connectance", "Max trophic level") ~ "Macro",
                            stat %in% c("Generality", "Vulnerability") ~ "Micro",
                            TRUE ~ "Meso"))

levs <- c("Macro", "Meso", "Micro")

plot_list_emm <- vector("list", length = 3)

for (i in seq_along(levs)) {
  
  plot_list_emm[[i]] <-
    ggplot(
      emm_df %>% filter(level == levs[i]),
      aes(x = model,
          y = emmean,
          colour = model)
    ) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = lower.CL,
                      ymax = upper.CL),
                  width = 0.2,
                  linewidth = 0.6) +
    geom_text(aes(label = .group,
                  y = upper.CL),
              vjust = -0.6,
              size = 3.5,
              colour = "black") +
    facet_wrap(
      vars(stat),
      scales = "free_y",
      ncol = 2
    ) +
    scale_colour_manual(values = pal_df$c,
                        breaks = pal_df$l) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      x = "Model",
      y = "Estimated marginal mean",
      title = levs[i]
    ) +
    coord_cartesian(clip = "off") +
    figure_theme
}

model_comparison_plot <-
  plot_list_emm[[1]] /
  plot_list_emm[[2]] /
  plot_list_emm[[3]] +
  plot_layout(
    guides = "collect",
    heights = c(1, 2, 1)
  )

model_comparison_plot
