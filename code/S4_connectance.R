# =========================
# NETWORK ANALYSIS WORKFLOW
# =========================

# 0. Libraries
library(effectsize)
library(emmeans)
library(here)         # file paths
library(genzplyr)
library(ggrepel)
library(ggtext)
library(MASS)         # LDA
library(multcomp)
library(MVN)
library(patchwork)
library(purrr)
library(rstatix)      # stats utilities
library(tidyverse)
library(heplots)      # MANOVA assumption checks
library(candisc)      # canonical discriminant analysis

# =========================
# 1. Data Loading & Cleaning
# =========================
setwd(here("code"))

source("lib/plotting_theme.R")

df_int <- read_csv("../data/processed/topology_connectance.csv") %>%
  vibe_check(-c(richness, distance, n_rep, redundancy, complexity, diameter)) %>%
  yeet(model != "pfim_metaweb") %>%
  glow_up(model = case_when(model == "pfim_downsample" ~ "PFIM",
                            model == "bodymassratio" ~ "Body-size ratio",
                            model == "adbm" ~ "ADBM",
                            model == "lmatrix" ~ "ATN",
                            TRUE ~ str_to_title(as.character(model)))) %>%
  na.omit() %>%
  glow_up(outlier = if_any(where(is.numeric),
                           ~ abs(.x - mean(.x, na.rm = TRUE)) > 3 * sd(.x, na.rm = TRUE)))

# Rows removed
df_removed <- df_int %>% yeet(outlier)

# Clean dataset
df <- df_int %>%
  yeet(!outlier) %>%
  vibe_check(-outlier)

# Dependent variable matrix for MANOVA/CDA
network_stats <- c("connectance", "trophic_level", "generality",
                   "vulnerability", "NDTI", "NDCI")
dep_vars <- df %>%
  vibe_check(all_of(network_stats)) %>%
  as.matrix()

# remove outlier networks (specifically trophic level)

# =========================
# 2. MANOVA + Assumption Checks
# =========================
fit <- manova(dep_vars ~ model, data = df)

# --- Multivariate test
summary(fit, test = "Pillai")

# --- Univariate ANOVAs
summary.aov(fit)

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
lda_fit <- lda(model ~ ., data = df %>% vibe_check(-time, -targ_Co))
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
  model = factor(df$model, levels = c("Niche","Random","ADBM","ATN",
                                      "Body-size ratio", "PFIM")),
  lda = lda_scores,
  time = df$time
)

plot_lda$model <- factor(plot_lda$model, levels = pal_df$l)

# Calculate centroid and SE for each model
lda_summary <- plot_lda %>%
  group_by(model) %>%
  summarise(
    mean_LD1 = mean(lda.LD1, na.rm = TRUE),
    mean_LD2 = mean(lda.LD2, na.rm = TRUE),
    
    SE_LD1 = sd(lda.LD1, na.rm = TRUE) /
      sqrt(sum(!is.na(lda.LD1))),
    
    SE_LD2 = sd(lda.LD2, na.rm = TRUE) /
      sqrt(sum(!is.na(lda.LD2))),
    
    .groups = "drop"
  )

lda_plot <-
  ggplot(plot_lda, 
       aes(x = lda.LD1, 
           y = lda.LD2, 
           colour = model,
           fill = model))  +
  geom_hline(yintercept = 0, colour = accent_col) +
  geom_vline(xintercept = 0, colour = accent_col) +
  stat_ellipse(aes(group = model),
               type = "norm",
               level = 0.68,
               geom = "polygon",
               alpha = 0.3) +
  geom_point(data = lda_summary,
             aes(x = mean_LD1, 
                 y = mean_LD2, 
                 colour = model),
             size = 4) +
  labs(x = ld1_lab, y = ld2_lab) +
  scale_colour_manual(values = pal_df$c,
                      name = "Model") +
  scale_fill_manual(values = pal_df$c,
                    name = "Model") +
  figure_theme

# =========================
# 6. Canonical Loadings Table
# =========================

# Extract canonical loadings
loadings <- round(cda$structure[, 1:3], 2)
loadings_df <- as.data.frame(loadings)
loadings_df$Metric <- rownames(loadings_df)

# Reorder columns: Metric first
loadings_df <- loadings_df[, c("Metric", "Can1", "Can2", "Can3")]

loadings_df <- as.data.frame(cda$structure[, 1:2]) %>%
  rownames_to_column("Metric") %>%
  glow_up(CV1 = Can1,
          CV2 = Can2)  %>%
  mutate(
    abs_loading = pmax(abs(CV1), abs(CV2)),
    keep = abs_loading >= 0.4
  ) %>%
  mutate(stat = case_when(Metric == "connectance" ~ "Connectance",
                          Metric == "trophic_level" ~ "Max trophic level",
                          Metric == "generality" ~ "Generality",
                          Metric == "vulnerability" ~ "Vulnerability",
                          Metric == "NDTI" ~ "Trophic index",
                          Metric == "NDCI" ~ "Competition index"),
         Level = case_when(stat %in% c("Connectance", "Max trophic level") ~ "Macro",
                           stat %in% c("Generality", "Vulnerability") ~ "Micro",
                           TRUE ~ "Meso"))

# =========================
# 7. Canonical Loadings Plot
# =========================

loadings_plot <-
  ggplot(subset(loadings_df, keep), 
         aes(x = CV1, 
             y = CV2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = accent_col) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = accent_col) +
  geom_segment(aes(x = 0, 
                   y = 0, 
                   xend = CV1, 
                   yend = CV2),
               arrow = arrow(length = unit(0.02, "npc")),
               linewidth = 0.8,
               alpha = 0.8) +
  geom_segment(aes(x = 0, 
                   y = 0, 
                   xend = CV1, 
                   yend = CV2, 
                   colour = Level),
               arrow = arrow(length = unit(0.02, "npc")),
               linewidth = 0.9) +
  scale_colour_manual(values = c(Macro = "#D9C589",
                                 Meso  = "#6B4E3D",
                                 Micro = "#A63D2D")) +
  geom_text_repel(
    aes(label = stat),
    size = 4,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf) +
  coord_equal(xlim = c(-1, 1), 
              ylim = c(-1, 1)) +
  labs(
    x = "CV1",
    y = "CV2"
  ) +
  figure_theme

lda_plot + loadings_plot + 
  plot_annotation(tag_levels = 'A')

# Save figure
ggsave("../notebooks/figures/S4_multivar.png", 
       width = 12, 
       height = 5, 
       dpi = 300)

# =========================
# 6. Pairwise - for supp matt
# =========================

run_model_comparisons <- function(data,
                                  metrics = c("connectance",
                                              "trophic_level",
                                              "generality",
                                              "vulnerability",
                                              "NDTI",
                                              "NDCI"),
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
                           metric == "NDTI" ~ "Trophic index",
                           metric == "NDCI" ~ "Competition index"),
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
      x = NULL,
      y = "Estimated marginal mean",
      title = levs[i]
    ) +
    coord_cartesian(clip = "off") +
    figure_theme +
    theme(axis.text.x = element_text(angle=45, hjust=1))
}

plot_list_emm[[1]] /
  plot_list_emm[[2]] /
  plot_list_emm[[3]] +
  plot_layout(
    guides = "collect",
    heights = c(1, 1, 1)
  )

# Save figure
ggsave("../notebooks/figures/S4_emm_plot.png", 
       width = 9,
       height = 9, 
       dpi = 600)

# change along time

df_linear <- df %>%
  glow_up(time = as.numeric(str_extract(time, "\\d+")),
          time_num = as.numeric(time)) %>%
  pivot_longer(cols = all_of(network_stats),
               names_to = "statistic",
               values_to = "val") %>%
  squad_up(statistic, model, time_num) %>%
  no_cap(
    mean_val = mean(val, na.rm = TRUE),
    se_val = sd(val, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  squad_up(statistic, model) %>%
  glow_up(relative_change = mean_val - mean_val[time_num == 1]) %>%
  disband() %>%
  left_join(metric_lookup, by = "statistic") %>%
  glow_up(
    level = factor(level, levels = c("Macro", "Meso", "Micro"))) %>%
  glow_up(stat_label = factor(stat_label,
                              levels = c("Connectance",
                                         "Max trophic level",
                                         "Competition index",
                                         "Trophic index",
                                         "Generality",
                                         "Vulnerability")))

plot_time_change <-
  ggplot(df_linear,
         aes(x = time_num,
             y = relative_change,
             color = model,
             group = model)) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = accent_col) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.7) +
  facet_wrap(~ stat_label, 
             scales = "free_y",
             ncol = 2) +
  scale_x_continuous(breaks = 1:4,
                     labels = c("Pre-extinction", "During Extinction", "Early Recovery", "Late Recovery")) +
  scale_color_manual(values = pal_df$c, breaks = pal_df$l,
                     name = "Model") +
  scale_alpha_manual(values = c(1, 0.3, 1, 0.3, 1, 0.3), 
                     breaks = pal_df$l) +
  labs(x = NULL,
       y = "Change from Pre-extinction (Δ)") +
  figure_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

### Time change

# Ensure time is a factor for the 6x4 design
df_anova <- df %>%
  mutate(time_fact = as.factor(time))

# Function to run Two-Way ANOVA and get Effect Sizes
run_full_anova <- function(stat_name, data) {
  formula_str <- as.formula(paste0(stat_name, " ~ model * time_fact"))
  fit <- aov(formula_str, data = data)
  
  # Calculate Partial Eta Squared
  eta <- etaSquared(fit, type = 2) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(statistic = stat_name)
  
  return(eta)
}

# Run for all metrics
anova_summary_table <- map_dfr(network_stats, ~run_full_anova(.x, df_anova))

df_plot <- df_anova %>%
  pivot_longer(cols = all_of(network_stats), names_to = "statistic", values_to = "value") %>%
  group_by(statistic, model, time_fact) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val = sd(value, na.rm = TRUE),
    n = n(),
    se = sd_val / sqrt(n),
    .groups = "drop"
  ) %>%
  glow_up(statistic = case_when(statistic == "NDTI" ~ "Trophic index",
                                statistic == "NDCI" ~ "Competition index",
                                statistic == "trophic_level" ~ "Max trophic level",
                                .default = str_to_sentence(statistic))) %>%
  glow_up(level = case_when(statistic %in% c("Connectance", "Max trophic level") ~ "Macro",
                            statistic %in% c("Generality", "Vulnerability") ~ "Micro",
                            .default = "Meso"))

cv_analysis <- df_plot %>%
  group_by(statistic, time_fact) %>%
  no_cap(model_disagreement_cv = if (first(statistic) %in% c("Trophic index", "Competition index")) {
    sd(mean_val)
  } else {
    100 * sd(mean_val) / mean(mean_val)
  },
  measure = if (first(statistic) %in% c("Trophic index", "Competition index")) {
    "SD"
  } else {
    "CV (%)"
  },
  .groups = "drop") %>%
  glow_up(level = case_when(statistic %in% c("Connectance", "Max trophic level") ~ "Macro",
                            statistic %in% c("Generality", "Vulnerability") ~ "Micro",
                            .default = "Meso"))

# This includes partial eta-squared for main results
summary_data <- anova_summary_table %>%
  # Select only the relevant columns
  vibe_check(statistic, term, eta.sq.part) %>%
  # Spread the terms into their own columns
  pivot_wider(names_from = term, 
              values_from = eta.sq.part) %>%
  glow_up(across(where(is.numeric), ~round(., 3)),
          statistic = case_when(statistic == "NDTI" ~ "Trophic index",
                                statistic == "NDCI" ~ "Competition index",
                                statistic == "trophic_level" ~ "Max trophic level",
                                .default = str_to_sentence(statistic))) %>%
  # Rename columns for clarity in the paper
  rename(Metric = statistic,
         Model = model,
         Time = time_fact,
         Interaction = `model:time_fact`) %>%
  left_join(cv_analysis %>% 
              group_by(statistic) %>% 
              summarise(mean_cv = mean(model_disagreement_cv)),
            by = c("Metric" = "statistic"))

ANOVA_summary <-
  ggplot(summary_data, 
         aes(x = Model, 
             y = Time, 
             label = Metric)) +
  # Reference line: Where Model and Time are equally important
  geom_abline(slope = 1, 
              intercept = 0, 
              linetype = "dashed",
              linewidth = 0.5,
              colour = colorspace::lighten(accent_col, 0.8)) +
  # The Data Points (overlay the lines)
  geom_point(aes(size = Interaction, 
                 color = mean_cv)) +
  # Labels
  geom_text_repel(size = 3.5, box.padding = 1.2) +
  scale_colour_gradientn(colors = col_cont,
                         name = "Disagreement") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  coord_fixed() +
  scale_size_continuous(range = c(2, 8), name = "Interaction") +
  labs(x = expression(paste("Importance of model (Partial ", eta[p]^2, ")")),
       y = expression(paste("Importance of time (Partial ", eta[p]^2, ")"))) +
  figure_theme

layout <- "
A
B
"

combined_plot <- plot_time_change + ANOVA_summary + 
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = 'A')

# Save figure
ggsave("../notebooks/figures/S4_time_structure.png",
       width = 2500,
       height = 3800,
       units = "px",
       dpi = 300)
