# libraries
library(car)
library(emmeans)
library(genzplyr)
library(here)
library(lsr)
library(mgcv)
library(paletteer)
library(patchwork)
library(reshape2)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

#load script that determines plotting aesthetics
source("lib/plotting_theme.R")

# import simulated data

df <- read_csv("../data/processed/topology.csv") %>%
  # remove metaweb pfims
  yeet(model != "pfim_metaweb") %>%
  # rename the remianing pfim col
  glow_up(model = case_when(model == "pfim_downsample" ~ "PFIM",
                            model == "bodymassratio" ~ "Body-size ratio",
                            model == "adbm" ~ "ADBM",
                            model == "lmatrix" ~ "ATN",
                            .default = as.character(model)),
          #make tiem numeric
          time = as.numeric(str_extract(time, "\\d+")))  %>%
  glow_up(model = as.factor(model)) %>%
  na.omit()

network_stats <- c("connectance", "trophic_level", "generality",
                   "vulnerability", "S1", "S2", "S4", "S5")

# Updated function to handle the 'unique' vector error
check_anova_assumptions <- function(stat_name, data) {
  
  # 1. Ensure the grouping variables are factors and not nested in a tibble
  # We use the interaction of model and time as the grouping variable for Levene's
  grouping_var <- interaction(data$model, as.factor(data$time))
  
  # 2. Extract the metric as a numeric vector
  metric_values <- data[[stat_name]]
  
  # 3. Fit the model
  formula_str <- as.formula(paste0(stat_name, " ~ model * as.factor(time)"))
  fit <- aov(formula_str, data = data)
  
  # 4. Run Levene's Test using the vectors directly
  # This avoids the "unique() applies only to vectors" error
  lev_test <- leveneTest(metric_values ~ grouping_var)
  
  return(lev_test)
}

levene_results <- map_dfr(network_stats, ~check_anova_assumptions(.x, df))
print(levene_results)

# 1. Prepare the diagnostic data
# We iterate through each metric, fit the model, and extract residuals
diag_data <- network_stats %>%
  map_dfr(function(stat) {
    # Fit the 6x4 ANOVA model
    formula_str <- as.formula(paste0(stat, " ~ model * as.factor(time)"))
    fit <- aov(formula_str, data = df)
    
    # Extract fitted values and residuals
    data.frame(
      metric = stat,
      fitted = fitted(fit),
      resid = residuals(fit)
    )
  })

# 2. Create the Residuals vs. Fitted Plot (Homogeneity Check)
res_vs_fit_plot <- ggplot(diag_data, aes(x = fitted, y = resid)) +
  # Use high transparency (alpha) because of the large N (2400)
  geom_point(alpha = 0.1, size = 0.5, color = "midnightblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  # Facet by metric with 'free' scales since metrics have different units
  facet_wrap(~metric, scales = "free", ncol = 4) +
  labs(
    title = "Residuals vs. Fitted Plots (All Metrics)",
    subtitle = "Check for constant variance: the vertical spread should be similar across the x-axis",
    x = "Fitted Values",
    y = "Residuals"
  ) +
  figure_theme +
  theme(strip.text = element_text(face = "bold"))

# 3. Create the Q-Q Plot (Normality Check)
qq_plot <- ggplot(diag_data, aes(sample = resid)) +
  stat_qq(alpha = 0.1, size = 0.5) +
  stat_qq_line(color = "red") +
  facet_wrap(~metric, scales = "free", ncol = 4) +
  labs(
    title = "Normal Q-Q Plots (All Metrics)",
    subtitle = "Check for normality: points should follow the red diagonal line",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  figure_theme +
  theme(strip.text = element_text(face = "bold"))

# 4. Display or Save
# You can stack them or look at them individually
res_vs_fit_plot
qq_plot

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

# Calculate means and confidence intervals for plotting
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
  glow_up(statistic = case_when(statistic == "S1" ~ "No. of linear chains",
                                statistic == "S2" ~ "No. of omnivory motifs",
                                statistic == "S4" ~ "No. of direct competition motifs",
                                statistic == "S5" ~ "No. of apparent competition motifs",
                                statistic == "trophic_level" ~ "Max trophic level",
                                .default = str_to_sentence(statistic))) %>%
  glow_up(level = case_when(statistic %in% c("Connectance", "Max trophic level") ~ "Macro",
                            statistic %in% c("Generality", "Vulnerability") ~ "Micro",
                            .default = "Meso"))

plot_list_raw <- vector("list", length = 3)
levs <- c("Macro", "Meso", "Micro")

for (i in seq_along(plot_list_raw)) {
  
  plot_list_raw[[i]] <-
    ggplot(
      df_plot %>% yeet(level == levs[i]),
      aes(x = time_fact, 
          y = mean_val, 
          colour = model, 
          group = model)
    ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1) +
    # Standard Error bars to show precision
    geom_errorbar(aes(ymin = mean_val - se, ymax = mean_val + se), width = 0.1) +
    facet_wrap(
      vars(statistic),
      scales = "free_y",
      ncol = 2
    ) +
    scale_x_discrete(breaks = c(1, 2, 3, 4),
                     labels = c("pre", "during", "early", "late")) +
    scale_colour_manual(values = pal_df$c, breaks = pal_df$l) +
    labs(
      x = NULL,
      y = NULL,
      title = levs[i]
    ) +
    coord_cartesian(clip = "off") +
    figure_theme
}

raw_lines_plot <- 
  plot_list_raw[[1]] /
  plot_list_raw[[2]] /
  plot_list_raw[[3]] +
  plot_layout(
    guides = "collect",
    heights = c(1, 2, 1)
  )

ggsave("../figures/raw_time_structure.png",
       plot = raw_lines_plot,
       width = 5000,
       height = 6500,
       units = "px",
       dpi = 600)

cv_analysis <- df_plot %>%
  group_by(statistic, time_fact) %>%
  summarise(
    # How much do the 6 models disagree?
    model_disagreement_cv = (sd(mean_val) / mean(mean_val)) * 100,
    .groups = "drop"
  ) %>%
  glow_up(level = case_when(statistic %in% c("Connectance", "Max trophic level") ~ "Macro",
                            statistic %in% c("Generality", "Vulnerability") ~ "Micro",
                            .default = "Meso"))

cv_plot_list <- vector("list", length = 3)
levs <- c("Macro", "Meso", "Micro")

#get max val to standardise y axis

max_y = max(cv_analysis$model_disagreement_cv)

# 2. Loop through the levels to create the 3 stacked panels
for (i in seq_along(cv_plot_list)) {
  
  # Filter data for the current level
  plot_data <- cv_analysis %>% filter(level == levs[i])
  
  cv_plot_list[[i]] <- ggplot(plot_data, aes(x = time_fact, 
                                             y = model_disagreement_cv, 
                                             group = statistic)) +
    # Adding a subtle area fill to emphasize the 'dip' in disagreement
    geom_area(fill = "grey90", alpha = 0.5) +
    geom_line(colour = "#36393B", linewidth = 0.8) +
    geom_point(colour = "#36393B", size = 2) +
    
    # Faceting by metric within the scale (Macro/Meso/Micro)
    facet_wrap(vars(statistic), scales = "free_y", ncol = 2) +
    
    # Matching your specific x-axis formatting
    scale_x_discrete(breaks = c(1, 2, 3, 4),
                     labels = c("pre", "during", "early", "late")) +
    ylim(c(0, max_y)) + 
    # Styling to match your previous figures
    labs(
      x = NULL,
      y = "CV (%)",
      title = paste(levs[i])
    ) +
    coord_cartesian(clip = "off") +
    figure_theme
}

# 3. Assemble the final multi-panel plot
# Adjust heights based on how many metrics are in each level (Meso usually has more)
cv_final_plot <- cv_plot_list[[1]] / 
  cv_plot_list[[2]] / 
  cv_plot_list[[3]] +
  plot_layout(heights = c(1, 2, 1))

ggsave("../figures/cv.png",
       plot = cv_final_plot,
       width = 5000,
       height = 6500,
       units = "px",
       dpi = 600)


# 1. Function to get pairwise model differences per time bin
run_posthoc <- function(stat_name, data) {
  # Fit the model
  formula_str <- as.formula(paste0(stat_name, " ~ model * as.factor(time)"))
  fit <- aov(formula_str, data = data)
  
  # Get estimated marginal means for model grouped by time
  # This compares models WITHIN each time step
  posthoc <- emmeans(fit, pairwise ~ model | time)
  
  # Convert to a tidy dataframe
  results <- as.data.frame(posthoc$contrasts) %>%
    mutate(statistic = stat_name)
  
  return(results)
}

# 2. Run for all metrics
all_posthocs <- map_dfr(network_stats, ~run_posthoc(.x, df))

# 3. Filter for significant differences only
significant_comparisons <- all_posthocs %>%
  filter(p.value < 0.05) %>%
  arrange(statistic, time, desc(abs(estimate)))

# Prepare Tukey data for plotting
plot_tukey_data <- all_posthocs %>%
  glow_up(
    # Clean up comparison names (e.g., "modelA - modelB")
    contrast = str_replace_all(contrast, "model", ""),
    contrast = str_replace_all(contrast, "\\(", ""),
    contrast = str_replace_all(contrast, "\\)", ""),
    # Define significance for alpha/coloring
    sig_label = ifelse(p.value < 0.05, "Significant", "Non-Significant"),
    # Create a cleaner time label
    time_label = case_when(time == "G1" ~ "Pre",
                           time == "G2" ~ "During",
                           time == "G3" ~ "Early",
                           time == "G4" ~ "Late"
    ),
    statistic = case_when(statistic == "S1" ~ "No. of linear chains",
                          statistic == "S2" ~ "No. of omnivory motifs",
                          statistic == "S4" ~ "No. of direct competition motifs",
                          statistic == "S5" ~ "No. of apparent competition motifs",
                          statistic == "trophic_level" ~ "Max trophic level",
                          .default = str_to_sentence(statistic))
  )

ggplot(plot_tukey_data, aes(x = time_label, y = contrast, fill = estimate)) +
  geom_tile(color = "white") +
  # Use a diverging scale (Red = Positive difference, Blue = Negative)
  scale_fill_gradient2(low = "#084594", mid = "white", high = "#8C2F02", midpoint = 0) +
  # Add a border or 'star' to significant cells
  geom_text(data = filter(plot_tukey_data, p.value < 0.05), label = "*", size = 5) +
  facet_wrap(~statistic, ncol = 4) +
  labs(
    subtitle = "Asterisks (*) denote p < 0.05 (Tukey HSD). Color intensity indicates magnitude of difference.",
    x = "Extinction Phase",
    fill = "Mean Diff"
  ) +
  figure_theme

ggsave("../figures/ANOVA_tukey.png",
       width = 6000,
       height = 5000,
       units = "px",
       dpi = 600)


# 1. Export ANOVA Summary
# This includes partial eta-squared for main results
manuscript_anova_table <- anova_summary_table %>%
  # Select only the relevant columns
  select(statistic, term, eta.sq.part) %>%
  # Spread the terms into their own columns
  pivot_wider(
    names_from = term, 
    values_from = eta.sq.part
  ) %>%
  glow_up(across(where(is.numeric), ~round(., 3)),
          statistic = case_when(statistic == "S1" ~ "No. of linear chains",
                                statistic == "S2" ~ "No. of omnivory motifs",
                                statistic == "S4" ~ "No. of direct competition motifs",
                                statistic == "S5" ~ "No. of apparent competition motifs",
                                statistic == "trophic_level" ~ "Max trophic level",
                                .default = str_to_sentence(statistic))) %>%
  # Rename columns for clarity in the paper
  rename(
    Metric = statistic,
    Model = model,
    Time = time_fact,
    Interaction = `model:time_fact`
  )

write_csv(manuscript_anova_table, "../tables/ANOVA_Results.csv")

# 2. Export CV Analysis (The Convergence Table)
# Formatted to show Time Bins as columns for easier reading
cv_analysis %>%
  glow_up(time_fact = case_when(time_fact == "1" ~ "Pre extinction",
                                time_fact == "2" ~ "During extinction",
                                time_fact == "3" ~ "Early extinction",
                                time_fact == "4" ~ "Late extinction")) %>%
  pivot_wider(names_from = time_fact,
              values_from = model_disagreement_cv) %>%
  write_csv("../tables/Model_Agreement_CV.csv")


# 1. Combine ANOVA importance with CV disagreement
# We take the mean CV across all time points for a general 'disagreement' score
summary_data <- manuscript_anova_table %>%
  left_join(
    cv_analysis %>% 
      group_by(statistic) %>% 
      summarise(mean_cv = mean(model_disagreement_cv)),
    by = c("Metric" = "statistic")
  )

summary_plot <- 
  ggplot(summary_data, aes(x = Model, y = Time, label = Metric)) +
  # Reference line: Where Model and Time are equally important
  geom_abline(slope = 1, 
              intercept = 0, 
              linetype = "dashed", 
              colour = "grey90") +
  # The Data Points
  geom_point(aes(size = Interaction, color = mean_cv), alpha = 0.8) +
  # Labels
  geom_text_repel(size = 3.5, box.padding = 0.5) +
  # Color scale (Viridis for clarity)
  scale_colour_gradientn(colors = c("#F2E8CF", "#6A994E", "#154734"),
                         name = "Disagreement (CV%)") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  coord_fixed() +
  scale_size_continuous(range = c(3, 12), name = "Interaction (CV%)") +
  labs(
    x = "Importance of model (Partial Eta-Squared)",
    y = "Importance of time (Partial Eta-Squared)"
  ) +
  figure_theme

ggsave("../figures/ANOVA_summary.png",
       plot = summary_plot,
       width = 6000,
       height = 5000,
       units = "px",
       dpi = 600)
