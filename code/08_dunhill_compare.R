# libraries
library(genzplyr)
library(here)
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
  glow_up(model = case_when(model == "pfim_downsample" ~ "pfim",
                            .default = as.character(model)),
          #make tiem numeric
          time = as.numeric(str_extract(time, "\\d+")))  %>%
  glow_up(model = as.factor(model)) %>%
  na.omit()

network_stats <- c("connectance", "trophic_level", "generality",
                   "vulnerability", "S1", "S2", "S4", "S5")

# ANOVA function (to make mapping easier)
fit_anova_per_time <- function(stat_name, time_bin, data) {
  
  # Filter data for the specific time bin
  data_filtered <- data %>% filter(time == time_bin)
  
  # Fit ANOVA: stat ~ model
  formula_anova <- as.formula(paste0(stat_name, " ~ model"))
  anova_fit <- aov(formula_anova, data = data_filtered)
  
  # Run Tukey Post-hoc
  tukey_res <- TukeyHSD(anova_fit)
  
  # Tidy results
  posthoc_df <- as.data.frame(tukey_res$model) %>%
    rownames_to_column("comparison") %>%
    as_tibble() %>%
    mutate(
      statistic = as.character(stat_name),
      time_bin = as.character(time_bin),
      significant = ifelse(`p adj` < 0.05, "*", "ns")
    ) %>%
    select(statistic, time_bin, comparison, diff, lwr, upr, p_adj = `p adj`, significant)
  
  return(posthoc_df)
}

# Define your time bins based on your labels (1=pre, 2=during, 3=early, 4=late)
time_bins <- unique(df$time)

# Run the analysis across all metrics AND all time bins
all_time_results <- expand_grid(stat = network_stats, time = time_bins) %>%
  mutate(results = map2(stat, time, ~fit_anova_per_time(.x, .y, df))) %>%
  unnest(results)

# Relabel time bins for the plot
all_time_results_plot <- all_time_results %>%
  mutate(time_label = case_when(
    time_bin == "1" ~ "Pre-extinction",
    time_bin == "2" ~ "During Extinction",
    time_bin == "3" ~ "Early Recovery",
    time_bin == "4" ~ "Late Recovery"
  ))

ggplot(all_time_results_plot, aes(x = comparison, y = diff, color = significant)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = lwr, ymax = upr)) +
  # Facet by Metric (rows) and Time Bin (columns)
  facet_grid(statistic ~ time_label, scales = "free") +
  coord_flip() +
  scale_colour_manual(values = c("*" = "#E41A1C", "ns" = "#377EB8"),
                      name = NULL,
                      labels = c("Significant", "Non-Significant")) +
  labs(
    x = "Comparison",
    y = "Mean Difference"
  ) +
  figure_theme +
  theme(axis.text.y = element_text(size = 7))


# 1. Calculate Means and SE for the Linear Plot
df_linear <- df %>%
  # Ensure time is numeric (1, 2, 3, 4)
  mutate(time_num = as.numeric(time)) %>%
  pivot_longer(cols = all_of(network_stats), names_to = "statistic", values_to = "val") %>%
  group_by(statistic, model, time_num) %>%
  summarise(
    mean_val = mean(val, na.rm = TRUE),
    se_val = sd(val, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  # Calculate Change Relative to Time 1 (Pre-extinction)
  group_by(statistic, model) %>%
  mutate(relative_change = mean_val - mean_val[time_num == 1]) %>%
  ungroup() %>%
  # Labeling for the plot
  glow_up(stat_label = case_when(statistic == "S1" ~ "No. of linear chains",
                                 statistic == "S2" ~ "No. of omnivory motifs",
                                 statistic == "S4" ~ "No. of direct competition motifs",
                                 statistic == "S5" ~ "No. of apparent competition motifs",
                                 statistic == "trophic_level" ~ "Max trophic level",
                                 .default = str_to_sentence(statistic)),
          model = case_when(
            model == "pfim" ~ "PFIM",
            model == "bodymassratio" ~ "Body-size ratio",
            model == "adbm" ~ "ADBM",
            model == "lmatrix" ~ "ATN",
            TRUE ~ as.character(model)
          ))

ggplot(df_linear, aes(x = time_num, y = mean_val, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), width = 0.1) +
  facet_wrap(~stat_label, scales = "free_y") +
  scale_x_continuous(breaks = 1:4, labels = c("Pre", "Dur", "Early", "Late")) +
  scale_color_manual(values = pal_df$c, breaks = pal_df$l) +
  labs(title = "Absolute Change Over Linear Time",
       x = "Extinction Phase", y = "Mean Metric Value") +
  figure_theme

ggplot(df_linear, aes(x = time_num, y = relative_change, color = model, group = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~stat_label, scales = "free_y") +
  scale_x_continuous(breaks = 1:4, labels = c("Pre", "Dur", "Early", "Late")) +
  scale_color_manual(values = pal_df$c, breaks = pal_df$l) +
  labs(title = "Relative Change from Pre-extinction Baseline",
       subtitle = "Spreading lines = Divergence in response; Parallel lines = Similarity of differences",
       x = "Extinction Phase", y = "Change in Value (Δ)") +
  figure_theme

df_raw_plot <- 
  df %>%
  glow_up(model = case_when(model == "pfim" ~ "PFIM",
                            model == "bodymassratio" ~ "Body-size ratio",
                            model == "adbm" ~ "ADBM",
                            model == "lmatrix" ~ "ATN",
                            .default = as.character(model))) %>%
  vibe_check(-c(distance, redundancy, complexity, diameter, n_rep, richness)) %>%
  pivot_longer(cols = -c(model, time),
               names_to = "stat",
               values_to = "stat_val") %>%
  squad_up(model, time, stat) %>%
  no_cap(mean = mean(stat_val, na.rm = TRUE)) %>%
  ungroup() %>%
  glow_up(stat = case_when(stat == "S1" ~ "No. of linear chains",
                           stat == "S2" ~ "No. of omnivory motifs",
                           stat == "S4" ~ "No. of direct competition motifs",
                           stat == "S5" ~ "No. of apparent competition motifs",
                           stat == "trophic_level" ~ "Max trophic level",
                           .default = str_to_sentence(stat)),
          model = factor(model, ordered = TRUE,
                         levels = c("niche", "random", "ADBM", "ATN", "Body-size ratio", "PFIM")),
          level = case_when(stat %in% c("Connectance", "Max trophic level") ~ "Macro",
                            stat %in% c("Generality", "Vulnerability") ~ "Micro",
                            .default = "Meso"))

plot_list_raw <- vector("list", length = 3)
levs <- c("Macro", "Meso", "Micro")

for (i in seq_along(plot_list_raw)) {
  
  plot_list_raw[[i]] <-
    ggplot(
      df_raw_plot %>% yeet(level == levs[i]),
      aes(x = time, y = mean, colour = model, group = model)
    ) +
    geom_line(linewidth = 0.8) +
    facet_wrap(
      vars(stat),
      scales = "free_y",
      ncol = 2
    ) +
    scale_x_continuous(breaks = c(1, 2, 3, 4),
                       labels = c("pre", "during", "early", "late")) +
    scale_colour_manual(values = pal_df$c, breaks = pal_df$l) +
    labs(
      x = "Time",
      y = "Value",
      title = levs[i]
    ) +
    coord_cartesian(clip = "off") +
    figure_theme
}

plot_list_raw[[1]] /
  plot_list_raw[[2]] /
  plot_list_raw[[3]] +
  plot_layout(
    guides = "collect",
    heights = c(1, 2, 1)
  )

ggsave("../notebooks/figures/Figure_S3_time_structure.png",
       width = 5000,
       height = 6500,
       units = "px",
       dpi = 600)

# other results export

# Copy your original all_results
si_table <- all_results %>%
  # 1. Rename columns
  rename(
    `Network metric` = statistic,
    `Term type` = type,
    `Estimate / EDF` = estimate_or_edf,
    `Std. Error / Ref.df` = se_or_refdf,
    `t / F value` = t_or_f_value,
    `p-value` = p_value
  ) %>%
  
  # 2. Re-label terms
  mutate(
    Term = case_when(
      `Term type` == "parametric" & str_detect(term, "model") ~
        str_replace(term, "model", "Model: ") %>%
        str_replace("bodymassratio", "Body-size ratio") %>%
        str_replace("lmatrix", "ATN") %>%
        str_replace("pfim", "PFIM") %>%
        str_replace("niche", "Niche") %>%
        str_replace("random", "Random") %>%
        str_replace("adbm", "ADBM"),
      
      `Term type` == "parametric" & term == "(Intercept)" ~ "(Intercept)",
      
      `Term type` == "smooth" ~
        str_replace(term, "s\\(time\\):model", "Smooth: time by ") %>%
        str_replace("bodymassratio", "Body-size ratio") %>%
        str_replace("lmatrix", "ATN") %>%
        str_replace("pfim", "PFIM") %>%
        str_replace("niche", "Niche") %>%
        str_replace("random", "Random") %>%
        str_replace("adbm", "ADBM"),
      
      `Term type` == "anova" ~ "ANOVA: shared vs model-specific smooths",
      
      TRUE ~ term
    )
  ) %>%
  
  # 3. Round numeric columns
  mutate(
    `Estimate / EDF` = round(`Estimate / EDF`, 3),
    `Std. Error / Ref.df` = round(`Std. Error / Ref.df`, 3),
    `t / F value` = round(`t / F value`, 3),
    `p-value` = ifelse(`p-value` < 0.001, "<0.001", round(`p-value`, 3))
  ) %>%
  
  # 4. Order metrics and rows for readability
  mutate(
    `Network metric` = factor(
      `Network metric`,
      levels = c("connectance", "trophic_level",
                 "generality", "vulnerability",
                 "S1", "S2", "S4", "S5")
    ),
    `Term type` = factor(`Term type`, levels = c("parametric", "smooth", "anova"))
  ) %>%
  arrange(`Network metric`, `Term type`)

# 5. Select final column order
si_table <- si_table %>%
  select(Term, `Network metric`, `Term type`, `Estimate / EDF`,
         `Std. Error / Ref.df`, `t / F value`, `p-value`)

# 6. Export CSV
write_csv(si_table,
          "../notebooks/tables/Table_S7_GAM_results.csv")

# model compare
anova_table <-
  all_results %>%
  filter(type == "anova") %>%
  select(
    statistic,
    t_or_f_value,
    p_value
  ) %>%
  # Rename columns for clarity
  rename(
    `Network metric` = statistic,
    `F value` = t_or_f_value,
    `p-value` = p_value
  ) %>%
  
  # Add descriptive Term column
  mutate(
    Term = "ANOVA: shared vs model-specific smooths",
    
    # Round numbers for readability
    `F value` = round(`F value`, 3),
    `p-value` = ifelse(`p-value` < 0.001, "<0.001", round(`p-value`, 3))
  ) %>%
  
  # Select relevant columns in order
  select(Term, `Network metric`, `F value`, `p-value`) %>%
  
  # Order metrics to match main figure
  mutate(
    `Network metric` = factor(
      `Network metric`,
      levels = c("connectance", "trophic_level",
                 "generality", "vulnerability",
                 "S1", "S2", "S4", "S5")
    )
  ) %>%
  arrange(`Network metric`)

write_csv(
  anova_table,
  "../notebooks/tables/Table_S8_GAM_model_comparison.csv")

# means abs differences between real and extinction simulations

mad_df <- read_csv("../data/processed/extinction_topology.csv") %>%
  # remove metaweb pfims
  yeet(model != "pfim_metaweb") %>%
  # rename the remianing pfim col
  glow_up(model = case_when(model == "pfim_downsample" ~ "pfim",
                            .default = as.character(model))) %>%
  vibe_check(-c(rep, distance, redundancy, diameter, resilience, richness, complexity)) %>%
  pivot_longer(
    cols = -c(model, extinction_mechanism, n_rep),
    names_to = "stat",
    values_to = "sim_val") %>%
  full_join(.,
            read_csv("../data/processed/topology.csv") %>%
              # remove metaweb pfims
              yeet(model != "pfim_metaweb") %>%
              # rename the remianing pfim col
              glow_up(model = case_when(model == "pfim_downsample" ~ "pfim",
                                        .default = as.character(model))) %>%
              yeet(time == "G2") %>%
              vibe_check(-c(time, distance, redundancy, diameter, richness, complexity)) %>%
              pivot_longer(
                cols = -c(model, n_rep),
                names_to = "stat",
                values_to = "real_val")) %>%
  glow_up(model = case_when(model == "pfim" ~ "PFIM",
                            model == "bodymassratio" ~ "Body-size ratio",
                            model == "adbm" ~ "ADBM",
                            model == "lmatrix" ~ "ATN",
                            .default = as.character(model))) %>%
  glow_up(diff = real_val - sim_val) %>%
  squad_up(model, extinction_mechanism, stat) %>%
  no_cap(MAD = abs(mean(diff, na.rm = TRUE))) %>%
  glow_up(model = factor(model, ordered = TRUE, 
                         levels = c("niche", "random", "ADBM", "ATN", "Body-size ratio", "PFIM"))) %>%
  lowkey(scenario = extinction_mechanism, metric = stat) %>%
  rbind(.,
        read_csv("../data/processed/extinction_tss.csv")  %>%
          # remove metaweb pfims
          yeet(model != "pfim_metaweb") %>%
          vibe_check(-n_rep) %>%
          pivot_longer(!c(model, extinction_mechanism)) %>%
          # get mean (for now might be worth looking at sd as well...)
          squad_up(model, extinction_mechanism, name) %>%
          no_cap(MAD = mean(value, na.rm = TRUE)*-1) %>%
          glow_up(metric = case_when(name == "tss_link" ~ "Link",
                                     name == "tss_node" ~ "Node"),
                  model = case_when(model == "pfim_downsample" ~ "PFIM",
                                    model == "bodymassratio" ~ "Body-size ratio",
                                    model == "adbm" ~ "ADBM",
                                    model == "lmatrix" ~ "ATN",
                                    .default = as.character(model)),
                  name = NULL) %>%
          lowkey(scenario = extinction_mechanism) %>%
          glow_up(model = factor(model, ordered = TRUE, 
                                 levels = c("niche", "random", "ADBM", "ATN", "Body-size ratio", "PFIM")))
  )

metrics <- unique(mad_df$metric)

results_list <- list()
kendal_results <- data.frame()

for (met in metrics) {
  
  # Filter for this metric
  df_met <- mad_df %>% filter(metric == met)
  
  # Convert to wide format: rows=scenarios, columns=models
  wide_met <- df_met %>%
    pivot_wider(names_from = model, values_from = MAD) %>%
    vibe_check(-metric)
  
  # Rank within each model (lower MAD = better match)
  ranked <- as.data.frame(apply(wide_met[,-1], 2, rank, ties.method = "average")) %>%
    glow_up(scenario = wide_met$scenario)
  
  # Compute Kendall & Spearman correlation matrices
  kendall_corr <- cor(ranked[,-ncol(ranked)], method = "kendall")
  
  # Print Kendall correlations nicely
  cat("\nKendall rank correlation matrix:\n")
  print(kendall_corr)
  
  # Pairwise significance tests
  cat("\nPairwise Kendall correlation tests:\n")
  pairs <- combn(colnames(ranked)[-ncol(ranked)], 2)
  
  for (i in 1:ncol(pairs)) {
    m1 <- pairs[1, i]
    m2 <- pairs[2, i]
    ct <- cor.test(ranked[[m1]], ranked[[m2]], method = "kendall")
    cat("\n", m1, "vs", m2, 
        "tau =", round(ct$estimate, 3),
        "p =", round(ct$p.value, 4), "\n")
  }
  
  # Heatmap for Kendall correlations
  melted <- melt(kendall_corr) %>%
    # record network metric
    glow_up(metric = met) %>%
    # rename cols
    lowkey(Model1 = Var1, Model2 = Var2, tau = value)
  
  kendal_results <- 
    rbind(kendal_results, melted)
  
}

kendal_results <-
  kendal_results %>%
  glow_up(metric = case_when(metric == "S1" ~ "No. of linear chains",
                             metric == "S2" ~ "No. of omnivory motifs",
                             metric == "S5" ~ "No. of apparent competition motifs",
                             metric == "S4" ~ "No. of direct competition motifs",
                             metric == "trophic_level" ~ "Max trophic level",
                             .default = str_to_sentence(metric)),
          level = case_when(
            metric %in% c("Complexity", "Connectance", "Max trophic level", "Richness", "Diameter") ~ "Macro",
            metric %in% c("Generality", "Vulnerability") ~ "Micro",
            metric %in% c("Node", "Link") ~ "TSS",
            .default = "Meso"
          ))

plot_list <- vector(mode = "list", length = 4)
levs = c("Macro", "Meso", "Micro", "TSS")

for (i in seq_along(plot_list)) {
  
  plot_list[[i]] <- ggplot(kendal_results %>%
                             yeet(level == levs[i]), 
                           aes(Model1, 
                               Model2, 
                               fill = tau)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(breaks = c(-1, 0, 1),
                         colours = c("#8C2F02", "#FFFFFF", "#084594"),
                         limits = c(-1, 1)) +
    facet_wrap(vars(metric),
               #scales = 'free',
               ncol = 2) +
    labs(
      fill = "Kendall τ",
      title = levs[i],
      x = NULL,
      y = NULL
    ) +
    figure_theme +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
}

plot_list[[1]] / plot_list[[2]] / plot_list[[3]] / plot_list[[4]] +
  plot_layout(guides = 'collect') +
  plot_layout(height = c(1, 2, 1, 1))

ggsave("../figures/kendal_tau.png",
       width = 5000,
       height = 7000,
       units = "px",
       dpi = 600)
