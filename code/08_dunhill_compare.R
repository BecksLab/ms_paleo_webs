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
