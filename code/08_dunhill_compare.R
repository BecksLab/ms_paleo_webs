# libraries
library(genzplyr)
library(here)
library(mgcv)
library(patchwork)
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

network_stats <- c("connectance", "complexity", "trophic_level", "generality",
                   "vulnerability", "S1", "S2", "S4", "S5")

fit_gam_with_anova <- function(stat_name, data) {
  
  formula_shared <- as.formula(paste0(stat_name, " ~ model + s(time, k = 4)"))
  formula_model <- as.formula(paste0(stat_name, " ~ model + s(time, by = model, k = 4)"))
  
  # Fit models
  gam_shared <- gam(formula_shared, data = data, method = "REML")
  gam_model <- gam(formula_model, data = data, method = "REML")
  
  # Parametric coefficients
  param <- summary(gam_model)$p.table %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(
      statistic = stat_name,
      type = "parametric",
      estimate_or_edf = Estimate,
      se_or_refdf = `Std. Error`,
      t_or_f_value = `t value`,
      p_value = `Pr(>|t|)`
    ) %>%
    select(term, statistic, type, estimate_or_edf, se_or_refdf, t_or_f_value, p_value)
  
  # Smooth terms
  smooths <- summary(gam_model)$s.table %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(
      statistic = stat_name,
      type = "smooth",
      estimate_or_edf = edf,
      se_or_refdf = Ref.df,
      t_or_f_value = F,
      p_value = `p-value`
    ) %>%
    select(term, statistic, type, estimate_or_edf, se_or_refdf, t_or_f_value, p_value)
  
  # ANOVA comparison
  anova_res <- anova(gam_shared, gam_model, test = "F")
  anova_df <- tibble(
    term = "anova_model_comparison",
    statistic = stat_name,
    type = "anova",
    estimate_or_edf = NA,
    se_or_refdf = NA,
    t_or_f_value = anova_res$F[2],
    p_value = anova_res$`Pr(>F)`[2]
  )
  
  # Combine all
  bind_rows(param, smooths, anova_df)
}

# Loop over all network statistics
all_results <- map_dfr(network_stats, ~fit_gam_with_anova(.x, df))

# View tidy results
print(all_results)

# prep data fro plotting

df_plot <- 
  df %>%
  glow_up(model = str_replace(model, "bodymassratio", "log ratio")) %>%
  #glow_up(across(matches("S[[:digit:]]"), log)) %>%
  vibe_check(-c(distance, redundancy, complexity, diameter, n_rep, richness)) %>%
  pivot_longer(
    cols = -c(model, time),
    names_to = "stat",
    values_to = "stat_val") %>%
  # get mean values
  squad_up(model, time, stat) %>%
  no_cap(mean = mean(stat_val)) %>%
  ungroup() %>%
  # standardise names
  glow_up(stat = case_when(stat == "S1" ~ "No. of linear chains",
                           stat == "S2" ~ "No. of omnivory motifs",
                           stat == "S5" ~ "No. of apparent competition motifs",
                           stat == "S4" ~ "No. of direct competition motifs",
                           .default = as.character(stat))) %>%
  glow_up(model = factor(model, ordered = TRUE, 
                         levels = c("niche", "random", "adbm", "lmatrix", "log ratio", "pfim")),
          time = str_extract(time, "\\d+")) %>%
  glow_up(level = case_when(
    stat %in% c("complexity", "connectance", "trophic_level", "redundancy", "diameter") ~ "Macro",
    stat %in% c("generality", "vulnerability") ~ "Micro",
    .default = "Meso"
  ))

plot_list <- vector(mode = "list", length = 3)
levs = c("Macro", "Meso", "Micro")

for (i in seq_along(plot_list)) {
  
  plot_list[[i]] <- ggplot(df_plot %>% 
                             yeet(level == levs[i]),
                           aes(x = time,
                               y = mean,
                               colour = model,
                               group = model)) +
    geom_line()  +
    facet_wrap(vars(stat),
               scales = 'free',
               ncol = 2) +
    scale_size(guide = 'none') +
    xlab("Time") +
    ylab("value") +
    coord_cartesian(clip = "off") +
    scale_colour_manual(values = pal_df$c,
                        breaks = pal_df$l) +
    scale_x_discrete(labels = c("pre", "during", "early", "late")) +
    labs(title = levs[i]) +
    figure_theme
}

plot_list[[1]] / plot_list[[2]] / plot_list[[3]] +
  plot_layout(guides = 'collect') +
  plot_layout(height = c(1, 2, 1))

ggsave("../figures/time_structure.png",
       width = 5000,
       height = 6500,
       units = "px",
       dpi = 600)

# means abs differences between real and extinction simulations

df <- read_csv("../data/processed/extinction_topology.csv") %>%
  # remove metaweb pfims
  yeet(model != "pfim_metaweb") %>%
  # rename the remianing pfim col
  glow_up(model = case_when(model == "pfim_downsample" ~ "pfim",
                            .default = as.character(model))) %>%
  vibe_check(-c(rep, distance, redundancy, complexity, diameter, S1, S2, S4, S5, resilience)) %>%
  pivot_longer(
    cols = -c(model, extinction_mechanism, n_rep),
    names_to = "stat",
    values_to = "sim_val") %>%
  everyone_in_the_groupchat(.,
                            read_csv("../data/processed/topology.csv") %>%
                              # remove metaweb pfims
                              yeet(model != "pfim_metaweb") %>%
                              # rename the remianing pfim col
                              glow_up(model = case_when(model == "pfim_downsample" ~ "pfim",
                                                        .default = as.character(model))) %>%
                              yeet(time == "G2") %>%
                              vibe_check(-c(distance, redundancy, complexity, diameter, time, S1, S2, S4, S5, trophic_level)) %>%
                              pivot_longer(
                                cols = -c(model, n_rep),
                                names_to = "stat",
                                values_to = "real_val")) %>%
  glow_up(model = str_replace(model, "bodymassratio", "log ratio")) %>%
  glow_up(diff = real_val - sim_val) %>%
  squad_up(model, extinction_mechanism) %>%
  no_cap(mean_diff = abs(mean(diff, na.rm = TRUE))) %>%
  glow_up(model = factor(model, ordered = TRUE, 
                         levels = c("niche", "random", "adbm", "lmatrix", "log ratio", "pfim")))

mean_diff <- ggplot(df,
                    aes(y = extinction_mechanism,
                        x = mean_diff,
                        fill = model)) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_cartesian(clip = "off") +
  ylab(NULL)  +
  scale_fill_manual(values = pal_df$c,
                    breaks = pal_df$l) +
  figure_theme

# tss scores

df <- read_csv("../data/processed/extinction_tss.csv") %>%
  # remove metaweb pfims
  yeet(model != "pfim_metaweb") %>%
  # rename the remianing pfim col
  glow_up(model = case_when(model == "pfim_downsample" ~ "pfim",
                            .default = as.character(model)))%>%
  glow_up(model = str_replace(model, "bodymassratio", "log ratio")) %>%
  squad_up(model, extinction_mechanism) %>%
  no_cap(mean_tss = mean(tss, na.rm = TRUE))

tss <-  ggplot(df,
               aes(y = extinction_mechanism,
                   x = mean_tss,
                   fill = model)) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_cartesian(clip = "off") +
  ylab(NULL)  +
  scale_fill_manual(values = pal_df$c,
                    breaks = pal_df$l) +
  figure_theme


(tss + struct) / (mean_diff + motif)  +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A')

ggsave("../figures/dunhill_comp.png",
       width = 7000,
       height = 7500,
       units = "px",
       dpi = 600)

